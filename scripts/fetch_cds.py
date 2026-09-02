"""
Download RefSeq mRNA GenBank records and extract Ensembl-keyed CDS features.
Designed for polite NCBI E-utilities use: a shared throttle, retry/backoff,
and provenance records make requests both respectful and reproducible.
"""

from __future__ import annotations

import argparse, json, logging, time, io
from pathlib import Path

import pandas as pd
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
import requests


class NCBIClient:
    def __init__(self, email: str, api_key: str | None = None):
        Entrez.email, Entrez.api_key = email, api_key  # record credentials for Biopython provenance
        self.interval = 0.11 if api_key else 0.34
        self.last = 0.0

    def call(self, endpoint: str, **params) -> str:
        params.update(email=Entrez.email, tool="homolog-cds-portfolio")
        if Entrez.api_key: params["api_key"] = Entrez.api_key
        for attempt in range(5):
            time.sleep(max(0, self.interval - (time.monotonic() - self.last)))
            try:
                self.last = time.monotonic()
                url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/{endpoint}"
                response = requests.get(url, params=params, timeout=30)
                if response.status_code == 429 or response.status_code >= 500:
                    raise requests.HTTPError(f"HTTP {response.status_code}")
                response.raise_for_status()
                return response.text
            except requests.RequestException as exc:  # 429/5xx and transient socket faults
                if attempt == 4:
                    raise RuntimeError(f"NCBI request failed after retries: {exc}") from exc
                wait = min(30, 2 ** attempt + 0.25)
                logging.warning("NCBI request failed (%s); retry in %.2fs", exc, wait)
                time.sleep(wait)


def extract_cds(record, gene: str, ensembl_gene_id: str) -> SeqRecord | None:
    for feature in record.features:
        if feature.type != "CDS":
            continue
        names = feature.qualifiers.get("gene", []) + feature.qualifiers.get("gene_synonym", [])
        if any(x.upper() == gene.upper() for x in names):
            seq = feature.extract(record.seq)
            product = feature.qualifiers.get("product", [""])[0]
            organism = record.annotations.get("organism", "unknown").replace(" ", "_")
            # Keep the Ensembl gene ID in FASTA field 1: qc_and_match uses this
            # field to join Ensembl BioMart homolog_pairs.tsv.  RefSeq remains a
            # separate provenance field rather than pretending to be a gene ID.
            header = f"{ensembl_gene_id}|{gene}|{organism}|{record.id}|{product}"
            return SeqRecord(seq, id=header, description="")
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--input", required=True,
        help="TSV: query_id, ensembl_gene_id, gene_symbol, organism, taxid",
    )
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--email", required=True)
    ap.add_argument("--api-key")
    a = ap.parse_args()

    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    client = NCBIClient(a.email, a.api_key)
    queries = pd.read_csv(a.input, sep="\t")
    rows = []
    sequences = []
    required = {"query_id", "ensembl_gene_id", "gene_symbol", "organism", "taxid"}
    missing = required - set(queries.columns)
    if missing:
        ap.error("--input is missing required column(s): " + ", ".join(sorted(missing)))
    blank_mask = queries["ensembl_gene_id"].isna() | (queries["ensembl_gene_id"].astype(str).str.strip() == "")
    if blank_mask.any():
        logging.warning(
            "%d / %d rows have empty ensembl_gene_id — "
            "FASTA headers will use 'NA:<symbol>' placeholder; "
            "qc_and_match --allow-symbol-fallback handles those downstream",
            blank_mask.sum(), len(queries),
        )
    for q in queries.itertuples(index=False):
        raw = q.ensembl_gene_id
        has_id = True
        if pd.isna(raw):
            has_id = False
        elif isinstance(raw, str) and not raw.strip():
            has_id = False

        if has_id:
            eid_truth = raw.strip() if isinstance(raw, str) else raw
        else:
            eid_truth = ""

        eid_safe = eid_truth if has_id else f"NA:{q.gene_symbol}"
        term = f'"{q.gene_symbol}"[Gene Name] AND txid{q.taxid}[Organism:exp] AND refseq[filter] AND biomol_mrna[PROP]'
        raw_xml = client.call("esearch.fcgi", db="nuccore", term=term, retmax=3, retmode="xml")
        search = Entrez.read(io.BytesIO(raw_xml.encode("utf-8")))
        accession = next(iter(search["IdList"]), None)
        row = {
            "query_id": q.query_id,
            "ensembl_gene_id": eid_truth,
            "gene_symbol": q.gene_symbol,
            "taxid": q.taxid,
            "ncbi_uid": accession,
            "status": "not_found",
        }
        if accession:
            gb = client.call("efetch.fcgi", db="nuccore", id=accession, rettype="gb", retmode="text")
            gb_path = out / f"{q.query_id}_{accession}.gb"
            gb_path.write_text(gb, encoding="utf-8")
            rec = SeqIO.read(io.StringIO(gb), "genbank")
            cds = extract_cds(rec, q.gene_symbol, eid_safe)
            if cds:
                sequences.append(cds)
                row.update(status="ok", accession=rec.id, length=len(cds.seq), raw_genbank=str(gb_path))
            else:
                row["status"] = "cds_feature_not_found"
        rows.append(row)
    SeqIO.write(sequences, out / "cds_raw.fasta", "fasta")
    pd.DataFrame(rows).to_csv(out / "fetch_manifest.tsv", sep="\t", index=False)
    provenance = {
        "email": a.email,
        "api_key_used": bool(a.api_key),
        "retrieved_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    (out / "provenance.json").write_text(json.dumps(provenance, indent=2))

if __name__ == "__main__":
    main()