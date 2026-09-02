"""FASTA QC and Ensembl stable-ID homolog matching."""

from __future__ import annotations

import argparse, re
from pathlib import Path

import pandas as pd
from Bio import SeqIO

VALID = set("ACGTN")

def parse_header(header: str):
    p = header.split("|")
    return {
        "gene_id": re.sub(r"\.\d+$", "", p[0]),
        "gene_symbol": p[1] if len(p) > 1 else "",
        "organism": p[2] if len(p) > 2 else "",
        "accession": p[3] if len(p) > 3 else "",
    }

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta",required=True)
    ap.add_argument("--homologs",required=True)
    ap.add_argument("--outdir",required=True)
    ap.add_argument(
        "--allow-symbol-fallback", action="store_true",
        help="Use symbols only when each symbol maps to exactly one QC record; disabled by default."
    )
    a = ap.parse_args()

    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    seen = set()
    rows = []

    for r in SeqIO.parse(a.fasta,"fasta"):
        seq = str(r.seq).upper().replace("U","T")
        meta = parse_header(r.id)
        issues = []
        if not seq:
            issues.append("empty")
        if set(seq)-VALID:
            issues.append("invalid_base")
        if len(seq)%3:
            issues.append("not_in_frame")
        if seq in seen:
            issues.append("duplicate_sequence")
        seen.add(seq)
        n_fraction = round(seq.count("N") / max(len(seq), 1), 4)
        rows.append({
            **meta,
            "sequence_length": len(seq),
            "n_fraction": n_fraction,
            "qc_status": "pass" if not issues else "flag",
            "qc_flags": ";".join(issues),
        })

    qc = pd.DataFrame(rows)
    qc.to_csv(out/"qc_report.tsv",sep="\t",index=False)
    h = pd.read_csv(a.homologs,sep="\t")
    # Field 1 of fetch_cds.py output is an Ensembl gene ID, so this is a
    # namespace-consistent Ensembl-to-Ensembl join.
    matched = h.merge(
        qc.add_prefix("mouse_"),
        left_on="mouse_gene_id",
        right_on="mouse_gene_id",
        how="left",
    )
    matched["match_method"] = "stable_id"

    missing = matched["mouse_sequence_length"].isna()
    if a.allow_symbol_fallback and missing.any():
        # Never choose an arbitrary duplicate symbol; retain an unmatched row.
        unique_symbols = qc[qc["gene_symbol"].ne("") & ~qc["gene_symbol"].duplicated(False)]
        fallback = h.loc[missing].merge(
            unique_symbols.add_prefix("mouse_"),
            left_on="mouse_symbol",
            right_on="mouse_gene_symbol",
            how="left",
        )
        for col in fallback.columns:
            if col in matched: matched.loc[missing,col] = fallback[col].to_numpy()
        used = missing & matched["mouse_sequence_length"].notna()
        matched.loc[used,"match_method"] = "symbol_fallback_unique"

    matched.loc[matched["mouse_sequence_length"].isna(), "match_method"] = "unmatched"
    matched.to_csv(out/"matched_homologs.tsv",sep="\t",index=False)
    print(f"QC records: {len(qc)}; matched homolog rows: {matched.mouse_sequence_length.notna().sum()}")

if __name__=="__main__":
    main()
