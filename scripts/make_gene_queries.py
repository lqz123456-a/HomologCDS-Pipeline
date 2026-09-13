"""Build NCBI gene queries and candidate homolog rows from ranked candidates."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Create candidate-only NCBI queries and homolog rows."
    )
    ap.add_argument("--pre-candidates", required=True)
    ap.add_argument("--queries-out", required=True)
    ap.add_argument("--homologs-out", required=True)
    ap.add_argument("--organism", default="Mus musculus")
    ap.add_argument("--taxid", type=int, default=10090)
    a = ap.parse_args()

    pre_candidates = pd.read_csv(a.pre_candidates, sep="\t")
    candidates = pre_candidates.loc[
        pre_candidates["is_candidate"] & pre_candidates["mouse_gene_id"].notna()
    ].copy()

    symbols = candidates["mouse_symbol"].astype("string").str.strip()
    symbol_keys = symbols.str.casefold()
    duplicate_symbols = symbol_keys[symbol_keys.duplicated(keep=False)].unique()
    query_ids = "mouse_" + symbols
    if len(duplicate_symbols):
        ambiguous = symbol_keys.isin(duplicate_symbols)
        query_ids.loc[ambiguous] = (
            "mouse_"
            + symbols.loc[ambiguous]
            + "_"
            + candidates.loc[ambiguous, "mouse_gene_id"]
            .astype("string")
            .str.strip()
        )

    if query_ids.str.casefold().duplicated().any():
        raise ValueError("Generated query_id values are not unique.")

    queries = pd.DataFrame(
        {
            "query_id": query_ids,
            "ensembl_gene_id": candidates["mouse_gene_id"].astype(str),
            "gene_symbol": candidates["mouse_symbol"].astype(str),
            "organism": a.organism,
            "taxid": a.taxid,
        }
    )
    candidate_homologs = candidates[
        ["human_gene_id", "mouse_gene_id", "human_symbol", "mouse_symbol"]
    ]

    for output in (Path(a.queries_out), Path(a.homologs_out)):
        output.parent.mkdir(parents=True, exist_ok=True)
    queries.to_csv(a.queries_out, sep="\t", index=False)
    candidate_homologs.to_csv(a.homologs_out, sep="\t", index=False)

    print(
        f"Wrote {len(queries)} gene queries and "
        f"{len(candidate_homologs)} candidate homolog rows."
    )


if __name__ == "__main__":
    main()
