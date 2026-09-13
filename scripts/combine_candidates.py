"""Combine expression candidates with matched homolog CDS evidence."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


OUTPUT_COLUMNS = [
    "human_gene_id",
    "human_symbol",
    "mouse_gene_id",
    "mouse_symbol",
    "log2FC",
    "p_value",
    "fdr_bh",
    "mouse_organism",
    "mouse_accession",
    "mouse_sequence_length",
    "mouse_n_fraction",
    "mouse_qc_status",
    "mouse_qc_flags",
    "match_method",
]


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Left-join expression candidates with CDS homolog matches."
    )
    ap.add_argument("--pre-candidates", required=True)
    ap.add_argument("--matched-homologs", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    pre_candidates = pd.read_csv(args.pre_candidates, sep="\t")
    candidates = pre_candidates.loc[pre_candidates["is_candidate"]].copy()

    matched = pd.read_csv(args.matched_homologs, sep="\t")
    duplicate_columns = [
        column
        for column in matched.columns
        if column in candidates.columns and column != "mouse_gene_id"
    ]
    matched = matched.drop(columns=duplicate_columns)
    result = candidates.merge(
        matched,
        on="mouse_gene_id",
        how="left",
        validate="many_to_one",
    )

    if "match_method" not in result.columns:
        result["match_method"] = "unmatched"
    else:
        result["match_method"] = result["match_method"].fillna("unmatched")
        blank_method = result["match_method"].astype("string").str.strip() == ""
        result.loc[blank_method, "match_method"] = "unmatched"

    result = result[OUTPUT_COLUMNS]
    output = Path(args.out)
    output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output, sep="\t", index=False)

    matched_count = result["match_method"].ne("unmatched").sum()
    print(
        f"Wrote {len(result)} candidates; "
        f"{matched_count} have matched homolog evidence."
    )


if __name__ == "__main__":
    main()
