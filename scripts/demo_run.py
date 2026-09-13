"""Run the offline pipeline demo using a normalized, pre-fetched CDS snapshot."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/"results"/"demo"


def run(script: str, *args: str | Path) -> None:
    subprocess.run(
        [
            sys.executable,
            str(ROOT/"scripts"/script),
            *(str(arg) for arg in args),
        ],
        cwd=ROOT,
        check=True,
    )


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    expression_matrix = ROOT/"data"/"expression_matrix_GSE18784.tsv"
    homologs = ROOT/"data"/"homolog_pairs.tsv"
    cds_fasta = ROOT/"data"/"demo_cds.fasta"
    pre_candidates = OUT/"pre_candidate_genes.tsv"
    gene_queries = OUT/"gene_queries.tsv"
    candidate_homologs = OUT/"candidate_homologs.tsv"
    qc_dir = OUT/"qc"
    matched_homologs = qc_dir/"matched_homologs.tsv"
    final_candidates = OUT/"candidate_genes.tsv"

    run(
        "rank_candidates.py",
        "--matrix", expression_matrix,
        "--homologs", homologs,
        "--out", pre_candidates,
    )
    run(
        "make_gene_queries.py",
        "--pre-candidates", pre_candidates,
        "--queries-out", gene_queries,
        "--homologs-out", candidate_homologs,
    )
    run(
        "qc_and_match.py",
        "--fasta", cds_fasta,
        "--homologs", candidate_homologs,
        "--outdir", qc_dir,
    )
    run(
        "combine_candidates.py",
        "--pre-candidates", pre_candidates,
        "--matched-homologs", matched_homologs,
        "--out", final_candidates,
    )

    print("Offline demo complete. See results/demo/.")


if __name__ == "__main__":
    main()