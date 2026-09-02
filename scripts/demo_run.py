from pathlib import Path
import subprocess, sys, json

ROOT=Path(__file__).resolve().parents[1]
out=ROOT/'results'
out.mkdir(exist_ok=True)

def run(*args):
    subprocess.run([sys.executable,*args],cwd=ROOT,check=True)

run(
    "scripts/qc_and_match.py",
    "--fasta", "data/demo_cds.fasta",
    "--homologs", "data/homolog_pairs.tsv",
    "--outdir", "results",
)

run(
    "scripts/rank_candidates.py",
    "--matrix", "data/expression_matrix.tsv",
    "--homologs", "data/homolog_pairs.tsv",
    "--out", "results/candidate_genes.tsv",
)

summary = {
    "mode": "offline demo",
    "note": "Expression values are illustrative log2-scale values, not a biological conclusion."
}

(out / "run_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

print("Done. See results/.")