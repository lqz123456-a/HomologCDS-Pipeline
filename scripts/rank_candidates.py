"""Rank genes from a log-scale public expression matrix; exploratory, not DESeq2."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import ttest_ind

def bh(p):
    """Benjamini-Hochberg FDR correction."""
    x = np.asarray(p, dtype=float)
    order = np.argsort(x)
    n = len(x)
    # BH formula: adjusted[i] = min over j>=i of (p[j] * n / j)
    ranked = x[order] * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(adjusted, 1.0)
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--matrix',required=True)
    ap.add_argument('--homologs',required=True)
    ap.add_argument('--out',required=True)
    ap.add_argument('--min-abs-log2fc',type=float,default=1)
    ap.add_argument('--fdr',type=float,default=0.05)
    a = ap.parse_args()

    df = pd.read_csv(a.matrix,sep='\t')
    controls = [c for c in df if c.startswith('control_')]
    cases = [c for c in df if c.startswith('case_')]

    if len(controls)<2 or len(cases)<2:
        raise ValueError('Each group needs at least 2 samples.')

    p_values = []
    for _, row in df.iterrows():
        t = ttest_ind(
            row[cases].to_numpy(dtype=float),
            row[controls].to_numpy(dtype=float),
            equal_var=False,
        )
        p_values.append(t.pvalue)
    df['log2FC'] = df[cases].mean(axis=1)-df[controls].mean(axis=1)
    df["p_value"] = p_values
    df['fdr_bh'] = bh(df.p_value.fillna(1))
    df['is_candidate'] = (df.log2FC.abs()>=a.min_abs_log2fc)&(df.fdr_bh<=a.fdr)

    h = pd.read_csv(a.homologs,sep='\t')[['human_gene_id','mouse_gene_id','human_symbol','mouse_symbol']]
    result = df.merge(
        h, left_on="gene_id", right_on="mouse_gene_id", how="left",
    ).sort_values(
        ["is_candidate", "fdr_bh", "log2FC"],
        ascending=[False, True, False],
    )
    Path(a.out).parent.mkdir(parents=True,exist_ok=True)
    result.to_csv(a.out,sep='\t',index=False)

if __name__=='__main__':
    main()
