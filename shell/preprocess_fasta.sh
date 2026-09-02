#!/usr/bin/env bash
set -euo pipefail
infile="$1"; outdir="$2"
mkdir -p "$outdir"
# 统一换行、序列大小写，并用 awk 折行到 60 字符；保留 FASTA 标题以便后续 ID 追踪。
tr -d '\r' < "$infile" | awk '
BEGIN { seq="" }
/^>/ {if (seq != "") wrap(seq); print; seq=""; next}
{gsub(/[[:space:]]/, ""); seq = seq toupper($0)}
END { if (seq != "") wrap(seq) }
function wrap(s,   i) {for (i = 1; i <= length(s); i += 60) print substr(s, i, 60)}
' > "$outdir/cds_normalized.fasta"

awk '
/^>/{if(n){print id"\t"n}; id=$0;n=0;next}
{n+=length($0)}
END{if(n)print id"\t"n}
' "$outdir/cds_normalized.fasta" > "$outdir/sequence_lengths.tsv"