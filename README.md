# 跨物种同源 CDS 与表达差异候选基因筛选

这是一个可复现的生物信息学小项目。它把 **NCBI E-utilities + Biopython** 的序列获取、序列质控与 ID 对齐、以及公开表达矩阵的快速候选基因筛选串成一条可审计流水线。

> 示例对象：人（`Homo sapiens`）与小鼠（`Mus musculus`）的同源基因。仓库内的 CDS 与表达矩阵为教学数据；`data/homolog_pairs.tsv` 是从 Ensembl BioMart 真实导出并筛选的 one-to-one 同源关系，脚本也可连接 NCBI 获取真实 CDS。

## 能回答的问题

“在目标物种中，与人类已知基因同源、且在病例/对照表达矩阵中显著变化的候选基因有哪些？”

## 项目结构

```text
data/                         输入数据与示例表达矩阵
results/                      每次运行生成的结果（已忽略）
scripts/
  fetch_cds.py                NCBI 获取 CDS，并保留输入的 Ensembl 基因 ID 与 RefSeq 溯源
  qc_and_match.py             FASTA 质控、Ensembl stable-ID 同源匹配
  rank_candidates.py          表达矩阵 log2FC / Welch t 检验 / BH-FDR
  demo_run.py                 使用仓库示例数据的一键演示
shell/
  preprocess_fasta.sh         Linux 命令行批处理示例
docs/project-guide.md         中文方法文档：每步“为什么这样做”
```

## 同源表数据来源

`data/homolog_pairs.tsv` 导自 Ensembl BioMart 的 `hsapiens_gene_ensembl` 数据集，返回字段为人和小鼠的 Ensembl Gene ID、基因符号及 `Mouse homology type`。文件仅保留 BioMart 标为 `ortholog_one2one` 的 17,146 条人—小鼠关系（下载日期：2026-08-30）；每一行的 `orthology_source` 列也保留了该来源与日期。

同源关系会随 Ensembl 注释版本更新而变化。需要更新时，应以相同字段重新从 BioMart 导出，并按 `Mouse homology type = ortholog_one2one` 筛选；不要把一对多或多对多关系混入此文件。

## 快速开始（示例数据）

```bash
python -m venv .venv
source .venv/bin/activate       # Windows PowerShell: .venv\\Scripts\\Activate.ps1
pip install -r requirements.txt
python scripts/demo_run.py
```

结果写入 `results/`：`qc_report.tsv`、`matched_homologs.tsv`、`candidate_genes.tsv` 和 `run_summary.json`。

## 获取真实 CDS（需要网络）

先复制并按你的研究对象填写 `data/gene_queries.example.tsv`。输入表必须包含 `ensembl_gene_id` 列；为本流程抓取小鼠 CDS 时，应填入 `homolog_pairs.tsv` 的 `mouse_gene_id`，而不是 RefSeq accession。该列的值可留空，但脚本会把 FASTA 首字段写为 `NA:<gene_symbol>`；这种记录只有在后续显式启用 `--allow-symbol-fallback`、且小鼠基因符号在 FASTA 中唯一时才可能匹配。`gene_symbol` 和 `taxid` 用于 NCBI 检索及 CDS feature 校验；FASTA 中的 `organism` 来自所下载 GenBank 记录。输入表仍要求 `organism` 列以保持统一格式，但当前脚本不使用其值作为检索或输出条件。NCBI 要求 Entrez 请求附带邮箱；建议另设 API key（可提高速率上限）。

```bash
python scripts/fetch_cds.py \
  --input data/gene_queries.example.tsv \
  --outdir results/ncbi_fetch \
  --email your_email@example.edu \
  --api-key "$NCBI_API_KEY"
```

`fetch_cds.py` 默认无 API key 时最多约 3 请求/秒；遇到 HTTP 429/5xx 或网络异常会在约 1.25、2.25、4.25、8.25 秒后重试，最多尝试 5 次。脚本保存原始 GenBank 记录、抓取清单和 provenance 信息。输出 FASTA 的字段顺序为 `Ensembl_gene_id|gene_symbol|organism|RefSeq_accession|product`：其中小鼠 Ensembl ID 用于与同源表的 `mouse_gene_id` 做稳定 ID 匹配，RefSeq accession 用于定位 NCBI 记录。这样不会再把 `NM_...` 当作 Ensembl ID。

## Linux 预处理（可选）

在真实抓取后，可用下面的命令检查格式，再进入 Python 分析：

```bash
bash shell/preprocess_fasta.sh results/ncbi_fetch/cds_raw.fasta results/preprocessed
python scripts/qc_and_match.py --fasta results/preprocessed/cds_normalized.fasta \
  --homologs data/homolog_pairs.tsv --outdir results/qc
```

当前流程以小鼠 CDS 为匹配目标：`qc_and_match.py` 将 FASTA 首字段的小鼠 Ensembl ID 连接到真实同源表的 `mouse_gene_id`，并在结果中保留对应的人类同源基因。默认只接受 stable-ID 匹配；未匹配记录会标为 `unmatched`。若必须兼容缺少 Ensembl ID 的历史 FASTA，可显式传入 `--allow-symbol-fallback`，且仅在小鼠基因符号在 FASTA 中唯一时才匹配，结果会标为 `symbol_fallback_unique`，应人工复核。

详见 [方法文档](docs/project-guide.md)。
