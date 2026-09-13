# 跨物种同源 CDS 与表达差异候选基因筛选

这是一个可重跑、可追溯的生物信息学小项目。它把公开表达矩阵的快速候选筛选、
**NCBI E-utilities + Biopython** 的 CDS 获取、序列质控与 Ensembl stable-ID
同源匹配串成一条可审计流水线。

> 示例对象：人（`Homo sapiens`）与小鼠（`Mus musculus`）。当前主流程使用
> `data/homolog_pairs.tsv` 和 GSE18784 表达矩阵；CDS 从 NCBI 实时获取。

## 能回答的问题

“在目标物种中，与人类已知基因同源、且在病例/对照表达矩阵中显著变化的候选
基因有哪些，哪些候选能够进一步匹配到可用 CDS？”

## 输入与总体流程

当前流程只需要两个原始输入：

1. `data/homolog_pairs.tsv`：Ensembl BioMart 导出的人—小鼠 one-to-one 同源表。
2. `data/expression_matrix_GSE18784.tsv`：已标准化的 log2 表达矩阵，样本列以
   `case_` 和 `control_` 开头。

调用顺序固定为：

```text
rank_candidates.py
  -> make_gene_queries.py
  -> fetch_cds.py
  -> preprocess_fasta.sh
  -> qc_and_match.py
  -> combine_candidates.py
```

`rank_candidates.py` 输出全部表达排序结果；`make_gene_queries.py` 只把
`is_candidate=True` 且已匹配完整同源表的候选传给 CDS 获取。最终合并脚本会
保留全部表达候选，包括同源匹配或 CDS 匹配失败的记录，并将失败项标记为
`unmatched`。

## 环境准备与离线演示

在项目根目录创建虚拟环境并安装依赖：

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
```

当前验证环境为 Python 3.14.6；`requirements.txt` 记录最低依赖版本，尚未提供
完全锁定的环境文件。

`demo_run.py` 使用仓库内已经规范化的 `data/demo_cds.fasta` 快照，执行表达
排序、查询生成、FASTA QC、同源匹配和最终合并。它不会访问 NCBI，也不会执行
`preprocess_fasta.sh`：

```bash
python scripts/demo_run.py
```

演示输出保存在 `results/demo/`，不会覆盖正式的 GSE18784 结果。

运行 `rank_candidates.py` 时，SciPy 可能输出
`Precision loss occurred in moment calculation` 警告。这是少数基因在样本间
表达值几乎相同时的数值精度提示，不是脚本错误，不会中断流程或造成结果文件
缺失。对这些近恒定基因，其单独 `p_value` 的数值稳定性相对较弱，应结合
`log2FC`、`fdr_bh` 和原始表达分布一起解读；这不会阻止流程完成。

## 项目结构

```text
data/                         原始输入数据
results/                      每次运行生成的结果（已忽略）
scripts/
  rank_candidates.py          表达矩阵 log2FC / Welch t 检验 / BH-FDR
  make_gene_queries.py        生成候选 NCBI 查询与候选同源子表
  fetch_cds.py                NCBI 获取 CDS，并保留 Ensembl ID 与 RefSeq 溯源
  qc_and_match.py             FASTA 质控与 Ensembl stable-ID 同源匹配
  combine_candidates.py       合并表达候选与 CDS/QC 证据
  demo_run.py                 预抓取 FASTA 离线演示：跳过 NCBI 和 shell 预处理
shell/
  preprocess_fasta.sh         Linux/Git Bash 批处理：换行、大小写和序列折行
docs/project-guide.md         中文方法文档：每一步的输入、输出与设计理由
```

## 同源表数据来源

`data/homolog_pairs.tsv` 导自 Ensembl BioMart 的
`hsapiens_gene_ensembl` 数据集，返回字段为人和小鼠的 Ensembl Gene ID、
基因符号及 `Mouse homology type`。文件仅保留 BioMart 标为
`ortholog_one2one` 的 17,146 条人—小鼠关系（下载日期：2026-08-30）；
每一行的 `orthology_source` 列也保留了该来源与日期。

同源关系会随 Ensembl 注释版本更新而变化。需要更新时，应以相同字段重新从
BioMart 导出，并按 `Mouse homology type = ortholog_one2one` 筛选；不要把
一对多或多对多关系混入此文件。

## GSE18784 表达矩阵来源

`data/expression_matrix_GSE18784.tsv` 来自 GEO 数据集
[GSE18784](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18784)，
平台为 Agilent Whole Mouse Genome G4122A（GPL2872）。当前矩阵使用 36 个
肿瘤样本和 13 个正常肺样本，样本列分别整理为 `case_*` 和 `control_*`。
GSE18784 是双通道实验，矩阵中的值是 GEO 已处理的相对 universal reference
的 log2 ratio；项目直接使用这些处理值，没有从原始荧光信号重新计算技术重复。

GEO 探针经过 GPL 平台注释映射到 Entrez Gene ID，再通过 MyGene.info 映射到
小鼠 Ensembl Gene ID，最后与 `data/homolog_pairs.tsv` 的 `mouse_gene_id`
取交集。映射到同一小鼠基因的多个探针按均值折叠，因此矩阵保留了 14,827 个
one-to-one 同源基因。正式复算应重新核对 GEO 样本元数据、批次和技术重复，
并考虑使用原始数据重新处理。

## GSE18784 主流程

以下命令在项目根目录运行。GSE18784 的肿瘤组包含各 18 个 KRAS 和 HER2
模型样本，正常组包含 13 个样本，已分别整理为 `case_*`、`control_*` 列。

### 1. 表达差异候选排序

```bash
python scripts/rank_candidates.py \
  --matrix data/expression_matrix_GSE18784.tsv \
  --homologs data/homolog_pairs.tsv \
  --out results/GSE18784/pre_candidate_genes.tsv
```

当前数据默认产生 14,827 行排序结果，其中 60 行的 `is_candidate=True`。
阈值可通过 `--min-abs-log2fc` 和 `--fdr` 调整。

### 2. 生成候选查询与候选同源子表

```bash
python scripts/make_gene_queries.py \
  --pre-candidates results/GSE18784/pre_candidate_genes.tsv \
  --queries-out results/GSE18784/gene_queries.tsv \
  --homologs-out results/GSE18784/candidate_homologs.tsv
```

`gene_queries.tsv` 只含当前 60 个候选，字段为
`query_id, ensembl_gene_id, gene_symbol, organism, taxid`。其中
`ensembl_gene_id` 使用小鼠 `mouse_gene_id`，`query_id` 默认格式为
`mouse_<mouse_symbol>`，例如 `mouse_Celf4`。如果多个候选的小鼠符号仅大小写
不同，则冲突组内所有 ID 都追加 Ensembl Gene ID，例如
`mouse_Actb_ENSMUSG00000000001`，以保证唯一。`candidate_homologs.tsv`
直接由候选排名结果提取 `human_gene_id`、`mouse_gene_id`、`human_symbol` 和
`mouse_symbol` 四列。脚本先保留 `is_candidate=True` 的行，再要求
`mouse_gene_id` 非空；表达差异候选如果在上一步没有匹配到完整同源表，会直接
被排除。CDS 是否成功匹配仍由下一阶段 `qc_and_match.py` 判断。

### 3. 从 NCBI 获取 CDS

```bash
python scripts/fetch_cds.py \
  --input results/GSE18784/gene_queries.tsv \
  --outdir results/GSE18784/ncbi_fetch \
  --email your_email@example.edu
```

如有 NCBI API key，可追加：

```bash
  --api-key "$NCBI_API_KEY"
```

NCBI 要求请求附带邮箱；API key 可提高速率上限。脚本无 key 时最多约
3 请求/秒，遇到 429/5xx 或网络异常时按约 1.25、2.25、4.25、8.25 秒重试，
最多尝试 5 次。输出 FASTA 字段顺序为
`Ensembl_gene_id|gene_symbol|organism|RefSeq_accession|product`。

当前脚本最多取得 3 个搜索结果，但只下载并检查第一个结果；如果该记录的 CDS
feature 无法与 `gene` 或 `gene_synonym` 匹配，会记为
`cds_feature_not_found`，不会继续尝试后续结果。因此这里的选择规则不等同于
MANE Select、最长 CDS 或人工确认的主转录本。

### 4. Linux 预处理与 FASTA QC

`preprocess_fasta.sh` 需要 Bash，在 Windows 上可使用 Git Bash 或 WSL：

```bash
bash shell/preprocess_fasta.sh \
  results/GSE18784/ncbi_fetch/cds_raw.fasta \
  results/GSE18784/preprocessed

python scripts/qc_and_match.py \
  --fasta results/GSE18784/preprocessed/cds_normalized.fasta \
  --homologs results/GSE18784/candidate_homologs.tsv \
  --outdir results/GSE18784/qc
```

预处理统一换行、大小写、空白和 60 字符折行。`qc_and_match.py` 检查空序列、
非法碱基、CDS 长度是否为 3 的倍数、重复序列和 N 比例，再以小鼠 Ensembl
Gene ID 连接 `candidate_homologs.tsv`。默认只接受 stable-ID 匹配，不使用
基因符号回退。

### 5. 合并最终候选证据

```bash
python scripts/combine_candidates.py \
  --pre-candidates results/GSE18784/pre_candidate_genes.tsv \
  --matched-homologs results/GSE18784/qc/matched_homologs.tsv \
  --out results/GSE18784/candidate_genes.tsv
```

最终文件只按 `is_candidate=True` 选择表达候选，再使用左连接补充 CDS/QC
证据。同源表匹配失败或 CDS 抓取/匹配失败的候选都会保留，相关字段为空并将
`match_method` 标为 `unmatched`。

## 输出结构

```text
results/GSE18784/
  pre_candidate_genes.tsv
  gene_queries.tsv
  candidate_homologs.tsv
  ncbi_fetch/
    genbank/
      <query_id>_<ncbi_uid>.gb
    cds_raw.fasta
    fetch_manifest.tsv
    provenance.json
  preprocessed/
    cds_normalized.fasta
    sequence_lengths.tsv
  qc/
    qc_report.tsv
    matched_homologs.tsv
  candidate_genes.tsv
```

`candidate_genes.tsv` 是 14 列证据汇总表，字段为 `human_gene_id`、
`human_symbol`、`mouse_gene_id`、`mouse_symbol`、`log2FC`、`p_value`、
`fdr_bh`、`mouse_organism`、`mouse_accession`、`mouse_sequence_length`、
`mouse_n_fraction`、`mouse_qc_status`、`mouse_qc_flags` 和 `match_method`。
所有行都是表达候选；逐样本表达值仍保留在 `pre_candidate_genes.tsv`，可用
`mouse_gene_id` 回查。

## 验收

- `pre_candidate_genes.tsv` 的每条记录都有 `is_candidate`，且全部输入基因
  都保留在排序结果中。
- `gene_queries.tsv` 与 `candidate_homologs.tsv` 只含候选，行数和
  `mouse_gene_id` 顺序一致，查询 ID 唯一。
- `qc_report.tsv` 每条 FASTA 记录都有 QC 状态。
- `matched_homologs.tsv` 每条候选都有 `match_method`，包括 `unmatched`。
- `candidate_genes.tsv` 使用左连接保留全部表达候选，抽查任一 accession
  可回查保存的 GenBank 记录，并确认 `gene` 或 `gene_synonym` 与 CDS
  feature 一致。

严格复现还受实时 NCBI 数据库、未锁定的依赖版本和 GEO 原始下载文件是否
保留影响。当前项目保存输入快照、访问日期和抓取记录，保证流程可重跑和结果
可追溯，但不承诺任意时间重跑都得到逐字节相同的结果。

详见 [方法文档](docs/project-guide.md)。
