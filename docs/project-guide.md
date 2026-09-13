# 项目方法文档：跨物种同源 CDS 与表达候选基因筛选

## 1. 项目定位与边界

本项目的目标不是证明某个基因导致疾病，而是把公开数据中值得进一步验证的基因
按证据强弱排队。最终交付物是候选列表；实验验证、批次效应评估和正式差异分析
仍是下一阶段工作。

当前主流程从两个原始输入开始：人—小鼠 one-to-one 同源表，以及 GSE18784
病例/对照表达矩阵。流程先生成表达候选，再仅对这些候选获取 CDS、执行质控和
同源匹配，最后合并表达与序列证据。

## 2. 输入、中间产物与输出

原始输入：

1. `data/homolog_pairs.tsv`：Ensembl BioMart 导出的同源关系。默认键为小鼠
   Ensembl Gene ID，基因符号不参与默认匹配。
2. `data/expression_matrix_GSE18784.tsv`：行为小鼠 Ensembl Gene ID、列为样本，
   数值为 GEO 已处理的双通道 log2 ratio，样本列前缀为 `control_` 和 `case_`。

GSE18784 使用 Agilent Whole Mouse Genome G4122A（GPL2872），当前矩阵包含
36 个肿瘤样本和 13 个正常肺样本。GEO 探针先经 GPL 平台注释映射到 Entrez
Gene ID，再通过 MyGene.info 映射到小鼠 Ensembl Gene ID，最后与 one-to-one
同源表取交集；多个探针映射到同一小鼠基因时按均值折叠。该矩阵直接采用 GEO
已处理值，没有从原始荧光信号重建 dye-swap 技术重复。

中间产物：

1. `pre_candidate_genes.tsv`：表达排序的完整结果，保留全部输入基因和
   `is_candidate` 标记。
2. `gene_queries.tsv`：只含 `is_candidate=True` 的小鼠查询，供 NCBI 抓取使用。
3. `candidate_homologs.tsv`：从候选排名结果提取的人—小鼠同源字段，供 QC
   匹配使用。
4. `cds_raw.fasta`、`cds_normalized.fasta`、`qc_report.tsv` 和
   `matched_homologs.tsv`：抓取、预处理、QC 与稳定 ID 匹配结果。

最终输出为 `candidate_genes.tsv`。它使用左连接保留全部表达候选；同源匹配
失败或没有可用 CDS 的行不会被删除，而是通过空字段和
`match_method=unmatched` 标记。该文件是 14 列证据汇总表，不重复保存
`pre_candidate_genes.tsv` 中的逐样本表达值。

推荐输出目录为 `results/GSE18784/`，其中分别使用 `ncbi_fetch/`、
`preprocessed/` 和 `qc/` 子目录。

## 3. 表达矩阵的快速候选排序

`scripts/rank_candidates.py` 假设输入是已经标准化的 log2 ratio：

```text
log2FC = mean(case) - mean(control)
```

脚本对每个基因做 Welch t 检验（不强制两组方差相等），再用
Benjamini-Hochberg 校正获得 FDR。默认候选条件为 `|log2FC| >= 1` 且
`FDR <= 0.05`，最后与同源表合并并按显著性排序。

当前 GSE18784 数据包含 14,827 个可匹配表达行，默认阈值下有 60 个候选。
`pre_candidate_genes.tsv` 仍保留全部排序行，下游通过
`is_candidate == True` 选择候选，因此修改阈值后无需改脚本接口。

本流程处理的是 GEO 已处理的微阵列 log2 ratio，不是 RNA-seq 原始 count。对于
原始 count，均值与方差关系很强，应使用 DESeq2、edgeR 或 limma-voom；本脚本
定位为“已归一化公开表达矩阵的快速筛选”。矩阵应在报告中注明数据集编号、
平台、样本纳排标准、归一化方式和批次信息。

SciPy 可能报告 `Precision loss occurred in moment calculation`。该提示来自
少数在样本间几乎恒定的表达行，表示其 Welch t 检验的矩计算接近数值精度
边界；它不是运行失败，也不会阻止后续输出。近恒定基因的单独 `p_value` 应结合
`log2FC`、`fdr_bh` 和原始表达分布谨慎解读，但无需因此停止整条流程。

## 4. 生成 NCBI 查询与候选同源子表

`scripts/make_gene_queries.py` 只读取 `pre_candidate_genes.tsv`：

- 只选择 `is_candidate=True`；
- 只保留 `mouse_gene_id` 非空的候选，即上一步同源匹配成功的行；
- 表达差异候选如果同源字段为空，直接排除，不进入 CDS 流程；
- 保留候选在表达排序中的顺序。

输出的 `gene_queries.tsv` 列固定为
`query_id, ensembl_gene_id, gene_symbol, organism, taxid`。对小鼠流程，
`ensembl_gene_id` 来自 `mouse_gene_id`，`gene_symbol` 来自 `mouse_symbol`，
默认 `organism=Mus musculus`、`taxid=10090`。`query_id` 默认使用
`mouse_<mouse_symbol>`。如果多个候选符号仅大小写不同，冲突组内所有
`query_id` 都追加各自的 Ensembl Gene ID，例如
`mouse_Actb_ENSMUSG00000000001`，避免 Windows 文件名冲突。

`candidate_homologs.tsv` 直接从候选行提取 `human_gene_id`、
`mouse_gene_id`、`human_symbol` 和 `mouse_symbol` 四列。后续
`qc_and_match.py` 只在这个候选子集内匹配。

## 5. 批量获取 CDS

`scripts/fetch_cds.py` 先用 `esearch` 在 `nuccore` 中按 gene symbol、
TaxID、RefSeq 和 mRNA 查找记录，再用 `efetch` 下载 GenBank 记录。选择
GenBank 而不是只下载 FASTA，是因为 CDS 坐标、gene/gene_synonym qualifier
和产物名都在 feature 注释中；脚本据此用 Biopython 的 `feature.extract()`
提取真正 CDS，而不是把整个 UTR+mRNA 当作编码区。

输出 FASTA 首字段固定为：

```text
Ensembl_gene_id|gene_symbol|organism|RefSeq_accession|product
```

其中小鼠 Ensembl ID 用于后续 stable-ID 匹配，RefSeq accession 只作为 NCBI
记录和转录本版本的溯源信息。原始 GenBank 文件统一保存在
`ncbi_fetch/genbank/`，抓取清单和 `provenance.json` 保留在
`ncbi_fetch/` 根目录。

当前实现会请求最多 3 个搜索结果，但只使用第一个结果。脚本随后检查该记录中
是否有一个 CDS feature 的 `gene` 或 `gene_synonym` 与查询符号一致；如果首个
记录不满足条件，则记为 `cds_feature_not_found`，不会尝试其余搜索结果。该
策略不是 MANE Select、最长 CDS 或人工确认的主转录本选择，正式项目应固定
转录本选择规则或人工指定 accession。

NCBI E-utilities 是共享服务。无 API key 时脚本最多约 3 请求/秒，有 key 时
约为 9 请求/秒；失败时按约 1.25、2.25、4.25、8.25 秒指数退避，最多尝试
5 次，以处理 429、短暂 5xx 和网络波动。

## 6. Linux 预处理与 FASTA 质控

`shell/preprocess_fasta.sh` 统一 CRLF/LF、去除序列内部空白、转为大写并按
60 字符折行。Linux 命令行适合流式处理这一类文本，也让不规范换行或大小写
问题在进入 Python 前被消除。

`scripts/qc_and_match.py` 再次将序列转为大写，并将 `U` 规范化为 `T`，然后
检查：

- 是否为空序列；
- 是否出现 A/C/G/T/N 之外的字符；
- CDS 长度是否为 3 的倍数；
- 是否存在完全相同的重复序列；
- N 的比例。

异常记录被标成 `flag` 而非直接删除，因为异常序列有时是生物学上真实的预测
转录本；保留报告让筛选阈值可审计。

## 7. ID 标准化与同源匹配

不同资源中的同一基因可能同时有 symbol、NCBI Gene ID、Ensembl gene/transcript
ID 和 RefSeq accession；这些标识符不能直接互换。当前 FASTA 输出同时携带
小鼠 Ensembl Gene ID 和 RefSeq accession，`qc_and_match.py` 只用前者连接
`candidate_homologs.tsv` 的 `mouse_gene_id`。

输出 `matched_homologs.tsv` 保留候选同源行和 QC 字段，并增加
`match_method`：

- `stable_id`：通过唯一 Ensembl Gene ID 匹配；
- `symbol_fallback_unique`：仅显式启用回退且符号唯一时使用；
- `unmatched`：没有可用 CDS 或无法唯一匹配。

当前主流程默认关闭 symbol fallback。若兼容缺少 Ensembl ID 的历史 FASTA，
可显式传入 `--allow-symbol-fallback`，但结果必须人工复核。

## 8. 合并表达与 CDS 证据

`scripts/combine_candidates.py` 按 `mouse_gene_id` 将
`pre_candidate_genes.tsv` 与 `matched_homologs.tsv` 左连接：

- 只按 `is_candidate=True` 选择表达候选，保持候选顺序；
- 同源匹配失败的候选同样保留，`mouse_gene_id` 为空且不会匹配到 CDS 证据；
- 加入 CDS accession、序列长度、N 比例、QC 状态/标志和 `match_method`；
- 缺少 CDS 匹配的候选保留在最终结果中，CDS/QC 字段为空，
  `match_method=unmatched`；
- 使用 `validate="many_to_one"` 要求匹配表键唯一，避免静默扩充行数。

最终 `candidate_genes.tsv` 固定输出 14 列：`human_gene_id`、
`human_symbol`、`mouse_gene_id`、`mouse_symbol`、`log2FC`、`p_value`、
`fdr_bh`、`mouse_organism`、`mouse_accession`、`mouse_sequence_length`、
`mouse_n_fraction`、`mouse_qc_status`、`mouse_qc_flags` 和 `match_method`。
所有行都是候选；逐样本表达值从 `pre_candidate_genes.tsv` 按
`mouse_gene_id` 回查。

## 9. 离线演示与完整调用

仓库中的 `data/demo_cds.fasta` 是当前抓取结果经过规范化后的固定 CDS 快照。
`scripts/demo_run.py` 使用该文件执行表达排序、查询生成、FASTA QC、同源
匹配和最终合并，但跳过 NCBI 下载与 shell 预处理：

```bash
python scripts/demo_run.py
```

结果写入 `results/demo/`，与正式 GSE18784 结果隔离。演示可以重复运行，
不会请求网络；网络抓取和无 API key/API key 限速仍需通过下面的正式流程验证。

```bash
python scripts/rank_candidates.py \
  --matrix data/expression_matrix_GSE18784.tsv \
  --homologs data/homolog_pairs.tsv \
  --out results/GSE18784/pre_candidate_genes.tsv

python scripts/make_gene_queries.py \
  --pre-candidates results/GSE18784/pre_candidate_genes.tsv \
  --queries-out results/GSE18784/gene_queries.tsv \
  --homologs-out results/GSE18784/candidate_homologs.tsv

python scripts/fetch_cds.py \
  --input results/GSE18784/gene_queries.tsv \
  --outdir results/GSE18784/ncbi_fetch \
  --email your_email@example.edu

bash shell/preprocess_fasta.sh \
  results/GSE18784/ncbi_fetch/cds_raw.fasta \
  results/GSE18784/preprocessed

python scripts/qc_and_match.py \
  --fasta results/GSE18784/preprocessed/cds_normalized.fasta \
  --homologs results/GSE18784/candidate_homologs.tsv \
  --outdir results/GSE18784/qc

python scripts/combine_candidates.py \
  --pre-candidates results/GSE18784/pre_candidate_genes.tsv \
  --matched-homologs results/GSE18784/qc/matched_homologs.tsv \
  --out results/GSE18784/candidate_genes.tsv
```

验收时检查：

- 表达排序是否保留全部输入基因并生成候选标记；
- `gene_queries.tsv` 和 `candidate_homologs.tsv` 是否只含候选；
- `qc_report.tsv` 是否每条 FASTA 记录都有状态；
- `matched_homologs.tsv` 是否每条记录都有 `match_method`；
- `candidate_genes.tsv` 是否保留全部表达候选，并恰好输出约定的 14 列；
- 抽查一个 accession：打开保存的 GenBank，确认 `gene` 或 `gene_synonym`
  与 CDS feature 一致。

## 10. 设计决策与可扩展方向

流程先筛表达候选再抓取 CDS，避免为全基因组下载不必要的序列。Ensembl 与
RefSeq 命名空间分开使用，避免 stable-ID 匹配率因 ID 混用而降为 0。最终
左连接保留未匹配候选，使抓取失败和序列缺失可审计。

实时 NCBI 记录、未锁定的依赖版本和外部 GEO 文件会让严格复现存在时间差异。
项目通过固定输入快照、记录抓取日期、保存原始 GenBank 与 manifest 来保证
流程可重跑和结果可追溯；如果要追求逐字节复现，还需要锁定全部依赖、保存
GEO 原始下载文件并固定每次 NCBI accession。

下一步可扩展为 Snakemake/Nextflow 工作流、requests 缓存、Docker 环境和
自动化测试，并对候选基因执行 GO/KEGG 富集或蛋白序列比对。
