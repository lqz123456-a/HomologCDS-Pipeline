# 项目方法文档：跨物种同源 CDS 与表达候选基因筛选

## 1. 项目定位与边界

本项目的目标不是证明某个基因导致疾病，而是把公开数据中值得进一步验证的基因按证据强弱排队。最终交付物是候选列表；实验验证、批次效应评估和正式差异分析仍是下一阶段工作。

选择“人—小鼠同源基因 + 病例/对照表达数据”有三个原因：同源关系使模型物种结果更容易解释；CDS 是可计算、可复核的序列实体；表达矩阵能为每个基因补充与表型相关的证据。两类证据合并后，比只看序列或只看表达更接近真实研发问题。

## 2. 输入与输出

输入分为三类：

1. `gene_queries.tsv`：待抓取基因、物种记录与 NCBI TaxID。表必须包含 `ensembl_gene_id` 列；当前小鼠流程应填写 `mouse_gene_id`。TaxID 可避免物种检索歧义，`gene_symbol` 与 `taxid` 是实际 NCBI 检索条件。输入表仍要求 `organism` 列以保持统一格式，但当前脚本不使用其值作为检索或输出条件；FASTA 中的物种名来自所下载 GenBank 记录。
2. `homolog_pairs.tsv`：从 Ensembl BioMart 真实导出并筛选的、人—小鼠 one-to-one 同源映射；默认键为小鼠 Ensembl Gene ID，基因符号不参与默认匹配。
3. `expression_matrix.tsv`：行是基因、列是样本、数值是已归一化的 log2 表达值，列名前缀为 `control_` 和 `case_`。

输出包括原始 GenBank、抓取清单、CDS FASTA、QC 报告、同源匹配表、候选基因表和运行摘要。结果可回查输入、accession 与 `provenance.json` 中记录的抓取时间和 API key 使用状态；当前版本不保存完整脚本参数快照。

## 3. 第一步：批量获取 CDS

`scripts/fetch_cds.py` 先用 `esearch` 在 `nuccore` 中按 gene symbol、TaxID、RefSeq 和 mRNA 查找记录，再用 `efetch` 下载 GenBank 记录。输入表必须包含来自同源表的 `ensembl_gene_id` 列；当前流程抓取小鼠记录时，它应为 `mouse_gene_id`，不能用 `NM_...` / `XM_...` 代替。单元格可留空，脚本会使用 `NA:<gene_symbol>` 作为 FASTA 首字段；该记录需在后续显式启用 `--allow-symbol-fallback` 且符号唯一，才可能匹配。选择 GenBank 而不是只下载 FASTA，是因为 CDS 坐标、gene qualifier 和产物名都在 feature 注释中；脚本据此用 Biopython 的 `feature.extract()` 提取真正 CDS，而不是误把整段 UTR+mRNA 当作编码区。

检索式限制 `refseq[filter]` 与 `biomol_mrna[PROP]`，目的是优先得到带稳定注释的转录本。真实项目还应明确“选主转录本”规则，例如 MANE Select、最长 CDS 或人工指定 accession；同一基因的多转录本不能悄悄混在一起。

### 为什么要控制请求频率

NCBI E-utilities 是共享服务。脚本无 API key 时将请求间隔设为约 0.34 秒（最多约 3 次/秒），有 key 时约为 0.11 秒；失败时按约 1.25、2.25、4.25、8.25 秒指数退避，最多尝试 5 次。这同时处理 429 限流、短暂 5xx 和网络波动。每次请求还携带邮箱，输出 `provenance.json` 记录抓取时间与是否使用 key。

这比“请求失败就立即重试”可靠：后者会放大拥堵，还容易让一批任务整体失败。当前版本会保存原始 GenBank、`fetch_manifest.tsv` 和 `provenance.json` 以供追溯，但尚未实现避免重复下载的本地缓存。

## 4. 第二步：Linux 预处理与 FASTA 质控

`shell/preprocess_fasta.sh` 统一 CRLF/LF、去除序列内部空白、转为大写并计算长度。Linux 命令行适合做这一步：流式处理文本，速度快、内存占用小，也让不规范换行或大小写问题在进入 Python 前被消除。

随后 `scripts/qc_and_match.py` 会再次将序列转为大写，并将 `U` 规范化为 `T`，然后检查：

- 是否为空序列；
- 是否出现 A/C/G/T/N 之外的字符；
- CDS 长度是否为 3 的倍数；
- 是否存在完全相同的重复序列；
- N 的比例。

这些检查能尽早暴露截断、拼接、字符编码和重复记录问题。异常记录被标成 `flag` 而非直接删除，因为异常序列有时是生物学上真实的预测转录本；保留报告让筛选阈值可审计。若做蛋白层面的严格分析，可再增加起始密码子、末端终止密码子和内部终止密码子检查。

## 5. 第三步：ID 标准化与同源匹配

不同资源中的同一基因可同时有基因符号、NCBI Gene ID、Ensembl gene/transcript ID 和 RefSeq accession；它们不是可直接互换的 ID。此前抓取结果把 RefSeq accession 放在 FASTA 首字段，而 `qc_and_match.py` 将首字段按 Ensembl Gene ID 连接 `homolog_pairs.tsv`，因此 stable-ID 匹配率为 0%。

现行格式固定为 `Ensembl_gene_id|gene_symbol|organism|RefSeq_accession|product`。抓取时通过输入表传入、并原样保留对应的小鼠 Ensembl Gene ID；质控脚本将该字段去版本号后与真实同源表的 `mouse_gene_id` 做 stable-ID 连接，再输出相应的 `human_gene_id` 与人类基因符号。RefSeq accession 仅作 NCBI 记录和转录本版本的可追溯信息。这样 Ensembl 对 Ensembl、RefSeq 对 RefSeq，命名空间不会混用。

默认不会以基因符号补齐 unmatched 记录。若处理没有 Ensembl ID 的历史 FASTA，可显式启用 `--allow-symbol-fallback`；脚本只接受在 FASTA 中唯一的小鼠符号，并标注 `symbol_fallback_unique`。有重名或歧义的符号保持 `unmatched`，必须补充权威 ID 映射后再分析。

匹配表保留 `orthology_source`，它回答“这个同源关系从哪里来”。实际展示时，建议固定一个版本化来源（Ensembl Compara、NCBI Gene 或 OMA），导出日期和版本一起写进项目日志。对一对多/多对多关系，不能任意取第一行，应根据 one-to-one 注释、蛋白覆盖度或研究问题制定规则。

## 6. 第四步：表达矩阵的快速候选排序

`scripts/rank_candidates.py` 假设输入是已经标准化的 log2 表达量：

`log2FC = mean(case) - mean(control)`。

脚本对每个基因做 Welch t 检验（不强制两组方差相等），再用 Benjamini-Hochberg 校正获得 FDR。默认候选条件为 `|log2FC| >= 1` 且 `FDR <= 0.05`，最后把结果与同源表合并并排序。

为什么不用原始计数直接做 t 检验？RNA-seq 原始 count 的均值—方差关系很强，应使用 DESeq2、edgeR 或 limma-voom；本脚本定位为“已归一化公开表达矩阵的快速筛选”。矩阵应在报告中注明数据集编号、平台、样本纳排标准、归一化方式和批次信息。示例数据只用于验证流程，不能作生物学结论。

## 7. 复现与验收

```bash
python scripts/demo_run.py
```

验收时检查：`qc_report.tsv` 是否每条记录都有状态；`matched_homologs.tsv` 是否显示匹配方法；`candidate_genes.tsv` 是否含 `log2FC`、`p_value`、`fdr_bh`、`is_candidate`；`run_summary.json` 是否声明数据属性。真实运行还应抽查一个 accession：打开保存的 GenBank，确认脚本提取的 gene qualifier 与 CDS feature 一致。

## 8. 设计决策与可扩展方向

NCBI 限流会导致批任务不稳定，因此流程实现了节流、指数退避、邮件标识与原始文件留存；RefSeq accession 与 Ensembl Gene ID 的命名空间不能混用，因此将 Ensembl ID 设为抓取输入和 FASTA 首字段，把 RefSeq 单列溯源，并默认禁用符号回退；表达数据的多重检验会带来假阳性，因此同时报告效应量和 BH-FDR，而不是只看 p 值。

下一步可扩展为 Snakemake/Nextflow 工作流、requests 缓存、Docker 环境、单元测试，以及对候选基因执行 GO/KEGG 富集或蛋白序列比对。这样能自然展示从一次性脚本到可维护分析流程的能力。
