---
title: RNA-seq DESeq2 Web App
colorFrom: green
colorTo: blue
sdk: docker
app_port: 8788
pinned: false
---

# RNA-seq Raw Counts Web App

## 公网开放应用

这个目录现在可以作为一个独立 Web App 部署到云平台。部署成功后，会得到一个 `https://...` 的公网网址，别人不需要连接你的电脑，也不需要你的电脑一直开着。

推荐用 Docker 部署，文件已经准备好：

- `Dockerfile`：安装 R / Bioconductor / Python 环境并启动网页。
- `render.yaml`：Render 一键部署配置。
- `fly.toml`：Fly.io 部署配置。
- `PUBLIC_DEPLOY.md`：公网部署步骤。

最简单的方式是把 `outputs/rna_seq_webapp` 这个目录上传到 GitHub，然后用 Render 或 Fly.io 创建 Web Service。详细步骤见 `PUBLIC_DEPLOY.md`。

本机运行：

```bash
/Users/zhangyu/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3 server.py
```

打开终端里显示的地址，例如 `http://127.0.0.1:8787`。

不要直接双击打开 `index.html`。如果地址栏是 `file:///.../index.html`，上传接口不可用；请使用 `http://127.0.0.1:8788` 或终端打印的局域网地址。

局域网运行：

```bash
HOST=0.0.0.0 PORT=8788 /Users/zhangyu/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3 server.py
```

也可以双击 `run_network.command`。终端会显示 `LAN access: http://你的电脑IP:8788`，同一个 Wi-Fi 或同一个内网里的其他设备可以打开这个地址。

注意：局域网模式会允许同一网络内的设备访问和上传文件；公网部署建议放在受控服务器上，并加登录/HTTPS。

## 当前支持的流程

- 支持拖拽 Excel/CSV 文件。
- 支持分批添加文件；每次拖拽或选择文件后点击 `添加文件`，只是加入文件队列。
- 点击 `读取队列` 后，可以在样本表里手动指定每个样本属于 `batch1`、`batch2` 或其他 batch。
- 同一个文件里的样本也可以被分到不同 batch，多个文件也可以归到同一个 batch。
- 可以在样本表里手动填写每个样本的 group，再选择 control 和 KO/treatment contrasts。
- 同一批数据里可以有多个 KO / treatment group。
- 选择一个 control group 后，网页会自动对每个选中的 KO 做 `KO vs control`。
- Volcano plot 使用显著 protein-coding genes 的 `gene_name` 做 label。
- Heatmap 使用 raw counts 经过 DESeq2 VST 后的显著 protein-coding DEGs union，并用 `gene_name` 作为行名；TPM 相关热图逻辑不在网页中执行。
- GO 和 KEGG 分开成独立视图，且各自按 `contrast + UP/DOWN` 分开显示。
- GO 和 KEGG 图使用 dotplot：横轴为 GeneRatio，点大小为 Count，颜色为 `-log10(padj)`。
- 可以导出 Excel：summary、recipe、`KO1_all` / `KO2_all` 等当前 gene_type filter 下的 DE 结果、all_DESeq_results、GO、KEGG 会分 sheet 保存。
- `代码` 页面显示网页实际采用的 R Markdown 风格分析代码；左侧参数变化后可刷新代码预览。

## Excel 输入格式

第一个 sheet 或名为 `counts` 的 sheet 放 raw counts：

- 必需列：`gene`
- 推荐列：`gene_name`、`gene_type`
- 样本列：raw counts 数值列

可选 sheet：

- `metadata`：列名包含 `sample` 和 `condition`，用于自动分组。

GO/KEGG 可在网页里选择 `Human` 或 `Mouse`：

- Human 使用 `org.Hs.eg.db` 和 KEGG `hsa`。
- Mouse 使用 `org.Mm.eg.db` 和 KEGG `mmu`。

建议提供：

- `gene_name`：对应物种的 gene symbol，例如 human 的 `TP53`、mouse 的 `Trp53`。
- 或 `entrezgene_id`：NCBI Entrez Gene ID。

如果没有 `metadata`，网页会从样本名推断 group，例如 `KO1_1` 会推断为 `KO1`。

## 示例文件

- `sample_batch1_counts.xlsx`
- `sample_batch2_counts.xlsx`

这两个文件模拟你的 R 代码结构：batch1 和 batch2 各有 `control`、`KO1`、`KO2`。

## Analysis recipe

`代码` 页面里的 JSON recipe 会真正参与分析和作图。可修改项包括：

- `deg.padj_cutoff`
- `deg.log2fc_cutoff`
- `deg.protein_coding_only_for_labels`
- `deg.protein_coding_only_for_heatmap`
- `deg.protein_coding_only_for_enrichment`
- `volcano.xlim`
- `volcano.max_labels`
- `volcano.colors`
- `heatmap.top_genes`：仅保留为兼容字段；当前 raw-count/VST combined heatmap 按你的代码使用所有显著 protein-coding DEGs。
- `heatmap.clip`
- `heatmap.group_colors`
- `enrichment.min_overlap`
- `enrichment.max_terms`

默认配色按你的 R 代码：

- UP: `#D62728`
- DOWN: `#1F77B4`
- control: `#4D4D4D`
- KO1: `#D62728`
- KO2: `#1F77B4`

## 统计说明

网页逻辑参考你的 Rmd：

- 多个 batch 取 common genes 后合并。
- 使用 raw counts 构建 `DESeqDataSetFromMatrix(countData = count_mat, colData = colData, design = ~ batch + group)`。
- 过滤 `rowSums(counts(dds)) >= 10` 后执行 `DESeq(dds)`。
- 每个 KO contrast 使用 `lfcShrink(type = "apeglm")`。
- 每个 KO 分别对 control 做 contrast。
- DEG 阈值默认是 `padj < 0.05` 且 `|log2FoldChange| > 1`。
- Volcano label 过滤 `regulation != NS`、`gene_type == protein_coding`、`gene_name` 非空。
- Heatmap 使用显著 protein-coding DEGs 的 union，并用 `vst(dds, blind = FALSE)` 和 `limma::removeBatchEffect()`；小数据集若 `vst()` 不适用，会自动退回 `varianceStabilizingTransformation()`。
- GO 使用所选物种的 OrgDb 执行 `enrichGO(..., ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)`。
- KEGG 使用所选物种执行 `enrichKEGG(..., organism = "hsa" 或 "mmu", keyType = "ncbi-geneid")`。

当前后端是真正的 R/Bioconductor 分析流程，由 `deseq_pipeline.R` 执行 `DESeq2`、`apeglm`、`limma`、`clusterProfiler`、`org.Hs.eg.db` 和 `org.Mm.eg.db`。
