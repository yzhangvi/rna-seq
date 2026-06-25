const DEFAULT_RECIPE = {
  deg: {
    padj_cutoff: 0.05,
    log2fc_cutoff: 1.0,
    gene_types: ["protein_coding"],
    protein_coding_only_for_labels: true,
    protein_coding_only_for_heatmap: true,
    protein_coding_only_for_enrichment: true,
  },
  volcano: {
    xlim: [-7, 7],
    max_labels: 20,
    colors: { UP: "#D62728", DOWN: "#1F77B4", NS: "#B3B3B3" },
  },
  heatmap: {
    top_genes: 50,
    clip: 2.0,
    colors: { low: "#1F77B4", mid: "#F8FAFC", high: "#D62728" },
    group_palette: ["#4D4D4D", "#D62728", "#1F77B4", "#047C7B", "#B7791F", "#28724F", "#8B4F9F"],
    group_colors: { control: "#4D4D4D", KO1: "#D62728", KO2: "#1F77B4" },
  },
  enrichment: {
    min_overlap: 2,
    max_terms: 80,
  },
};

const state = {
  batches: [],
  pendingFiles: [],
  inspection: null,
  results: null,
  geneTypes: ["protein_coding"],
  activeContrast: null,
  volcanoPoints: [],
  recipe: structuredClone(DEFAULT_RECIPE),
};

const $ = (id) => document.getElementById(id);
const API_BASE = String(window.RNASEQ_API_BASE || localStorage.getItem("rnaseqApiBase") || "").replace(/\/$/, "");

function apiUrl(path) {
  return `${API_BASE}${path}`;
}

function makeId() {
  if (window.crypto?.randomUUID) return window.crypto.randomUUID();
  return `id_${Date.now()}_${Math.random().toString(36).slice(2, 10)}`;
}

const fmt = (value, digits = 3) => {
  if (value === null || value === undefined || Number.isNaN(Number(value))) return "-";
  const number = Number(value);
  if (number !== 0 && Math.abs(number) < 0.001) return number.toExponential(2);
  return number.toLocaleString("zh-CN", { maximumSignificantDigits: digits });
};

function setStatus(message, tone = "neutral") {
  const el = $("status");
  el.textContent = message;
  el.style.borderColor = tone === "error" ? "#c43c4a" : "#d9e0ea";
  el.style.color = tone === "error" ? "#a62f3b" : "#637083";
}

function analysisCodeText(recipe = state.recipe) {
  const padj = recipe.deg?.padj_cutoff ?? 0.05;
  const lfc = recipe.deg?.log2fc_cutoff ?? 1;
  const species = $("speciesSelect")?.value || "human";
  const orgDb = species === "mouse" ? "org.Mm.eg.db" : "org.Hs.eg.db";
  const kegg = species === "mouse" ? "mmu" : "hsa";
  const speciesLabel = species === "mouse" ? "Mouse" : "Human";
  return `---
title: "DEseq volcano, heatmap, KEGG & GO"
output: html_document
---

\`\`\`{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
\`\`\`

\`\`\`{r}
library(data.table)
library(DESeq2)
library(apeglm)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(openxlsx)
library(limma)
library(clusterProfiler)
library(enrichplot)
library(${orgDb})
library(AnnotationDbi)
\`\`\`

\`\`\`{r}
# 网页中这一步由上传文件生成：
# counts1 <- fread("batch1/annotated_counts.tsv")
# counts2 <- fread("batch2/annotated_counts.tsv")

gene_anno <- unique(
  rbind(
    counts1[, .(gene, gene_name, gene_type)],
    counts2[, .(gene, gene_name, gene_type)]
  )
)

count_cols1 <- setdiff(colnames(counts1), c("gene","gene_name","gene_type"))
count_cols2 <- setdiff(colnames(counts2), c("gene","gene_name","gene_type"))

count_mat1 <- as.matrix(counts1[, ..count_cols1])
rownames(count_mat1) <- counts1$gene

count_mat2 <- as.matrix(counts2[, ..count_cols2])
rownames(count_mat2) <- counts2$gene

common_genes <- intersect(rownames(count_mat1), rownames(count_mat2))
count_mat <- cbind(
  count_mat1[common_genes, ],
  count_mat2[common_genes, ]
)
\`\`\`

\`\`\`{r}
# 网页中这一步来自用户在 sample table 中选择的 batch/group。
group <- factor(c("control","KO1","KO2","control","KO1","KO2"),
                levels = c("control","KO1","KO2"))
batch <- factor(c(rep("batch1",3), rep("batch2",3)))

colData <- data.frame(
  row.names = colnames(count_mat),
  group = group,
  batch = batch
)
\`\`\`

\`\`\`{r}
dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData   = colData,
  design    = ~ batch + group
)

dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)
\`\`\`

\`\`\`{r}
run_deseq_contrast <- function(dds, coef_name, label, gene_anno) {
  res <- lfcShrink(dds, coef = coef_name, type = "apeglm")
  df <- as.data.frame(res)
  df$gene <- rownames(df)
  df <- df[!is.na(df$padj), ]
  df <- left_join(df, gene_anno, by = "gene")

  df$regulation <- "NS"
  df$regulation[df$padj < ${padj} & df$log2FoldChange >  ${lfc}] <- "UP"
  df$regulation[df$padj < ${padj} & df$log2FoldChange < -${lfc}] <- "DOWN"

  list(
    label = label,
    results = df,
    deg_up   = df$gene_name[df$regulation == "UP"],
    deg_down = df$gene_name[df$regulation == "DOWN"]
  )
}

KO1 <- run_deseq_contrast(dds, "group_KO1_vs_control", "KO1", gene_anno)
KO2 <- run_deseq_contrast(dds, "group_KO2_vs_control", "KO2", gene_anno)
\`\`\`

\`\`\`{r}
plot_volcano <- function(df, title) {
  ggplot(df, aes(log2FoldChange, -log10(padj))) +
    geom_point(color="grey70", size=1) +
    geom_point(data = subset(df, regulation=="UP"), color="#D62728", size=1.5) +
    geom_point(data = subset(df, regulation=="DOWN"), color="#1F77B4", size=1.5) +
    geom_text_repel(
      data = subset(df, regulation!="NS" & gene_name!="" & gene_type == "protein_coding"),
      aes(label = gene_name),
      size=3, max.overlaps=20
    ) +
    geom_vline(xintercept=c(-${lfc},${lfc}), linetype="dashed") +
    geom_hline(yintercept=-log10(${padj}), linetype="dashed") +
    coord_cartesian(xlim = c(-7, 7)) +
    scale_x_continuous(expand = expansion(mult = 0.05)) +
    theme_classic() +
    labs(title=title, x="log2 Fold Change", y="-log10(FDR)")
}

plot_volcano(KO1$results, "KO1 vs Control (DESeq2)")
plot_volcano(KO2$results, "KO2 vs Control (DESeq2)")
\`\`\`

\`\`\`{r}
# raw-count / VST heatmap only; TPM-related code is intentionally omitted.
deg_genes <- bind_rows(KO1$results, KO2$results) %>%
  dplyr::filter(regulation != "NS", gene_type == "protein_coding") %>%
  dplyr::pull(gene) %>%
  unique()

mat <- assay(vsd)[deg_genes, ]
mat_corrected <- limma::removeBatchEffect(mat, batch = colData(dds)$batch)

gene_name_map <- gene_anno$gene_name[match(rownames(mat), gene_anno$gene)]
gene_name_map[is.na(gene_name_map) | gene_name_map == ""] <- rownames(mat)
rownames(mat_corrected) <- make.unique(gene_name_map)

mat_scaled <- t(scale(t(mat_corrected)))
mat_scaled[mat_scaled > 2] <- 2
mat_scaled[mat_scaled < -2] <- -2

annotation_col <- as.data.frame(colData(dds)[, "group", drop = FALSE])
colnames(annotation_col) <- "Group"
rownames(annotation_col) <- colnames(mat_scaled)

sample_order <- order(colData(dds)$group)
mat_scaled <- mat_scaled[, sample_order]
annotation_col <- annotation_col[sample_order, , drop = FALSE]

heat_colors <- colorRampPalette(rev(brewer.pal(n = 11, "RdBu")))(100)
ann_colors <- list(Group = c(control="#4D4d4D", KO1="#D62728", KO2="#1F77B4"))

pheatmap(
  mat_scaled,
  color = heat_colors,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  show_rownames = TRUE,
  fontsize_row = 6,
  border_color = NA,
  main = "KO1 & KO2 vs Control (protein-coding DEGs)"
)
\`\`\`

\`\`\`{r}
# Export DESeq2 results used by GO / KEGG.
write.xlsx(
  list(
    KO1_all = KO1$results %>% filter(gene_type == "protein_coding"),
    KO2_all = KO2$results %>% filter(gene_type == "protein_coding")
  ),
  file = "DESeq2_all_results.xlsx",
  rowNames = FALSE
)
\`\`\`

\`\`\`{r}
map_entrez <- function(df) {
  df$entrezgene_id <- mapIds(
    ${orgDb},
    keys = df$gene_name,
    column = "ENTREZID",
    keytype = "SYMBOL",
    multiVals = "first"
  )
  df %>%
    tidyr::drop_na(entrezgene_id) %>%
    dplyr::distinct(entrezgene_id, .keep_all = TRUE)
}

KO1_data <- map_entrez(KO1$results)
KO2_data <- map_entrez(KO2$results)

KO1_up <- KO1_data %>% filter(log2FoldChange > 0, padj < 0.05, gene_type == "protein_coding")
KO1_down <- KO1_data %>% filter(log2FoldChange < 0, padj < 0.05, gene_type == "protein_coding")
KO2_up <- KO2_data %>% filter(log2FoldChange > 0, padj < 0.05, gene_type == "protein_coding")
KO2_down <- KO2_data %>% filter(log2FoldChange < 0, padj < 0.1, gene_type == "protein_coding")
\`\`\`

\`\`\`{r}
plot_kegg <- function(genes, title) {
  kegg <- enrichKEGG(
    gene = as.character(genes),
    organism = "${kegg}",
    keyType = "ncbi-geneid",
    pvalueCutoff = 0.05
  )
  dotplot(kegg, showCategory = 10) + ggtitle(title)
}

plot_kegg(KO1_up$entrezgene_id, "KO1 upregulated KEGG (${speciesLabel})")
plot_kegg(KO1_down$entrezgene_id, "KO1 downregulated KEGG (${speciesLabel})")
plot_kegg(KO2_up$entrezgene_id, "KO2 upregulated KEGG (${speciesLabel})")
plot_kegg(KO2_down$entrezgene_id, "KO2 downregulated KEGG (${speciesLabel})")
\`\`\`

\`\`\`{r}
plot_go <- function(genes, title) {
  go <- enrichGO(
    gene = genes,
    OrgDb = ${orgDb},
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
  )
  dotplot(go, showCategory = 10) + ggtitle(title)
}

plot_go(KO1_up$entrezgene_id, "KO1 upregulated GO (BP)")
plot_go(KO1_down$entrezgene_id, "KO1 downregulated GO (BP)")
plot_go(KO2_up$entrezgene_id, "KO2 upregulated GO (BP)")
plot_go(KO2_down$entrezgene_id, "KO2 downregulated GO (BP)")
\`\`\`
`;
}

function selectedCodeModules() {
  const checked = Array.from(document.querySelectorAll('input[name="codeModule"]:checked')).map((input) => input.value);
  return new Set(checked);
}

function rQuote(value) {
  return `"${String(value).replaceAll("\\", "\\\\").replaceAll('"', '\\"')}"`;
}

function rVector(values) {
  return `c(${values.map(rQuote).join(", ")})`;
}

function rColor(value, fallback) {
  return rQuote(value || fallback);
}

function selectedGeneTypes() {
  return Array.from(document.querySelectorAll('input[name="geneType"]:checked')).map((input) => input.value);
}

function geneTypeFilterCode(dataName = "df") {
  return `filter_gene_type(${dataName})`;
}

function recipeColors(recipe = state.recipe) {
  return {
    up: recipe.volcano?.colors?.UP || DEFAULT_RECIPE.volcano.colors.UP,
    down: recipe.volcano?.colors?.DOWN || DEFAULT_RECIPE.volcano.colors.DOWN,
    ns: recipe.volcano?.colors?.NS || DEFAULT_RECIPE.volcano.colors.NS,
    heatLow: recipe.heatmap?.colors?.low || DEFAULT_RECIPE.heatmap.colors.low,
    heatMid: recipe.heatmap?.colors?.mid || DEFAULT_RECIPE.heatmap.colors.mid,
    heatHigh: recipe.heatmap?.colors?.high || DEFAULT_RECIPE.heatmap.colors.high,
    group1: recipe.heatmap?.group_palette?.[0] || DEFAULT_RECIPE.heatmap.group_palette[0],
    group2: recipe.heatmap?.group_palette?.[1] || DEFAULT_RECIPE.heatmap.group_palette[1],
    group3: recipe.heatmap?.group_palette?.[2] || DEFAULT_RECIPE.heatmap.group_palette[2],
  };
}

function makeSafeName(value) {
  const cleaned = String(value)
    .replace(/[^A-Za-z0-9_]+/g, "_")
    .replace(/^_+|_+$/g, "");
  if (!cleaned) return "group_x";
  return /^[A-Za-z]/.test(cleaned) ? cleaned : `group_${cleaned}`;
}

function currentGroupPlan() {
  const meta = currentSampleMeta().filter((row) => row.sample && row.group);
  if (meta.length) {
    const control = $("controlGroup")?.value || meta[0].group;
    const ko = selectedKoGroups();
    return {
      groups: meta.map((row) => row.group),
      batches: meta.map((row) => row.batch || "batch1"),
      control,
      ko: ko.length ? ko : Array.from(new Set(meta.map((row) => row.group))).filter((group) => group !== control),
    };
  }
  return {
    groups: ["control", "KO1", "KO2", "control", "KO1", "KO2"],
    batches: ["batch1", "batch1", "batch1", "batch2", "batch2", "batch2"],
    control: "control",
    ko: ["KO1", "KO2"],
  };
}

function analysisCodeText(recipe = state.recipe) {
  const padj = recipe.deg?.padj_cutoff ?? 0.05;
  const lfc = recipe.deg?.log2fc_cutoff ?? 1;
  const species = $("speciesSelect")?.value || "human";
  const orgDb = species === "mouse" ? "org.Mm.eg.db" : "org.Hs.eg.db";
  const kegg = species === "mouse" ? "mmu" : "hsa";
  const speciesLabel = species === "mouse" ? "Mouse" : "Human";
  const modules = selectedCodeModules();
  if (modules.size === 0) {
    return "# 目前没有选择要显示的代码模块。点击“显示全部”，或勾选 DESeq2 / Volcano / Heatmap / KEGG / GO。";
  }
  const plan = currentGroupPlan();
  const geneTypes = recipe.deg?.gene_types || [];
  const geneTypeLabel = geneTypes.length ? geneTypes.join(" / ") : "all gene types";
  const colors = recipeColors(recipe);
  const koObjects = plan.ko.map((group) => makeSafeName(group));
  const contrastCalls = plan.ko
    .map((group, index) => `${koObjects[index]} <- run_deseq_contrast(dds, ${rQuote(`group_${group}_vs_${plan.control}`)}, ${rQuote(group)}, gene_anno)`)
    .join("\n");
  const resultObjects = koObjects.map((name) => `${name}$results`).join(", ");
  const base = [`---
title: "DEseq volcano, heatmap, KEGG & GO"
output: html_document
---

\`\`\`{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
\`\`\`

\`\`\`{r}
library(data.table)
library(DESeq2)
library(apeglm)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(openxlsx)
library(limma)
library(clusterProfiler)
library(enrichplot)
library(${orgDb})
library(AnnotationDbi)
\`\`\`

\`\`\`{r}
# 网页中这一步由上传文件生成：
# counts1 <- fread("batch1/annotated_counts.tsv")
# counts2 <- fread("batch2/annotated_counts.tsv")

gene_anno <- unique(rbind(
  counts1[, .(gene, gene_name, gene_type)],
  counts2[, .(gene, gene_name, gene_type)]
))

count_cols1 <- setdiff(colnames(counts1), c("gene","gene_name","gene_type"))
count_cols2 <- setdiff(colnames(counts2), c("gene","gene_name","gene_type"))

count_mat1 <- as.matrix(counts1[, ..count_cols1])
rownames(count_mat1) <- counts1$gene
count_mat2 <- as.matrix(counts2[, ..count_cols2])
rownames(count_mat2) <- counts2$gene

common_genes <- intersect(rownames(count_mat1), rownames(count_mat2))
count_mat <- cbind(count_mat1[common_genes, ], count_mat2[common_genes, ])
\`\`\`

\`\`\`{r}
# 这里反映当前网页 sample table 中的分组；例如用户输入 1/2/3 时，
# 1 就是当前 control，2/3 就是当前 KO/treatment，不会再额外保留 control/KO1/KO2。
group <- factor(${rVector(plan.groups)}, levels = ${rVector([plan.control, ...plan.ko])})
batch <- factor(${rVector(plan.batches)})

colData <- data.frame(row.names = colnames(count_mat), group = group, batch = batch)
\`\`\`

\`\`\`{r}
dds <- DESeqDataSetFromMatrix(countData = count_mat, colData = colData, design = ~ batch + group)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)
\`\`\`

\`\`\`{r}
selected_gene_types <- ${geneTypes.length ? rVector(geneTypes) : "character(0)"}
filter_gene_type <- function(df) {
  if (!length(selected_gene_types) || !"gene_type" %in% colnames(df)) return(df)
  df %>% dplyr::filter(gene_type %in% selected_gene_types)
}

run_deseq_contrast <- function(dds, coef_name, label, gene_anno) {
  res <- lfcShrink(dds, coef = coef_name, type = "apeglm")
  df <- as.data.frame(res)
  df$gene <- rownames(df)
  df <- df[!is.na(df$padj), ]
  df <- left_join(df, gene_anno, by = "gene")
  df$regulation <- "NS"
  df$regulation[df$padj < ${padj} & df$log2FoldChange >  ${lfc}] <- "UP"
  df$regulation[df$padj < ${padj} & df$log2FoldChange < -${lfc}] <- "DOWN"
  list(label = label, results = df)
}

${contrastCalls}
\`\`\``];

  if (modules.has("volcano")) {
    base.push(`\`\`\`{r}
plot_volcano <- function(df, title) {
  ggplot(df, aes(log2FoldChange, -log10(padj))) +
    geom_point(color=${rColor(colors.ns, "#B3B3B3")}, size=1) +
    geom_point(data = subset(df, regulation=="UP"), color=${rColor(colors.up, "#D62728")}, size=1.5) +
    geom_point(data = subset(df, regulation=="DOWN"), color=${rColor(colors.down, "#1F77B4")}, size=1.5) +
    geom_text_repel(
      data = filter_gene_type(subset(df, regulation!="NS" & gene_name!="")),
      aes(label = gene_name),
      size=3, max.overlaps=20
    ) +
    geom_vline(xintercept=c(-${lfc}, ${lfc}), linetype="dashed") +
    geom_hline(yintercept=-log10(${padj}), linetype="dashed") +
    coord_cartesian(xlim = c(-7, 7)) +
    scale_x_continuous(expand = expansion(mult = 0.05)) +
    theme_classic() +
    labs(title=title, x="log2 Fold Change", y="-log10(FDR)")
}

${koObjects.map((name, index) => `plot_volcano(${name}$results, "${plan.ko[index]} vs ${plan.control} (DESeq2)")`).join("\n")}
\`\`\``);
  }

  if (modules.has("heatmap")) {
    base.push(`\`\`\`{r}
deg_genes <- bind_rows(${resultObjects}) %>%
  dplyr::filter(regulation != "NS") %>%
  filter_gene_type() %>%
  dplyr::pull(gene) %>%
  unique()

mat <- assay(vsd)[deg_genes, ]
mat_corrected <- limma::removeBatchEffect(mat, batch = colData(dds)$batch)

gene_name_map <- gene_anno$gene_name[match(rownames(mat), gene_anno$gene)]
gene_name_map[is.na(gene_name_map) | gene_name_map == ""] <- rownames(mat)
rownames(mat_corrected) <- make.unique(gene_name_map)

mat_scaled <- t(scale(t(mat_corrected)))
mat_scaled[mat_scaled > 2] <- 2
mat_scaled[mat_scaled < -2] <- -2

annotation_col <- as.data.frame(colData(dds)[, "group", drop = FALSE])
colnames(annotation_col) <- "Group"
rownames(annotation_col) <- colnames(mat_scaled)

sample_order <- order(colData(dds)$group)
mat_scaled <- mat_scaled[, sample_order]
annotation_col <- annotation_col[sample_order, , drop = FALSE]

heat_colors <- colorRampPalette(c(${rColor(colors.heatLow, "#1F77B4")}, ${rColor(colors.heatMid, "#F8FAFC")}, ${rColor(colors.heatHigh, "#D62728")}))(100)
palette_values <- c(${rVector([colors.group1, colors.group2, colors.group3, "#047C7B", "#B7791F", "#28724F", "#8B4F9F"])})
ann_colors <- list(Group = setNames(palette_values[seq_along(levels(group))], levels(group)))

pheatmap(
  mat_scaled,
  color = heat_colors,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  show_rownames = TRUE,
  fontsize_row = 6,
  border_color = NA,
  main = "${plan.ko.join(" & ")} vs ${plan.control} (${geneTypeLabel} DEGs)"
)
\`\`\``);
  }

  if (modules.has("kegg") || modules.has("go")) {
    base.push(`\`\`\`{r}
map_entrez <- function(df) {
  df$entrezgene_id <- mapIds(
    ${orgDb},
    keys = df$gene_name,
    column = "ENTREZID",
    keytype = "SYMBOL",
    multiVals = "first"
  )
  df %>% tidyr::drop_na(entrezgene_id) %>% dplyr::distinct(entrezgene_id, .keep_all = TRUE)
}

${koObjects.map((name) => `${name}_data <- map_entrez(${name}$results)`).join("\n")}
\`\`\``);
  }

  if (modules.has("kegg")) {
    base.push(`\`\`\`{r}
plot_kegg <- function(genes, title) {
  kegg <- enrichKEGG(gene = as.character(genes), organism = "${kegg}", keyType = "ncbi-geneid", pvalueCutoff = 0.05)
  dotplot(kegg, showCategory = 10) + ggtitle(title)
}

${koObjects.map((name, index) => {
  const group = plan.ko[index];
  const cutoff = group === "KO2" ? 0.1 : 0.05;
  return `${name}_up <- ${name}_data %>% filter(log2FoldChange > 0, padj < 0.05) %>% filter_gene_type()
${name}_down <- ${name}_data %>% filter(log2FoldChange < 0, padj < ${cutoff}) %>% filter_gene_type()
plot_kegg(${name}_up$entrezgene_id, "${group} upregulated KEGG (${speciesLabel})")
plot_kegg(${name}_down$entrezgene_id, "${group} downregulated KEGG (${speciesLabel})")`;
}).join("\n\n")}
\`\`\``);
  }

  if (modules.has("go")) {
    base.push(`\`\`\`{r}
plot_go <- function(genes, title) {
  go <- enrichGO(
    gene = genes,
    OrgDb = ${orgDb},
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
  )
  dotplot(go, showCategory = 10) + ggtitle(title)
}

${koObjects.map((name, index) => {
  const group = plan.ko[index];
  return `if (!exists("${name}_up")) ${name}_up <- ${name}_data %>% filter(log2FoldChange > 0, padj < 0.05) %>% filter_gene_type()
if (!exists("${name}_down")) ${name}_down <- ${name}_data %>% filter(log2FoldChange < 0, padj < 0.05) %>% filter_gene_type()
plot_go(${name}_up$entrezgene_id, "${group} upregulated GO (BP)")
plot_go(${name}_down$entrezgene_id, "${group} downregulated GO (BP)")`;
}).join("\n\n")}
\`\`\``);
  }

  return base.join("\n\n");
}

function resetRecipe() {
  state.recipe = structuredClone(DEFAULT_RECIPE);
  renderGeneTypeControls(state.geneTypes, state.recipe.deg.gene_types);
  syncRecipeInputs();
  $("analysisCode").value = analysisCodeText();
  setStatus("已恢复默认参数并刷新 R 代码");
}

function syncRecipeInputs() {
  $("fcCutoff").value = state.recipe.deg?.log2fc_cutoff ?? 1;
  $("padjCutoff").value = state.recipe.deg?.padj_cutoff ?? 0.05;
  const colors = recipeColors(state.recipe);
  if ($("colorUp")) $("colorUp").value = colors.up;
  if ($("colorDown")) $("colorDown").value = colors.down;
  if ($("colorNs")) $("colorNs").value = colors.ns;
  if ($("colorHeatLow")) $("colorHeatLow").value = colors.heatLow;
  if ($("colorHeatMid")) $("colorHeatMid").value = colors.heatMid;
  if ($("colorHeatHigh")) $("colorHeatHigh").value = colors.heatHigh;
  if ($("colorGroup1")) $("colorGroup1").value = colors.group1;
  if ($("colorGroup2")) $("colorGroup2").value = colors.group2;
  if ($("colorGroup3")) $("colorGroup3").value = colors.group3;
}

function syncInputsToRecipe() {
  state.recipe.deg = state.recipe.deg || {};
  state.recipe.volcano = state.recipe.volcano || {};
  state.recipe.heatmap = state.recipe.heatmap || {};
  state.recipe.deg.log2fc_cutoff = Number($("fcCutoff").value || 1);
  state.recipe.deg.padj_cutoff = Number($("padjCutoff").value || 0.05);
  state.recipe.deg.gene_types = selectedGeneTypes();
  const onlyProteinCoding = state.recipe.deg.gene_types.length === 1 && state.recipe.deg.gene_types[0] === "protein_coding";
  state.recipe.deg.protein_coding_only_for_labels = onlyProteinCoding;
  state.recipe.deg.protein_coding_only_for_heatmap = onlyProteinCoding;
  state.recipe.deg.protein_coding_only_for_enrichment = onlyProteinCoding;
  state.recipe.volcano.colors = {
    UP: $("colorUp")?.value || DEFAULT_RECIPE.volcano.colors.UP,
    DOWN: $("colorDown")?.value || DEFAULT_RECIPE.volcano.colors.DOWN,
    NS: $("colorNs")?.value || DEFAULT_RECIPE.volcano.colors.NS,
  };
  state.recipe.heatmap.colors = {
    low: $("colorHeatLow")?.value || DEFAULT_RECIPE.heatmap.colors.low,
    mid: $("colorHeatMid")?.value || DEFAULT_RECIPE.heatmap.colors.mid,
    high: $("colorHeatHigh")?.value || DEFAULT_RECIPE.heatmap.colors.high,
  };
  const groupPalette = [
    $("colorGroup1")?.value || DEFAULT_RECIPE.heatmap.group_colors.control,
    $("colorGroup2")?.value || DEFAULT_RECIPE.heatmap.group_colors.KO1,
    $("colorGroup3")?.value || DEFAULT_RECIPE.heatmap.group_colors.KO2,
    "#047C7B",
    "#B7791F",
    "#28724F",
    "#8B4F9F",
  ];
  state.recipe.heatmap.group_palette = groupPalette;
  const plan = currentGroupPlan();
  state.recipe.heatmap.group_colors = Object.fromEntries([plan.control, ...plan.ko].map((group, index) => [group, groupPalette[index % groupPalette.length]]));
  $("analysisCode").value = analysisCodeText();
}

function flattenBatches() {
  const files = [];
  const overrides = [];
  state.batches.forEach((batch) => {
    batch.files.forEach((file) => {
      files.push(file);
      overrides.push({ filename: file.name, batch: batch.name });
    });
  });
  return { files, overrides };
}

function formDataWithBatches(extra = {}) {
  const { files, overrides } = flattenBatches();
  if (!files.length) throw new Error("请先添加一个或多个文件。");
  const form = new FormData();
  files.forEach((file) => form.append("file", file));
  const params = { ...(extra.params || {}), batch_overrides: overrides };
  form.append("params", JSON.stringify(params));
  return form;
}

async function apiPost(path, form) {
  let response;
  try {
    response = await fetch(apiUrl(path), { method: "POST", body: form });
  } catch (error) {
    const hint = API_BASE
      ? `无法连接后端 ${API_BASE}。请确认后端服务在线，且已经允许跨域访问。`
      : "无法连接后端 API。请确认公网部署不是纯静态网页，并且前端和 Python/R 后端部署在同一个服务里。";
    throw new Error(`${hint} 浏览器原始错误：${error.message}`);
  }
  const contentType = response.headers.get("content-type") || "";
  if (!contentType.includes("application/json")) {
    const text = await response.text();
    const shortText = text.replace(/\s+/g, " ").slice(0, 180);
    throw new Error(`服务器没有返回可读取的数据（HTTP ${response.status}）。${shortText || "请查看云平台日志。"}`);
  }
  const payload = await response.json();
  if (!payload.ok) throw new Error(payload.error || "请求失败。");
  return payload;
}

function defaultBatchName() {
  return $("batchName").value.trim() || `set${state.batches.length + 1}`;
}

function addBatch(files, name = defaultBatchName()) {
  const cleanFiles = Array.from(files || []).filter((file) => file.name);
  if (!cleanFiles.length) {
    setStatus("没有可添加的文件", "error");
    return;
  }
  const start = state.batches.length;
  cleanFiles.forEach((file, index) => {
    const batchName = cleanFiles.length === 1 ? name : `${name}_${index + 1}`;
    state.batches.push({ id: makeId(), name: batchName || `set${start + index + 1}`, files: [file] });
  });
  state.inspection = null;
  state.pendingFiles = [];
  $("fileInput").value = "";
  $("sampleGrid").innerHTML = "";
  $("batchName").value = `set${state.batches.length + 1}`;
  renderBatchList();
  setStatus(`已加入队列：${cleanFiles.length} 个文件`);
}

function renderBatchList() {
  const list = $("batchList");
  list.innerHTML = "";
  state.batches.forEach((batch) => {
    const item = document.createElement("div");
    item.className = "batch-item";
    const header = document.createElement("div");
    header.className = "batch-header";
    const text = document.createElement("div");
    const title = document.createElement("strong");
    title.textContent = batch.name;
    const meta = document.createElement("span");
    meta.textContent = `${batch.files.length} 个文件`;
    text.append(title, meta);
    const remove = document.createElement("button");
    remove.className = "icon-button";
    remove.type = "button";
    remove.textContent = "×";
    remove.title = "移除整个文件组";
    remove.addEventListener("click", () => {
      state.batches = state.batches.filter((entry) => entry.id !== batch.id);
      state.inspection = null;
      $("sampleGrid").innerHTML = "";
      renderBatchList();
      setStatus("已移除文件组");
    });
    header.append(text, remove);
    const files = document.createElement("div");
    files.className = "batch-files";
    batch.files.forEach((file, fileIndex) => {
      const row = document.createElement("div");
      row.className = "batch-file";
      const name = document.createElement("span");
      name.textContent = file.name;
      const deleteFile = document.createElement("button");
      deleteFile.className = "file-remove";
      deleteFile.type = "button";
      deleteFile.textContent = "删除";
      deleteFile.title = `删除 ${file.name}`;
      deleteFile.addEventListener("click", () => {
        batch.files.splice(fileIndex, 1);
        state.batches = state.batches.filter((entry) => entry.files.length > 0);
        state.inspection = null;
        $("sampleGrid").innerHTML = "";
        renderBatchList();
        setStatus(`已删除文件：${file.name}`);
      });
      row.append(name, deleteFile);
      files.appendChild(row);
    });
    item.append(header, files);
    list.appendChild(item);
  });
  const totalFiles = flattenBatches().files.length;
  $("fileLabel").textContent = totalFiles ? `待读取：${totalFiles} 个文件` : "拖拽文件到这里，或点击选择";
}

function handlePickedFiles(files) {
  state.pendingFiles = Array.from(files || []);
  if (state.pendingFiles.length) {
    addBatch(state.pendingFiles, defaultBatchName());
  }
}

function queuePendingFilesIfNeeded() {
  if (state.pendingFiles.length) {
    addBatch(state.pendingFiles, defaultBatchName());
  }
}

function uniqueGroups() {
  const groups = new Set();
  document.querySelectorAll(".sample-group-input").forEach((input) => {
    if (input.value.trim()) groups.add(input.value.trim());
  });
  return Array.from(groups).sort((a, b) => {
    if (a.toLowerCase() === "control") return -1;
    if (b.toLowerCase() === "control") return 1;
    return a.localeCompare(b);
  });
}

function renderGroupControls() {
  const groups = uniqueGroups();
  const controlSelect = $("controlGroup");
  const detectedControl = groups.find((g) => g.toLowerCase() === "control") || groups[0] || "";
  const previousControl = groups.includes(controlSelect.value) ? controlSelect.value : detectedControl;
  controlSelect.innerHTML = "";
  groups.forEach((group) => {
    const option = document.createElement("option");
    option.value = group;
    option.textContent = group;
    option.selected = group === previousControl;
    controlSelect.appendChild(option);
  });
  controlSelect.value = previousControl;

  const holder = $("contrastChecks");
  holder.innerHTML = "";
  groups
    .filter((group) => group !== controlSelect.value)
    .forEach((group) => {
      const label = document.createElement("label");
      label.className = "sample-chip";
      const input = document.createElement("input");
      input.type = "checkbox";
      input.name = "koGroup";
      input.value = group;
      input.checked = true;
      const span = document.createElement("span");
      span.textContent = group;
      label.append(input, span);
      holder.appendChild(label);
    });
}

function renderGeneTypeControls(types = state.geneTypes, selected = state.recipe.deg?.gene_types || DEFAULT_RECIPE.deg.gene_types) {
  const holder = $("geneTypeChecks");
  if (!holder) return;
  const cleanTypes = Array.from(new Set((types || []).map((type) => String(type || "").trim()).filter(Boolean))).sort();
  const list = cleanTypes.length ? cleanTypes : ["protein_coding"];
  state.geneTypes = list;
  const selectedSet = new Set(selected && selected.length ? selected : []);
  holder.innerHTML = "";
  list.forEach((type) => {
    const label = document.createElement("label");
    label.className = "gene-type-option";
    const input = document.createElement("input");
    input.type = "checkbox";
    input.name = "geneType";
    input.value = type;
    input.checked = selectedSet.size ? selectedSet.has(type) : type === "protein_coding";
    input.addEventListener("change", () => {
      updateGeneTypeSummary();
      syncInputsToRecipe();
    });
    const span = document.createElement("span");
    span.textContent = type;
    label.append(input, span);
    holder.appendChild(label);
  });
  updateGeneTypeSummary();
}

function updateGeneTypeSummary() {
  const summary = $("geneTypeSummary");
  if (!summary) return;
  const selected = selectedGeneTypes();
  if (!selected.length) {
    summary.textContent = "All gene types";
  } else if (selected.length <= 2) {
    summary.textContent = selected.join(", ");
  } else {
    summary.textContent = `${selected.length} selected`;
  }
}

function renderSampleTable() {
  const grid = $("sampleGrid");
  grid.innerHTML = "";
  const meta = state.inspection?.sample_meta || [];
  const batchChoices = Array.from(new Set(["batch1", "batch2", "batch3", ...meta.map((row) => row.batch)])).filter(Boolean);
  meta.forEach((row) => {
    const name = document.createElement("div");
    name.className = "sample-name";
    const strong = document.createElement("strong");
    strong.textContent = row.sample;
    const small = document.createElement("span");
    small.textContent = `${row.group} · ${row.batch}`;
    name.append(strong, small);

    const batchSelect = document.createElement("select");
    batchSelect.className = "sample-batch-select";
    batchSelect.dataset.sample = row.sample;
    batchSelect.dataset.file = row.file || "";
    batchSelect.dataset.originalSample = row.original_sample || row.sample;
    batchChoices.forEach((batch) => {
      const option = document.createElement("option");
      option.value = batch;
      option.textContent = batch;
      option.selected = batch === row.batch;
      batchSelect.appendChild(option);
    });

    const groupInput = document.createElement("input");
    groupInput.type = "text";
    groupInput.className = "sample-group-input";
    groupInput.dataset.sample = row.sample;
    groupInput.value = row.group;
    groupInput.addEventListener("input", () => {
      small.textContent = `${groupInput.value} · ${batchSelect.value}`;
      renderGroupControls();
    });
    batchSelect.addEventListener("change", () => {
      small.textContent = `${groupInput.value} · ${batchSelect.value}`;
    });
    grid.append(name, batchSelect, groupInput);
  });
  renderGroupControls();
}

async function inspectFiles() {
  try {
    queuePendingFilesIfNeeded();
    setStatus("读取并合并文件队列中...");
    const payload = await apiPost("/api/inspect", formDataWithBatches());
    state.inspection = payload;
    state.geneTypes = payload.gene_types || state.geneTypes;
    renderGeneTypeControls(state.geneTypes, state.recipe.deg?.gene_types || DEFAULT_RECIPE.deg.gene_types);
    renderSampleTable();
    const fileCount = payload.info.files.length;
    const flags = [
      payload.has_gene_name ? "gene_name" : "无 gene_name",
      payload.has_gene_type ? "gene_type" : "无 gene_type",
      payload.has_annotation ? "GO/KEGG 注释" : "无 GO/KEGG 注释",
    ].join(" · ");
    setStatus(`已读取 ${fileCount} 个文件，common genes ${payload.info.common_gene_count.toLocaleString("zh-CN")}；${flags}`);
  } catch (error) {
    setStatus(error.message, "error");
  }
}

function currentSampleMeta() {
  return Array.from(document.querySelectorAll(".sample-group-input")).map((input) => {
    const batchSelect = document.querySelector(`.sample-batch-select[data-sample="${CSS.escape(input.dataset.sample)}"]`);
    return {
      sample: input.dataset.sample,
      original_sample: batchSelect?.dataset.originalSample || input.dataset.sample,
      batch: batchSelect?.value || "batch1",
      file: batchSelect?.dataset.file || "",
      group: input.value.trim(),
    };
  });
}

function selectedKoGroups() {
  return Array.from(document.querySelectorAll('input[name="koGroup"]:checked')).map((input) => input.value);
}

async function analyzeFiles() {
  try {
    queuePendingFilesIfNeeded();
    syncInputsToRecipe();
    const params = {
      control_group: $("controlGroup").value,
      ko_groups: selectedKoGroups(),
      species: $("speciesSelect").value,
      sample_meta: currentSampleMeta(),
      analysis_recipe: state.recipe,
    };
    setStatus("多 KO 分析中...");
    const payload = await apiPost("/api/analyze", formDataWithBatches({ params }));
    renderResults(payload);
    setStatus(`完成：${payload.summary.contrasts} 个 contrast，${payload.summary.batches} 个 batch`);
  } catch (error) {
    setStatus(error.message, "error");
  }
}

async function loadDemo() {
  try {
    setStatus("加载多 KO 示例...");
    const response = await fetch(apiUrl("/api/demo"));
    const payload = await response.json();
    if (!payload.ok) throw new Error(payload.error || "示例数据加载失败。");
    renderResults(payload);
    setStatus("多 KO 示例数据已完成");
  } catch (error) {
    setStatus(error.message, "error");
  }
}

function renderContrastSelect(contrasts) {
  const select = $("contrastSelect");
  select.innerHTML = "";
  contrasts.forEach((contrast) => {
    const option = document.createElement("option");
    option.value = contrast.label;
    option.textContent = `${contrast.label}  ↑${contrast.summary.up} ↓${contrast.summary.down}`;
    option.selected = contrast.label === state.activeContrast;
    select.appendChild(option);
  });
}

function renderContrastMetrics(contrasts = []) {
  const holder = $("contrastMetrics");
  if (!holder) return;
  holder.innerHTML = "";
  if (!contrasts.length) {
    const empty = document.createElement("div");
    empty.className = "contrast-metric";
    empty.innerHTML = '<div class="contrast-metric-title">Contrasts</div><div class="contrast-metric-counts">-</div>';
    holder.appendChild(empty);
    return;
  }
  contrasts.forEach((contrast) => {
    const card = document.createElement("div");
    card.className = "contrast-metric";
    const title = document.createElement("div");
    title.className = "contrast-metric-title";
    title.textContent = contrast.label;
    const counts = document.createElement("div");
    counts.className = "contrast-metric-counts";
    const up = document.createElement("span");
    up.className = "up";
    up.textContent = `Up ${Number(contrast.summary?.up || 0).toLocaleString("zh-CN")}`;
    const down = document.createElement("span");
    down.className = "down";
    down.textContent = `Down ${Number(contrast.summary?.down || 0).toLocaleString("zh-CN")}`;
    counts.append(up, down);
    card.append(title, counts);
    holder.appendChild(card);
  });
}

function activeContrast() {
  return (state.results?.contrasts || []).find((contrast) => contrast.label === state.activeContrast) || state.results?.contrasts?.[0];
}

function renderResults(payload) {
  state.results = payload;
  state.recipe = payload.recipe || state.recipe;
  state.geneTypes = payload.gene_types || state.geneTypes;
  renderGeneTypeControls(state.geneTypes, state.recipe.deg?.gene_types || DEFAULT_RECIPE.deg.gene_types);
  syncRecipeInputs();
  if (payload.params?.species && $("speciesSelect")) $("speciesSelect").value = payload.params.species;
  $("analysisCode").value = analysisCodeText();
  state.activeContrast = payload.active_contrast || payload.contrasts?.[0]?.label || null;
  $("genesMetric").textContent = payload.summary.genes.toLocaleString("zh-CN");
  $("methodMetric").textContent = payload.summary.method;
  $("downloadBtn").disabled = false;
  renderContrastSelect(payload.contrasts || []);
  renderContrastMetrics(payload.contrasts || []);
  renderVolcano();
  renderHeatmap();
  renderEnrichment("go");
  renderEnrichment("kegg");
  renderResultTable();
}

function clearCanvas(canvas, title) {
  const ctx = canvas.getContext("2d");
  ctx.clearRect(0, 0, canvas.width, canvas.height);
  ctx.fillStyle = "#f7f9fc";
  ctx.fillRect(0, 0, canvas.width, canvas.height);
  ctx.fillStyle = "#637083";
  ctx.font = "16px system-ui, sans-serif";
  if (title) ctx.fillText(title, 32, 42);
  return ctx;
}

function drawAxes(ctx, plot, xLabel, yLabel) {
  ctx.strokeStyle = "#8793a5";
  ctx.lineWidth = 1;
  ctx.beginPath();
  ctx.moveTo(plot.x, plot.y + plot.h);
  ctx.lineTo(plot.x + plot.w, plot.y + plot.h);
  ctx.moveTo(plot.x, plot.y);
  ctx.lineTo(plot.x, plot.y + plot.h);
  ctx.stroke();
  ctx.fillStyle = "#4d5a6c";
  ctx.font = "13px system-ui, sans-serif";
  ctx.fillText(xLabel, plot.x + plot.w / 2 - 50, plot.y + plot.h + 42);
  ctx.save();
  ctx.translate(plot.x - 50, plot.y + plot.h / 2 + 50);
  ctx.rotate(-Math.PI / 2);
  ctx.fillText(yLabel, 0, 0);
  ctx.restore();
}

function renderPlotGallery(kind, emptyMessage) {
  const gallery = $(`${kind}Gallery`);
  const canvas = $(`${kind}Canvas`);
  if (!gallery) return false;
  const plots = state.results?.plots?.[kind] || [];
  gallery.innerHTML = "";
  if (canvas) canvas.style.display = "none";
  if (!plots.length) {
    const empty = document.createElement("div");
    empty.className = "empty-plot";
    empty.textContent = emptyMessage;
    gallery.appendChild(empty);
    return false;
  }
  plots.forEach((plot, index) => {
    const card = document.createElement("article");
    card.className = "plot-card";
    const header = document.createElement("div");
    header.className = "plot-card-header";
    const title = document.createElement("div");
    title.className = "plot-card-title";
    title.textContent = plot.title || `${kind} plot ${index + 1}`;
    const link = document.createElement("a");
    link.href = plot.src;
    link.download = plot.filename || `${kind}_${index + 1}.png`;
    link.textContent = "下载 PNG";
    const img = document.createElement("img");
    img.src = plot.src;
    img.alt = plot.title || `${kind} plot ${index + 1}`;
    header.append(title, link);
    card.append(header, img);
    gallery.appendChild(card);
  });
  return true;
}

function renderVolcano() {
  state.volcanoPoints = [];
  if (renderPlotGallery("volcano", "等待 R 生成 volcano plot")) return;
  const canvas = $("volcanoCanvas");
  if (!canvas) return;
  canvas.style.display = "block";
  const ctx = clearCanvas(canvas, "");
  const contrast = activeContrast();
  const rows = contrast?.results || [];
  if (!rows.length) {
    ctx.fillText("等待分析结果", 32, 42);
    return;
  }
  const recipe = state.results.recipe || state.recipe;
  const fcCut = Number(recipe.deg.log2fc_cutoff || 1);
  const padjCut = Number(recipe.deg.padj_cutoff || 0.05);
  const xlim = recipe.volcano?.xlim || [-7, 7];
  const colors = recipe.volcano?.colors || DEFAULT_RECIPE.volcano.colors;
  const plot = { x: 82, y: 42, w: canvas.width - 128, h: canvas.height - 116 };
  const ys = rows.map((d) => -Math.log10(Math.max(Number(d.padj || d.pvalue || 1), 1e-300))).filter(Number.isFinite);
  const xMin = Number(xlim[0] ?? -7);
  const xMax = Number(xlim[1] ?? 7);
  const yMax = Math.max(2, Math.ceil(Math.max(...ys, 2)));
  const sx = (x) => plot.x + ((Math.max(xMin, Math.min(xMax, x)) - xMin) / (xMax - xMin)) * plot.w;
  const sy = (y) => plot.y + plot.h - (y / yMax) * plot.h;

  ctx.strokeStyle = "#d5dce7";
  ctx.lineWidth = 1;
  [-fcCut, fcCut].forEach((value) => {
    ctx.beginPath();
    ctx.moveTo(sx(value), plot.y);
    ctx.lineTo(sx(value), plot.y + plot.h);
    ctx.stroke();
  });
  const padjLine = -Math.log10(padjCut);
  ctx.beginPath();
  ctx.moveTo(plot.x, sy(padjLine));
  ctx.lineTo(plot.x + plot.w, sy(padjLine));
  ctx.stroke();

  state.volcanoPoints = rows.map((row) => {
    const x = Number(row.log2FoldChange);
    const y = -Math.log10(Math.max(Number(row.padj || row.pvalue || 1), 1e-300));
    const regulation = row.regulation || "NS";
    const color = colors[regulation] || colors.NS || "#B3B3B3";
    return { row, x: sx(x), y: sy(y), color, regulation };
  });

  state.volcanoPoints.forEach((point) => {
    ctx.fillStyle = point.color;
    ctx.globalAlpha = point.regulation === "NS" ? 0.5 : 0.9;
    ctx.beginPath();
    ctx.arc(point.x, point.y, point.regulation === "NS" ? 2.1 : 3.2, 0, Math.PI * 2);
    ctx.fill();
  });
  ctx.globalAlpha = 1;

  const maxLabels = Number(recipe.volcano?.max_labels || 20);
  const selectedTypes = recipe.deg?.gene_types || [];
  const labelRows = rows
    .filter((row) => row.regulation !== "NS" && (!selectedTypes.length || selectedTypes.includes(row.gene_type)) && row.display_gene)
    .sort((a, b) => Number(a.padj) - Number(b.padj))
    .slice(0, maxLabels);
  ctx.font = "11px system-ui, sans-serif";
  ctx.fillStyle = "#16202a";
  ctx.strokeStyle = "rgba(255,255,255,0.88)";
  ctx.lineWidth = 3;
  labelRows.forEach((row, index) => {
    const x = sx(Number(row.log2FoldChange));
    const y = sy(-Math.log10(Math.max(Number(row.padj || row.pvalue || 1), 1e-300)));
    const dx = Number(row.log2FoldChange) >= 0 ? 8 : -68;
    const dy = index % 2 === 0 ? -7 : 12;
    const label = String(row.display_gene).slice(0, 14);
    ctx.strokeText(label, x + dx, y + dy);
    ctx.fillText(label, x + dx, y + dy);
  });

  drawAxes(ctx, plot, "log2 Fold Change", "-log10(FDR)");
  ctx.fillStyle = "#16202a";
  ctx.font = "700 17px system-ui, sans-serif";
  ctx.fillText(`${contrast.label} volcano`, plot.x, 25);
  ctx.font = "12px system-ui, sans-serif";
  ctx.fillStyle = "#637083";
  ctx.fillText(`label: ${selectedTypes.length ? selectedTypes.join(" / ") : "all gene types"} gene_name`, plot.x + 210, 25);
}

function renderHeatmap() {
  if (renderPlotGallery("heatmap", "等待 R 生成 heatmap")) return;
  const canvas = $("heatmapCanvas");
  if (!canvas) return;
  canvas.style.display = "block";
  const ctx = clearCanvas(canvas, "");
  const heatmap = state.results?.heatmap;
  if (!heatmap || !heatmap.genes.length) {
    ctx.fillText("等待热图数据", 32, 42);
    return;
  }
  const recipe = state.results.recipe || state.recipe;
  const genes = heatmap.genes;
  const samples = heatmap.samples;
  const groups = heatmap.groups || [];
  const matrix = heatmap.matrix;
  const left = 170;
  const top = 72;
  const cellW = Math.max(42, Math.floor((canvas.width - left - 36) / samples.length));
  const cellH = Math.max(10, Math.min(18, Math.floor((canvas.height - top - 36) / genes.length)));

  ctx.fillStyle = "#16202a";
  ctx.font = "700 17px system-ui, sans-serif";
  ctx.fillText("DEseq heatmap", 32, 28);
  ctx.fillStyle = "#637083";
  ctx.font = "12px system-ui, sans-serif";
  const geneTypeText = (recipe.deg?.gene_types || []).length ? recipe.deg.gene_types.join(" / ") : "all gene types";
  ctx.fillText(`${geneTypeText} DEGs, row-scaled, batch adjusted`, 170, 28);

  samples.forEach((sample, col) => {
    ctx.fillStyle = groupColor(groups[col], recipe);
    ctx.fillRect(left + col * cellW, top - 18, cellW, 8);
    ctx.save();
    ctx.translate(left + col * cellW + cellW / 2, top - 24);
    ctx.rotate(-Math.PI / 5);
    ctx.fillStyle = "#4d5a6c";
    ctx.font = "12px system-ui, sans-serif";
    ctx.fillText(sample, 0, 0);
    ctx.restore();
  });

  genes.forEach((gene, row) => {
    ctx.fillStyle = "#4d5a6c";
    ctx.font = "11px system-ui, sans-serif";
    const label = gene.length > 22 ? `${gene.slice(0, 22)}...` : gene;
    ctx.fillText(label, 20, top + row * cellH + cellH - 3);
    matrix[row].forEach((value, col) => {
      ctx.fillStyle = heatColor(Number(value), recipe);
      ctx.fillRect(left + col * cellW, top + row * cellH, cellW, cellH);
    });
  });
}

function groupColor(group = "", recipe = state.recipe) {
  const palette = recipe.heatmap?.group_colors || DEFAULT_RECIPE.heatmap.group_colors;
  return palette[group] || stringColor(group);
}

function stringColor(text) {
  let hash = 0;
  for (let i = 0; i < text.length; i += 1) hash = (hash * 31 + text.charCodeAt(i)) >>> 0;
  const colors = ["#047c7b", "#b7791f", "#28724f", "#8b4f9f", "#5367a6", "#b44b62"];
  return colors[hash % colors.length];
}

function heatColor(value, recipe = state.recipe) {
  const clip = Number(recipe.heatmap?.clip || 2);
  const colors = recipeColors(recipe);
  const mid = hexToRgb(colors.heatMid) || [248, 250, 252];
  const high = hexToRgb(colors.heatHigh) || [196, 60, 74];
  const low = hexToRgb(colors.heatLow) || [49, 93, 186];
  const clipped = Math.max(-clip, Math.min(clip, value));
  if (clipped >= 0) {
    const t = clipped / clip;
    return mixColor(mid, high, t);
  }
  const t = Math.abs(clipped) / clip;
  return mixColor(mid, low, t);
}

function mixColor(a, b, t) {
  const rgb = a.map((value, index) => Math.round(value + (b[index] - value) * t));
  return `rgb(${rgb[0]}, ${rgb[1]}, ${rgb[2]})`;
}

function hexToRgb(hex) {
  const match = String(hex || "").match(/^#?([a-f\d]{2})([a-f\d]{2})([a-f\d]{2})$/i);
  if (!match) return null;
  return [parseInt(match[1], 16), parseInt(match[2], 16), parseInt(match[3], 16)];
}

function renderEnrichment(kind) {
  renderPlotGallery(kind, `等待 R 生成 ${kind.toUpperCase()} dotplot`);
  const canvas = $(`${kind}Canvas`);
  const ctx = canvas ? clearCanvas(canvas, "") : null;
  if (canvas) canvas.style.display = "none";
  const tbody = $(`${kind}Body`);
  tbody.innerHTML = "";
  const enrichment = state.results?.enrichment;
  const rows = enrichment?.[kind] || [];
  const title = kind === "go" ? "GO enrichment" : "KEGG enrichment";
  if (!rows.length) {
    if (ctx) ctx.fillText(enrichment?.message || `等待 ${title} 结果`, 32, 42);
    return;
  }
  if (state.results?.plots?.[kind]?.length) {
    rows.slice(0, 120).forEach((term) => {
      const tr = document.createElement("tr");
      [term.contrast, term.direction, term.category, term.term_name || term.term_id, `${term.overlap}/${term.term_size}`, fmt(term.padj, 3)].forEach((value) => {
        const td = document.createElement("td");
        td.textContent = value;
        tr.appendChild(td);
      });
      tbody.appendChild(tr);
    });
    return;
  }
  if (canvas) canvas.style.display = "block";
  const terms = rows.slice(0, 18).reverse();
  const plot = { x: 300, y: 46, w: canvas.width - 380, h: canvas.height - 112 };
  const maxRatio = Math.max(...terms.map((d) => Number(d.gene_ratio)), 0.01);
  const maxOverlap = Math.max(...terms.map((d) => Number(d.overlap)), 1);
  const maxScore = Math.max(...terms.map((d) => -Math.log10(Math.max(Number(d.padj || d.pvalue || 1), 1e-300))), 1);
  const step = plot.h / Math.max(terms.length - 1, 1);
  const xScale = (ratio) => plot.x + (Number(ratio) / maxRatio) * plot.w;
  const yScale = (index) => plot.y + index * step;

  ctx.fillStyle = "#16202a";
  ctx.font = "700 17px system-ui, sans-serif";
  ctx.fillText(title, 32, 25);
  ctx.strokeStyle = "#8793a5";
  ctx.lineWidth = 1;
  ctx.beginPath();
  ctx.moveTo(plot.x, plot.y + plot.h + 20);
  ctx.lineTo(plot.x + plot.w, plot.y + plot.h + 20);
  ctx.moveTo(plot.x, plot.y - 12);
  ctx.lineTo(plot.x, plot.y + plot.h + 20);
  ctx.stroke();

  terms.forEach((term, index) => {
    const y = yScale(index);
    const x = xScale(Number(term.gene_ratio));
    const score = -Math.log10(Math.max(Number(term.padj || term.pvalue || 1), 1e-300));
    const radius = 4 + (Number(term.overlap) / maxOverlap) * 10;
    ctx.fillStyle = "#4d5a6c";
    ctx.font = "11px system-ui, sans-serif";
    const label = `${term.contrast} ${term.direction} · ${term.term_name || term.term_id}`;
    ctx.fillText(label.length > 42 ? `${label.slice(0, 42)}...` : label, 18, y + 4);
    ctx.strokeStyle = "#e2e8f0";
    ctx.beginPath();
    ctx.moveTo(plot.x, y);
    ctx.lineTo(plot.x + plot.w, y);
    ctx.stroke();
    ctx.fillStyle = dotColor(score / maxScore, kind);
    ctx.beginPath();
    ctx.arc(x, y, radius, 0, Math.PI * 2);
    ctx.fill();
    ctx.strokeStyle = "rgba(22,32,42,0.25)";
    ctx.stroke();
  });

  ctx.fillStyle = "#4d5a6c";
  ctx.font = "12px system-ui, sans-serif";
  ctx.fillText("GeneRatio", plot.x + plot.w / 2 - 28, canvas.height - 22);
  [0, 0.5, 1].forEach((tick) => {
    const x = plot.x + tick * plot.w;
    ctx.fillText(fmt(tick * maxRatio, 2), x - 10, plot.y + plot.h + 38);
    ctx.beginPath();
    ctx.moveTo(x, plot.y + plot.h + 14);
    ctx.lineTo(x, plot.y + plot.h + 22);
    ctx.stroke();
  });

  ctx.fillStyle = "#637083";
  ctx.fillText("dot size = Count, color = -log10(padj)", plot.x, 25);

  rows.slice(0, 120).forEach((term) => {
    const tr = document.createElement("tr");
    [term.contrast, term.direction, term.category, term.term_name || term.term_id, `${term.overlap}/${term.term_size}`, fmt(term.padj, 3)].forEach((value) => {
      const td = document.createElement("td");
      td.textContent = value;
      tr.appendChild(td);
    });
    tbody.appendChild(tr);
  });
}

function dotColor(t, kind) {
  const low = kind === "kegg" ? [246, 224, 181] : [199, 236, 229];
  const high = kind === "kegg" ? [183, 121, 31] : [4, 124, 123];
  return mixColor(low, high, Math.max(0, Math.min(1, t)));
}

function renderResultTable() {
  const tbody = $("resultBody");
  tbody.innerHTML = "";
  const contrast = activeContrast();
  (contrast?.results || []).slice(0, 700).forEach((row) => {
    const tr = document.createElement("tr");
    [
      row.gene,
      row.display_gene || row.gene_name || "",
      contrast.label,
      row.regulation,
      row.baseMean,
      row.log2FoldChange,
      row.pvalue,
      row.padj,
      row.dispersion,
    ].forEach((value, index) => {
      const td = document.createElement("td");
      td.textContent = index <= 3 ? value : fmt(value, 4);
      tr.appendChild(td);
    });
    tbody.appendChild(tr);
  });
}

async function downloadExcel() {
  if (!state.results?.contrasts?.length) return;
  const response = await fetch(apiUrl("/api/export_deseq"), {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(state.results),
  });
  if (!response.ok) {
    setStatus("Excel 导出失败", "error");
    return;
  }
  const blob = await response.blob();
  const url = URL.createObjectURL(blob);
  const link = document.createElement("a");
  link.href = url;
  link.download = "DESeq2_all_results.xlsx";
  link.click();
  URL.revokeObjectURL(url);
  setStatus("DESeq 结果 Excel 已导出");
}

function switchTab(event) {
  const tab = event.currentTarget.dataset.tab;
  document.querySelectorAll(".tab").forEach((button) => button.classList.toggle("active", button.dataset.tab === tab));
  document.querySelectorAll(".view").forEach((view) => view.classList.toggle("active", view.id === tab));
}

function handleVolcanoMove(event) {
  const tip = $("volcanoTip");
  if (!state.volcanoPoints.length) return;
  const canvas = $("volcanoCanvas");
  const rect = canvas.getBoundingClientRect();
  const scaleX = canvas.width / rect.width;
  const scaleY = canvas.height / rect.height;
  const x = (event.clientX - rect.left) * scaleX;
  const y = (event.clientY - rect.top) * scaleY;
  let nearest = null;
  let dist = Infinity;
  state.volcanoPoints.forEach((point) => {
    const d = Math.hypot(point.x - x, point.y - y);
    if (d < dist) {
      dist = d;
      nearest = point;
    }
  });
  if (!nearest || dist > 10) {
    tip.style.display = "none";
    return;
  }
  tip.innerHTML = `<strong>${nearest.row.display_gene || nearest.row.gene}</strong><br>${nearest.row.gene}<br>log2FC ${fmt(nearest.row.log2FoldChange, 4)}<br>padj ${fmt(nearest.row.padj, 4)}`;
  tip.style.left = `${event.offsetX + 18}px`;
  tip.style.top = `${event.offsetY + 18}px`;
  tip.style.display = "block";
}

let uploadEventsBound = false;

function bindUploadEvents() {
  if (uploadEventsBound) return;
  const dropZone = $("dropZone");
  const fileInput = $("fileInput");
  const addBatchBtn = $("addBatchBtn");
  const clearFilesBtn = $("clearFilesBtn");
  if (!dropZone || !fileInput || !addBatchBtn || !clearFilesBtn) return;
  uploadEventsBound = true;

  dropZone.addEventListener("click", () => fileInput.click());
  const handleDroppedFiles = (event) => {
    event.preventDefault();
    event.stopPropagation();
    dropZone.classList.remove("drag-over");
    const files = event.dataTransfer?.files || [];
    if (files.length) addBatch(files, defaultBatchName());
  };
  dropZone.addEventListener("dragover", (event) => {
    event.preventDefault();
    dropZone.classList.add("drag-over");
  });
  dropZone.addEventListener("dragleave", () => dropZone.classList.remove("drag-over"));
  dropZone.addEventListener("drop", handleDroppedFiles);
  document.addEventListener("dragover", (event) => {
    event.preventDefault();
    if (Array.from(event.dataTransfer?.types || []).includes("Files")) dropZone.classList.add("drag-over");
  });
  document.addEventListener("drop", handleDroppedFiles);
  fileInput.addEventListener("change", (event) => handlePickedFiles(event.target.files));
  addBatchBtn.addEventListener("click", () => addBatch(state.pendingFiles, defaultBatchName()));
  clearFilesBtn.addEventListener("click", () => {
    state.batches = [];
    state.pendingFiles = [];
    state.inspection = null;
    fileInput.value = "";
    $("sampleGrid").innerHTML = "";
    renderBatchList();
    setStatus("已清空文件队列");
  });
}

function wireEvents() {
  $("inspectBtn").addEventListener("click", inspectFiles);
  $("analyzeBtn").addEventListener("click", analyzeFiles);
  $("demoBtn").addEventListener("click", loadDemo);
  $("downloadBtn").addEventListener("click", downloadExcel);
  $("contrastSelect").addEventListener("change", (event) => {
    state.activeContrast = event.target.value;
    renderVolcano();
    renderResultTable();
  });
  $("controlGroup").addEventListener("change", renderGroupControls);
  $("fcCutoff").addEventListener("change", syncInputsToRecipe);
  $("padjCutoff").addEventListener("change", syncInputsToRecipe);
  $("speciesSelect").addEventListener("change", syncInputsToRecipe);
  $("geneTypeAllBtn").addEventListener("click", () => {
    document.querySelectorAll('input[name="geneType"]').forEach((input) => {
      input.checked = true;
    });
    updateGeneTypeSummary();
    syncInputsToRecipe();
  });
  $("geneTypeDefaultBtn").addEventListener("click", () => {
    document.querySelectorAll('input[name="geneType"]').forEach((input) => {
      input.checked = input.value === "protein_coding";
    });
    updateGeneTypeSummary();
    syncInputsToRecipe();
  });
  $("codeAllBtn").addEventListener("click", () => {
    document.querySelectorAll('input[name="codeModule"]').forEach((input) => {
      input.checked = true;
    });
    syncInputsToRecipe();
  });
  $("codeNoneBtn").addEventListener("click", () => {
    document.querySelectorAll('input[name="codeModule"]').forEach((input) => {
      input.checked = false;
    });
    syncInputsToRecipe();
  });
  [
    "colorUp",
    "colorDown",
    "colorNs",
    "colorHeatHigh",
    "colorHeatMid",
    "colorHeatLow",
    "colorGroup1",
    "colorGroup2",
    "colorGroup3",
  ].forEach((id) => $(id).addEventListener("input", syncInputsToRecipe));
  document.querySelectorAll('input[name="codeModule"]').forEach((input) => {
    input.addEventListener("change", syncInputsToRecipe);
  });
  $("applyCodeBtn").addEventListener("click", () => {
    syncInputsToRecipe();
    setStatus("R 代码已按当前参数刷新");
  });
  $("resetCodeBtn").addEventListener("click", resetRecipe);
  $("volcanoCanvas").addEventListener("mousemove", handleVolcanoMove);
  $("volcanoCanvas").addEventListener("mouseleave", () => ($("volcanoTip").style.display = "none"));
  document.querySelectorAll(".tab").forEach((button) => button.addEventListener("click", switchTab));
}

async function checkBackendHealth() {
  try {
    const response = await fetch(apiUrl("/api/health"), { method: "GET" });
    const payload = await response.json();
    if (!payload.ok) throw new Error("health check failed");
  } catch (error) {
    const message = API_BASE
      ? `后端 ${API_BASE} 暂时无法连接，读取队列会失败。`
      : "后端 API 暂时无法连接。如果这是公网部署，请确认不是只部署了静态网页，而是部署了包含 Python/R 的 Docker Web Service。";
    setStatus(message, "error");
  }
}

bindUploadEvents();
try {
  renderGeneTypeControls();
  syncRecipeInputs();
  $("analysisCode").value = analysisCodeText();
  wireEvents();
  checkBackendHealth();
} catch (error) {
  console.error(error);
  setStatus(`页面部分控件初始化失败，但文件上传仍可用：${error.message}`, "error");
}
