suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(limma)
  library(dplyr)
  library(jsonlite)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
  library(enrichplot)
  library(base64enc)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript deseq_pipeline.R counts.csv metadata.csv params.json output.json")
}

counts_path <- args[[1]]
metadata_path <- args[[2]]
params_path <- args[[3]]
output_path <- args[[4]]

params <- fromJSON(params_path, simplifyVector = FALSE)
recipe <- params$analysis_recipe
control_group <- params$control_group
ko_groups <- unlist(params$ko_groups)
species <- params$species
if (is.null(species) || !species %in% c("human", "mouse")) species <- "human"
org_db <- if (species == "mouse") org.Mm.eg.db else org.Hs.eg.db
kegg_organism <- if (species == "mouse") "mmu" else "hsa"
species_label <- if (species == "mouse") "Mouse" else "Human"

`%||%` <- function(x, y) if (is.null(x)) y else x

padj_cutoff <- as.numeric(recipe$deg$padj_cutoff %||% 0.05)
log2fc_cutoff <- as.numeric(recipe$deg$log2fc_cutoff %||% 1)
top_genes_n <- as.integer(recipe$heatmap$top_genes %||% 50)
heat_clip <- as.numeric(recipe$heatmap$clip %||% 2)
max_terms <- as.integer(recipe$enrichment$max_terms %||% 80)
selected_gene_types <- unlist(recipe$deg$gene_types %||% list("protein_coding"))
selected_gene_types <- selected_gene_types[!is.na(selected_gene_types) & selected_gene_types != ""]
volcano_colors <- recipe$volcano$colors %||% list(UP = "#D62728", DOWN = "#1F77B4", NS = "#B3B3B3")
heatmap_colors_recipe <- recipe$heatmap$colors %||% list(low = "#1F77B4", mid = "#F8FAFC", high = "#D62728")
group_palette_recipe <- unlist(recipe$heatmap$group_palette %||% list("#4D4D4D", "#D62728", "#1F77B4", "#047C7B", "#B7791F", "#28724F", "#8B4F9F"))

filter_gene_type <- function(df) {
  if (!length(selected_gene_types) || !"gene_type" %in% colnames(df)) return(df)
  df[df$gene_type %in% selected_gene_types, , drop = FALSE]
}

gene_type_label <- if (length(selected_gene_types)) paste(selected_gene_types, collapse = " / ") else "all gene types"

counts_df <- read.csv(counts_path, check.names = FALSE, stringsAsFactors = FALSE)
metadata <- read.csv(metadata_path, check.names = FALSE, stringsAsFactors = FALSE)

if (!"gene" %in% colnames(counts_df)) stop("counts table must include gene column")
if (!all(c("sample", "batch", "group") %in% colnames(metadata))) {
  stop("metadata must include sample, batch, group columns")
}

sample_cols <- metadata$sample
missing_samples <- setdiff(sample_cols, colnames(counts_df))
if (length(missing_samples) > 0) {
  stop(paste("counts table missing samples:", paste(missing_samples, collapse = ", ")))
}

count_mat <- as.matrix(counts_df[, sample_cols, drop = FALSE])
mode(count_mat) <- "numeric"
count_mat[is.na(count_mat)] <- 0
count_mat <- round(count_mat)
rownames(count_mat) <- counts_df$gene

gene_anno_cols <- intersect(c("gene", "gene_name", "gene_type", "entrezgene_id"), colnames(counts_df))
gene_anno <- unique(counts_df[, gene_anno_cols, drop = FALSE])
if (!"gene_name" %in% colnames(gene_anno)) gene_anno$gene_name <- gene_anno$gene
if (!"gene_type" %in% colnames(gene_anno)) gene_anno$gene_type <- ""

metadata$group <- factor(metadata$group, levels = c(control_group, ko_groups))
metadata$batch <- factor(metadata$batch)
rownames(metadata) <- metadata$sample
metadata <- metadata[sample_cols, , drop = FALSE]

dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData = metadata,
  design = ~ batch + group
)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)
vsd <- tryCatch(
  vst(dds, blind = FALSE),
  error = function(e) varianceStabilizingTransformation(dds, blind = FALSE)
)

clean_df <- function(df) {
  df <- as.data.frame(df)
  df[] <- lapply(df, function(col) {
    if (is.factor(col)) as.character(col) else col
  })
  df
}

records <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(list())
  df <- clean_df(df)
  rows <- lapply(seq_len(nrow(df)), function(i) {
    row <- as.list(df[i, , drop = FALSE])
    names(row) <- colnames(df)
    row
  })
  unname(rows)
}

plot_uri <- function(filename, draw, width = 1344, height = 960, res = 144) {
  path <- file.path(dirname(output_path), filename)
  png(path, width = width, height = height, res = res)
  ok <- tryCatch({
    draw()
    TRUE
  }, error = function(e) {
    message("Plot failed: ", conditionMessage(e))
    FALSE
  })
  dev.off()
  if (!ok || !file.exists(path)) return(NULL)
  paste0("data:image/png;base64,", base64enc::base64encode(path))
}

safe_name <- function(x) {
  gsub("[^A-Za-z0-9_-]+", "_", x)
}

plot_record <- function(kind, title, filename, draw, width = 1344, height = 960, res = 144) {
  src <- plot_uri(filename, draw, width = width, height = height, res = res)
  if (is.null(src)) return(NULL)
  list(kind = kind, title = title, filename = filename, src = src)
}

find_coef <- function(group, control) {
  exact <- paste0("group_", group, "_vs_", control)
  rn <- resultsNames(dds)
  if (exact %in% rn) return(exact)
  group_safe <- make.names(group)
  control_safe <- make.names(control)
  candidates <- rn[grepl("^group_.*_vs_.*$", rn)]
  hit <- candidates[grepl(group_safe, candidates, fixed = TRUE) & grepl(control_safe, candidates, fixed = TRUE)]
  if (length(hit) > 0) hit[[1]] else NA_character_
}

run_deseq_contrast <- function(group) {
  coef_name <- find_coef(group, control_group)
  if (!is.na(coef_name)) {
    res <- lfcShrink(dds, coef = coef_name, type = "apeglm")
  } else {
    res <- results(dds, contrast = c("group", group, control_group))
  }
  df <- as.data.frame(res)
  df$gene <- rownames(df)
  df <- df[!is.na(df$padj), , drop = FALSE]
  df <- left_join(df, gene_anno, by = "gene")
  df$gene_name[is.na(df$gene_name)] <- ""
  df$gene_type[is.na(df$gene_type)] <- ""
  df$display_gene <- ifelse(df$gene_name != "", df$gene_name, df$gene)
  df$regulation <- "NS"
  df$regulation[df$padj < padj_cutoff & df$log2FoldChange > log2fc_cutoff] <- "UP"
  df$regulation[df$padj < padj_cutoff & df$log2FoldChange < -log2fc_cutoff] <- "DOWN"
  df
}

contrast_results <- list()
flat_results <- list()
total_up <- 0
total_down <- 0

for (group in ko_groups) {
  df <- run_deseq_contrast(group)
  label <- paste(group, "vs", control_group)
  up <- sum(df$regulation == "UP", na.rm = TRUE)
  down <- sum(df$regulation == "DOWN", na.rm = TRUE)
  total_up <- total_up + up
  total_down <- total_down + down
  flat <- df
  flat$contrast <- label
  flat_results[[label]] <- flat
  contrast_results[[length(contrast_results) + 1]] <- list(
    label = label,
    group = group,
    control = control_group,
    samples = list(
      control = metadata$sample[metadata$group == control_group],
      ko = metadata$sample[metadata$group == group]
    ),
    summary = list(up = up, down = down, genes = nrow(df)),
    results = records(df)
  )
}

all_results <- bind_rows(flat_results)

volcano_plots <- list()
for (contrast in contrast_results) {
  df <- flat_results[[contrast$label]]
  plot_title <- paste0(contrast$group, " vs ", control_group, " (DESeq2)")
  volcano <- ggplot(df, aes(log2FoldChange, -log10(padj))) +
    geom_point(color = volcano_colors$NS %||% "#B3B3B3", size = 1) +
    geom_point(data = subset(df, regulation == "UP"), color = volcano_colors$UP %||% "#D62728", size = 1.5) +
    geom_point(data = subset(df, regulation == "DOWN"), color = volcano_colors$DOWN %||% "#1F77B4", size = 1.5) +
    geom_text_repel(
      data = filter_gene_type(subset(df, regulation != "NS" & gene_name != "")),
      aes(label = gene_name),
      size = 3,
      max.overlaps = 20
    ) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    coord_cartesian(xlim = c(-7, 7)) +
    scale_x_continuous(expand = expansion(mult = 0.05)) +
    theme_classic() +
    labs(title = plot_title, x = "log2 Fold Change", y = "-log10(FDR)")
  rec <- plot_record(
    kind = "volcano",
    title = plot_title,
    filename = paste0("volcano_", safe_name(contrast$group), ".png"),
    draw = function() print(volcano)
  )
  if (!is.null(rec)) volcano_plots[[length(volcano_plots) + 1]] <- rec
}

heatmap_payload <- list(genes = list(), gene_ids = list(), samples = sample_cols, groups = as.character(metadata$group), batches = as.character(metadata$batch), matrix = list())
deg_genes <- all_results %>%
  filter(regulation != "NS") %>%
  filter_gene_type() %>%
  pull(gene) %>%
  unique()

if (length(deg_genes) > 0) {
  mat <- assay(vsd)[deg_genes, , drop = FALSE]
  mat_corrected <- limma::removeBatchEffect(mat, batch = metadata$batch)
  gene_name_map <- gene_anno$gene_name[match(rownames(mat_corrected), gene_anno$gene)]
  gene_name_map[is.na(gene_name_map) | gene_name_map == ""] <- rownames(mat_corrected)[is.na(gene_name_map) | gene_name_map == ""]
  rownames(mat_corrected) <- make.unique(gene_name_map)
  mat_scaled <- t(scale(t(mat_corrected)))
  mat_scaled[is.na(mat_scaled)] <- 0
  mat_scaled[mat_scaled > heat_clip] <- heat_clip
  mat_scaled[mat_scaled < -heat_clip] <- -heat_clip
  sample_order <- order(metadata$group)
  mat_scaled <- mat_scaled[, sample_order, drop = FALSE]
  ordered_meta <- metadata[sample_order, , drop = FALSE]
  heatmap_payload <- list(
    genes = rownames(mat_scaled),
    gene_ids = deg_genes,
    samples = colnames(mat_scaled),
    groups = as.character(ordered_meta$group),
    batches = as.character(ordered_meta$batch),
    matrix = unname(lapply(seq_len(nrow(mat_scaled)), function(i) unname(as.numeric(round(mat_scaled[i, ], 4)))))
  )
}

heat_colors <- colorRampPalette(c(
  heatmap_colors_recipe$low %||% "#1F77B4",
  heatmap_colors_recipe$mid %||% "#F8FAFC",
  heatmap_colors_recipe$high %||% "#D62728"
))(100)
active_levels <- levels(droplevels(metadata$group))
palette_values <- rep(c(group_palette_recipe, "#047C7B", "#B7791F", "#28724F", "#8B4F9F"), length.out = length(active_levels))
all_group_colors <- setNames(palette_values, active_levels)
heatmap_plots <- list()

make_heatmap_plot <- function(gene_ids, selected_samples, title, cluster_cols = FALSE, ann_colors = all_group_colors, fontsize_row = 6) {
  gene_ids <- intersect(gene_ids, rownames(vsd))
  if (length(gene_ids) == 0) return(NULL)
  mat <- assay(vsd)[gene_ids, , drop = FALSE]
  mat <- limma::removeBatchEffect(mat, batch = metadata$batch)
  mat <- mat[, selected_samples, drop = FALSE]
  gene_name_map <- gene_anno$gene_name[match(rownames(mat), gene_anno$gene)]
  gene_name_map[is.na(gene_name_map) | gene_name_map == ""] <- rownames(mat)[is.na(gene_name_map) | gene_name_map == ""]
  rownames(mat) <- make.unique(gene_name_map)
  ann_col <- data.frame(Group = metadata[colnames(mat), "group"])
  rownames(ann_col) <- colnames(mat)
  ann_col$Group <- droplevels(factor(as.character(ann_col$Group), levels = levels(metadata$group)))
  pheatmap::pheatmap(
    mat,
    scale = "row",
    color = heat_colors,
    cluster_cols = cluster_cols,
    cluster_rows = TRUE,
    annotation_col = ann_col,
    annotation_colors = list(Group = ann_colors),
    show_rownames = TRUE,
    fontsize_row = fontsize_row,
    border_color = NA,
    main = title
  )
}

if (length(deg_genes) > 0) {
  ordered_samples <- rownames(metadata)[order(metadata$group)]
  all_heatmap_title <- paste(paste(ko_groups, collapse = " & "), "vs", control_group, paste0("(", gene_type_label, " DEGs)"))
  rec <- plot_record(
    kind = "heatmap",
    title = all_heatmap_title,
    filename = "heatmap_all_contrasts_ordered.png",
    draw = function() make_heatmap_plot(deg_genes, ordered_samples, all_heatmap_title, cluster_cols = FALSE, fontsize_row = 6)
  )
  if (!is.null(rec)) heatmap_plots[[length(heatmap_plots) + 1]] <- rec
}

map_entrez <- function(df) {
  if ("entrezgene_id" %in% colnames(df) && any(!is.na(df$entrezgene_id) & df$entrezgene_id != "")) {
    df$entrezgene_id <- as.character(df$entrezgene_id)
  } else {
    mapped <- tryCatch(
      AnnotationDbi::mapIds(
        org_db,
        keys = df$gene_name,
        column = "ENTREZID",
        keytype = "SYMBOL",
        multiVals = "first"
      ),
      error = function(e) character()
    )
    if (length(mapped) == 0) return(df[0, , drop = FALSE])
    df$entrezgene_id <- unname(as.character(mapped[df$gene_name]))
  }
  df %>%
    filter(!is.na(entrezgene_id), entrezgene_id != "") %>%
    distinct(entrezgene_id, .keep_all = TRUE)
}

mapped_results <- lapply(flat_results, map_entrez)

format_enrich <- function(enrich_obj, contrast, direction, category, significant_n) {
  if (is.null(enrich_obj)) return(data.frame())
  df <- as.data.frame(enrich_obj)
  if (nrow(df) == 0) return(data.frame())
  ratio <- vapply(strsplit(df$GeneRatio, "/", fixed = TRUE), function(x) as.numeric(x[[1]]) / as.numeric(x[[2]]), numeric(1))
  data.frame(
    contrast = contrast,
    direction = direction,
    category = category,
    term_id = df$ID,
    term_name = df$Description,
    overlap = df$Count,
    term_size = NA_integer_,
    significant_genes = significant_n,
    pvalue = df$pvalue,
    padj = df$p.adjust,
    gene_ratio = ratio,
    genes = df$geneID,
    stringsAsFactors = FALSE
  )
}

enrichment_padj_cutoff <- function(group, direction) {
  if (group == "KO2" && direction == "DOWN") return(0.1)
  padj_cutoff
}

direction_word <- function(direction) {
  if (direction == "UP") "upregulated" else "downregulated"
}

empty_enrich_plot <- function(title) {
  ggplot() +
    theme_void() +
    annotate("text", x = 0, y = 0, label = "No enriched terms", size = 5, color = "grey40") +
    ggtitle(title)
}

run_enrichments <- function(label, group, df, direction) {
  sig_padj <- enrichment_padj_cutoff(group, direction)
  go_title <- paste(group, direction_word(direction), "GO (BP)")
  kegg_title <- paste0(group, " ", direction_word(direction), " KEGG (", species_label, ")")
  if (direction == "UP") {
    sig <- df %>% filter(log2FoldChange > 0, padj < sig_padj) %>% filter_gene_type()
  } else {
    sig <- df %>% filter(log2FoldChange < 0, padj < sig_padj) %>% filter_gene_type()
  }
  genes <- as.character(sig$entrezgene_id)
  if (length(genes) == 0) {
    go_plot <- plot_record(
      kind = "go",
      title = go_title,
      filename = paste0("go_", safe_name(group), "_", tolower(direction), ".png"),
      draw = function() print(empty_enrich_plot(go_title))
    )
    kegg_plot <- plot_record(
      kind = "kegg",
      title = kegg_title,
      filename = paste0("kegg_", safe_name(group), "_", tolower(direction), ".png"),
      draw = function() print(empty_enrich_plot(kegg_title))
    )
    return(list(go = data.frame(), kegg = data.frame(), go_plot = go_plot, kegg_plot = kegg_plot))
  }
  go_obj <- tryCatch(
    enrichGO(
      gene = genes,
      OrgDb = org_db,
      keyType = "ENTREZID",
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05,
      readable = TRUE
    ),
    error = function(e) NULL
  )
  kegg_obj <- tryCatch(
    enrichKEGG(
      gene = genes,
      organism = kegg_organism,
      keyType = "ncbi-geneid",
      pvalueCutoff = 0.05
    ),
    error = function(e) NULL
  )
  go_plot <- NULL
  if (!is.null(go_obj) && nrow(as.data.frame(go_obj)) > 0) {
    go_plot <- plot_record(
      kind = "go",
      title = go_title,
      filename = paste0("go_", safe_name(group), "_", tolower(direction), ".png"),
      draw = function() print(dotplot(go_obj, showCategory = 10) + ggtitle(go_title))
    )
  } else {
    go_plot <- plot_record(
      kind = "go",
      title = go_title,
      filename = paste0("go_", safe_name(group), "_", tolower(direction), ".png"),
      draw = function() print(empty_enrich_plot(go_title))
    )
  }
  kegg_plot <- NULL
  if (!is.null(kegg_obj) && nrow(as.data.frame(kegg_obj)) > 0) {
    kegg_plot <- plot_record(
      kind = "kegg",
      title = kegg_title,
      filename = paste0("kegg_", safe_name(group), "_", tolower(direction), ".png"),
      draw = function() print(dotplot(kegg_obj, showCategory = 10) + ggtitle(kegg_title))
    )
  } else {
    kegg_plot <- plot_record(
      kind = "kegg",
      title = kegg_title,
      filename = paste0("kegg_", safe_name(group), "_", tolower(direction), ".png"),
      draw = function() print(empty_enrich_plot(kegg_title))
    )
  }
  list(
    go = format_enrich(go_obj, label, direction, "GO", length(genes)),
    kegg = format_enrich(kegg_obj, label, direction, "KEGG", length(genes)),
    go_plot = go_plot,
    kegg_plot = kegg_plot
  )
}

go_frames <- list()
kegg_frames <- list()
go_plots <- list()
kegg_plots <- list()
contrast_groups <- setNames(vapply(contrast_results, function(x) x$group, character(1)), vapply(contrast_results, function(x) x$label, character(1)))
for (label in names(flat_results)) {
  group <- contrast_groups[[label]]
  for (direction in c("UP", "DOWN")) {
    enr <- run_enrichments(label, group, mapped_results[[label]], direction)
    if (nrow(enr$go) > 0) go_frames[[length(go_frames) + 1]] <- enr$go
    if (nrow(enr$kegg) > 0) kegg_frames[[length(kegg_frames) + 1]] <- enr$kegg
    if (!is.null(enr$go_plot)) go_plots[[length(go_plots) + 1]] <- enr$go_plot
    if (!is.null(enr$kegg_plot)) kegg_plots[[length(kegg_plots) + 1]] <- enr$kegg_plot
  }
}

go_df <- if (length(go_frames)) bind_rows(go_frames) %>% arrange(padj, pvalue) %>% slice_head(n = max_terms) else data.frame()
kegg_df <- if (length(kegg_frames)) bind_rows(kegg_frames) %>% arrange(padj, pvalue) %>% slice_head(n = max_terms) else data.frame()
terms_df <- bind_rows(go_df, kegg_df)

payload <- list(
  ok = TRUE,
  params = params,
  recipe = recipe,
  gene_types = sort(unique(as.character(gene_anno$gene_type[!is.na(gene_anno$gene_type) & gene_anno$gene_type != ""]))),
  sample_meta = records(metadata),
  summary = list(
    genes = nrow(dds),
    samples = length(sample_cols),
    batches = length(unique(metadata$batch)),
    contrasts = length(contrast_results),
    up = total_up,
    down = total_down,
        method = paste0("R/Bioconductor DESeq2 + apeglm (", species_label, " enrichment)")
  ),
  contrasts = contrast_results,
  active_contrast = if (length(contrast_results)) contrast_results[[1]]$label else "",
  de = records(all_results),
  heatmap = heatmap_payload,
  plots = list(
    volcano = volcano_plots,
    heatmap = heatmap_plots,
    go = go_plots,
    kegg = kegg_plots
  ),
  enrichment = list(
    available = TRUE,
    message = "",
    terms = records(terms_df),
    go = records(go_df),
    kegg = records(kegg_df)
  )
)

write_json(payload, output_path, auto_unbox = TRUE, dataframe = "rows", null = "null", na = "null", digits = NA)
