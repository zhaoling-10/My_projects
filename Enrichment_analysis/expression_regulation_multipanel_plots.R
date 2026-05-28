#!/usr/bin/env Rscript
# =============================================================================
# Multi-panel expression plots: LT vs LE with KOBAS enrichment context
#
# Data (same 8 RNA-seq samples in every count file: LE_1–4 + LT_1–4):
#   LE_counts_*  → mapped to Lepus europaeus reference (NC gene IDs)
#   LT_counts_*  → mapped to Lepus timidus reference (CM gene IDs)
#   Filenames refer to REFERENCE GENOME, not “LE-only” vs “LT-only” samples.
#   DESeq2: LT_vs_LE_adj.csv (europaeus ref) + LE_vs_LT_adj.csv (timidus ref; LFC flipped)
#   KOBAS:  kobas_all_pathway_genes_with_direction.csv (regulation already harmonized)
#
# Outputs: mmu_analysis/figures/regulation_multipanel/
#
# Note on log10(padj): see README section at end of script / console on run.
# Volcano Y-axis uses -log10(padj) (standard); values increase upward for smaller padj.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(readr)
  library(DESeq2)
})

has_pheatmap <- requireNamespace("pheatmap", quietly = TRUE)
has_ggraph   <- requireNamespace("ggraph", quietly = TRUE) &&
  requireNamespace("igraph", quietly = TRUE)

# ---- Paths (symmetric: one folder per reference genome) ----
base_dir <- path.expand("~/Documents/Enrichment_analysis")
deseq_base <- path.expand("~/Documents/DESeq2/DESEq2_analysis")

ref_eur <- list(
  label = "Lepus_europaeus (NC IDs)",
  dir = file.path(deseq_base, "1.Lepus_europaeus_as_reference"),
  counts_adj = "in/LE_counts_adj.txt",
  counts_raw = "in/LE_counts_raw.txt",
  metadata = "in/LE_metadata.txt",
  deseq = "LT_vs_LE_adj.csv",
  flip_lfc = FALSE
)
ref_tim <- list(
  label = "Lepus_timidus (CM IDs)",
  dir = file.path(deseq_base, "2.Lepus timidus_as_reference"),
  counts_adj = "in/LT_counts_adj.txt",
  counts_raw = "in/LT_counts_raw.txt",
  metadata = "in/LT_metadata.txt",
  deseq = "LE_vs_LT_adj.csv",
  flip_lfc = TRUE
)

out_fig  <- file.path(base_dir, "mmu_analysis", "figures", "regulation_multipanel")
out_tab  <- file.path(base_dir, "mmu_analysis", "tables", "regulation_multipanel")
pathway_file <- file.path(base_dir, "mmu_analysis/tables/kobas_all_pathway_genes_with_direction.csv")
dir.create(out_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(out_tab, recursive = TRUE, showWarnings = FALSE)

resolve_first_existing <- function(paths) {
  for (p in paths) {
    p <- path.expand(p)
    if (file.exists(p)) return(p)
  }
  NULL
}

padj_thresh  <- 0.05
logfc_thresh <- 1.0
top_genes_heatmap <- 50
top_path_network  <- 8
top_path_mosaic   <- 10

sample_levels <- c("LE_1", "LE_2", "LE_3", "LE_4", "LT_1", "LT_2", "LT_3", "LT_4")

save_plot <- function(p, name, w, h) {
  if (is.null(p)) return(invisible(NULL))
  ggsave(file.path(out_fig, paste0(name, ".pdf")), p, width = w, height = h, dpi = 300)
  ggsave(file.path(out_fig, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
  cat("Saved:", name, "\n")
}

# =============================================================================
# Load counts & metadata
# =============================================================================

read_counts_adj <- function(path) {
  m <- as.matrix(read.table(path, header = TRUE, row.names = 1, check.names = FALSE))
  storage.mode(m) <- "integer"
  m[, sample_levels, drop = FALSE]
}

read_counts_raw <- function(path) {
  df <- read.delim(path, check.names = FALSE)
  cn <- grep("^aln/", colnames(df), value = TRUE)
  if (length(cn) != 8) stop("Expected 8 count columns in ", path)
  m <- as.matrix(df[, cn])
  colnames(m) <- sample_levels
  rownames(m) <- df$Geneid
  storage.mode(m) <- "integer"
  m
}

load_counts_for_ref <- function(ref) {
  adj <- file.path(ref$dir, ref$counts_adj)
  raw <- file.path(ref$dir, ref$counts_raw)
  if (file.exists(adj)) {
    cat("Counts [", ref$label, "]:", adj, "\n", sep = "")
    read_counts_adj(adj)
  } else if (file.exists(raw)) {
    cat("Counts [", ref$label, " raw]:", raw, "\n", sep = "")
    read_counts_raw(raw)
  } else {
    stop("No count matrix for ", ref$label, " under ", ref$dir)
  }
}

cat("=== Expression + regulation multi-panel plots ===\n")
cat("Note: LE_counts_* and LT_counts_* both contain the same 8 samples;\n")
cat("      they differ only by reference genome (NC vs CM gene IDs).\n\n")

counts_eur <- load_counts_for_ref(ref_eur)
counts_tim <- load_counts_for_ref(ref_tim)

if (!identical(colnames(counts_eur), colnames(counts_tim))) {
  stop("Sample column names differ between europaeus and timidus count matrices.")
}

meta_eur <- read.table(file.path(ref_eur$dir, ref_eur$metadata), header = TRUE, stringsAsFactors = FALSE)
meta_tim <- read.table(file.path(ref_tim$dir, ref_tim$metadata), header = TRUE, stringsAsFactors = FALSE)
if (!identical(meta_eur, meta_tim)) {
  warning("Metadata differs between reference folders; using europaeus metadata.")
}
meta <- meta_eur
rownames(meta) <- meta$sampleName
meta <- meta[colnames(counts_eur), , drop = FALSE]
meta$species_group <- ifelse(grepl("^LE", rownames(meta)), "LE", "LT")

# =============================================================================
# DESeq2 (or load existing results)
# =============================================================================

run_deseq <- function(counts, meta) {
  meta$group_genotype <- paste0(meta$group, "_", meta$Type)
  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = meta,
    design = ~ group_genotype
  )
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("group_genotype", "LT_parental_Parental", "LE_parental_Parental"))
  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)
  list(
    dds = dds,
    res = res,
    df = res_df %>%
      mutate(
        padj = replace_na(padj, 1),
        neg_log10_padj = -log10(pmax(padj, 1e-300)),
        log10_padj = log10(pmax(padj, 1e-300)),
        regulation_LT_vs_LE = case_when(
          padj < padj_thresh & log2FoldChange >=  logfc_thresh ~ "Upregulated",
          padj < padj_thresh & log2FoldChange <= -logfc_thresh ~ "Downregulated",
          TRUE ~ "Not significant"
        ),
        regulation = case_when(
          regulation_LT_vs_LE == "Upregulated" ~ "Up in LT",
          regulation_LT_vs_LE == "Downregulated" ~ "Up in LE",
          TRUE ~ "Not significant"
        )
      )
  )
}

load_deseq_csv <- function(path, flip_lfc = FALSE, deseq_source = path) {
  deseq_df <- read.csv(path, check.names = FALSE)
  if (colnames(deseq_df)[1] == "" || is.na(colnames(deseq_df)[1])) {
    colnames(deseq_df)[1] <- "gene_id"
  }
  names(deseq_df) <- tolower(names(deseq_df))
  if (names(deseq_df)[1] %in% c("", "x")) {
    names(deseq_df)[1] <- "gene_id"
  }
  if (!"gene_id" %in% names(deseq_df)) {
    deseq_df$gene_id <- deseq_df[[1]]
  }
  if (!"log2foldchange" %in% names(deseq_df)) {
    stop("log2FoldChange column not found in ", path)
  }
  deseq_df$log2FoldChange <- deseq_df$log2foldchange
  if (flip_lfc) deseq_df$log2FoldChange <- -deseq_df$log2FoldChange
  deseq_df$deseq_source <- deseq_source
  deseq_df %>%
    mutate(
      padj = replace_na(padj, 1),
      neg_log10_padj = -log10(pmax(padj, 1e-300)),
      log10_padj = log10(pmax(padj, 1e-300)),
      regulation_LT_vs_LE = case_when(
        padj < padj_thresh & log2FoldChange >=  logfc_thresh ~ "Upregulated",
        padj < padj_thresh & log2FoldChange <= -logfc_thresh ~ "Downregulated",
        TRUE ~ "Not significant"
      ),
      regulation = case_when(
        regulation_LT_vs_LE == "Upregulated" ~ "Up in LT",
        regulation_LT_vs_LE == "Downregulated" ~ "Up in LE",
        TRUE ~ "Not significant"
      )
    )
}

make_dds_for_vst <- function(counts, meta) {
  meta$group_genotype <- paste0(meta$group, "_", meta$Type)
  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = meta,
    design = ~ group_genotype
  )
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds)
  dds
}

deseq_eur_path <- file.path(ref_eur$dir, ref_eur$deseq)
deseq_tim_path <- resolve_first_existing(c(
  file.path(ref_tim$dir, ref_tim$deseq),
  file.path(deseq_base, "2.Lepus_timidus_as_reference", ref_tim$deseq)
))
if (!file.exists(deseq_eur_path)) stop("DESeq2 results not found: ", deseq_eur_path)
if (is.null(deseq_tim_path)) stop("DESeq2 results not found for timidus reference.")

cat("DESeq2 [", ref_eur$label, "]:", deseq_eur_path, "\n", sep = "")
cat("DESeq2 [", ref_tim$label, "]:", deseq_tim_path, "\n", sep = "")

de_eur <- load_deseq_csv(deseq_eur_path, flip_lfc = ref_eur$flip_lfc, deseq_source = ref_eur$label)
de_tim <- load_deseq_csv(deseq_tim_path, flip_lfc = ref_tim$flip_lfc, deseq_source = ref_tim$label)

# Genome-wide volcano: europaeus reference (NC IDs)
de_all <- de_eur

cat("Fitting dispersions for VST [europaeus ref]...\n")
dds_eur <- make_dds_for_vst(counts_eur, meta)
cat("Fitting dispersions for VST [timidus ref]...\n")
dds_tim <- make_dds_for_vst(counts_tim, meta)

vst_eur <- assay(varianceStabilizingTransformation(dds_eur, blind = TRUE))
vst_tim <- assay(varianceStabilizingTransformation(dds_tim, blind = TRUE))
vst_mat_all <- rbind(vst_eur, vst_tim[setdiff(rownames(vst_tim), rownames(vst_eur)), , drop = FALSE])

# KOBAS pathway genes
if (!file.exists(pathway_file)) {
  stop("Run kobas_enrichment_expression_direction.R first.")
}
pg <- read.csv(pathway_file, stringsAsFactors = FALSE)

# Per-gene DE for KOBAS (NC + CM), harmonized in pg — same logic as kobas_enrichment_expression_direction.R
gene_de <- pg %>%
  dplyr::distinct(gene_id, gene_symbol, log2FC, padj, regulation_LT_vs_LE, deseq_source) %>%
  mutate(
    log2FoldChange = log2FC,
    neg_log10_padj = -log10(pmax(padj, 1e-300)),
    log10_padj = log10(pmax(padj, 1e-300))
  )
de_kobas <- gene_de

write.csv(de_kobas, file.path(out_tab, "kobas_genes_expression_regulation.csv"), row.names = FALSE)

# =============================================================================
# 1. Volcano plot (all genes + highlight KOBAS)
# =============================================================================

vol <- de_all %>%
  mutate(
    highlight = ifelse(gene_id %in% unique(pg$gene_id), "KOBAS enriched (NC)", "Other"),
    label_gene = gene_id %in% (
      gene_de %>%
        filter(regulation_LT_vs_LE != "Not significant", grepl("^NC_", gene_id)) %>%
        arrange(padj) %>%
        head(15) %>%
        pull(gene_id)
    )
  )

vol_kobas_nc <- filter(vol, highlight == "KOBAS enriched (NC)") %>%
  dplyr::select(-dplyr::any_of("regulation_LT_vs_LE")) %>%
  left_join(gene_de %>% dplyr::select(gene_id, regulation_LT_vs_LE), by = "gene_id")

p_volcano <- ggplot(vol, aes(x = log2FoldChange, y = neg_log10_padj)) +
  geom_point(
    data = filter(vol, highlight == "Other"),
    color = "grey80", size = 0.6, alpha = 0.4
  ) +
  geom_point(
    data = vol_kobas_nc,
    aes(color = regulation_LT_vs_LE),
    size = 1.1, alpha = 0.75
  ) +
  scale_color_manual(
    values = c(
      "Upregulated" = "#D62728",
      "Downregulated" = "#1F77B4",
      "Not significant" = "grey50"
    ),
    name = "LT vs LE"
  ) +
  geom_vline(xintercept = c(-logfc_thresh, logfc_thresh), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(padj_thresh), linetype = "dashed", color = "grey40") +
  labs(
    title = "Volcano plot: LT vs LE (L. europaeus reference, NC genes)",
    subtitle = "Y-axis: -log10(padj) — larger values = more significant (see script note on log10(padj))",
    x = "log2 fold change (LT / LE)",
    y = expression(-log[10](padj))
  ) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

save_plot(p_volcano, "01_volcano_all_genes_KOBAS_highlight", 10, 7)

# =============================================================================
# 2. Bar plot: up / down counts (genome-wide + KOBAS subset)
# =============================================================================

count_bar <- de_all %>%
  dplyr::count(regulation_LT_vs_LE, name = "n") %>%
  mutate(dataset = "All genes (NC ref)") %>%
  bind_rows(
    gene_de %>%
      dplyr::distinct(gene_id, regulation_LT_vs_LE) %>%
      dplyr::count(regulation_LT_vs_LE, name = "n") %>%
      mutate(dataset = "KOBAS pathway genes (NC + CM)")
  ) %>%
  mutate(
    regulation_LT_vs_LE = factor(
      regulation_LT_vs_LE,
      levels = c("Upregulated", "Downregulated", "Not significant")
    )
  )

p_bar <- ggplot(count_bar, aes(x = dataset, y = n, fill = regulation_LT_vs_LE)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_text(aes(label = n), position = position_dodge(width = 0.7), vjust = -0.3, size = 3) +
  scale_fill_manual(
    values = c("Upregulated" = "#D62728", "Downregulated" = "#1F77B4", "Not significant" = "grey75"),
    name = "LT vs LE"
  ) +
  labs(title = "Upregulated vs downregulated gene counts",
       subtitle = sprintf("|log2FC| >= %.1f, padj < %.2f", logfc_thresh, padj_thresh),
       x = NULL, y = "Gene count") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

save_plot(p_bar, "02_barplot_up_down_counts", 9, 5.5)

# By KOBAS gene set
if ("GeneSet_label" %in% names(pg)) {
  bar_set <- pg %>%
    distinct(gene_id, GeneSet_label, regulation_LT_vs_LE) %>%
    dplyr::count(GeneSet_label, regulation_LT_vs_LE) %>%
    mutate(regulation_LT_vs_LE = factor(regulation_LT_vs_LE,
      levels = c("Upregulated", "Downregulated", "Not significant")))

  p_bar_set <- ggplot(bar_set, aes(x = GeneSet_label, y = n, fill = regulation_LT_vs_LE)) +
    geom_col(position = "stack") +
    scale_fill_manual(
      values = c("Upregulated" = "#D62728", "Downregulated" = "#1F77B4", "Not significant" = "grey75"),
      name = "LT vs LE"
    ) +
    labs(title = "Regulation in KOBAS input gene sets", x = NULL, y = "Unique genes") +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))

  save_plot(p_bar_set, "02b_barplot_up_down_by_KOBAS_geneset", 8, 5)
}

# =============================================================================
# 3. Heatmap (VST) — top variable KOBAS genes
# =============================================================================

if (has_pheatmap) {
  kobas_ids <- intersect(rownames(vst_mat_all), unique(pg$gene_id))
  vst_k <- vst_mat_all[kobas_ids, , drop = FALSE]

  gene_var <- apply(vst_k, 1, var)
  top_ids <- names(sort(gene_var, decreasing = TRUE))[seq_len(min(top_genes_heatmap, length(gene_var)))]
  mat <- vst_k[top_ids, , drop = FALSE]

  ann_col <- data.frame(
    Species = meta[colnames(mat), "species_group"],
    row.names = colnames(mat)
  )

  ann_row_df <- gene_de %>%
    filter(gene_id %in% top_ids) %>%
    dplyr::select(gene_id, regulation_LT_vs_LE) %>%
    dplyr::distinct(gene_id, .keep_all = TRUE)
  ann_row <- data.frame(
    Regulation = ann_row_df$regulation_LT_vs_LE[match(rownames(mat), ann_row_df$gene_id)],
    row.names = rownames(mat)
  )

  pdf(file.path(out_fig, "03_heatmap_VST_top_KOBAS_genes.pdf"), width = 9, height = 11)
  pheatmap::pheatmap(
    mat,
    scale = "row",
    annotation_col = ann_col,
    annotation_row = ann_row,
    show_rownames = nrow(mat) <= 60,
    fontsize_row = 6,
    main = "VST expression: top variable KOBAS genes (NC + CM refs)"
  )
  dev.off()
  cat("Saved: 03_heatmap_VST_top_KOBAS_genes.pdf\n")
}

# =============================================================================
# 4. Stacked bar — pathway direction (KEGG, top terms)
# =============================================================================

path_sum <- pg %>%
  filter(DB_type == "KEGG") %>%
  group_by(GeneSet_label, Term, ID) %>%
  summarise(
    n_up = sum(regulation_LT_vs_LE == "Upregulated", na.rm = TRUE),
    n_down = sum(regulation_LT_vs_LE == "Downregulated", na.rm = TRUE),
    n_ns = sum(regulation_LT_vs_LE == "Not significant" | is.na(regulation_LT_vs_LE), na.rm = TRUE),
    min_padj = min(Corrected_PValue, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(GeneSet_label) %>%
  slice_min(min_padj, n = top_path_mosaic, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(Term_short = str_trunc(Term, 40))

path_long <- path_sum %>%
  pivot_longer(c(n_up, n_down, n_ns), names_to = "dir", values_to = "n") %>%
  mutate(
    direction = recode(dir,
      n_up = "Upregulated",
      n_down = "Downregulated",
      n_ns = "Not significant"
    ),
    direction = factor(direction, levels = c("Upregulated", "Downregulated", "Not significant"))
  )

p_stack <- ggplot(path_long, aes(x = reorder(Term_short, -min_padj), y = n, fill = direction)) +
  geom_col(position = "stack", width = 0.75) +
  facet_wrap(~GeneSet_label, scales = "free_y", ncol = 1) +
  coord_flip() +
  scale_fill_manual(
    values = c("Upregulated" = "#D62728", "Downregulated" = "#1F77B4", "Not significant" = "grey75"),
    name = "LT vs LE"
  ) +
  labs(
    title = "Stacked bar: regulation within top KEGG pathways",
    x = NULL, y = "Gene count"
  ) +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold"), strip.text = element_text(face = "bold"))

save_plot(p_stack, "04_stacked_bar_pathway_direction_KEGG", 10, max(8, top_path_mosaic * 0.35))

# =============================================================================
# 5. Boxplot — VST by species (KOBAS DE genes, sample-level)
# =============================================================================

vst_df <- as.data.frame(vst_mat_all)
vst_df$gene_id <- rownames(vst_df)
vst_long <- vst_df %>%
  pivot_longer(-gene_id, names_to = "sample", values_to = "vst") %>%
  left_join({
    m2 <- meta
    m2$sample <- rownames(m2)
    m2
  }, by = "sample") %>%
  left_join(
    gene_de %>% dplyr::select(gene_id, regulation_LT_vs_LE, gene_symbol),
    by = "gene_id"
  ) %>%
  filter(regulation_LT_vs_LE != "Not significant")

top_box_genes <- gene_de %>%
  filter(regulation_LT_vs_LE != "Not significant") %>%
  arrange(padj) %>%
  head(12) %>%
  pull(gene_id)

vst_box <- vst_long %>%
  filter(gene_id %in% top_box_genes) %>%
  mutate(
    gene_label = ifelse(!is.na(gene_symbol) & gene_symbol != "", gene_symbol, gene_id),
    gene_label = paste0(gene_label, "\n(", regulation_LT_vs_LE, ")")
  )

p_box <- ggplot(vst_box, aes(x = species_group, y = vst, fill = species_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.55) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.8) +
  facet_wrap(~gene_label, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c(LE = "#1F77B4", LT = "#D62728")) +
  labs(
    title = "VST expression by species (top significant KOBAS genes)",
    x = NULL, y = "VST normalized expression"
  ) +
  theme_bw(base_size = 9) +
  theme(legend.position = "none", strip.text = element_text(size = 7))

save_plot(p_box, "05_boxplot_VST_top_KOBAS_genes", 12, 8)

# =============================================================================
# 6. Mosaic-style plot (pathway × regulation)
# =============================================================================

mosaic_df <- pg %>%
  filter(DB_type == "KEGG") %>%
  dplyr::count(GeneSet_label, Term, regulation_LT_vs_LE, name = "n") %>%
  group_by(GeneSet_label) %>%
  slice_max(n, n = top_path_mosaic, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    Term_short = str_trunc(Term, 35),
    regulation_LT_vs_LE = factor(
      regulation_LT_vs_LE,
      levels = c("Upregulated", "Downregulated", "Not significant")
    )
  )

p_mosaic <- ggplot(mosaic_df, aes(x = Term_short, y = GeneSet_label, fill = n)) +
  geom_tile(color = "white") +
  geom_text(aes(label = n), size = 2.8, color = "white") +
  facet_grid(regulation_LT_vs_LE ~ ., scales = "free", space = "free") +
  scale_fill_gradient(low = "white", high = "#7F000D", name = "Gene count") +
  labs(
    title = "Mosaic-style view: pathway × gene set × regulation",
    subtitle = "Tile area/color = number of genes (KEGG, top pathways per set)",
    x = NULL, y = NULL
  ) +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1), strip.text = element_text(face = "bold"))

save_plot(p_mosaic, "06_mosaic_pathway_regulation_KEGG", 12, 7)

# =============================================================================
# 7. Network plot — pathways ↔ genes (top KEGG terms)
# =============================================================================

if (has_ggraph) {
  net_terms <- pg %>%
    filter(DB_type == "KEGG") %>%
    group_by(GeneSet_label) %>%
    slice_min(Corrected_PValue, n = top_path_network, with_ties = FALSE) %>%
    ungroup()

  net_edges <- pg %>%
    filter(DB_type == "KEGG", Term %in% net_terms$Term) %>%
    transmute(
      from = str_trunc(Term, 35),
      to = gene_id,
      regulation_LT_vs_LE = regulation_LT_vs_LE
    )

  nodes <- tibble(
    name = unique(c(net_edges$from, net_edges$to)),
    type = ifelse(name %in% net_edges$from, "Pathway", "Gene")
  ) %>%
    left_join(
      gene_de %>%
        dplyr::select(gene_id, regulation_LT_vs_LE, gene_symbol) %>%
        mutate(name = gene_id),
      by = "name"
    )

  g <- igraph::graph_from_data_frame(net_edges %>% select(from, to), vertices = nodes, directed = FALSE)

  p_net <- ggraph::ggraph(g, layout = "fr") +
    ggraph::geom_edge_link(alpha = 0.15, colour = "grey50") +
    ggraph::geom_node_point(
      aes(color = ifelse(type == "Pathway", "Pathway", regulation_LT_vs_LE)),
      size = ifelse(nodes$type == "Pathway", 4, 1.8), alpha = 0.85
    ) +
    scale_color_manual(
      values = c(
        "Pathway" = "#333333",
        "Upregulated" = "#D62728",
        "Downregulated" = "#1F77B4",
        "Not significant" = "grey70"
      ),
      name = NULL, na.value = "grey70"
    ) +
    ggraph::geom_node_text(
      aes(label = ifelse(type == "Pathway", name, gene_symbol)),
      size = 2.5, repel = TRUE, max.overlaps = 25
    ) +
    labs(title = "Gene–pathway network (top KEGG terms)") +
    theme_void(base_size = 10)

  save_plot(p_net, "07_network_pathway_gene_KEGG", 12, 10)
} else {
  cat("Skip network plot (install igraph + ggraph).\n")
}

# =============================================================================
# log10(padj) explanation
# =============================================================================

cat("\n")
cat(strrep("=", 72), "\n", sep = "")
cat("WHY -log10(padj) ON VOLCANO PLOTS?\n", sep = "")
cat(strrep("=", 72), "\n", sep = "")
cat("
padj is an adjusted p-value between 0 and 1 (smaller = more significant).

1) -log10(padj)  [used in 01_volcano_all_genes_KOBAS_highlight.pdf — STANDARD]
   - Example: padj = 0.05  -> -log10(padj) = 1.30
   - Example: padj = 1e-10 -> -log10(padj) = 10
   - MORE significant genes are plotted HIGHER on the Y-axis.
   - This is what most RNA-seq tools use (DESeq2 plotMA, EnhancedVolcano, etc.).
For presentations, prefer -log10(padj).
", sep = "")
cat(strrep("=", 72), "\n", sep = "")
cat("\nDone. Figures:", out_fig, "\n")
