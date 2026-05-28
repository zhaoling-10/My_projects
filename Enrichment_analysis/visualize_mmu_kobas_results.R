#!/usr/bin/env Rscript
# =============================================================================
# Visualize Mus musculus (mmu) KOBAS enrichment results
#
# Reads KOBAS *_identify.txt outputs (GO + KEGG × 3 gene sets).
# Produces dot/bubble plots and bar charts (publication-style).
#
# Run after KOBAS; independent of mmu_enrichment_expression_analysis.R
# For up/down regulation of these same enriched genes, run:
#   kobas_enrichment_expression_direction.R
# Outputs: mmu_analysis/figures/
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(readr)
})

if (requireNamespace("pheatmap", quietly = TRUE)) {
  has_pheatmap <- TRUE
} else {
  has_pheatmap <- FALSE
  message("Note: install 'pheatmap' for enrichment heatmaps (optional).")
}

# ---- Configuration ----
base_dir   <- path.expand("~/Documents/Enrichment_analysis")
out_dir    <- file.path(base_dir, "mmu_analysis", "figures")
fdr_thresh <- 0.05
top_n_dot  <- 20   # terms per gene set in dot plots
top_n_bar  <- 15
top_n_heat <- 25

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

identify_files <- list(
  list(path = file.path(base_dir, "KOBAS_output_mmu_KEGG_1",
                        "KOBAS_output_mmu_KEGG_1.tsv_identify.txt"),
       gs = "Shared_508",    db = "KEGG"),
  list(path = file.path(base_dir, "KOBAS_output_mmu_KEGG_2",
                        "KOBAS_output_mmu_KEGG_2.tsv_identify.txt"),
       gs = "LE_unique_158", db = "KEGG"),
  list(path = file.path(base_dir, "KOBAS_output_mmu_KEGG_3",
                        "KOBAS_output_mmu_KEGG_3.tsv_identify.txt"),
       gs = "LT_unique_133", db = "KEGG"),
  list(path = file.path(base_dir, "KOBAS_output_mmu_GO_1",
                        "KOBAS_output_mmu_GO_1.tsv_identify.txt"),
       gs = "Shared_508",    db = "GO"),
  list(path = file.path(base_dir, "KOBAS_output_mmu_GO_2",
                        "KOBAS_output_mmu_GO_2.tsv_identify.txt"),
       gs = "LE_unique_158", db = "GO"),
  list(path = file.path(base_dir, "KOBAS_output_mmu_GO_3",
                        "KOBAS_output_mmu_GO_3.tsv_identify.txt"),
       gs = "LT_unique_133", db = "GO")
)

gene_set_levels <- c("Shared_508", "LE_unique_158", "LT_unique_133")
gene_set_labels <- c(
  Shared_508    = "508 shared",
  LE_unique_158 = "158 LE-unique",
  LT_unique_133 = "133 LT-unique"
)

# =============================================================================
# Parse KOBAS identify files
# =============================================================================

parse_kobas_identify <- function(filepath, gene_set, db_type) {
  if (!file.exists(filepath)) {
    warning("Not found: ", filepath)
    return(NULL)
  }
  lines <- readLines(filepath)
  data_lines <- lines[!grepl("^#", lines) & nchar(trimws(lines)) > 0]
  if (length(data_lines) == 0) return(NULL)

  df <- read.table(
    text = paste(data_lines, collapse = "\n"),
    sep = "\t", header = FALSE, quote = "", fill = TRUE,
    stringsAsFactors = FALSE
  )
  n <- min(ncol(df), 8)
  colnames(df)[seq_len(n)] <- c(
    "Term", "Database", "ID", "Input_n", "Background_n",
    "PValue", "Corrected_PValue", "Input_genes"
  )[seq_len(n)]

  df %>%
    mutate(
      GeneSet = gene_set,
      DB_type = db_type,
      RefSpecies = "mmu",
      Input_n = as.numeric(Input_n),
      Background_n = as.numeric(Background_n),
      PValue = as.numeric(PValue),
      Corrected_PValue = as.numeric(Corrected_PValue),
      GeneRatio = Input_n / pmax(Background_n, 1),
      neg_log10_padj = -log10(pmax(Corrected_PValue, .Machine$double.xmin)),
      GeneSet_label = gene_set_labels[gene_set]
    )
}

cat("=== mmu KOBAS enrichment visualization ===\n")
cat("Output:", out_dir, "\n\n")

all_results <- bind_rows(lapply(identify_files, function(x) {
  parse_kobas_identify(x$path, x$gs, x$db)
}))

if (is.null(all_results) || nrow(all_results) == 0) {
  stop("No KOBAS identify data loaded. Check mmu KOBAS output folders.")
}

sig <- all_results %>%
  filter(!is.na(Corrected_PValue), Corrected_PValue <= fdr_thresh) %>%
  mutate(
    GeneSet = factor(GeneSet, levels = gene_set_levels),
    GeneSet_label = factor(GeneSet_label, levels = gene_set_labels[gene_set_levels]),
    Term_short = str_trunc(Term, 50)
  )

write.csv(
  sig %>% select(GeneSet, GeneSet_label, DB_type, Term, ID, Input_n,
                 Background_n, GeneRatio, PValue, Corrected_PValue),
  file.path(dirname(out_dir), "mmu_significant_enrichment_FDR05.csv"),
  row.names = FALSE
)

cat("Significant terms (FDR <=", fdr_thresh, "):\n")
print(table(sig$DB_type, sig$GeneSet_label))
cat("\n")

save_plot <- function(p, name, width, height) {
  if (is.null(p)) return(invisible(NULL))
  ggsave(file.path(out_dir, paste0(name, ".pdf")), p,
         width = width, height = height, dpi = 300)
  ggsave(file.path(out_dir, paste0(name, ".png")), p,
         width = width, height = height, dpi = 300)
  cat("Saved:", name, ".pdf / .png\n")
}

# =============================================================================
# 1. Dot plot (clusterProfiler-style): GeneRatio vs pathway, one panel per gene set
# =============================================================================

plot_dotplot <- function(data, db_type, top_n = top_n_dot) {
  df <- data %>%
    filter(DB_type == db_type) %>%
    group_by(GeneSet) %>%
    slice_min(Corrected_PValue, n = top_n, with_ties = FALSE) %>%
    ungroup()

  if (nrow(df) == 0) {
    message("No significant ", db_type, " terms for dot plot.")
    return(NULL)
  }

  df <- df %>%
    group_by(GeneSet) %>%
    mutate(
      Term_ordered = factor(Term_short,
                            levels = rev(unique(Term_short[order(Corrected_PValue)])))
    ) %>%
    ungroup()

  ggplot(df, aes(x = GeneRatio, y = Term_ordered)) +
    geom_point(aes(size = Input_n, color = neg_log10_padj), alpha = 0.9) +
    facet_wrap(~GeneSet_label, ncol = 1, scales = "free_y") +
    scale_color_gradient(
      low = "#56B1F7", high = "#132B43",
      name = expression(-log[10](FDR))
    ) +
    scale_size_continuous(name = "Gene count", range = c(2.5, 10)) +
    labs(
      title = paste0(db_type, " enrichment — Mus musculus reference"),
      subtitle = paste0("Top ", top_n, " terms per gene set (FDR <= ", fdr_thresh, ")"),
      x = "Gene ratio (input / background)",
      y = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
}

# =============================================================================
# 2. Bubble plot: pathways × gene sets (size = count, color = FDR)
# =============================================================================

plot_bubble <- function(data, db_type, top_n = top_n_dot) {
  top_terms <- data %>%
    filter(DB_type == db_type) %>%
    group_by(Term) %>%
    slice_min(Corrected_PValue, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(Corrected_PValue) %>%
    slice_head(n = top_n) %>%
    pull(Term)

  if (length(top_terms) == 0) return(NULL)

  df_all <- data %>%
    filter(DB_type == db_type, Term %in% top_terms) %>%
    mutate(Term_short = factor(Term_short, levels = rev(unique(Term_short))))

  ggplot(df_all, aes(x = GeneSet_label, y = Term_short)) +
    geom_point(aes(size = Input_n, color = neg_log10_padj), alpha = 0.88) +
    scale_color_gradient(
      low = "#56B1F7", high = "#132B43",
      name = expression(-log[10](FDR))
    ) +
    scale_size_continuous(name = "Gene count", range = c(3, 12)) +
    labs(
      title = paste0(db_type, " pathway bubble plot — mmu reference"),
      subtitle = paste0("Top ", length(top_terms), " pathways across gene sets"),
      x = "Gene set", y = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 30, hjust = 1),
      axis.text.y = element_text(size = 9)
    )
}

# =============================================================================
# 3. Horizontal bar chart: -log10(FDR) per term
# =============================================================================

plot_bar <- function(data, db_type, top_n = top_n_bar) {
  df <- data %>%
    filter(DB_type == db_type) %>%
    group_by(GeneSet) %>%
    slice_min(Corrected_PValue, n = top_n, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(Term_short = factor(Term_short, levels = unique(Term_short)))

  if (nrow(df) == 0) return(NULL)

  ggplot(df, aes(x = neg_log10_padj, y = reorder(Term_short, neg_log10_padj),
                 fill = GeneSet_label)) +
    geom_col(position = "dodge") +
    facet_wrap(~GeneSet_label, scales = "free_y", ncol = 1) +
    scale_fill_brewer(palette = "Set2", name = "Gene set") +
    labs(
      title = paste0("Top ", db_type, " terms — mmu reference"),
      x = expression(-log[10](FDR)),
      y = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "none",
      strip.text = element_text(face = "bold")
    )
}

# =============================================================================
# 4. Combined publication dot plot (KEGG, all gene sets)
# =============================================================================

plot_publication_dot <- function(data, db_type = "KEGG", top_n = 8) {
  df <- data %>%
    filter(DB_type == db_type) %>%
    group_by(GeneSet) %>%
    slice_min(Corrected_PValue, n = top_n, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(
      Term_short = factor(Term_short,
                          levels = rev(unique(Term_short[order(Corrected_PValue)])))
    )

  if (nrow(df) == 0) return(NULL)

  ggplot(df, aes(x = GeneSet_label, y = Term_short)) +
    geom_point(aes(size = Input_n, color = neg_log10_padj), alpha = 0.92) +
    scale_color_gradient2(
      low = "steelblue", mid = "#FFC107", high = "#D62728",
      midpoint = median(df$neg_log10_padj, na.rm = TRUE),
      name = expression(-log[10](FDR))
    ) +
    scale_size_continuous(range = c(3, 14), name = "Gene count") +
    labs(
      title = paste0("Top enriched ", db_type, " pathways (mmu)"),
      subtitle = "Lepus europaeus vs Lepus timidus (skin transcriptome)",
      x = NULL, y = NULL
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.text.y = element_text(size = 9),
      axis.text.x = element_text(size = 11, face = "bold"),
      panel.grid.major.x = element_blank()
    )
}

# =============================================================================
# Generate and save figures
# =============================================================================

for (db in c("KEGG", "GO")) {
  db_lower <- tolower(db)

  p_dot <- plot_dotplot(sig, db)
  save_plot(p_dot, paste0(db_lower, "_dotplot_mmu_by_geneset"),
            width = 10, height = max(8, top_n_dot * 0.35 * 3))

  p_bub <- plot_bubble(sig, db)
  save_plot(p_bub, paste0(db_lower, "_bubble_mmu_combined"),
            width = 11, height = max(7, top_n_dot * 0.28))

  p_bar <- plot_bar(sig, db)
  save_plot(p_bar, paste0(db_lower, "_barplot_mmu_by_geneset"),
            width = 10, height = max(8, top_n_bar * 0.4 * 3))
}

p_pub_kegg <- plot_publication_dot(sig, "KEGG", top_n = 8)
save_plot(p_pub_kegg, "KEGG_publication_dotplot_mmu", width = 13, height = 9)

p_pub_go <- plot_publication_dot(sig, "GO", top_n = 8)
save_plot(p_pub_go, "GO_publication_dotplot_mmu", width = 13, height = 10)

# Optional heatmap (KEGG)
if (has_pheatmap) {
  mat_df <- sig %>%
    filter(DB_type == "KEGG") %>%
    group_by(Term) %>%
    slice_min(Corrected_PValue, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(Corrected_PValue) %>%
    slice_head(n = top_n_heat)

  if (nrow(mat_df) >= 2) {
    mat <- sig %>%
      filter(DB_type == "KEGG", Term %in% mat_df$Term) %>%
      select(Term_short, GeneSet_label, neg_log10_padj) %>%
      pivot_wider(names_from = GeneSet_label, values_from = neg_log10_padj,
                  values_fill = 0) %>%
      as.data.frame()
    rownames(mat) <- mat$Term_short
    mat$Term_short <- NULL
    mat <- as.matrix(mat)

    pdf(file.path(out_dir, "KEGG_heatmap_mmu.pdf"), width = 8, height = 10)
    pheatmap::pheatmap(
      mat,
      color = colorRampPalette(c("white", "#56B1F7", "#132B43"))(100),
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      main = "KEGG enrichment (mmu, -log10 FDR)",
      fontsize_row = 8
    )
    dev.off()
    cat("Saved: KEGG_heatmap_mmu.pdf\n")
  }
}

cat("\nDone. KOBAS figures in:", out_dir, "\n")
cat("  Dot plots:    *_dotplot_mmu_by_geneset.pdf\n")
cat("  Bubble plots: *_bubble_mmu_combined.pdf\n")
cat("  Bar charts:   *_barplot_mmu_by_geneset.pdf\n")
cat("  Publication:  KEGG/GO_publication_dotplot_mmu.pdf\n")
