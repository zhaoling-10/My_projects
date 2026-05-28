#!/usr/bin/env Rscript
# =============================================================================
# KOBAS mmu enrichment plots — UPREGULATED genes only (LT vs LE)
#
# Uses the same significant pathways as visualize_mmu_kobas_results.R, but
# counts and plots only genes with regulation_LT_vs_LE == "Upregulated"
# (higher in Lepus timidus; log2FC >= threshold, padj < threshold).
#
# Also writes comparison figures: all enriched genes vs upregulated-only.
#
# Prerequisite: run kobas_enrichment_expression_direction.R once, OR this
# script will load DESeq2 + KOBAS identify files itself.
#
# Outputs:
#   mmu_analysis/tables/upregulated_only/
#   mmu_analysis/figures/upregulated_only/
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(readr)
})

base_dir   <- path.expand("~/Documents/Enrichment_analysis")
table_dir  <- file.path(base_dir, "mmu_analysis", "tables", "upregulated_only")
fig_dir    <- file.path(base_dir, "mmu_analysis", "figures", "upregulated_only")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

fdr_thresh   <- 0.05
padj_thresh  <- 0.05
logfc_thresh <- 1.0
top_n_dot    <- 20
top_n_bar    <- 15
top_n_path   <- 12

gene_set_levels <- c("Shared_508", "LE_unique_158", "LT_unique_133")
gene_set_labels <- c(
  Shared_508    = "508 shared",
  LE_unique_158 = "158 LE-unique",
  LT_unique_133 = "133 LT-unique"
)

save_plot <- function(p, name, width, height, subdir = fig_dir) {
  if (is.null(p)) return(invisible(NULL))
  ggsave(file.path(subdir, paste0(name, ".pdf")), p, width = width, height = height, dpi = 300)
  ggsave(file.path(subdir, paste0(name, ".png")), p, width = width, height = height, dpi = 300)
  cat("Saved:", name, ".pdf / .png\n")
}

# =============================================================================
# Load pathway-gene table (prefer cached output from direction script)
# =============================================================================

pathway_file <- file.path(
  base_dir, "mmu_analysis", "tables", "kobas_all_pathway_genes_with_direction.csv"
)

if (file.exists(pathway_file)) {
  cat("Loading:", pathway_file, "\n")
  pathway_genes <- read.csv(pathway_file, stringsAsFactors = FALSE)
  if (!"regulation_LT_vs_LE" %in% names(pathway_genes)) {
    stop("Run kobas_enrichment_expression_direction.R first (missing regulation_LT_vs_LE).")
  }
} else {
  cat("Cached table not found — run:\n",
      "  Rscript kobas_enrichment_expression_direction.R\n", sep = "")
  stop("Required file missing: ", pathway_file, call. = FALSE)
}

pathway_genes <- pathway_genes %>%
  mutate(
    GeneSet = factor(GeneSet, levels = gene_set_levels),
    GeneSet_label = factor(
      gene_set_labels[as.character(GeneSet)],
      levels = gene_set_labels[gene_set_levels]
    ),
    Term_short = str_trunc(Term, 50),
    neg_log10_padj = -log10(pmax(Corrected_PValue, .Machine$double.xmin))
  )

pathway_up <- pathway_genes %>%
  filter(regulation_LT_vs_LE == "Upregulated")

cat("\n=== Upregulated-only KOBAS plots ===\n")
cat("Total gene-pathway rows:", nrow(pathway_genes), "\n")
cat("Upregulated rows:", nrow(pathway_up), "\n")
cat("Unique upregulated genes:", n_distinct(pathway_up$gene_id), "\n\n")

# =============================================================================
# Tables: term stats (upregulated count vs original Input_n)
# =============================================================================

term_stats <- pathway_genes %>%
  group_by(GeneSet, GeneSet_label, DB_type, Term, ID, Corrected_PValue, Input_n) %>%
  summarise(
    n_genes_total_in_term = first(Input_n),
    n_upregulated = sum(regulation_LT_vs_LE == "Upregulated", na.rm = TRUE),
    n_downregulated = sum(regulation_LT_vs_LE == "Downregulated", na.rm = TRUE),
    n_not_significant = sum(regulation_LT_vs_LE == "Not significant", na.rm = TRUE),
    pct_upregulated = round(100 * n_upregulated / n(), 1),
    mean_log2FC_up = round(mean(log2FC[regulation_LT_vs_LE == "Upregulated"], na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  mutate(
    neg_log10_padj = -log10(pmax(Corrected_PValue, .Machine$double.xmin)),
    GeneRatio_up = n_upregulated / pmax(n_genes_total_in_term, 1),
    Term_short = str_trunc(Term, 50)
  ) %>%
  filter(n_upregulated > 0) %>%
  arrange(DB_type, GeneSet, Corrected_PValue)

write.csv(term_stats,
          file.path(table_dir, "kobas_terms_upregulated_only.csv"),
          row.names = FALSE)
cat("Saved: kobas_terms_upregulated_only.csv\n")

write.csv(
  pathway_up %>%
    select(GeneSet, GeneSet_label, DB_type, Term, ID, gene_id, gene_symbol,
           log2FC, padj, regulation, regulation_LT_vs_LE),
  file.path(table_dir, "kobas_pathway_genes_upregulated_only.csv"),
  row.names = FALSE
)
cat("Saved: kobas_pathway_genes_upregulated_only.csv\n")

# =============================================================================
# 1. Bar: unique upregulated genes per gene set
# =============================================================================

by_set_up <- pathway_up %>%
  distinct(gene_id, GeneSet, GeneSet_label) %>%
  count(GeneSet, GeneSet_label, name = "n_upregulated")

by_set_all <- pathway_genes %>%
  distinct(gene_id, GeneSet, GeneSet_label, regulation_LT_vs_LE) %>%
  count(GeneSet, GeneSet_label, regulation_LT_vs_LE, name = "n") %>%
  pivot_wider(names_from = regulation_LT_vs_LE, values_from = n, values_fill = 0)

write.csv(by_set_up, file.path(table_dir, "upregulated_count_by_geneset.csv"), row.names = FALSE)

p_counts <- ggplot(by_set_up, aes(x = GeneSet_label, y = n_upregulated, fill = GeneSet_label)) +
  geom_col(width = 0.65, show.legend = FALSE) +
  geom_text(aes(label = n_upregulated), vjust = -0.4, size = 3.5) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Upregulated genes in significant KOBAS enrichments (mmu)",
    subtitle = sprintf(
      "LT vs LE: higher in mountain hare | |log2FC| >= %.1f, padj < %.2f",
      logfc_thresh, padj_thresh
    ),
    x = NULL, y = "Unique upregulated genes"
  ) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 25, hjust = 1))

save_plot(p_counts, "01_upregulated_count_by_geneset", 8, 5)

# =============================================================================
# 2. Comparison: All regulation types vs upregulated-only (side-by-side bars)
# =============================================================================

cmp_long <- pathway_genes %>%
  distinct(gene_id, GeneSet_label, regulation_LT_vs_LE) %>%
  count(GeneSet_label, regulation_LT_vs_LE, name = "n") %>%
  mutate(
    plot_type = "All genes in enrichment",
    regulation_LT_vs_LE = factor(
      regulation_LT_vs_LE,
      levels = c("Upregulated", "Downregulated", "Not significant")
    )
  )

cmp_up <- by_set_up %>%
  mutate(
    plot_type = "Upregulated only",
    regulation_LT_vs_LE = factor("Upregulated", levels = c("Upregulated", "Downregulated", "Not significant")),
    n = n_upregulated
  ) %>%
  select(GeneSet_label, regulation_LT_vs_LE, n, plot_type)

cmp_df <- bind_rows(cmp_long, cmp_up) %>%
  mutate(
    plot_type = factor(plot_type, levels = c("All genes in enrichment", "Upregulated only"))
  )

p_cmp <- ggplot(cmp_df, aes(x = GeneSet_label, y = n, fill = regulation_LT_vs_LE)) +
  geom_col(position = "stack", width = 0.7) +
  facet_wrap(~plot_type, ncol = 1) +
  scale_fill_manual(
    values = c(
      "Upregulated" = "#D62728",
      "Downregulated" = "#1F77B4",
      "Not significant" = "grey75"
    ),
    name = "LT vs LE"
  ) +
  labs(
    title = "Comparison: all enriched genes vs upregulated-only",
    subtitle = "Same KOBAS pathways (FDR <= 0.05); right panel is upregulated subset only",
    x = NULL, y = "Unique gene count"
  ) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 25, hjust = 1),
        strip.text = element_text(face = "bold"))

save_plot(p_cmp, "02_comparison_all_vs_upregulated_only", 9, 8)

# =============================================================================
# 3. Dot plot — upregulated gene count per pathway (KEGG + GO)
# =============================================================================

plot_dotplot_up <- function(data, db_type, top_n = top_n_dot) {
  df <- data %>%
    filter(DB_type == db_type) %>%
    group_by(GeneSet) %>%
    slice_max(order_by = n_upregulated, n = top_n, with_ties = FALSE) %>%
    ungroup()

  if (nrow(df) == 0) return(NULL)

  df <- df %>%
    group_by(GeneSet) %>%
    mutate(
      Term_ordered = factor(
        Term_short,
        levels = rev(Term_short[order(n_upregulated)])
      )
    ) %>%
    ungroup()

  ggplot(df, aes(x = GeneRatio_up, y = Term_ordered)) +
    geom_point(aes(size = n_upregulated, color = neg_log10_padj), alpha = 0.9) +
    facet_wrap(~GeneSet_label, ncol = 1, scales = "free_y") +
    scale_color_gradient(low = "#FEC44F", high = "#7F000D", name = expression(-log[10](FDR))) +
    scale_size_continuous(name = "Upregulated\ngene count", range = c(2.5, 11)) +
    labs(
      title = paste0(db_type, " enrichment — upregulated genes only"),
      subtitle = paste0("Top ", top_n, " pathways by # upregulated genes (FDR <= ", fdr_thresh, ")"),
      x = "Fraction upregulated (n_up / KOBAS Input_n)",
      y = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold"), strip.text = element_text(face = "bold"))
}

for (db in c("KEGG", "GO")) {
  p <- plot_dotplot_up(term_stats, db)
  save_plot(p, paste0("03_dotplot_upregulated_", tolower(db)),
            width = 10, height = max(8, top_n_dot * 0.35 * 3))
}

# =============================================================================
# 4. Bubble plot — upregulated counts across gene sets
# =============================================================================

plot_bubble_up <- function(data, db_type, top_n = top_n_dot) {
  top_terms <- data %>%
    filter(DB_type == db_type) %>%
    group_by(Term) %>%
    slice_max(n_upregulated, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(desc(n_upregulated)) %>%
    slice_head(n = top_n) %>%
    pull(Term)

  if (length(top_terms) == 0) return(NULL)

  df_all <- data %>%
    filter(DB_type == db_type, Term %in% top_terms) %>%
    mutate(
      Term_label = ifelse(
        duplicated(Term_short) | duplicated(Term_short, fromLast = TRUE),
        paste0(Term_short, " (", ID, ")"),
        Term_short
      ),
      Term_label = factor(
        Term_label,
        levels = rev(unique(Term_label[order(n_upregulated)]))
      )
    )

  ggplot(df_all, aes(x = GeneSet_label, y = Term_label)) +
    geom_point(aes(size = n_upregulated, color = neg_log10_padj), alpha = 0.88) +
    scale_color_gradient(low = "#FEC44F", high = "#7F000D", name = expression(-log[10](FDR))) +
    scale_size_continuous(name = "Upregulated\ngenes", range = c(2, 14)) +
    labs(
      title = paste0(db_type, " bubble plot — upregulated genes only"),
      subtitle = paste0("Top ", length(top_terms), " pathways (by upregulated gene count)"),
      x = "Gene set", y = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold"),
          axis.text.x = element_text(angle = 30, hjust = 1),
          axis.text.y = element_text(size = 9))
}

for (db in c("KEGG", "GO")) {
  p <- plot_bubble_up(term_stats, db)
  save_plot(p, paste0("04_bubble_upregulated_", tolower(db)), 11, max(7, top_n_dot * 0.28))
}

# =============================================================================
# 5. Horizontal bar — top pathways by upregulated gene count (KEGG)
# =============================================================================

plot_bar_pathways_up <- function(data, db_type, top_n = top_n_bar) {
  df <- data %>%
    filter(DB_type == db_type) %>%
    group_by(GeneSet) %>%
    slice_max(n_upregulated, n = top_n, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(
      Term_short = factor(Term_short, levels = unique(Term_short[order(n_upregulated)]))
    )

  if (nrow(df) == 0) return(NULL)

  ggplot(df, aes(x = n_upregulated, y = Term_short, fill = GeneSet_label)) +
    geom_col() +
    facet_wrap(~GeneSet_label, ncol = 1, scales = "free_y") +
    scale_fill_brewer(palette = "Set2", guide = "none") +
    labs(
      title = paste0("Top ", db_type, " pathways — upregulated gene count"),
      x = "Upregulated genes in pathway",
      y = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "bold"), strip.text = element_text(face = "bold"))
}

for (db in c("KEGG", "GO")) {
  p <- plot_bar_pathways_up(term_stats, db)
  save_plot(p, paste0("05_barplot_pathways_upregulated_", tolower(db)),
            10, max(8, top_n_bar * 0.4 * 3))
}

# =============================================================================
# 6. Top KEGG pathways — upregulated genes per pathway (faceted bars)
# =============================================================================

top_kegg <- term_stats %>%
  filter(DB_type == "KEGG") %>%
  group_by(GeneSet_label) %>%
  slice_min(Corrected_PValue, n = top_n_path, with_ties = FALSE) %>%
  ungroup()

if (nrow(top_kegg) > 0) {
  plot_kegg <- pathway_up %>%
    filter(
      DB_type == "KEGG",
      Term %in% top_kegg$Term,
      GeneSet_label %in% top_kegg$GeneSet_label
    ) %>%
    count(GeneSet_label, Term, Term_short = str_trunc(Term, 42), name = "n_upregulated") %>%
    right_join(
      top_kegg %>% select(GeneSet_label, Term, Term_short, Corrected_PValue),
      by = c("GeneSet_label", "Term", "Term_short")
    ) %>%
    mutate(n_upregulated = replace_na(n_upregulated, 0L))

  p_kegg_facets <- ggplot(plot_kegg, aes(x = reorder(Term_short, -Corrected_PValue), y = n_upregulated)) +
    geom_col(fill = "#D62728", width = 0.75) +
    facet_wrap(~GeneSet_label, ncol = 1, scales = "free_y") +
    coord_flip() +
    labs(
      title = "Top KEGG pathways: upregulated genes only",
      subtitle = "One bar = count of upregulated genes assigned to that pathway",
      x = NULL, y = "Upregulated gene count"
    ) +
    theme_bw(base_size = 9) +
    theme(plot.title = element_text(face = "bold"), strip.text = element_text(face = "bold"))

  save_plot(p_kegg_facets, "06_top_KEGG_pathways_upregulated_bars",
            10, max(8, nrow(top_kegg) * 0.22))
}

# =============================================================================
# 7. Pathway-level comparison: total Input_n vs n_upregulated (dot)
# =============================================================================

cmp_terms <- term_stats %>%
  filter(DB_type == "KEGG") %>%
  group_by(GeneSet_label) %>%
  slice_min(Corrected_PValue, n = 10, with_ties = FALSE) %>%
  ungroup()

if (nrow(cmp_terms) > 0) {
  p_scatter <- ggplot(cmp_terms, aes(x = n_genes_total_in_term, y = n_upregulated)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(aes(color = GeneSet_label, size = neg_log10_padj), alpha = 0.85) +
    scale_color_brewer(palette = "Set2", name = "Gene set") +
    scale_size_continuous(name = expression(-log[10](FDR)), range = c(2, 8)) +
    facet_wrap(~GeneSet_label, ncol = 3) +
    labs(
      title = "KEGG pathways: total enriched genes vs upregulated subset",
      subtitle = "Points below diagonal = pathway has downregulated or NS genes",
      x = "KOBAS Input_n (all genes in pathway)",
      y = "Upregulated genes in pathway"
    ) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

  save_plot(p_scatter, "07_KEGG_total_vs_upregulated_scatter", 12, 5)
}

cat("\n--- Upregulated genes per gene set ---\n")
print(by_set_up)

cat("\nDone.\n")
cat("  Tables:", table_dir, "\n")
cat("  Figures:", fig_dir, "\n")
cat("  Key comparison: 02_comparison_all_vs_upregulated_only.pdf\n")
