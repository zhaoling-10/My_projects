# =============================================================================
# KOBAS Enrichment Analysis - Lepus europaeus vs Lepus timidus
# Gene sets: 508 shared, 158 LE-unique, 133 LT-unique
# Reference species: mmu (mouse), ocu (rabbit)
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(ggrepel)
library(RColorBrewer)
library(pheatmap)
library(stringr)
library(gridExtra)
library(reshape2)

# =============================================================================
# SECTION 1: HELPER FUNCTIONS
# =============================================================================

# Parse KOBAS _identify.txt enrichment result
parse_kobas_identify <- function(filepath, gene_set_label, ref_species) {
  if (!file.exists(filepath)) {
    warning("File not found: ", filepath)
    return(NULL)
  }
  lines <- readLines(filepath)
  # Skip comment lines starting with #
  data_lines <- lines[!grepl("^#", lines) & nchar(trimws(lines)) > 0]
  if (length(data_lines) == 0) return(NULL)

  df <- read.table(
    text   = paste(data_lines, collapse = "\n"),
    sep    = "\t",
    header = FALSE,
    quote  = "",
    fill   = TRUE,
    stringsAsFactors = FALSE
  )
  # KOBAS identify columns:
  # Term | Database | ID | Input number | Background number | P-Value | Corrected P-Value | Input | Hyperlink
  cols_needed <- min(ncol(df), 8)
  colnames(df)[1:cols_needed] <- c("Term","Database","ID","Input_n","Background_n",
                                    "PValue","Corrected_PValue","Input_genes")[1:cols_needed]
  df$GeneSet      <- gene_set_label
  df$RefSpecies   <- ref_species
  df$Input_n      <- as.numeric(df$Input_n)
  df$Background_n <- as.numeric(df$Background_n)
  df$PValue       <- as.numeric(df$PValue)
  df$Corrected_PValue <- as.numeric(df$Corrected_PValue)
  df$GeneRatio    <- df$Input_n / df$Background_n
  df
}

# =============================================================================
# SECTION 2: LOAD DATA
# Build the file path vectors to match your folder structure.
# Adjust base_dir to wherever your KOBAS output folders live.
# =============================================================================

base_dir <- "~/Documents/Enrichment_analysis"   # <-- CHANGE if needed

# ------ define all 12 result files (3 gene sets × 2 species × 2 databases) ------
identify_files <- list(
  # mmu – KEGG
  list(path = file.path(base_dir, "KOBAS_output_mmu_KEGG_1",
                        "KOBAS_output_mmu_KEGG_1.tsv_identify.txt"),
       gene_set = "Shared_508", ref = "mmu", db = "KEGG"),
  list(path = file.path(base_dir, "KOBAS_output_mmu_KEGG_2",
                        "KOBAS_output_mmu_KEGG_2.tsv_identify.txt"),
       gene_set = "LE_unique_158", ref = "mmu", db = "KEGG"),
  list(path = file.path(base_dir, "KOBAS_output_mmu_KEGG_3",
                        "KOBAS_output_mmu_KEGG_3.tsv_identify.txt"),
       gene_set = "LT_unique_133", ref = "mmu", db = "KEGG"),
  # mmu – GO
  list(path = file.path(base_dir, "KOBAS_output_mmu_GO_1",
                        "KOBAS_output_mmu_GO_1.tsv_identify.txt"),
       gene_set = "Shared_508", ref = "mmu", db = "GO"),
  list(path = file.path(base_dir, "KOBAS_output_mmu_GO_2",
                        "KOBAS_output_mmu_GO_2.tsv_identify.txt"),
       gene_set = "LE_unique_158", ref = "mmu", db = "GO"),
  list(path = file.path(base_dir, "KOBAS_output_mmu_GO_3",
                        "KOBAS_output_mmu_GO_3.tsv_identify.txt"),
       gene_set = "LT_unique_133", ref = "mmu", db = "GO"),
  # ocu – KEGG
  list(path = file.path(base_dir, "KOBAS_output_ocu_KEGG_1",
                        "KOBAS_output_ocu_KEGG_1.tsv_identify.txt"),
       gene_set = "Shared_508", ref = "ocu", db = "KEGG"),
  list(path = file.path(base_dir, "KOBAS_output_ocu_KEGG_2",
                        "KOBAS_output_ocu_KEGG_2.tsv_identify.txt"),
       gene_set = "LE_unique_158", ref = "ocu", db = "KEGG"),
  list(path = file.path(base_dir, "KOBAS_output_ocu_KEGG_3",
                        "KOBAS_output_ocu_KEGG_3.tsv_identify.txt"),
       gene_set = "LT_unique_133", ref = "ocu", db = "KEGG"),
  # ocu – GO
  list(path = file.path(base_dir, "KOBAS_output_ocu_GO_1",
                        "KOBAS_output_ocu_GO_1.tsv_identify.txt"),
       gene_set = "Shared_508", ref = "ocu", db = "GO"),
  list(path = file.path(base_dir, "KOBAS_output_ocu_GO_2",
                        "KOBAS_output_ocu_GO_2.tsv_identify.txt"),
       gene_set = "LE_unique_158", ref = "ocu", db = "GO"),
  list(path = file.path(base_dir, "KOBAS_output_ocu_GO_3",
                        "KOBAS_output_ocu_GO_3.tsv_identify.txt"),
       gene_set = "LT_unique_133", ref = "ocu", db = "GO")
)

# Load and combine
all_results <- bind_rows(lapply(identify_files, function(x) {
  df <- parse_kobas_identify(x$path, x$gene_set, x$ref)
  if (!is.null(df)) df$DB_type <- x$db
  df
}))

# Filter significant (FDR ≤ 0.05)
sig_results <- all_results %>%
  filter(!is.na(Corrected_PValue) & Corrected_PValue <= 0.05) %>%
  mutate(
    neg_log10_padj = -log10(Corrected_PValue),
    GeneSet = factor(GeneSet,
                     levels = c("Shared_508", "LE_unique_158", "LT_unique_133"),
                     labels = c("Shared (508)", "LE unique (158)", "LT unique (133)"))
  )

cat("Significant enriched terms loaded:\n")
print(table(sig_results$GeneSet, sig_results$RefSpecies, sig_results$DB_type))

# =============================================================================
# SECTION 3: BUBBLE PLOT — KEGG pathways, mmu reference
# =============================================================================

plot_bubble_kegg <- function(data, ref_sp, top_n = 20,
                              title_suffix = "") {
  df <- data %>%
    filter(DB_type == "KEGG", RefSpecies == ref_sp) %>%
    group_by(Term) %>%
    slice_min(order_by = Corrected_PValue, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(Corrected_PValue) %>%
    slice_head(n = top_n)

  if (nrow(df) == 0) {
    message("No significant KEGG terms for ref: ", ref_sp)
    return(NULL)
  }

  # Re-attach all gene sets for those top terms
  df_all <- data %>%
    filter(DB_type == "KEGG", RefSpecies == ref_sp,
           Term %in% df$Term,
           !is.na(Corrected_PValue), Corrected_PValue <= 0.05) %>%
    mutate(Term = factor(Term, levels = rev(unique(df$Term))))

  ggplot(df_all,
         aes(x = GeneSet, y = Term,
             size  = Input_n,
             color = neg_log10_padj)) +
    geom_point(alpha = 0.85) +
    scale_color_gradient(low = "#56B1F7", high = "#132B43",
                         name = expression(-log[10](FDR))) +
    scale_size_continuous(name = "Gene count", range = c(3, 12)) +
    labs(
      title    = paste0("KEGG Pathway Enrichment — ", ref_sp, title_suffix),
      subtitle = "Lepus europaeus vs Lepus timidus (skin, adult)",
      x = "Gene set", y = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.y    = element_text(size = 9),
      legend.position = "right",
      plot.title     = element_text(face = "bold")
    )
}

p_kegg_mmu <- plot_bubble_kegg(sig_results, "mmu", top_n = 20)
p_kegg_ocu <- plot_bubble_kegg(sig_results, "ocu", top_n = 20)

if (!is.null(p_kegg_mmu)) {
  ggsave("KEGG_bubble_mmu.pdf", p_kegg_mmu, width = 14, height = 9, dpi = 300)
  ggsave("KEGG_bubble_mmu.png", p_kegg_mmu, width = 14, height = 9, dpi = 300)
  cat("Saved: KEGG_bubble_mmu.pdf\n")
}
if (!is.null(p_kegg_ocu)) {
  ggsave("KEGG_bubble_ocu.pdf", p_kegg_ocu, width = 14, height = 9, dpi = 300)
  cat("Saved: KEGG_bubble_ocu.pdf\n")
}

# =============================================================================
# SECTION 4: BAR CHART — top enriched KEGG terms per gene set
# =============================================================================

plot_bar_kegg <- function(data, ref_sp, top_n = 15) {
  df <- data %>%
    filter(DB_type == "KEGG", RefSpecies == ref_sp) %>%
    group_by(GeneSet) %>%
    slice_min(order_by = Corrected_PValue, n = top_n, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(GeneSet, Corrected_PValue) %>%
    mutate(Term = str_trunc(Term, 50))

  if (nrow(df) == 0) return(NULL)

  ggplot(df, aes(x = reorder(Term, neg_log10_padj),
                 y = neg_log10_padj,
                 fill = GeneSet)) +
    geom_bar(stat = "identity", position = "dodge") +
    coord_flip() +
    scale_fill_brewer(palette = "Set2") +
    labs(
      title = paste0("Top KEGG Enrichment — ", ref_sp),
      x = NULL,
      y = expression(-log[10](FDR)),
      fill = "Gene set"
    ) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))
}

p_bar_mmu <- plot_bar_kegg(sig_results, "mmu", top_n = 15)
p_bar_ocu <- plot_bar_kegg(sig_results, "ocu", top_n = 15)

if (!is.null(p_bar_mmu)) {
  ggsave("KEGG_barplot_mmu.pdf", p_bar_mmu, width = 14, height = 9, dpi = 300)
  cat("Saved: KEGG_barplot_mmu.pdf\n")
}
if (!is.null(p_bar_ocu)) {
  ggsave("KEGG_barplot_ocu.pdf", p_bar_ocu, width = 14, height = 9, dpi = 300)
  cat("Saved: KEGG_barplot_ocu.pdf\n")
}

# =============================================================================
# SECTION 5: GO ENRICHMENT BAR CHART (top GO terms per gene set)
# =============================================================================

plot_bar_go <- function(data, ref_sp, top_n = 15) {
  df <- data %>%
    filter(DB_type == "GO", RefSpecies == ref_sp) %>%
    group_by(GeneSet) %>%
    slice_min(order_by = Corrected_PValue, n = top_n, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(Term = str_trunc(Term, 55))

  if (nrow(df) == 0) return(NULL)

  ggplot(df, aes(x = reorder(Term, neg_log10_padj),
                 y = neg_log10_padj,
                 fill = GeneSet)) +
    geom_bar(stat = "identity", position = "dodge") +
    coord_flip() +
    scale_fill_brewer(palette = "Set1") +
    labs(
      title = paste0("Top GO Enrichment — ", ref_sp),
      x = NULL,
      y = expression(-log[10](FDR)),
      fill = "Gene set"
    ) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))
}

p_go_mmu <- plot_bar_go(sig_results, "mmu", top_n = 15)
p_go_ocu <- plot_bar_go(sig_results, "ocu", top_n = 15)

if (!is.null(p_go_mmu)) {
  ggsave("GO_barplot_mmu.pdf", p_go_mmu, width = 14, height = 9, dpi = 300)
  cat("Saved: GO_barplot_mmu.pdf\n")
}
if (!is.null(p_go_ocu)) {
  ggsave("GO_barplot_ocu.pdf", p_go_ocu, width = 14, height = 9, dpi = 300)
  cat("Saved: GO_barplot_ocu.pdf\n")
}

# =============================================================================
# SECTION 6: HEATMAP — enrichment matrix
# Shows which pathways are enriched in which gene set (mmu KEGG shown)
# =============================================================================

make_enrichment_heatmap <- function(data, ref_sp, db_type = "KEGG",
                                    top_n = 30) {
  df <- data %>%
    filter(DB_type == db_type, RefSpecies == ref_sp) %>%
    group_by(Term) %>%
    slice_min(Corrected_PValue, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(Corrected_PValue) %>%
    slice_head(n = top_n)

  if (nrow(df) < 2) return(NULL)

  top_terms <- df$Term

  mat <- data %>%
    filter(DB_type == db_type, RefSpecies == ref_sp,
           Term %in% top_terms) %>%
    select(Term, GeneSet, neg_log10_padj) %>%
    pivot_wider(names_from = GeneSet, values_from = neg_log10_padj,
                values_fill = 0) %>%
    as.data.frame()

  rownames(mat) <- str_trunc(mat$Term, 55)
  mat$Term <- NULL
  mat <- as.matrix(mat)

  pheatmap(
    mat,
    color           = colorRampPalette(c("white", "#56B1F7", "#132B43"))(100),
    cluster_rows    = TRUE,
    cluster_cols    = FALSE,
    fontsize_row    = 8,
    fontsize_col    = 11,
    main            = paste0(db_type, " Enrichment Heatmap (", ref_sp, ")"),
    na_col          = "grey95",
    filename        = paste0(db_type, "_heatmap_", ref_sp, ".pdf"),
    width           = 12,
    height          = max(6, round(nrow(mat) * 0.28))
  )
  cat("Saved:", paste0(db_type, "_heatmap_", ref_sp, ".pdf\n"))
}

make_enrichment_heatmap(sig_results, "mmu", "KEGG", top_n = 30)
make_enrichment_heatmap(sig_results, "ocu", "KEGG", top_n = 30)
make_enrichment_heatmap(sig_results, "mmu", "GO",   top_n = 30)
make_enrichment_heatmap(sig_results, "ocu", "GO",   top_n = 30)

# =============================================================================
# SECTION 7: COMPARE mmu vs ocu (consistency check)
# Scatter: -log10(FDR) in mmu vs ocu for KEGG pathways (shared gene set)
# =============================================================================

compare_species <- function(data, gene_set_label, db_type = "KEGG") {
  df <- data %>%
    filter(DB_type == db_type, GeneSet == gene_set_label) %>%
    select(Term, RefSpecies, neg_log10_padj) %>%
    pivot_wider(names_from = RefSpecies,
                values_from = neg_log10_padj,
                values_fill = 0)

  if (!all(c("mmu","ocu") %in% colnames(df))) return(NULL)
  if (nrow(df) < 3) return(NULL)

  ggplot(df, aes(x = mmu, y = ocu, label = str_trunc(Term, 40))) +
    geom_point(color = "#2171B5", size = 2.5, alpha = 0.7) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    geom_text_repel(size = 2.8, max.overlaps = 12) +
    labs(
      title    = paste0(db_type, " enrichment: mmu vs ocu (", gene_set_label, ")"),
      x        = expression(-log[10](FDR)~"[mouse reference]"),
      y        = expression(-log[10](FDR)~"[rabbit reference]")
    ) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))
}

p_comp_shared <- compare_species(sig_results, "Shared (508)", "KEGG")
p_comp_le     <- compare_species(sig_results, "LE unique (158)", "KEGG")
p_comp_lt     <- compare_species(sig_results, "LT unique (133)", "KEGG")

for (nm in c("Shared_508", "LE_unique_158", "LT_unique_133")) {
  label_for_plot <- switch(nm,
    "Shared_508"    = "Shared (508)",
    "LE_unique_158" = "LE unique (158)",
    "LT_unique_133" = "LT unique (133)"
  )
  p <- compare_species(sig_results, label_for_plot, "KEGG")
  if (!is.null(p)) {
    fname <- paste0("KEGG_mmu_vs_ocu_", nm, ".pdf")
    ggsave(fname, p, width = 10, height = 8, dpi = 300)
    cat("Saved:", fname, "\n")
  }
}

# =============================================================================
# SECTION 8: UNIQUE PATHWAYS PER GENE SET
# Which pathways are enriched ONLY in shared genes, ONLY in LE, ONLY in LT?
# =============================================================================

find_unique_pathways <- function(data, ref_sp, db_type = "KEGG") {
  df <- data %>%
    filter(RefSpecies == ref_sp, DB_type == db_type) %>%
    select(Term, GeneSet)

  if (nrow(df) == 0) return(NULL)

  # Which gene sets contain each term?
  term_sets <- df %>%
    group_by(Term) %>%
    summarise(gene_sets = paste(sort(unique(as.character(GeneSet))),
                                collapse = " | ")) %>%
    ungroup()

  # Exclusively enriched in one gene set
  exclusive <- term_sets %>%
    filter(!grepl("\\|", gene_sets)) %>%
    arrange(gene_sets, Term)

  # Enriched in all three
  shared_all <- term_sets %>%
    filter(str_count(gene_sets, "\\|") == 2)

  list(exclusive = exclusive, shared_all = shared_all, term_sets = term_sets)
}

for (ref in c("mmu","ocu")) {
  res <- find_unique_pathways(sig_results, ref, "KEGG")
  if (!is.null(res)) {
    write.csv(res$exclusive,
              paste0("KEGG_", ref, "_exclusive_enrichment.csv"),
              row.names = FALSE)
    write.csv(res$shared_all,
              paste0("KEGG_", ref, "_enriched_in_all3_genesets.csv"),
              row.names = FALSE)
    cat("Saved exclusive/shared pathway tables for ref:", ref, "\n")
  }
}

# =============================================================================
# SECTION 9: UPSET-STYLE SUMMARY TABLE
# Count of significant terms per gene set × database × species
# =============================================================================

summary_counts <- sig_results %>%
  group_by(GeneSet, RefSpecies, DB_type) %>%
  summarise(n_significant = n(), .groups = "drop") %>%
  pivot_wider(names_from = c(RefSpecies, DB_type),
              values_from = n_significant,
              values_fill = 0)

print(summary_counts)
write.csv(summary_counts, "enrichment_summary_counts.csv", row.names = FALSE)
cat("Saved: enrichment_summary_counts.csv\n")

# =============================================================================
# SECTION 10: EXPORT FULL SIGNIFICANT RESULTS TABLE
# =============================================================================

write.csv(sig_results,
          "all_significant_enrichment_FDR05.csv",
          row.names = FALSE)
cat("Saved: all_significant_enrichment_FDR05.csv\n")

cat("\n===== Analysis complete =====\n")
cat("Output files in working directory:\n")
cat("  Bubble plots:  KEGG_bubble_mmu.pdf, KEGG_bubble_ocu.pdf\n")
cat("  Bar charts:    KEGG_barplot_mmu.pdf, KEGG_barplot_ocu.pdf\n")
cat("                 GO_barplot_mmu.pdf, GO_barplot_ocu.pdf\n")
cat("  Heatmaps:      KEGG_heatmap_mmu.pdf, KEGG_heatmap_ocu.pdf\n")
cat("                 GO_heatmap_mmu.pdf, GO_heatmap_ocu.pdf\n")
cat("  mmu vs ocu:    KEGG_mmu_vs_ocu_*.pdf\n")
cat("  Tables:        KEGG_*_exclusive_enrichment.csv\n")
cat("                 KEGG_*_enriched_in_all3_genesets.csv\n")
cat("                 enrichment_summary_counts.csv\n")
cat("                 all_significant_enrichment_FDR05.csv\n")
