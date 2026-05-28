#!/usr/bin/env Rscript
# =============================================================================
# Mus musculus KOBAS enrichment + DESeq2 up/down regulation
#
# Integrates:
#   - LT vs LE differential expression (Lepus europaeus reference by default)
#   - Three gene sets: Shared_508, LE_unique_158, LT_unique_133
#   - mmu KOBAS GO and KEGG enrichment (_identify.txt, FDR <= 0.05)
#
# KOBAS enrichment dot/bubble plots: run visualize_mmu_kobas_results.R
#
# Outputs (in mmu_analysis/):
#   mmu_gene_expression_classification.csv
#   mmu_pathway_gene_with_direction.csv
#   mmu_pathway_direction_summary.csv
#   mmu_gene_set_direction_counts.csv
#   volcano_mmu_gene_sets.pdf
#   pathway_direction_top_KEGG.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(readr)
})

# ---- Configuration (edit paths/thresholds if needed) ----
base_dir <- path.expand("~/Documents/Enrichment_analysis")
out_dir  <- file.path(base_dir, "mmu_analysis")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

deseq_file <- path.expand(
  "~/Documents/DESeq2/DESEq2_analysis/1.Lepus_europaeus_as_reference/LT_vs_LE_adj.csv"
)
# Alternative: mountain-hare reference
# deseq_file <- path.expand("~/Documents/DESeq2/DESEq2_analysis/2.Lepus timidus_as_reference/LT_vs_LE_adj.csv")

padj_thresh  <- 0.05
logfc_thresh <- 1.0
fdr_pathway  <- 0.05

gene_sets <- list(
  Shared_508    = file.path(base_dir, "1.int_ref_Genes508_adj.fa"),
  LE_unique_158 = file.path(base_dir, "2.int_ref_Genes158_adj.fa"),
  LT_unique_133 = file.path(base_dir, "3.int_ref_Genes133_adj.fa")
)

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

# =============================================================================
# Helpers
# =============================================================================

read_fasta_ids <- function(fasta_path) {
  if (!file.exists(fasta_path)) stop("FASTA not found: ", fasta_path)
  lines <- readLines(fasta_path)
  ids <- lines[grepl("^>", lines)]
  sub("^>", "", trimws(ids))
}

load_deseq <- function(filepath) {
  if (!file.exists(filepath)) stop("DESeq2 file not found: ", filepath)
  df <- read.csv(filepath, stringsAsFactors = FALSE, check.names = FALSE)
  # First column is often unnamed row names from write.csv(DESeq2 results)
  if (ncol(df) > 0 && (colnames(df)[1] == "" || is.na(colnames(df)[1]))) {
    colnames(df)[1] <- "gene_id"
  }
  if (!"gene_id" %in% colnames(df)) {
    if (ncol(df) > 0 && grepl("^X", colnames(df)[1])) {
      df$gene_id <- df[[1]]
    } else {
      df$gene_id <- rownames(df)
    }
  }
  names(df) <- tolower(names(df))
  if (is.na(names(df)["gene_id"]) || names(df)["gene_id"] == "") {
    names(df)[1] <- "gene_id"
  }
  lfc_col  <- grep("log2fold", names(df), value = TRUE)[1]
  padj_col <- grep("^padj$|^fdr$", names(df), value = TRUE)[1]
  if (is.na(lfc_col) || is.na(padj_col)) {
    stop("Cannot find log2FoldChange / padj columns in ", filepath)
  }
  df %>%
    transmute(
      gene_id = .data$gene_id,
      baseMean = if ("basemean" %in% names(df)) .data$basemean else NA_real_,
      log2FC   = .data[[lfc_col]],
      padj     = .data[[padj_col]]
    ) %>%
    mutate(
      regulation = case_when(
        padj < padj_thresh & log2FC >=  logfc_thresh ~ "Up in LT",
        padj < padj_thresh & log2FC <= -logfc_thresh ~ "Up in LE",
        TRUE ~ "Not significant"
      ),
      regulation_short = case_when(
        regulation == "Up in LT" ~ "Upregulated",
        regulation == "Up in LE" ~ "Downregulated",
        TRUE ~ "NS"
      )
    )
}

parse_identify <- function(filepath, gene_set_label, db_type) {
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
  n <- min(ncol(df), 9)
  cnames <- c("Term", "Database", "ID", "Input_n", "Background_n",
              "PValue", "FDR", "Genes", "Hyperlink")[seq_len(n)]
  colnames(df)[seq_len(n)] <- cnames
  df %>%
    mutate(
      GeneSet = gene_set_label,
      DB_type = db_type,
      RefSpecies = "mmu",
      FDR = as.numeric(FDR),
      Input_n = as.numeric(Input_n)
    )
}

cat("=== Mus musculus enrichment + expression analysis ===\n")
cat("DESeq2:", deseq_file, "\n")
cat("Thresholds: |log2FC| >=", logfc_thresh, ", padj <", padj_thresh, "\n")
cat("Output:", out_dir, "\n\n")

# =============================================================================
# 1. Gene sets + expression classification
# =============================================================================

deseq <- load_deseq(deseq_file)

gene_set_tbl <- bind_rows(lapply(names(gene_sets), function(gs) {
  tibble(gene_id = read_fasta_ids(gene_sets[[gs]]), GeneSet = gs)
}))

gene_expr <- gene_set_tbl %>%
  left_join(deseq, by = "gene_id") %>%
  mutate(
    in_deseq = !is.na(log2FC),
    GeneSet_label = recode(
      GeneSet,
      Shared_508 = "508 shared (both species)",
      LE_unique_158 = "158 LE-unique (brown hare)",
      LT_unique_133 = "133 LT-unique (mountain hare)"
    )
  )

write.csv(gene_expr,
          file.path(out_dir, "mmu_gene_expression_classification.csv"),
          row.names = FALSE)
cat("Saved: mmu_gene_expression_classification.csv\n")

cat("\nRegulation counts per gene set:\n")
gs_counts <- gene_expr %>%
  filter(in_deseq) %>%
  count(GeneSet_label, regulation) %>%
  arrange(GeneSet_label, regulation)
print(gs_counts)
write.csv(gs_counts,
          file.path(out_dir, "mmu_gene_set_direction_counts.csv"),
          row.names = FALSE)

# =============================================================================
# 2. Parse mmu KOBAS enrichment and attach direction
# =============================================================================

kobas <- bind_rows(lapply(identify_files, function(x) {
  parse_identify(x$path, x$gs, x$db)
}))

if (is.null(kobas) || nrow(kobas) == 0) {
  stop("No KOBAS identify data parsed. Check mmu KOBAS output folders.")
}

sig <- kobas %>% filter(!is.na(FDR), FDR <= fdr_pathway)

pathway_genes <- sig %>%
  filter("Genes" %in% names(.)) %>%
  mutate(gene_list = str_split(Genes, "\\|")) %>%
  unnest(gene_list) %>%
  mutate(gene_list = trimws(gene_list)) %>%
  filter(nchar(gene_list) > 0) %>%
  left_join(
    deseq %>% select(gene_id, log2FC, padj, regulation, regulation_short),
    by = c("gene_list" = "gene_id")
  )

write.csv(
  pathway_genes %>%
    select(GeneSet, DB_type, Term, ID, FDR, Input_n, gene_list,
           log2FC, padj, regulation, regulation_short),
  file.path(out_dir, "mmu_pathway_gene_with_direction.csv"),
  row.names = FALSE
)
cat("Saved: mmu_pathway_gene_with_direction.csv\n")

pathway_summary <- pathway_genes %>%
  group_by(GeneSet, DB_type, Term, ID, FDR, Input_n) %>%
  summarise(
    n_genes = n(),
    n_up_LT = sum(regulation == "Up in LT", na.rm = TRUE),
    n_up_LE = sum(regulation == "Up in LE", na.rm = TRUE),
    n_NS = sum(regulation == "Not significant" | is.na(regulation), na.rm = TRUE),
    mean_log2FC = mean(log2FC, na.rm = TRUE),
    dominant_direction = case_when(
      sum(regulation == "Up in LT", na.rm = TRUE) >
        sum(regulation == "Up in LE", na.rm = TRUE) &
        sum(regulation == "Up in LT", na.rm = TRUE) >=
        sum(regulation == "Not significant" | is.na(regulation), na.rm = TRUE) ~
        "Mostly up in LT",
      sum(regulation == "Up in LE", na.rm = TRUE) >
        sum(regulation == "Up in LT", na.rm = TRUE) &
        sum(regulation == "Up in LE", na.rm = TRUE) >=
        sum(regulation == "Not significant" | is.na(regulation), na.rm = TRUE) ~
        "Mostly up in LE",
      TRUE ~ "Mixed / NS"
    ),
    .groups = "drop"
  ) %>%
  arrange(DB_type, FDR)

write.csv(pathway_summary,
          file.path(out_dir, "mmu_pathway_direction_summary.csv"),
          row.names = FALSE)
cat("Saved: mmu_pathway_direction_summary.csv\n")

# =============================================================================
# 3. Figures
# =============================================================================

volcano_df <- deseq %>%
  inner_join(gene_set_tbl %>% distinct(gene_id, GeneSet), by = "gene_id") %>%
  mutate(
    GeneSet_label = recode(
      GeneSet,
      Shared_508 = "508 shared",
      LE_unique_158 = "158 LE-unique",
      LT_unique_133 = "133 LT-unique"
    ),
    neg_log10_padj = pmin(-log10(padj), 50),
    highlight = regulation != "Not significant"
  )

p_volcano <- ggplot(volcano_df, aes(x = log2FC, y = neg_log10_padj)) +
  geom_point(data = filter(volcano_df, !highlight),
             color = "grey80", size = 0.8, alpha = 0.5) +
  geom_point(data = filter(volcano_df, highlight),
             aes(color = regulation), size = 1.2, alpha = 0.85) +
  facet_wrap(~GeneSet_label, ncol = 1, scales = "free_y") +
  scale_color_manual(
    values = c("Up in LT" = "#D62728", "Up in LE" = "#1F77B4", "Not significant" = "grey50"),
    name = "Regulation (LT vs LE)"
  ) +
  geom_vline(xintercept = c(-logfc_thresh, logfc_thresh), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(padj_thresh), linetype = "dashed", color = "grey40") +
  labs(
    title = "Differential expression in KOBAS gene sets (mmu enrichment input)",
    subtitle = sprintf("LT vs LE | |log2FC| >= %.1f, padj < %.2f", logfc_thresh, padj_thresh),
    x = "log2 Fold Change (positive = higher in mountain hare LT)",
    y = expression(-log[10](padj))
  ) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

ggsave(file.path(out_dir, "volcano_mmu_gene_sets.pdf"), p_volcano,
       width = 9, height = 11, dpi = 300)
cat("Saved: volcano_mmu_gene_sets.pdf\n")

top_kegg <- pathway_summary %>%
  filter(DB_type == "KEGG") %>%
  group_by(GeneSet) %>%
  slice_min(FDR, n = 8, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(Term_short = str_trunc(Term, 45))

if (nrow(top_kegg) > 0) {
  kegg_long <- pathway_genes %>%
    filter(DB_type == "KEGG", Term %in% top_kegg$Term, GeneSet %in% top_kegg$GeneSet) %>%
    mutate(
      regulation_plot = factor(
        ifelse(is.na(regulation) | regulation == "Not significant", "Not significant", regulation),
        levels = c("Up in LT", "Up in LE", "Not significant")
      ),
      Term_short = str_trunc(Term, 45)
    )

  p_path <- ggplot(kegg_long, aes(x = regulation_plot, fill = regulation_plot)) +
    geom_bar(position = "stack") +
    facet_grid(GeneSet ~ Term_short, scales = "free", space = "free_y") +
    scale_fill_manual(
      values = c("Up in LT" = "#D62728", "Up in LE" = "#1F77B4", "Not significant" = "grey70"),
      name = NULL
    ) +
    labs(
      title = "Regulation of genes in top mmu KEGG pathways (FDR <= 0.05)",
      x = NULL, y = "Gene count"
    ) +
    theme_bw(base_size = 9) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text.y = element_text(size = 7)
    )

  ggsave(file.path(out_dir, "pathway_direction_top_KEGG.pdf"), p_path,
         width = 12, height = max(8, nrow(top_kegg) * 0.4), dpi = 300,
         limitsize = FALSE)
  cat("Saved: pathway_direction_top_KEGG.pdf\n")
}

# =============================================================================
# 4. Console summary
# =============================================================================

cat("\n--- Quick interpretation ---\n")
cat("Positive log2FC  = higher expression in Lepus timidus (mountain hare, LT)\n")
cat("Negative log2FC  = higher expression in Lepus europaeus (brown hare, LE)\n")
cat("LE_unique_158 genes are enriched for brown-hare-specific biology;\n")
cat("LT_unique_133 genes are enriched for mountain-hare-specific biology.\n")
cat("\nDone. All outputs in:", out_dir, "\n")
