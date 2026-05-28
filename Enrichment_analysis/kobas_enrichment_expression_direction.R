#!/usr/bin/env Rscript
# =============================================================================
# Expression direction for ALL genes in significant mmu KOBAS enrichments
# (same gene lists as visualize_mmu_kobas_results.R — not limited to seasonal module)
#
# Joins KOBAS pathway membership with DESeq2:
#   - Shared_508 + LE_unique_158: LT_vs_LE (europaeus ref, NC IDs)
#   - LT_unique_133: LE_vs_LT (timidus ref, CM IDs; LFC flipped to LT-vs-LE)
#
# Outputs: mmu_analysis/tables/
#   kobas_all_pathway_genes_with_direction.csv
#   kobas_pathway_direction_summary_all.csv
#   kobas_gene_direction_unique.csv
#   kobas_regulation_by_geneset.csv
#   kobas_regulation_LT_vs_LE_by_geneset.csv
# Figures: mmu_analysis/figures/
#   regulation_summary_by_geneset.pdf
#   top_pathways_regulation_KEGG.pdf (optional stacked bars)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(readr)
})

base_dir   <- path.expand("~/Documents/Enrichment_analysis")
table_dir  <- file.path(base_dir, "mmu_analysis", "tables")
fig_dir    <- file.path(base_dir, "mmu_analysis", "figures")
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

fdr_thresh   <- 0.05
padj_thresh  <- 0.05
logfc_thresh <- 1.0
top_n_path   <- 12

deseq_file_le_ref <- "~/Documents/DESeq2/DESEq2_analysis/1.Lepus_europaeus_as_reference/LT_vs_LE_adj.csv"
deseq_file_lt_ref <- c(
  "~/Documents/DESeq2/DESEq2_analysis/2.Lepus_timidus_as_reference/LE_vs_LT_adj.csv",
  "~/Documents/DESeq2/DESEq2_analysis/2.Lepus timidus_as_reference/LE_vs_LT_adj.csv"
)

kobas_tsv_files <- list(
  Shared_508    = file.path(base_dir, "KOBAS_output_mmu_KEGG_1", "KOBAS_output_mmu_KEGG_1.tsv"),
  LE_unique_158 = file.path(base_dir, "KOBAS_output_mmu_KEGG_2", "KOBAS_output_mmu_KEGG_2.tsv"),
  LT_unique_133 = file.path(base_dir, "KOBAS_output_mmu_KEGG_3", "KOBAS_output_mmu_KEGG_3.tsv")
)

identify_files <- list(
  list(file.path(base_dir, "KOBAS_output_mmu_KEGG_1", "KOBAS_output_mmu_KEGG_1.tsv_identify.txt"), "Shared_508",    "KEGG"),
  list(file.path(base_dir, "KOBAS_output_mmu_KEGG_2", "KOBAS_output_mmu_KEGG_2.tsv_identify.txt"), "LE_unique_158", "KEGG"),
  list(file.path(base_dir, "KOBAS_output_mmu_KEGG_3", "KOBAS_output_mmu_KEGG_3.tsv_identify.txt"), "LT_unique_133", "KEGG"),
  list(file.path(base_dir, "KOBAS_output_mmu_GO_1",   "KOBAS_output_mmu_GO_1.tsv_identify.txt"),   "Shared_508",    "GO"),
  list(file.path(base_dir, "KOBAS_output_mmu_GO_2",   "KOBAS_output_mmu_GO_2.tsv_identify.txt"),   "LE_unique_158", "GO"),
  list(file.path(base_dir, "KOBAS_output_mmu_GO_3",   "KOBAS_output_mmu_GO_3.tsv_identify.txt"),   "LT_unique_133", "GO")
)

gene_set_labels <- c(
  Shared_508    = "508 shared",
  LE_unique_158 = "158 LE-unique",
  LT_unique_133 = "133 LT-unique"
)

# ---- Shared helpers (aligned with rank_candidates_and_module_tables.R) ----

resolve_deseq_file <- function(paths) {
  for (p in paths) {
    p <- path.expand(p)
    if (file.exists(p)) return(p)
  }
  stop("DESeq2 file not found: ", paste(paths, collapse = ", "), call. = FALSE)
}

normalize_gene_id <- function(x) sub("\\.t[0-9]+$", "", as.character(x))

clean_gene_symbol <- function(s) {
  s <- as.character(s)
  bad <- is.na(s) | s == "" | grepl("http|www|\\.|/", s, ignore.case = TRUE) |
    nchar(s) > 25
  ifelse(bad, NA_character_, s)
}

load_deseq <- function(filepath, deseq_source, flip_lfc = FALSE) {
  df <- read.csv(filepath, stringsAsFactors = FALSE, check.names = FALSE)
  if (ncol(df) > 0 && (colnames(df)[1] == "" || is.na(colnames(df)[1]))) {
    colnames(df)[1] <- "gene_id"
  }
  names(df) <- tolower(names(df))
  if (is.na(names(df)["gene_id"]) || names(df)["gene_id"] == "") names(df)[1] <- "gene_id"
  lfc_col  <- grep("log2fold", names(df), value = TRUE)[1]
  padj_col <- grep("^padj$|^fdr$", names(df), value = TRUE)[1]
  if (is.na(lfc_col) || is.na(padj_col)) {
    stop("Cannot find log2FoldChange / padj in ", filepath)
  }
  out <- df %>%
    transmute(
      gene_id = normalize_gene_id(.data$gene_id),
      baseMean = if ("basemean" %in% names(df)) .data$basemean else NA_real_,
      log2FC_raw = .data[[lfc_col]],
      padj = .data[[padj_col]],
      deseq_source = deseq_source
    )
  out$log2FC <- if (flip_lfc) -out$log2FC_raw else out$log2FC_raw
  out %>%
    mutate(
      regulation = case_when(
        padj < padj_thresh & log2FC >=  logfc_thresh ~ "Up in LT",
        padj < padj_thresh & log2FC <= -logfc_thresh ~ "Up in LE",
        TRUE ~ "Not significant"
      ),
      # LT vs LE contrast: positive log2FC = upregulated in LT
      regulation_LT_vs_LE = case_when(
        regulation == "Up in LT" ~ "Upregulated",
        regulation == "Up in LE" ~ "Downregulated",
        TRUE ~ "Not significant"
      ),
      is_DE = regulation != "Not significant"
    )
}

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
      GeneSet_label = gene_set_labels[gene_set],
      Input_n = as.numeric(Input_n),
      Corrected_PValue = as.numeric(Corrected_PValue)
    )
}

parse_kobas_blast <- function(filepath) {
  if (!file.exists(filepath)) return(NULL)
  lines <- readLines(filepath)
  data_lines <- lines[!grepl("^#", lines) & nchar(trimws(lines)) > 0]
  df <- read.table(
    text = paste(data_lines, collapse = "\n"),
    sep = "\t", header = FALSE, quote = "", fill = TRUE,
    stringsAsFactors = FALSE
  )
  colnames(df)[1:2] <- c("query_gene", "ref_info")
  df %>%
    filter(grepl("^(NC_|CM)", query_gene)) %>%
    mutate(
      ref_gene_id = str_extract(ref_info, "^[^|]+"),
      gene_symbol = clean_gene_symbol(str_extract(ref_info, "(?<=\\|)[^|,]+"))
    ) %>%
    group_by(query_gene) %>%
    slice(1) %>%
    ungroup() %>%
    transmute(
      gene_id = as.character(query_gene),
      ref_gene_id = as.character(ref_gene_id),
      gene_symbol = as.character(gene_symbol)
    )
}

# =============================================================================
# Load data
# =============================================================================

cat("=== KOBAS enrichment genes + expression direction ===\n")

deseq_le_path <- path.expand(deseq_file_le_ref)
deseq_lt_path <- resolve_deseq_file(deseq_file_lt_ref)
cat("DESeq2 (NC):", deseq_le_path, "\n")
cat("DESeq2 (CM):", deseq_lt_path, "\n\n")

deseq <- bind_rows(
  load_deseq(deseq_le_path, "LT_vs_LE (europaeus ref)", flip_lfc = FALSE),
  load_deseq(deseq_lt_path, "LE_vs_LT (timidus ref; LFC flipped)", flip_lfc = TRUE)
)

annotations <- bind_rows(lapply(kobas_tsv_files, parse_kobas_blast))

all_terms <- bind_rows(lapply(identify_files, function(x) {
  parse_kobas_identify(x[[1]], x[[2]], x[[3]])
}))

sig <- all_terms %>%
  filter(!is.na(Corrected_PValue), Corrected_PValue <= fdr_thresh)

pathway_genes <- sig %>%
  mutate(gene_list = str_split(Input_genes, "\\|")) %>%
  unnest(gene_list) %>%
  mutate(gene_id = trimws(gene_list)) %>%
  filter(nchar(gene_id) > 0, grepl("^(NC_|CM)", gene_id)) %>%
  left_join(deseq, by = "gene_id") %>%
  left_join(annotations, by = "gene_id") %>%
  select(
    GeneSet, GeneSet_label, DB_type, Term, ID, Input_n, Corrected_PValue,
    gene_id, gene_symbol, ref_gene_id,
    deseq_source, log2FC, log2FC_raw, padj, baseMean,
    regulation, regulation_LT_vs_LE, is_DE
  )

write.csv(
  pathway_genes,
  file.path(table_dir, "kobas_all_pathway_genes_with_direction.csv"),
  row.names = FALSE
)
cat("Saved: kobas_all_pathway_genes_with_direction.csv (",
    nrow(pathway_genes), " gene-pathway rows)\n", sep = "")

# Per-pathway direction summary (all terms from visualization script)
pathway_summary <- pathway_genes %>%
  group_by(GeneSet, GeneSet_label, DB_type, Term, ID, Input_n, Corrected_PValue) %>%
  summarise(
    n_genes = n(),
    n_up_LT = sum(regulation == "Up in LT", na.rm = TRUE),
    n_up_LE = sum(regulation == "Up in LE", na.rm = TRUE),
    n_upregulated = sum(regulation_LT_vs_LE == "Upregulated", na.rm = TRUE),
    n_downregulated = sum(regulation_LT_vs_LE == "Downregulated", na.rm = TRUE),
    n_NS = sum(regulation_LT_vs_LE == "Not significant" | is.na(regulation_LT_vs_LE), na.rm = TRUE),
    pct_up_LT = round(100 * n_up_LT / n_genes, 1),
    pct_up_LE = round(100 * n_up_LE / n_genes, 1),
    pct_upregulated = round(100 * n_upregulated / n_genes, 1),
    pct_downregulated = round(100 * n_downregulated / n_genes, 1),
    mean_log2FC = round(mean(log2FC, na.rm = TRUE), 3),
    dominant_direction = case_when(
      n_up_LT > n_up_LE & n_up_LT >= n_NS ~ "Mostly up in LT",
      n_up_LE > n_up_LT & n_up_LE >= n_NS ~ "Mostly up in LE",
      TRUE ~ "Mixed / not significant"
    ),
    dominant_regulation_LT_vs_LE = case_when(
      n_upregulated > n_downregulated & n_upregulated >= n_NS ~ "Mostly upregulated",
      n_downregulated > n_upregulated & n_downregulated >= n_NS ~ "Mostly downregulated",
      TRUE ~ "Mixed / not significant"
    ),
    example_upregulated = paste(
      head(gene_symbol[regulation_LT_vs_LE == "Upregulated" & !is.na(gene_symbol)], 3),
      collapse = ", "
    ),
    example_downregulated = paste(
      head(gene_symbol[regulation_LT_vs_LE == "Downregulated" & !is.na(gene_symbol)], 3),
      collapse = ", "
    ),
    example_up_LT = paste(
      head(gene_symbol[regulation == "Up in LT" & !is.na(gene_symbol)], 3),
      collapse = ", "
    ),
    example_up_LE = paste(
      head(gene_symbol[regulation == "Up in LE" & !is.na(gene_symbol)], 3),
      collapse = ", "
    ),
    .groups = "drop"
  ) %>%
  arrange(DB_type, GeneSet, Corrected_PValue)

write.csv(
  pathway_summary,
  file.path(table_dir, "kobas_pathway_direction_summary_all.csv"),
  row.names = FALSE
)
cat("Saved: kobas_pathway_direction_summary_all.csv\n")

# One row per unique enriched gene (may appear in multiple pathways)
gene_unique <- pathway_genes %>%
  group_by(gene_id, GeneSet, GeneSet_label) %>%
  summarise(
    gene_symbol = first(na.omit(gene_symbol)),
    n_pathways = n_distinct(Term),
    pathways = paste(sort(unique(Term)), collapse = "; "),
    log2FC = first(log2FC),
    padj = first(padj),
    regulation = first(regulation),
    regulation_LT_vs_LE = first(regulation_LT_vs_LE),
    is_DE = first(is_DE),
    DB_types = paste(sort(unique(DB_type)), collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(GeneSet, desc(is_DE), desc(abs(log2FC)))

write.csv(
  gene_unique,
  file.path(table_dir, "kobas_gene_direction_unique.csv"),
  row.names = FALSE
)
cat("Saved: kobas_gene_direction_unique.csv (",
    nrow(gene_unique), " unique genes)\n", sep = "")

# Counts per gene set (species labels + LT vs LE labels)
by_set <- pathway_genes %>%
  distinct(gene_id, GeneSet, GeneSet_label, regulation, regulation_LT_vs_LE, is_DE) %>%
  group_by(GeneSet, GeneSet_label, regulation) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = regulation, values_from = n, values_fill = 0)

by_set_lt_vs_le <- pathway_genes %>%
  distinct(gene_id, GeneSet, GeneSet_label, regulation_LT_vs_LE) %>%
  group_by(GeneSet, GeneSet_label, regulation_LT_vs_LE) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = regulation_LT_vs_LE, values_from = n, values_fill = 0)

write.csv(by_set, file.path(table_dir, "kobas_regulation_by_geneset.csv"), row.names = FALSE)
write.csv(by_set_lt_vs_le,
          file.path(table_dir, "kobas_regulation_LT_vs_LE_by_geneset.csv"),
          row.names = FALSE)
cat("Saved: kobas_regulation_by_geneset.csv\n")
cat("Saved: kobas_regulation_LT_vs_LE_by_geneset.csv\n")

# =============================================================================
# Figures
# =============================================================================

reg_cols <- setdiff(names(by_set_lt_vs_le), c("GeneSet", "GeneSet_label"))
p_set <- by_set_lt_vs_le %>%
  pivot_longer(all_of(reg_cols), names_to = "regulation_LT_vs_LE", values_to = "n") %>%
  mutate(
    regulation_LT_vs_LE = factor(
      regulation_LT_vs_LE,
      levels = c("Upregulated", "Downregulated", "Not significant")
    )
  )

p1 <- ggplot(p_set, aes(x = GeneSet_label, y = n, fill = regulation_LT_vs_LE)) +
  geom_col(position = "stack") +
  scale_fill_manual(
    values = c(
      "Upregulated" = "#D62728",
      "Downregulated" = "#1F77B4",
      "Not significant" = "grey75"
    ),
    name = "LT vs LE"
  ) +
  labs(
    title = "Regulation of genes in significant KOBAS enrichments (mmu)",
    subtitle = sprintf(
      "Upregulated = higher in LT; Downregulated = higher in LE | |log2FC| >= %.1f, padj < %.2f",
      logfc_thresh, padj_thresh
    ),
    x = NULL, y = "Unique genes"
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), plot.title = element_text(face = "bold"))

ggsave(file.path(fig_dir, "regulation_summary_by_geneset.pdf"), p1,
       width = 9, height = 5.5, dpi = 300)
ggsave(file.path(fig_dir, "regulation_summary_by_geneset.png"), p1,
       width = 9, height = 5.5, dpi = 300)
cat("Saved: regulation_summary_by_geneset.pdf / .png\n")

# Top KEGG pathways: stacked regulation
top_kegg_terms <- pathway_summary %>%
  filter(DB_type == "KEGG") %>%
  group_by(GeneSet_label) %>%
  slice_min(Corrected_PValue, n = top_n_path, with_ties = FALSE) %>%
  ungroup()

if (nrow(top_kegg_terms) > 0) {
  plot_df <- pathway_genes %>%
    filter(DB_type == "KEGG", Term %in% top_kegg_terms$Term,
           GeneSet_label %in% top_kegg_terms$GeneSet_label) %>%
    mutate(
      Term_short = str_trunc(Term, 42),
      regulation_LT_vs_LE = factor(
        regulation_LT_vs_LE,
        levels = c("Upregulated", "Downregulated", "Not significant")
      )
    )

  p2 <- ggplot(plot_df, aes(x = regulation_LT_vs_LE, fill = regulation_LT_vs_LE)) +
    geom_bar() +
    facet_grid(GeneSet_label ~ Term_short, scales = "free_y", space = "free_y") +
    scale_fill_manual(
      values = c(
        "Upregulated" = "#D62728",
        "Downregulated" = "#1F77B4",
        "Not significant" = "grey75"
      ),
      guide = "none"
    ) +
    labs(
      title = "Expression direction in top KEGG pathways (all enriched genes)",
      x = "LT vs LE",
      y = "Gene count"
    ) +
    theme_bw(base_size = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.text.y = element_text(size = 8))

  ggsave(file.path(fig_dir, "top_pathways_regulation_KEGG.pdf"), p2,
         width = 12, height = max(8, nrow(top_kegg_terms) * 0.35), dpi = 300,
         limitsize = FALSE)
  cat("Saved: top_pathways_regulation_KEGG.pdf\n")
}

# =============================================================================
# Console summary
# =============================================================================

cat("\n--- Regulation counts (unique genes per gene set) ---\n")
cat("By species (regulation):\n")
print(
  pathway_genes %>%
    distinct(gene_id, GeneSet_label, regulation) %>%
    count(GeneSet_label, regulation) %>%
    pivot_wider(names_from = regulation, values_from = n, values_fill = 0)
)
cat("\nBy LT vs LE (regulation_LT_vs_LE):\n")
print(
  pathway_genes %>%
    distinct(gene_id, GeneSet_label, regulation_LT_vs_LE) %>%
    count(GeneSet_label, regulation_LT_vs_LE) %>%
    pivot_wider(names_from = regulation_LT_vs_LE, values_from = n, values_fill = 0)
)

cat("\nInterpretation:\n")
cat("  regulation column:\n")
cat("    Up in LT  = higher in Lepus timidus (mountain hare)\n")
cat("    Up in LE  = higher in Lepus europaeus (brown hare)\n")
cat("  regulation_LT_vs_LE column (LT vs LE contrast):\n")
cat("    Upregulated   = same as Up in LT (log2FC > 0)\n")
cat("    Downregulated = same as Up in LE (log2FC < 0)\n")
cat("  Positive log2FC always means upregulated in LT (LE_vs_LT file is sign-flipped).\n")
cat("\nDone. Tables:", table_dir, "| Figures:", fig_dir, "\n")
