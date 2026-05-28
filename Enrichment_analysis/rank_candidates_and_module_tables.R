#!/usr/bin/env Rscript
# =============================================================================
# High-impact tables for supervisor meetings (no extra enrichment plots)
#
# 1. ranked_candidate_genes.csv — all genes ranked by evidence score
# 2. top_candidates_by_gene_set.csv — top 15 per set for quick review
# 3. seasonal_module_genes_with_direction.csv — melanogenesis/circadian/hormone module
# 4. seasonal_module_summary.txt — one-page text summary for meetings
#
# Run after KOBAS + DESeq2. Outputs: mmu_analysis/tables/
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
})

# ---- Configuration ----
base_dir <- path.expand("~/Documents/Enrichment_analysis")
out_dir  <- file.path(base_dir, "mmu_analysis", "tables")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# L. europaeus reference (NC gene IDs) — Shared_508 + LE_unique_158
deseq_file_le_ref <- "~/Documents/DESeq2/DESEq2_analysis/1.Lepus_europaeus_as_reference/LT_vs_LE_adj.csv"
# L. timidus reference (CM gene IDs) — LT_unique_133
# Folder name may use space or underscore on disk
deseq_file_lt_ref <- c(
  "~/Documents/DESeq2/DESEq2_analysis/2.Lepus_timidus_as_reference/LE_vs_LT_adj.csv",
  "~/Documents/DESeq2/DESEq2_analysis/2.Lepus timidus_as_reference/LE_vs_LT_adj.csv"
)

padj_thresh  <- 0.05
logfc_thresh <- 1.0
fdr_pathway  <- 0.05
top_n_report <- 15

gene_sets <- list(
  Shared_508    = file.path(base_dir, "1.int_ref_Genes508_adj.fa"),
  LE_unique_158 = file.path(base_dir, "2.int_ref_Genes158_adj.fa"),
  LT_unique_133 = file.path(base_dir, "3.int_ref_Genes133_adj.fa")
)

gene_set_labels <- c(
  Shared_508    = "508 shared (both species)",
  LE_unique_158 = "158 LE-unique (brown hare)",
  LT_unique_133 = "133 LT-unique (mountain hare)"
)

kobas_tsv_files <- list(
  Shared_508    = file.path(base_dir, "KOBAS_output_mmu_KEGG_1", "KOBAS_output_mmu_KEGG_1.tsv"),
  LE_unique_158 = file.path(base_dir, "KOBAS_output_mmu_KEGG_2", "KOBAS_output_mmu_KEGG_2.tsv"),
  LT_unique_133 = file.path(base_dir, "KOBAS_output_mmu_KEGG_3", "KOBAS_output_mmu_KEGG_3.tsv")
)

identify_files <- list(
  list(file.path(base_dir, "KOBAS_output_mmu_KEGG_1", "KOBAS_output_mmu_KEGG_1.tsv_identify.txt"), "Shared_508"),
  list(file.path(base_dir, "KOBAS_output_mmu_KEGG_2", "KOBAS_output_mmu_KEGG_2.tsv_identify.txt"), "LE_unique_158"),
  list(file.path(base_dir, "KOBAS_output_mmu_KEGG_3", "KOBAS_output_mmu_KEGG_3.tsv_identify.txt"), "LT_unique_133")
)

# Coat colour / photoperiod / endocrine module (KEGG terms from your enrichment)
seasonal_pathways <- c(
  "Melanogenesis",
  "Circadian entrainment",
  "Thyroid hormone synthesis",
  "Cortisol synthesis and secretion",
  "Aldosterone synthesis and secretion",
  "Cushing syndrome"
)

generic_pathways <- c(
  "Pathways in cancer",
  "Metabolic pathways",
  "Human papillomavirus infection",
  "Proteoglycans in cancer"
)

clean_gene_symbol <- function(s) {
  s <- as.character(s)
  bad <- is.na(s) |
    s == "" |
    grepl("http|www|\\.|/", s, ignore.case = TRUE) |
    nchar(s) > 25 |
    grepl("^\\s", s)
  ifelse(bad, NA_character_, s)
}

# =============================================================================
# Helpers
# =============================================================================

read_fasta_ids <- function(path) {
  lines <- readLines(path)
  sub("^>", "", trimws(lines[grepl("^>", lines)]))
}

resolve_deseq_file <- function(paths, require_padj = TRUE) {
  for (p in paths) {
    p <- path.expand(p)
    if (!file.exists(p)) next
    hdr <- tolower(names(read.csv(p, nrows = 0, check.names = FALSE)))
    if (!any(grepl("log2fold", hdr))) next
    if (require_padj && !any(grepl("^padj", hdr))) next
    return(p)
  }
  stop(
    "DESeq2 file not found (or missing log2FoldChange/padj). Tried:\n",
    paste0("  - ", path.expand(paths), collapse = "\n"),
    call. = FALSE
  )
}

normalize_gene_id <- function(x) {
  x <- as.character(x)
  sub("\\.t[0-9]+$", "", x)
}

# Load DESeq2 and optionally flip LFC so positive = higher in LT (mountain hare)
# LE_vs_LT raw: positive LFC = higher in LE; flip for harmonized LT-vs-LE view
load_deseq <- function(filepath, deseq_source, flip_lfc = FALSE) {
  df <- read.csv(filepath, stringsAsFactors = FALSE, check.names = FALSE)
  if (ncol(df) > 0 && (colnames(df)[1] == "" || is.na(colnames(df)[1]))) {
    colnames(df)[1] <- "gene_id"
  }
  names(df) <- tolower(names(df))
  if (is.na(names(df)["gene_id"]) || names(df)["gene_id"] == "") {
    names(df)[1] <- "gene_id"
  }
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
      lfc_harmonized = flip_lfc,
      regulation = case_when(
        padj < padj_thresh & log2FC >=  logfc_thresh ~ "Up in LT",
        padj < padj_thresh & log2FC <= -logfc_thresh ~ "Up in LE",
        TRUE ~ "Not significant"
      ),
      is_DE = regulation != "Not significant"
    )
}

parse_kobas_blast <- function(filepath) {
  if (!file.exists(filepath)) return(NULL)
  lines <- readLines(filepath)
  data_lines <- lines[!grepl("^#", lines) & nchar(trimws(lines)) > 0]
  if (length(data_lines) == 0) return(NULL)

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
      gene_symbol = as.character(str_extract(ref_info, "(?<=\\|)[^|,]+")),
      gene_symbol = ifelse(is.na(gene_symbol) | gene_symbol == "", NA_character_, gene_symbol)
    ) %>%
    group_by(query_gene) %>%
    slice(1) %>%
    ungroup() %>%
    transmute(
      gene_id = as.character(query_gene),
      ref_gene_id = as.character(ref_gene_id),
      gene_symbol = clean_gene_symbol(gene_symbol)
    )
}

parse_identify_kegg <- function(filepath, gene_set) {
  if (!file.exists(filepath)) return(NULL)
  lines <- readLines(filepath)
  data_lines <- lines[!grepl("^#", lines) & nchar(trimws(lines)) > 0]
  df <- read.table(
    text = paste(data_lines, collapse = "\n"),
    sep = "\t", header = FALSE, quote = "", fill = TRUE,
    stringsAsFactors = FALSE
  )
  n <- min(ncol(df), 8)
  colnames(df)[seq_len(n)] <- c(
    "Term", "Database", "ID", "Input_n", "Background_n",
    "PValue", "FDR", "Genes"
  )[seq_len(n)]
  df %>%
    mutate(
      GeneSet = gene_set,
      FDR = as.numeric(FDR)
    ) %>%
    filter(!is.na(FDR), FDR <= fdr_pathway) %>%
    mutate(gene_list = str_split(Genes, "\\|")) %>%
    unnest(gene_list) %>%
    mutate(gene_list = trimws(gene_list)) %>%
    filter(nchar(gene_list) > 0, grepl("^(NC_|CM)", gene_list)) %>%
    select(gene_id = gene_list, GeneSet, Term, ID, FDR)
}

cat("=== Candidate ranking + seasonal module tables ===\n")
cat("Output:", out_dir, "\n\n")

# ---- Load core data ----
deseq_le_path <- path.expand(deseq_file_le_ref)
if (!file.exists(deseq_le_path)) {
  stop("DESeq2 file not found: ", deseq_le_path, call. = FALSE)
}
deseq_lt_path <- resolve_deseq_file(deseq_file_lt_ref)

cat("DESeq2 (NC IDs, Shared + LE-unique):", deseq_le_path, "\n")
cat("DESeq2 (CM IDs, LT-unique):", deseq_lt_path, "\n\n")

deseq_le <- load_deseq(
  deseq_le_path,
  deseq_source = "LT_vs_LE (L. europaeus reference, NC IDs)",
  flip_lfc = FALSE
)
deseq_lt <- load_deseq(
  deseq_lt_path,
  deseq_source = "LE_vs_LT (L. timidus reference, CM IDs; LFC sign flipped to LT-vs-LE)",
  flip_lfc = TRUE
)
deseq <- bind_rows(deseq_le, deseq_lt)

gene_set_tbl <- bind_rows(lapply(names(gene_sets), function(gs) {
  tibble(gene_id = read_fasta_ids(gene_sets[[gs]]), GeneSet = gs)
}))

annotations <- bind_rows(lapply(kobas_tsv_files, function(f) parse_kobas_blast(f)))
if (is.null(annotations) || nrow(annotations) == 0) {
  stop("No KOBAS BLAST annotations parsed.")
}

pathway_long <- bind_rows(lapply(identify_files, function(x) {
  parse_identify_kegg(x[[1]], x[[2]])
}))

# Per-gene pathway summary (KEGG, mmu, sig)
pathway_by_gene <- pathway_long %>%
  group_by(gene_id, GeneSet) %>%
  summarise(
    n_pathways = n_distinct(Term),
    best_pathway = Term[which.min(FDR)][1],
    best_pathway_FDR = min(FDR, na.rm = TRUE),
    all_pathways = paste(sort(unique(Term)), collapse = "; "),
    in_seasonal_module = any(Term %in% seasonal_pathways),
    seasonal_pathways = paste(intersect(unique(Term), seasonal_pathways), collapse = "; "),
    only_generic = all(Term %in% generic_pathways),
    .groups = "drop"
  )

# =============================================================================
# PART 1: Rank all candidate genes
# =============================================================================

candidates <- gene_set_tbl %>%
  left_join(deseq, by = "gene_id") %>%
  left_join(annotations, by = "gene_id") %>%
  left_join(pathway_by_gene, by = c("gene_id", "GeneSet")) %>%
  mutate(
    gene_symbol = clean_gene_symbol(gene_symbol),
    GeneSet_label = gene_set_labels[GeneSet],
    id_assembly = ifelse(grepl("^CM", gene_id), "L. timidus (CM)",
                         ifelse(grepl("^NC", gene_id), "L. europaeus (NC)", "other")),
    expression_note = case_when(
      is.na(log2FC) ~ "Gene ID not found in DESeq2 results",
      TRUE ~ NA_character_
    ),
    has_symbol = !is.na(gene_symbol) & gene_symbol != "",
    n_pathways = replace_na(n_pathways, 0L),
    in_seasonal_module = replace_na(in_seasonal_module, FALSE),
    only_generic = replace_na(only_generic, TRUE),
    # Scoring components (higher = stronger candidate)
    score_DE = ifelse(is_DE, pmin(-log10(pmax(padj, 1e-300)), 50), 0),
    score_LFC = pmin(abs(replace_na(log2FC, 0)), 5),
    score_symbol = ifelse(has_symbol, 3, 0),
    score_pathway = pmin(n_pathways, 5) * 0.5,
    score_seasonal = ifelse(in_seasonal_module, 8, 0),
    score_generic_penalty = ifelse(only_generic & n_pathways > 0, -3, 0),
    # Bonus if DE direction matches species-specific expectation
    score_direction_match = case_when(
      GeneSet == "LE_unique_158" & regulation == "Up in LE" ~ 5,
      GeneSet == "LT_unique_133" & regulation == "Up in LT" ~ 5,
      GeneSet == "Shared_508" & is_DE ~ 2,
      TRUE ~ 0
    ),
    priority_score = score_DE + score_LFC + score_symbol + score_pathway +
      score_seasonal + score_generic_penalty + score_direction_match,
    rank_within_set = NA_integer_
  ) %>%
  group_by(GeneSet) %>%
  mutate(rank_within_set = rank(-priority_score, ties.method = "min")) %>%
  ungroup() %>%
  arrange(GeneSet, rank_within_set) %>%
  select(
    rank_within_set, GeneSet, GeneSet_label, gene_id, id_assembly, gene_symbol, ref_gene_id,
    deseq_source, log2FC, log2FC_raw, padj, regulation, is_DE, baseMean, expression_note,
    n_pathways, best_pathway, best_pathway_FDR,
    in_seasonal_module, seasonal_pathways,
    priority_score,
    score_DE, score_LFC, score_symbol, score_seasonal, score_direction_match,
    all_pathways
  )

write.csv(candidates,
          file.path(out_dir, "ranked_candidate_genes.csv"),
          row.names = FALSE)

top_candidates <- candidates %>%
  filter(rank_within_set <= top_n_report) %>%
  mutate(
    one_line_note = case_when(
      in_seasonal_module ~ paste0("Seasonal module (", seasonal_pathways, ")"),
      GeneSet == "LE_unique_158" & regulation == "Up in LE" ~ "LE-unique & higher in brown hare",
      GeneSet == "LT_unique_133" & regulation == "Up in LT" ~ "LT-unique & higher in mountain hare",
      is_DE ~ paste0("DE (", regulation, ")"),
      TRUE ~ "In gene set; not DE at current thresholds"
    )
  ) %>%
  select(rank_within_set, GeneSet_label, gene_id, gene_symbol,
         log2FC, padj, regulation, best_pathway, priority_score, one_line_note)

write.csv(top_candidates,
          file.path(out_dir, "top_candidates_by_gene_set.csv"),
          row.names = FALSE)

cat("Saved: ranked_candidate_genes.csv (", nrow(candidates), " genes)\n", sep = "")
cat("Saved: top_candidates_by_gene_set.csv (top ", top_n_report, " per set)\n", sep = "")

# =============================================================================
# PART 2: Seasonal adaptation module + expression direction
# =============================================================================

module_genes <- pathway_long %>%
  filter(GeneSet == "Shared_508", Term %in% seasonal_pathways) %>%
  distinct(gene_id, Term, FDR) %>%
  left_join(deseq, by = "gene_id") %>%
  left_join(annotations, by = "gene_id") %>%
  group_by(gene_id) %>%
  summarise(
    gene_symbol = first(na.omit(gene_symbol)),
    ref_gene_id = first(na.omit(ref_gene_id)),
    log2FC = first(log2FC),
    padj = first(padj),
    baseMean = first(baseMean),
    regulation = first(regulation),
    is_DE = first(is_DE),
    module_pathways = paste(sort(unique(Term)), collapse = "; "),
    n_module_pathways = n_distinct(Term),
    best_module_FDR = min(FDR),
    .groups = "drop"
  ) %>%
  mutate(
    expression_interpretation = case_when(
      regulation == "Up in LT" ~ "Higher in mountain hare (LT)",
      regulation == "Up in LE" ~ "Higher in brown hare (LE)",
      log2FC > 0 ~ "Trend higher in LT (NS)",
      log2FC < 0 ~ "Trend higher in LE (NS)",
      TRUE ~ "No expression data / flat"
    ),
    module_rank = rank(-(n_module_pathways * 2 +
                           ifelse(is_DE, -log10(pmax(padj, 1e-10)), 0) +
                           abs(replace_na(log2FC, 0)))),
    validation_priority = case_when(
      gene_symbol %in% c("Mc1r", "Mtnr1a", "Tyr", "Mitf", "Asip", "Pomc") ~
        "High (coat-colour / pigment literature)",
      is_DE & n_module_pathways >= 2 ~ "High (multi-pathway module + DE)",
      is_DE ~ "Medium-high (module + DE)",
      !is.na(gene_symbol) & gene_symbol != "" ~ "Medium (annotated)",
      TRUE ~ "Lower (weak evidence)"
    )
  ) %>%
  arrange(module_rank)

write.csv(module_genes,
          file.path(out_dir, "seasonal_module_genes_with_direction.csv"),
          row.names = FALSE)
cat("Saved: seasonal_module_genes_with_direction.csv (",
    nrow(module_genes), " unique genes)\n", sep = "")

# Pathway-level summary: how many module genes go each direction per pathway
module_pathway_summary <- pathway_long %>%
  filter(GeneSet == "Shared_508", Term %in% seasonal_pathways) %>%
  left_join(deseq %>% select(gene_id, regulation, log2FC, padj, is_DE), by = "gene_id") %>%
  group_by(Term, ID) %>%
  summarise(
    n_genes = n_distinct(gene_id),
    n_up_LT = sum(regulation == "Up in LT", na.rm = TRUE),
    n_up_LE = sum(regulation == "Up in LE", na.rm = TRUE),
    n_NS = sum(regulation == "Not significant" | is.na(regulation), na.rm = TRUE),
    mean_log2FC = round(mean(log2FC, na.rm = TRUE), 3),
    dominant_direction = case_when(
      n_up_LT > n_up_LE & n_up_LT >= n_NS ~ "Mostly up in LT (mountain hare)",
      n_up_LE > n_up_LT & n_up_LE >= n_NS ~ "Mostly up in LE (brown hare)",
      TRUE ~ "Mixed / not significant"
    ),
    example_genes = paste(
      head(unique(gene_id[regulation != "Not significant"]), 5),
      collapse = ", "
    ),
    .groups = "drop"
  ) %>%
  arrange(Term)

write.csv(module_pathway_summary,
          file.path(out_dir, "seasonal_module_pathway_direction_summary.csv"),
          row.names = FALSE)
cat("Saved: seasonal_module_pathway_direction_summary.csv\n")

# =============================================================================
# PART 3: Plain-text meeting summary
# =============================================================================

n_de <- candidates %>% group_by(GeneSet_label) %>%
  summarise(
    n = n(),
    n_DE = sum(is_DE),
    n_up_LT = sum(regulation == "Up in LT"),
    n_up_LE = sum(regulation == "Up in LE"),
    .groups = "drop"
  )

display_symbol <- function(sym, gid) {
  ifelse(!is.na(sym) & sym != "", sym, gid)
}

top3 <- top_candidates %>%
  filter(rank_within_set <= 3) %>%
  mutate(line = sprintf(
    "  %d. %s (%s) — %s | LFC=%.2f padj=%.2e | %s",
    rank_within_set, display_symbol(gene_symbol, gene_id), gene_id,
    regulation, log2FC, padj, one_line_note
  ))

summary_lines <- c(
  "======================================================================",
  "  MEETING SUMMARY — Lepus europaeus vs Lepus timidus (skin, mmu KOBAS)",
  paste("  Generated:", format(Sys.time(), "%Y-%m-%d %H:%M")),
  "======================================================================",
  "",
  "DESeq2 sources:",
  "  - Shared_508 + LE_unique_158: LT_vs_LE (europaeus ref, NC IDs)",
  "  - LT_unique_133: LE_vs_LT (timidus ref, CM IDs); LFC flipped so + = higher in LT",
  "",
  "GENE SETS & DIFFERENTIAL EXPRESSION (|log2FC|>=1, padj<0.05):",
  paste(sprintf("  - %s: %d genes, %d DE (%d up LT, %d up LE)",
                n_de$GeneSet_label, n_de$n, n_de$n_DE, n_de$n_up_LT, n_de$n_up_LE),
        collapse = "\n"),
  "",
  "SEASONAL / COAT-RELATED MODULE (Shared 508, KEGG mmu):",
  sprintf("  - %d unique genes in melanogenesis/circadian/hormone pathways",
          nrow(module_genes)),
  sprintf("  - %d with significant DE; %d higher in LT; %d higher in LE",
          sum(module_genes$is_DE),
          sum(module_genes$regulation == "Up in LT", na.rm = TRUE),
          sum(module_genes$regulation == "Up in LE", na.rm = TRUE)),
  "",
  "KEY NAMED GENES IN MODULE (for discussion):",
  paste("  -",
        paste(
          module_genes %>%
            filter(!is.na(gene_symbol), gene_symbol != "") %>%
            arrange(module_rank) %>%
            head(12) %>%
            mutate(s = sprintf("%s (%s, %s)", gene_symbol, regulation,
                               str_trunc(module_pathways, 35))) %>%
            pull(s),
          collapse = "\n  - "
        )),
  "",
  "TOP CANDIDATES PER GENE SET (see top_candidates_by_gene_set.csv):",
  paste(top3$line, collapse = "\n"),
  "",
  "SUGGESTED VALIDATION TARGETS (qPCR / in situ):",
  paste("  -",
        module_genes %>%
          filter(validation_priority %in% c(
            "High (coat-colour / pigment literature)",
            "High (multi-pathway module + DE)",
            "High (module + DE)"
          )) %>%
          filter(!is.na(gene_symbol), gene_symbol != "") %>%
          head(8) %>%
          mutate(s = sprintf("%s — %s", gene_symbol, expression_interpretation)) %>%
          pull(s),
        collapse = "\n  - "),
  "",
  "FILES IN mmu_analysis/tables/:",
  "  ranked_candidate_genes.csv",
  "  top_candidates_by_gene_set.csv",
  "  seasonal_module_genes_with_direction.csv",
  "  seasonal_module_pathway_direction_summary.csv",
  "",
  "======================================================================"
)

writeLines(summary_lines, file.path(out_dir, "seasonal_module_summary.txt"))
cat("Saved: seasonal_module_summary.txt\n\n")

cat("--- Top 5 candidates per gene set ---\n")
print(
  top_candidates %>%
    filter(rank_within_set <= 5) %>%
    select(rank_within_set, GeneSet_label, gene_symbol, log2FC, padj, one_line_note)
)

cat("\n--- Seasonal module (top 10 by rank) ---\n")
print(
  module_genes %>%
    head(10) %>%
    select(gene_symbol, regulation, log2FC, padj, module_pathways, validation_priority)
)

cat("\nDone. Tables in:", out_dir, "\n")
