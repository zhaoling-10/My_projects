# =============================================================================
# Gene-level drill-down for significant KOBAS pathways
# Run after KOBAS_enrichment_analysis.R
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(readr)

base_dir <- "~/Documents/Enrichment_analysis"   # <-- adjust if needed

# =============================================================================
# 1. Parse the _identify.txt files to extract genes-per-pathway
# =============================================================================

parse_identify_with_genes <- function(filepath, gene_set_label, ref_species, db_type) {
  if (!file.exists(filepath)) { warning("Not found: ", filepath); return(NULL) }
  lines <- readLines(filepath)
  data_lines <- lines[!grepl("^#", lines) & nchar(trimws(lines)) > 0]
  if (length(data_lines) == 0) return(NULL)

  df <- read.table(text = paste(data_lines, collapse = "\n"),
                   sep = "\t", header = FALSE, quote = "", fill = TRUE,
                   stringsAsFactors = FALSE)
  # Cols: Term | Database | ID | Input_n | Background_n | PValue | FDR | Genes | Hyperlink
  n <- min(ncol(df), 9)
  cnames <- c("Term","Database","ID","Input_n","Background_n",
               "PValue","FDR","Genes","Hyperlink")[1:n]
  colnames(df)[1:n] <- cnames
  df$GeneSet    <- gene_set_label
  df$RefSpecies <- ref_species
  df$DB_type    <- db_type
  df$Input_n    <- as.numeric(df$Input_n)
  df$FDR        <- as.numeric(df$FDR)
  df
}

identify_files <- list(
  list(path = file.path(base_dir,"KOBAS_output_mmu_KEGG_1","KOBAS_output_mmu_KEGG_1.tsv_identify.txt"),
       gs="Shared_508",   ref="mmu", db="KEGG"),
  list(path = file.path(base_dir,"KOBAS_output_mmu_KEGG_2","KOBAS_output_mmu_KEGG_2.tsv_identify.txt"),
       gs="LE_unique_158",ref="mmu", db="KEGG"),
  list(path = file.path(base_dir,"KOBAS_output_mmu_KEGG_3","KOBAS_output_mmu_KEGG_3.tsv_identify.txt"),
       gs="LT_unique_133",ref="mmu", db="KEGG"),
  list(path = file.path(base_dir,"KOBAS_output_mmu_GO_1","KOBAS_output_mmu_GO_1.tsv_identify.txt"),
       gs="Shared_508",   ref="mmu", db="GO"),
  list(path = file.path(base_dir,"KOBAS_output_mmu_GO_2","KOBAS_output_mmu_GO_2.tsv_identify.txt"),
       gs="LE_unique_158",ref="mmu", db="GO"),
  list(path = file.path(base_dir,"KOBAS_output_mmu_GO_3","KOBAS_output_mmu_GO_3.tsv_identify.txt"),
       gs="LT_unique_133",ref="mmu", db="GO")
)

all_data <- bind_rows(lapply(identify_files, function(x)
  parse_identify_with_genes(x$path, x$gs, x$ref, x$db)))

sig <- all_data %>% filter(!is.na(FDR), FDR <= 0.05)

# =============================================================================
# 2. Expand genes: one row per gene-pathway pair
# =============================================================================

gene_pathway <- sig %>%
  filter("Genes" %in% colnames(sig) | ncol(sig) >= 8) %>%
  mutate(gene_list = str_split(Genes, "\\|")) %>%
  unnest(gene_list) %>%
  mutate(gene_list = trimws(gene_list)) %>%
  filter(nchar(gene_list) > 0)

# Save full gene-pathway table
write.csv(gene_pathway %>%
            select(GeneSet, DB_type, RefSpecies, Term, ID, FDR, Input_n, gene_list),
          "gene_pathway_membership.csv", row.names = FALSE)
cat("Saved: gene_pathway_membership.csv\n")

# =============================================================================
# 3. Count how many pathways each gene appears in (multifunctionality score)
# =============================================================================

gene_freq <- gene_pathway %>%
  group_by(GeneSet, gene_list) %>%
  summarise(n_pathways = n_distinct(Term), .groups = "drop") %>%
  arrange(desc(n_pathways))

write.csv(gene_freq, "gene_pathway_frequency.csv", row.names = FALSE)
cat("Saved: gene_pathway_frequency.csv\n")

# Top hub genes per gene set
cat("\n--- Top 10 hub genes (most pathway memberships) ---\n")
gene_freq %>%
  group_by(GeneSet) %>%
  slice_max(n_pathways, n = 10) %>%
  print(n = 50)

# =============================================================================
# 4. Bar plot: top hub genes for LE-unique (most biologically interesting)
# =============================================================================

plot_hub_genes <- function(data, gene_set, top_n = 20, db = "KEGG") {
  df <- data %>%
    filter(GeneSet == gene_set, DB_type == db) %>%
    group_by(gene_list) %>%
    summarise(n_pathways = n_distinct(Term), .groups = "drop") %>%
    arrange(desc(n_pathways)) %>%
    slice_head(n = top_n)

  if (nrow(df) == 0) return(NULL)

  ggplot(df, aes(x = reorder(gene_list, n_pathways), y = n_pathways)) +
    geom_bar(stat = "identity", fill = "#2171B5") +
    coord_flip() +
    labs(title = paste0("Hub genes: ", gene_set, " (", db, ", mmu)"),
         subtitle = "Number of significant KEGG pathways each gene belongs to",
         x = "Gene ID", y = "Number of enriched pathways") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))
}

for (gs in c("Shared_508","LE_unique_158","LT_unique_133")) {
  p <- plot_hub_genes(gene_pathway, gs, top_n = 20, db = "KEGG")
  if (!is.null(p)) {
    fname <- paste0("hub_genes_KEGG_", gs, "_mmu.pdf")
    ggsave(fname, p, width = 10, height = 7, dpi = 300)
    cat("Saved:", fname, "\n")
  }
}

# =============================================================================
# 5. Pathway-gene membership heatmap for the top priority pathways
# =============================================================================

# Focus on the biologically most interesting: LE-unique KEGG
priority_pathways <- c(
  "Herpes simplex virus 1 infection",
  "Ribosome",
  "Spliceosome",
  "Steroid biosynthesis",
  "Metabolism of xenobiotics by cytochrome P450"
)

for (gs in c("LE_unique_158","LT_unique_133","Shared_508")) {
  sub <- gene_pathway %>%
    filter(GeneSet == gs, DB_type == "KEGG", RefSpecies == "mmu") %>%
    select(Term, gene_list)

  if (nrow(sub) < 2) next

  # Binary membership matrix
  mat <- sub %>%
    mutate(present = 1) %>%
    pivot_wider(names_from = gene_list, values_from = present, values_fill = 0) %>%
    as.data.frame()
  rownames(mat) <- str_trunc(mat$Term, 50)
  mat$Term <- NULL
  mat <- as.matrix(mat)

  if (nrow(mat) < 2 || ncol(mat) < 2) next

  # Cluster and save
  fname <- paste0("pathway_gene_heatmap_", gs, "_mmu.pdf")
  pheatmap::pheatmap(
    mat,
    color       = c("white","#2171B5"),
    cluster_rows = TRUE, cluster_cols = TRUE,
    show_colnames = (ncol(mat) <= 60),
    fontsize_row = 7, fontsize_col = 6,
    main = paste0("Pathway × Gene membership: ", gs, " (KEGG, mmu)"),
    filename = fname,
    width = max(10, ncol(mat)*0.15),
    height = max(6, nrow(mat)*0.3)
  )
  cat("Saved:", fname, "\n")
}

# =============================================================================
# 6. Cross-species Venn of enriched KEGG pathway IDs (mmu)
# =============================================================================
# Which pathway TERMS are enriched in exactly one vs multiple gene sets?

term_membership <- sig %>%
  filter(RefSpecies == "mmu", DB_type == "KEGG") %>%
  select(Term, GeneSet) %>%
  distinct() %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = GeneSet, values_from = present, values_fill = 0)

# Categorize
term_membership <- term_membership %>%
  mutate(
    category = case_when(
      Shared_508 == 1 & LE_unique_158 == 0 & LT_unique_133 == 0 ~ "Shared only",
      Shared_508 == 0 & LE_unique_158 == 1 & LT_unique_133 == 0 ~ "LE only",
      Shared_508 == 0 & LE_unique_158 == 0 & LT_unique_133 == 1 ~ "LT only",
      Shared_508 == 1 & LE_unique_158 == 1 & LT_unique_133 == 0 ~ "Shared + LE",
      Shared_508 == 1 & LE_unique_158 == 0 & LT_unique_133 == 1 ~ "Shared + LT",
      Shared_508 == 0 & LE_unique_158 == 1 & LT_unique_133 == 1 ~ "LE + LT",
      TRUE ~ "All three"
    )
  )

write.csv(term_membership, "KEGG_pathway_overlap_categories.csv", row.names = FALSE)
cat("Saved: KEGG_pathway_overlap_categories.csv\n")

cat("\n--- KEGG pathway overlap summary (mmu) ---\n")
print(table(term_membership$category))

# Bar chart version of the Venn overlap
overlap_summary <- term_membership %>%
  count(category) %>%
  mutate(category = factor(category, levels = c(
    "Shared only","LE only","LT only",
    "Shared + LE","Shared + LT","LE + LT","All three")))

p_overlap <- ggplot(overlap_summary, aes(x = category, y = n, fill = category)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = n), vjust = -0.4, size = 4) +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  labs(title = "KEGG pathway overlap across gene sets (mmu, FDR ≤ 0.05)",
       x = "Enrichment category", y = "Number of pathways") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        plot.title = element_text(face = "bold"))

ggsave("KEGG_pathway_overlap_barplot.pdf", p_overlap, width = 9, height = 6, dpi = 300)
cat("Saved: KEGG_pathway_overlap_barplot.pdf\n")

cat("\n===== Gene drill-down analysis complete =====\n")
