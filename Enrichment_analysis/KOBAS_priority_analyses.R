# =============================================================================
# Priority next analyses after KOBAS enrichment
# 1. Annotate gene IDs with functional descriptions from KOBAS blast output
# 2. Integrate DESeq2 fold-change data with enriched pathways
# 3. Highlight the seasonal adaptation gene module
# 4. Produce a publication-ready combined figure
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(readr)
library(pheatmap)
library(ggrepel)
library(RColorBrewer)

base_dir <- "~/Documents/Enrichment_analysis"

# =============================================================================
# PART 1: Annotate gene IDs using KOBAS BLAST output (.tsv files)
# These map your gene IDs (NC_084827.1-g413) to known gene symbols
# =============================================================================

parse_kobas_blast <- function(filepath, gene_set_label) {
  # The main .tsv file from KOBAS has format:
  # QueryGeneID | RefGeneID | Gene name + aliases | URL
  if (!file.exists(filepath)) { warning("Not found: ", filepath); return(NULL) }
  lines <- readLines(filepath)
  data_lines <- lines[!grepl("^#", lines) & nchar(trimws(lines)) > 0]
  if (length(data_lines) == 0) return(NULL)

  df <- read.table(text = paste(data_lines, collapse="\n"),
                   sep="\t", header=FALSE, quote="", fill=TRUE,
                   stringsAsFactors=FALSE)
  # Cols: query_gene | ref_gene_id|symbol|aliases|url
  colnames(df)[1:min(2,ncol(df))] <- c("query_gene","ref_info")[1:min(2,ncol(df))]
  df$GeneSet <- gene_set_label

  # Extract gene symbol from ref_info (format: "mmu:12345|Gapdh, alias1, alias2|url")
  if ("ref_info" %in% colnames(df)) {
    df <- df %>%
      mutate(
        ref_gene_id  = str_extract(ref_info, "^[^|]+"),
        gene_symbol  = str_extract(ref_info, "(?<=\\|)[^|,]+"),
        gene_aliases = str_extract(ref_info, "(?<=, ).*(?=\\|)")
      )
  }
  df
}

blast_files <- list(
  list(path = file.path(base_dir,"KOBAS_output_mmu_KEGG_1",
                        "KOBAS_output_mmu_KEGG_1.tsv"),
       gs = "Shared_508"),
  list(path = file.path(base_dir,"KOBAS_output_mmu_KEGG_2",
                        "KOBAS_output_mmu_KEGG_2.tsv"),
       gs = "LE_unique_158"),
  list(path = file.path(base_dir,"KOBAS_output_mmu_KEGG_3",
                        "KOBAS_output_mmu_KEGG_3.tsv"),
       gs = "LT_unique_133")
)

annotations <- bind_rows(lapply(blast_files, function(x)
  parse_kobas_blast(x$path, x$gs)))

# Clean up: keep best hit per query gene (first line = best BLAST hit)
annotation_best <- annotations %>%
  group_by(query_gene, GeneSet) %>%
  slice(1) %>%
  ungroup() %>%
  select(query_gene, GeneSet, ref_gene_id, gene_symbol, gene_aliases)

write.csv(annotation_best, "gene_annotations_best_hit.csv", row.names=FALSE)
cat("Saved: gene_annotations_best_hit.csv\n")

cat("\nAnnotation coverage per gene set:\n")
print(annotation_best %>% count(GeneSet))

# =============================================================================
# PART 2: Integrate DESeq2 results with enriched pathways
# Requires your DESeq2 output CSVs. Adjust file paths to match your setup.
# The expected format: gene_id, log2FoldChange, padj
# =============================================================================

# Load DESeq2 results (LT vs LE, brown hare reference)
deseq_file_le <- "LT_vs_LE_adj.csv"   # from your DESeq2 analysis
deseq_file_lt <- "LT_vs_LE_adj.csv"   # adjust if you have LT-ref results separately

load_deseq <- function(filepath, label) {
  if (!file.exists(filepath)) {
    message("DESeq2 file not found: ", filepath,
            "\nSkipping pathway-expression integration for ", label)
    return(NULL)
  }
  df <- read.csv(filepath, stringsAsFactors=FALSE)
  # Normalise column names
  colnames(df) <- tolower(colnames(df))
  if (!"gene_id" %in% colnames(df) && !is.null(rownames(df))) {
    df$gene_id <- rownames(df)
  }
  # Accept log2foldchange or log2FoldChange etc.
  lfc_col  <- grep("log2fold", colnames(df), ignore.case=TRUE, value=TRUE)[1]
  padj_col <- grep("^padj$|^fdr$|adj.*p|p.*adj", colnames(df), ignore.case=TRUE, value=TRUE)[1]
  if (is.na(lfc_col) || is.na(padj_col)) {
    message("Cannot identify LFC/padj columns in ", filepath); return(NULL)
  }
  df <- df %>%
    rename(log2FC = !!lfc_col, padj = !!padj_col) %>%
    mutate(
      direction   = case_when(log2FC > 1  & padj < 0.05 ~ "Up in LT",
                              log2FC < -1 & padj < 0.05 ~ "Up in LE",
                              TRUE ~ "NS"),
      sig_label   = label
    )
  df
}

deseq_le <- load_deseq(deseq_file_le, "LE_ref")

if (!is.null(deseq_le) && nrow(annotation_best) > 0) {
  # Join DESeq2 LFC to KOBAS annotations
  deseq_annotated <- deseq_le %>%
    left_join(annotation_best, by = c("gene_id" = "query_gene"))

  # Then join to pathway membership
  # Load gene_pathway_membership.csv from previous script
  if (file.exists("gene_pathway_membership.csv")) {
    gpm <- read.csv("gene_pathway_membership.csv", stringsAsFactors=FALSE)

    pathway_deseq <- gpm %>%
      left_join(deseq_le %>% select(gene_id, log2FC, padj, direction),
                by = c("gene_list" = "gene_id")) %>%
      left_join(annotation_best %>% select(query_gene, gene_symbol),
                by = c("gene_list" = "query_gene"))

    write.csv(pathway_deseq, "pathway_gene_with_LFC.csv", row.names=FALSE)
    cat("Saved: pathway_gene_with_LFC.csv\n")

    # Volcano coloured by pathway membership
    top_pathways <- c("Herpes simplex virus 1 infection",
                      "Ribosome", "Glutamatergic synapse",
                      "Metabolic pathways", "Glycolysis / Gluconeogenesis",
                      "Steroid biosynthesis")

    volcano_df <- deseq_le %>%
      mutate(gene_id_clean = gene_id) %>%
      left_join(
        pathway_deseq %>%
          filter(Term %in% top_pathways) %>%
          select(gene_list, Term) %>%
          distinct() %>%
          group_by(gene_list) %>%
          summarise(top_pathway = first(Term), .groups="drop"),
        by = c("gene_id" = "gene_list")
      ) %>%
      mutate(top_pathway = ifelse(is.na(top_pathway), "Other / not enriched", top_pathway),
             neg_log10_padj = pmin(-log10(padj), 50))

    p_volcano_pathway <- ggplot(volcano_df,
                                aes(x=log2FC, y=neg_log10_padj,
                                    color=top_pathway, alpha=top_pathway)) +
      geom_point(size=1.2) +
      scale_color_manual(
        values = c(
          "Herpes simplex virus 1 infection" = "#D62728",
          "Ribosome"                         = "#1F77B4",
          "Glutamatergic synapse"            = "#2CA02C",
          "Metabolic pathways"               = "#FF7F0E",
          "Glycolysis / Gluconeogenesis"     = "#9467BD",
          "Steroid biosynthesis"             = "#8C564B",
          "Other / not enriched"             = "grey75"
        ),
        name = "Top KEGG pathway"
      ) +
      scale_alpha_manual(
        values = c(setNames(rep(0.9, 6), top_pathways),
                   "Other / not enriched" = 0.2),
        guide = "none"
      ) +
      geom_hline(yintercept = -log10(0.05), linetype="dashed", color="grey40") +
      geom_vline(xintercept = c(-1, 1),     linetype="dashed", color="grey40") +
      labs(title   = "Volcano plot coloured by KEGG pathway membership",
           subtitle = "LT vs LE (brown hare reference); dashed lines: |LFC|=1, FDR=0.05",
           x = "log2 Fold Change (LT / LE)",
           y = expression(-log[10](FDR))) +
      theme_bw(base_size=12) +
      theme(plot.title=element_text(face="bold"),
            legend.position="right")

    ggsave("volcano_pathway_coloured.pdf", p_volcano_pathway,
           width=12, height=8, dpi=300)
    cat("Saved: volcano_pathway_coloured.pdf\n")
  }
}

# =============================================================================
# PART 3: Seasonal adaptation gene module
# Focus on the Melanogenesis + Circadian + Thyroid/Cortisol/Aldosterone cluster
# visible in the Shared_508 pathway-gene heatmap
# =============================================================================

seasonal_pathways <- c(
  "Melanogenesis",
  "Circadian entrainment",
  "Thyroid hormone synthesis",
  "Cortisol synthesis and secretion",
  "Aldosterone synthesis and secretion",
  "Cushing syndrome"
)

if (file.exists("gene_pathway_membership.csv")) {
  gpm <- read.csv("gene_pathway_membership.csv", stringsAsFactors=FALSE)

  seasonal_genes <- gpm %>%
    filter(GeneSet == "Shared_508",
           DB_type == "KEGG",
           RefSpecies == "mmu",
           Term %in% seasonal_pathways) %>%
    select(Term, gene_list) %>%
    distinct()

  # Add annotations
  seasonal_annotated <- seasonal_genes %>%
    left_join(annotation_best %>% select(query_gene, gene_symbol, ref_gene_id),
              by = c("gene_list" = "query_gene")) %>%
    mutate(display_name = ifelse(!is.na(gene_symbol) & gene_symbol != "",
                                 paste0(gene_symbol, "\n(", gene_list, ")"),
                                 gene_list))

  write.csv(seasonal_annotated,
            "seasonal_adaptation_genes.csv", row.names=FALSE)
  cat("Saved: seasonal_adaptation_genes.csv\n")
  cat("\nSeasonal adaptation genes found:\n")
  print(seasonal_annotated)

  # Heatmap of seasonal module
  if (nrow(seasonal_genes) >= 2) {
    mat_s <- seasonal_genes %>%
      mutate(present=1) %>%
      pivot_wider(names_from=gene_list, values_from=present, values_fill=0) %>%
      as.data.frame()
    rownames(mat_s) <- mat_s$Term
    mat_s$Term <- NULL
    mat_s <- as.matrix(mat_s)

    pheatmap(
      mat_s,
      color        = c("white","#B22222"),
      cluster_rows = FALSE,
      cluster_cols = TRUE,
      fontsize_row = 10,
      fontsize_col = 7,
      main         = "Seasonal adaptation gene module (Shared genes, KEGG mmu)",
      filename     = "seasonal_adaptation_heatmap.pdf",
      width        = max(8, ncol(mat_s)*0.2),
      height       = 5
    )
    cat("Saved: seasonal_adaptation_heatmap.pdf\n")
  }
}

# =============================================================================
# PART 4: Publication-ready summary figure
# Dot plot showing top 6 pathways per gene set with LFC direction overlay
# =============================================================================

if (file.exists("gene_pathway_membership.csv")) {
  gpm <- read.csv("gene_pathway_membership.csv", stringsAsFactors=FALSE)

  # Pick top 6 most significant pathways per gene set (mmu, KEGG)
  if (file.exists("all_significant_enrichment_FDR05.csv")) {
    all_sig <- read.csv("all_significant_enrichment_FDR05.csv",
                        stringsAsFactors=FALSE)
    top_paths <- all_sig %>%
      filter(RefSpecies=="mmu", DB_type=="KEGG") %>%
      group_by(GeneSet) %>%
      slice_min(Corrected_PValue, n=6) %>%
      ungroup() %>%
      select(GeneSet, Term, Input_n, Corrected_PValue, GeneRatio) %>%
      mutate(neg_log10_fdr = -log10(Corrected_PValue),
             GeneSet = factor(GeneSet,
               levels=c("Shared (508)","LE unique (158)","LT unique (133)")))

    p_pub <- ggplot(top_paths,
                    aes(x=GeneSet, y=reorder(Term, neg_log10_fdr),
                        size=Input_n, color=neg_log10_fdr)) +
      geom_point(alpha=0.9) +
      scale_color_gradient2(
        low="steelblue", mid="#FFC107", high="#D62728",
        midpoint=median(top_paths$neg_log10_fdr, na.rm=TRUE),
        name=expression(-log[10](FDR))
      ) +
      scale_size_continuous(range=c(3,14), name="Gene count") +
      scale_x_discrete(drop=FALSE) +
      labs(
        title    = "Top enriched KEGG pathways per gene set",
        subtitle = "Brown hare (LE) vs Mountain hare (LT) — skin transcriptome",
        x=NULL, y=NULL
      ) +
      theme_bw(base_size=13) +
      theme(
        plot.title    = element_text(face="bold", size=14),
        plot.subtitle = element_text(size=11, color="grey30"),
        axis.text.y   = element_text(size=10),
        axis.text.x   = element_text(size=11, face="bold"),
        legend.position="right",
        panel.grid.major.x = element_blank()
      )

    ggsave("publication_dotplot_KEGG.pdf", p_pub, width=13, height=9, dpi=300)
    ggsave("publication_dotplot_KEGG.png", p_pub, width=13, height=9, dpi=300)
    cat("Saved: publication_dotplot_KEGG.pdf / .png\n")
  }
}

cat("\n===== All priority analyses complete =====\n")
cat("Key outputs:\n")
cat("  gene_annotations_best_hit.csv       — gene IDs mapped to gene symbols\n")
cat("  pathway_gene_with_LFC.csv           — pathways + fold changes\n")
cat("  volcano_pathway_coloured.pdf        — volcano coloured by pathway\n")
cat("  seasonal_adaptation_genes.csv       — the melanogenesis/circadian module\n")
cat("  seasonal_adaptation_heatmap.pdf     — heatmap of seasonal module\n")
cat("  publication_dotplot_KEGG.pdf        — publication-ready summary figure\n")
