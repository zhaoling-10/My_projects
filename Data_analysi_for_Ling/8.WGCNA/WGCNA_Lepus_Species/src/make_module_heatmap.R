library("ComplexHeatmap")

make_module_heatmap <- function(module_name,
                                expression_mat = normalized_counts,
                                metadata_df = metadata,
                                groups = condition,
                                gene_module_key_df = gene_module_key,
                                module_eigengenes_df = module_eigengenes) {
  # Create a summary heatmap of a given module.
  #
  # Args:
  # module_name: a character indicating what module should be plotted, e.g. "ME19"
  # expression_mat: The full gene expression matrix. Default is `normalized_counts`.
  # metadata_df: a data frame with sample_code and conditions (factor) as columns. Default is `metadata`.
  # gene_module_key: a data.frame indicating what genes are a part of what modules. Default is `gene_module_key`.
  # module_eigengenes: a sample x eigengene data.frame with samples as row names. Default is `module_eigengenes`.
  #
  # Returns:
  # A heatmap of expression matrix for a module's genes, with a barplot of the
  # eigengene expression for that module.
  
  #Ensure 'cross' is a factor with correct levels corresponding to the color map
  metadata_df$condition <- as.factor(metadata_df$condition)
  
  # Set up the module eigengene with its sample code
  module_eigengene <- module_eigengenes_df %>%
    dplyr::select(all_of(module_name)) %>%
    tibble::rownames_to_column("sample")
  
  # Set up column annotation from metadata
  col_annot_df <- metadata_df %>%
    # Only select the treatment and sample ID columns
    dplyr::select(sample, group) %>%
    # Add on the eigengene expression by joining with sample IDs
    dplyr::inner_join(module_eigengene, by = "sample") %>%
    # Arrange by sample and condition
    dplyr::arrange(cross) %>%
    # Store sample
    tibble::column_to_rownames("sample")
  
  # Define a fixed color palette for the 'cross' levels
  cross_colors <- c("#004777","#a30000",  "#ff7700", "#efd28d", "#00afb5")
  names(cross_colors) <- unique(metadata_df$cross)
  
  # Create the ComplexHeatmap column annotation object
  col_annot <- ComplexHeatmap::HeatmapAnnotation(
    # Supply treatment labels
    condition = col_annot_df$cross,
    # Add annotation barplot
    module_eigengene = ComplexHeatmap::anno_barplot(dplyr::select(col_annot_df, module_name)),
    # Pick colors for each experimental group in conditions
    col = list(cross = cross_colors)
  )
  
  # Get a vector of the Ensembl gene IDs that correspond to this module
  module_genes <- gene_module_key_df %>%
    dplyr::filter(module == module_name) %>%
    dplyr::pull(gene)
  
  # Set up the gene expression data frame
  mod_mat <- expression_mat %>%
    t() %>%
    as.data.frame() %>%
    # Only keep genes from this module
    dplyr::filter(rownames(.) %in% module_genes) %>%
    # Order the samples to match col_annot_df
    dplyr::select(rownames(col_annot_df)) %>%
    # Data needs to be a matrix
    as.matrix()
  
  # scale the gene expression values
  mod_mat <- mod_mat %>%
    # Scale can work on matrices, but it does it by column so we will need to
    # transpose first
    t() %>%
    scale() %>%
    # And now we need to transpose back
    t()
  
  # Create a color function based on standardized scale
  color_func <- circlize::colorRamp2(
    c(-2, 0, 2),
    c("#67a9cf", "#f7f7f7", "#ef8a62")
  )
  
  # Plot on a heatmap
  heatmap <- ComplexHeatmap::Heatmap(mod_mat,
                                     name = module_name,
                                     # Supply color function
                                     col = color_func,
                                     # Supply column annotation
                                     bottom_annotation = col_annot,
                                     # We don't want to cluster samples
                                     cluster_columns = FALSE,
                                     # We don't need to show sample or gene labels
                                     show_row_names = FALSE,
                                     show_column_names = FALSE
  )
  
  # Return heatmap
  return(heatmap)
}