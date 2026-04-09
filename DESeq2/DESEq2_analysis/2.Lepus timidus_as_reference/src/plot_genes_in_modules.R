#plotting number of genes in modules

library(dplyr)
library(ggplot2)

plot_genes_in_modules <- function(network = bwnet){
  # Initialize an empty data frame
  plot_data <- data.frame()
  
  # Loop through each unique module color and count the number of genes
  for (i in unique(network$colors)) {
    plot_data <- rbind(plot_data, data.frame(mod_num = i, Module = paste0("ME", i), n_genes = sum(network$colors == i)))
  }
  
  # Plot using ggplot2
  p <- ggplot(plot_data, aes(x = Module, y = n_genes, fill = Module)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = labels2colors(plot_data$mod_num)) + # Assigns a color per module
    labs(x = "Module", y = "Number of Genes", title = "Number of Genes in Each Module") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p)
}

plot_genes_in_modules_km <- function(network = bwnet){
  # Initialize an empty data frame
  plot_data <- data.frame()
  
  # Loop through each unique module color and count the number of genes
  for (i in unique(network$colors)) {
    plot_data <- rbind(plot_data, data.frame(Module = i, n_genes = sum(network$colors == i)))
  }
  
  plot_data$Module <- factor(plot_data$Module, levels = plot_data$Module)
  
  # Plot using ggplot2
  p <- ggplot(plot_data, aes(x = Module, y = n_genes, fill = Module)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = levels(plot_data$Module)) + # Assigns a color per module
    labs(x = "Module", y = "Number of Genes", title = "Number of Genes in Each Module") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p)
}
