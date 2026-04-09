

##############################################################################

##                      KEGG enrichment analysis

##############################################################################


# ----------------------------------------
# Load required libraries for KEGG graph
# ----------------------------------------
# Uncomment the following lines if the packages are not installed
# install.packages("ggplot2")
# install.packages("readxl")
# install.packages("tidyr")


library(ggplot2)
library(readxl)
library(tidyr)


# ----------------------------------------
# Load the KEGG enrichment data
# ----------------------------------------
data <- read.table("Mouse.txt", 
                   header = TRUE, 
                   sep = "\t", 
                   strip.white = TRUE)

# ----------------------------------------


# Extract columns for clarity (optional)
# ----------------------------------------
Description <- data$Description
GeneRatio   <- data$GeneRatio
Count       <- data$Count
padj        <- data$padj
Comparison  <- data$Comparison

# ----------------------------------------


# KEGG enrichment bubble plot
# ----------------------------------------
kegg_plot <- ggplot(data, aes(
  y = Description,
  x = GeneRatio,
  size = Count,
  fill = padj,
  color = Comparison,
  shape = Comparison)) +
  
  # Draw points with sizes corresponding to Count
  geom_point() + 
  
  # Set x-axis to start from 0 and limit to 0.65
  scale_x_continuous(expand = c(0, 0), limits = c(0, 0.65)) +
  
  # Set bubble sizes and define size legend breaks
  scale_size_continuous(range = c(4, 10), breaks = c(2, 5, 10, 15, 20)) +
  
  # Set color gradient for padj values
  scale_fill_gradientn(colours = rainbow(5), name = "Adjusted p-value (padj)", limits = c(0, 1)) +
  
  # Define shape for each comparison group (all same shape: 16 = solid circle)
  scale_shape_manual(values = c(
    "Genes508" = 16,
    #"2.Unique_571 DEGs" = 16,
    "Genes158" = 16,
    "Genes133" = 16)) +
  
  # Customize the plot theme
  theme(
    panel.background = element_rect(fill = 'white', colour = "gray"),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),       # Border around plot panel
    plot.background = element_rect(fill = "white", color = "black", size = 1),# Outer border of the whole plot
    panel.grid = element_line(colour = "grey"),
    
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.text.x  = element_text(size = 17, face = "bold"),
    axis.text.y  = element_text(size = 20, face = "bold"),
    
    plot.title   = element_text(size = 20, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text  = element_text(size = 18, face = "bold"),
    legend.position = "right",
    legend.box = "vertical",
    strip.text = element_text(size = 18, face = "bold")) +
  
  # Title of the plot
  ggtitle("KEGG Pathway Enrichment")

# ----------------------------------------
# Display the plot
# ----------------------------------------
print(kegg_plot)

# ----------------------------------------
# Save the plot as a PDF
# ----------------------------------------
ggsave("Rabbit3.pdf", kegg_plot, width = 20, height = 20, units = "in", dpi = 600)

