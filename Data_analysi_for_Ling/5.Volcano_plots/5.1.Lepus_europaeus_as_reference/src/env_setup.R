#Environment setup for WGCNA

#Bioconductor version '3.18' requires R version '4.3', Bioconductor version '3.18' is necessary to install CoExpNets
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.18")

# if (!("DESeq2" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("DESeq2", update = FALSE)
}

if (!("devtools" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  install.packages("devtools")
}

if (!("impute" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("impute")
}

if (!("WGCNA" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  install.packages("WGCNA")
}

if (!("ggforce" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  install.packages("ggforce")
}

if (!("ComplexHeatmap" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("ComplexHeatmap")
}

if (!("dplyr" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  install.packages("dplyr")
}

if (!("tidyr" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  install.packages("tidyr")
}

if (!("readr" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  install.packages("readr")
}

if (!("limma" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("limma")
}

if (!("sva" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("sva")
}

if (!("CorLevelPlot" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  devtools::install_github("kevinblighe/CorLevelPlot")
}
  
if (!requireNamespace(c("graph", "RBGL", "topGO", "GOSim"), quietly = TRUE)){
  BiocManager::install(c("graph", "RBGL", "topGO", "GOSim"))
}

if (!("CoExpNets" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  devtools::install_github('juanbot/CoExpNets', force = T)
}

library(DESeq2)
library(magrittr)
library(WGCNA)
library(ggplot2)
library(sva)
library(gridExtra)
library(CorLevelPlot)
library(CoExpNets)
source("src/make_module_heatmap.R")
source("src/export_to_cytoscape.R")
source("src/plot_genes_in_modules.R")

# Check if the "plots" folder exists; if not, create it
if (!dir.exists("plots")) {
  dir.create("plots")
  cat("Created 'plots' folder.\n")
} else {
  cat("'plots' folder already exists.\n")
}

# Check if the "out" folder exists; if not, create it
if (!dir.exists("out")) {
  dir.create("out")
  cat("Created 'out' folder.\n")
} else {
  cat("'out' folder already exists.\n")
}

# Check if the "out/cytoscape" subdirectory exists; if not, create it
if (!dir.exists("out/cytoscape")) {
  dir.create("out/cytoscape", recursive = TRUE)
  cat("Created 'out/cytoscape' subdirectory.\n")
} else {
  cat("'out/cytoscape' subdirectory already exists.\n")
}