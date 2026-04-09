rm(list=ls())

options(stringsAsFactors = F)

setwd("C:/Users/danil/OneDrive/Documenti/Didattica/code")
######################################
source("src/script/getLibrary.R")
source("src/script/getSource.R")
######################################
getLibrary()
getSource()
input_parameter <- config()
input_file <- inputFiles()
output_file <- outputFiles()
######################################

# Module 1: Exploratory Data Analysis -> filtering and removing genes containing zero values.Indeed, we obtained a matrix of differentially expressed genes.
data <- ExploratoryDataAnalysis()

# Module 3: Network Analysis
network <- NetworkAnalysis(data,checkNetIntegrity = T, screePlot = T)

# Module 4: Switch Mining
switch <- SwitchMining()

saveParameters()

# Stage 5: Stability diagnostics (LOO + 80% subsampling; NO permutation)
source("src/script/stabilityInf_vs_Leaf.R")
run_stability(data = data, input_parameter = input_parameter, input_file = input_file)

# Stage 6: Resilience analysis
if(input_parameter$removal_node == "yes") resilience <- ResilienceAnalysis()
