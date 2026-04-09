rm(list=ls())

options(stringsAsFactors = F)

setwd("C:/Users/danil/OneDrive/Documenti/Didattica/code")
######################################
source("src/script/getLibrary.R")
source("src/script/getSource.R")
#source("src/script/lib/ResilienceAnalysis/ResilienceAnalysis.R")
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

# >>> NEW: run stability diagnostics (subsampling + permutation)
normalizePath("src/script/stability.R")
source("src/script/stability.R")
txt <- deparse(run_stability)
any(grepl("Could not find adjacency in 'network' object", txt))

cat(paste(grep("run_once <- function", txt):length(txt) |> 
            (\(i) txt[i])(), collapse = "\n"))

run_stability(data = data, input_parameter = input_parameter, input_file = input_file)


if(input_parameter$removal_node == "yes") resilience <- ResilienceAnalysis()

