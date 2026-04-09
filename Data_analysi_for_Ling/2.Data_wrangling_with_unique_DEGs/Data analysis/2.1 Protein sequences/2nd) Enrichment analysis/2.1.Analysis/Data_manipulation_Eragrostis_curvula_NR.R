#####
#packages
rm(list=ls())
library(readxl)
library(tidyverse)
library(dplyr)
#set directory
#setwd("E:/Transcriptomic project on Eragrostis curvula/Processed data/1.2 DEGs/Unique DEGs/Augusto_Unique_DEGs/1. Protein sequences/2nd) Enrichment analysis")
#setwd("E:/Transcriptomic project on Eragrostis curvula/Processed data/1.2 DEGs/3. Unique DEGs_adjusted_Emidio&Jose/3.3 Data analysis/3.3.1 Protein sequences/2nd) Enrichment analysis")
#setwd("E:/Transcriptomic project on Eragrostis curvula/Processed data/1.5 DEGs adjusted according to the Josè's recommendations/out/45. Data analysis/3.3.1 Protein sequences/2nd) Enrichment analysis")
#setwd("E:/Transcriptomic project on Eragrostis curvula/Processed data/1.5 DEGs adjusted according to the Josè's recommendations/out/45. Data analysis/3.3.2 Nucleotide sequences/2nd) Enrichment analysis")
#load data
DEGs_508<- read_excel("LT-LE_shared.xlsx")
DEGs_158 <- read_excel("LE_unique.xlsx")
DEGs_133 <- read_excel("LT_unique.xlsx")

#Reference <- read.table("UnigeneLepus_europaeus_shared_proteins.fa", sep = "", header = F, row.names = NULL)
Reference <- read.table("UnigeneLepus_timidus_unique_proteins.txt", sep = "", header = F, row.names = NULL)

# Creiamo un vettore per i codici dei geni e uno per le sequenze
# Il commando di seguito riportato mi ha permesso di capire dove stava l'errore. Difatti, abbiamo appurato una gene la cui sequenza non veniva riportata per intero ma, invece, veniva "frammentata" tale da costituire due sequenze diverse: in questo modo, le righe non erano uguali, alla luce di questo gene in più.

gene_codes <- Reference$V1[seq(1, nrow(Reference), by = 2)]
sequences <- Reference$V1[seq(2, nrow(Reference), by = 2)]
# Controlla se il numero di righe è diverso e, in tal caso, aggiungi un valore NA a sequences
if (length(gene_codes) != length(sequences)) {
  sequences <- c(sequences, NA)
}

# data wrangling
Genes508 <- DEGs_508[,1]
Genes158 <- DEGs_158[,1]
Genes133 <- DEGs_133[,1]

# Manipulation ####
ref <- data.frame(Column1 = Reference$V1[c(TRUE, FALSE)], 
           Column2 = Reference$V1[c(FALSE, TRUE)])

# renaming column name in order to perform right join later on
colnames(Genes508) <- "Column1"
colnames(Genes158) <- "Column1"
colnames(Genes133) <- "Column1"

# series of joins between different df against a reference (ref)
int_ref_Genes508 <- right_join(ref, Genes508)
int_ref_Genes158 <- right_join(ref, Genes158)
int_ref_Genes133 <- right_join(ref, Genes133)

###########################
int_ref_Genes508_adj <- data.frame(data=c(t(int_ref_Genes508)))
int_ref_Genes158_adj <- data.frame(data=c(t(int_ref_Genes158)))
int_ref_Genes133_adj <- data.frame(data=c(t(int_ref_Genes133)))

# saving
write.table(int_ref_Genes508_adj, file = "1.int_ref_Genes508_adj.txt", append = FALSE, sep="\t", row.names=FALSE, col.names = FALSE)
write.table(int_ref_Genes158_adj, file = "1.int_ref_Genes158_adj.txt", append = FALSE, sep="\t", row.names=FALSE, col.names = FALSE)
write.table(int_ref_Genes133_adj, file = "1.int_ref_Genes133_adj.txt", append = FALSE, sep="\t", row.names=FALSE, col.names = FALSE)

###########################