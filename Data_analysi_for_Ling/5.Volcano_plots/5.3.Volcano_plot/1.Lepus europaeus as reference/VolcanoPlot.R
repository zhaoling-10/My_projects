# Activate the libraries we need

install.packages("ggpubr")
install.packages("reshape2")
install.packages("clValid")
install.packages("corrplot")
install.packages("factoextra")

library(ggplot2)
library(ggpubr)
library(reshape2)
library(datasets)
library(factoextra)
library(clValid)
library(corrplot)

##Volcanoplot con funzione plot

### LT vs LE
LT_vs_LE=read.table("LT_vs_LE_Novel.txt", header=T, row.names = 1)

tiff('LT vs LE.tiff',3000,3000,res=400)
plot(LT_vs_LE, pch=20,
     col=ifelse(test = LT_vs_LE$Log10padj<1.3,
                yes = "darkgray",no = ifelse(test = LT_vs_LE$log2FoldChange>0,
                                             yes = "red",no = "green")), 
     main="LT vs LE", xlab("Log2FoldChange"), ylab("-Log10padj"),xlim=c(-20,20),ylim=c(0,90),cex.axis=1.0,cex.lab = 1.0, font.lab=2,font.axis = 2)
legend(x = 7.7, y=90, title = "DEG:754",legend = c("Up:406","Down:322"), 
       col = c("red", "green"), pch = 19, cex = 1.1, box.lty = 0, text.font = 2)
abline(h=1.3, lty=6, lwd= 2)
abline(v=c(0.58,-0.58), lty=6, lwd= 2)
dev.off()

help(ggsave)
help(legend) 
help(abline)

#### Save the graph in PDF format
pdf(file = "LT vs LE.pdf", width = 10, height = 10)
plot(LT_vs_LE, pch=20,
     col=ifelse(test = LT_vs_LE$Log10padj<1.3,
                yes = "darkgray",no = ifelse(test = LT_vs_LE$log2FoldChange>0,
                                             yes = "red",no = "green")), 
     main="LT vs LE", xlab("Log2FoldChange"), ylab("-Log10padj"),xlim=c(-20,20),ylim=c(0,90),cex.axis=1.0,cex.lab = 1.0, font.lab=2,font.axis = 2)
legend(x = 13.0, y=90, title = "DEG:754",legend = c("Up:406","Down:322"), 
       col = c("red", "green"), pch = 19, cex = 1.1, box.lty = 0, text.font = 2)
abline(h=1.3, lty=6, lwd= 2)
abline(v=c(0.58,-0.58), lty=6, lwd= 2)
dev.off()



#### Novel

pdf(file = "LT vs LE.pdf", width = 10, height = 10)
plot(LT_vs_LE, pch=20,
     col=ifelse(test = LT_vs_LE$Log10padj<1.3,
                yes = "darkgray",no = ifelse(test = LT_vs_LE$log2FoldChange>0,
                                             yes = "red",no = "green")), 
     main="LT vs LE", xlab("Log2FoldChange"), ylab("-Log10padj"),xlim=c(-20,20),ylim=c(0,90),cex.axis=1.0,cex.lab = 1.0, font.lab=2,font.axis = 2)
legend(x = 13.0, y=90, title = "DEG:754",legend = c("Up:451","Down:303"), 
       col = c("red", "green"), pch = 19, cex = 1.1, box.lty = 0, text.font = 2)
abline(h=1.3, lty=6, lwd= 2)
abline(v=c(0.58,-0.58), lty=6, lwd= 2)
dev.off()
