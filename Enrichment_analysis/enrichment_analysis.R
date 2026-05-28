# =============================================================================
# GO & KEGG Enrichment Analysis
# Species: Lepus europaeus & Lepus timidus
# Gene sets: 508 shared DEGs, 158 unique europaeus, 133 unique timidus
# Author: [Your Name]
# Date: 2026-04-22
# =============================================================================
# 
# 脚本用途 (Script Purpose):
#   本脚本对三组基因集（共享508条、europaeus特有、timidus特有）分别进行：
#   1. GO富集分析（BP/MF/CC三个本体）
#   2. KEGG通路富集分析
#   3. 自动生成可视化图表（dotplot, barplot, upset plot, pathway map）
#   4. 输出结果表格供论文使用
#
# =============================================================================


# ─────────────────────────────────────────────
# SECTION 1: 安装 & 加载包
# ─────────────────────────────────────────────
.libPaths()

install.packages("BiocManager")

install.packages("fgsea", type = "source", verbose = TRUE)

# 首次运行时取消注释以下代码安装包：
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    BiocManager::install(c("clusterProfiler", "org.Mm.eg.db", "enrichplot",
                       "pathview", "ggplot2", "dplyr", "openxlsx",
                       "AnnotationDbi", "DOSE"))
install.packages(c("ggplot2", "dplyr", "openxlsx", "ggupset", "patchwork"))

# 先单独安装 fgsea  使用BiocManager来安装fgsea（这会自动解决所有依赖）
BiocManager::install("fgsea")

# 如果 fgsea 成功，再安装 treeio
BiocManager::install("treeio")

if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")
devtools::install_github("YuLab-SMU/ggtree")
BiocManager::install("GOSemSim", force = TRUE)

# 安装 BiocManager
install.packages("BiocManager")

# **关键**：明确指定使用 Bioconductor 3.18（与 R 4.3 匹配）
BiocManager::install(version = "3.18", update = FALSE, ask = FALSE)

# 安装 clusterProfiler（会自动安装所有依赖）
BiocManager::install("clusterProfiler")

packageVersion("GOSemSim")
# 2. 从 Bioconductor 安装稳定版的 enrichplot
BiocManager::install("enrichplot")
# 1. 安装 enrichplot 开发版
devtools::install_github("YuLab-SMU/enrichplot")

# 最后安装 ggtree, DOSE, enrichplot, clusterProfiler
BiocManager::install(c("ggtree", "DOSE", "enrichplot", "clusterProfiler"))

BiocManager::install("clusterProfiler")

library(clusterProfiler)
library(org.Mm.eg.db)       # 小鼠注释数据库（因best hits多为Mus musculus）
library(enrichplot)
library(pathview)
library(ggplot2)
library(dplyr)
library(openxlsx)
library(patchwork)

# ─────────────────────────────────────────────
# SECTION 2: 读取基因列表
# ─────────────────────────────────────────────

# !! 请修改为你的实际文件路径 !!
setwd("your/working/directory")   # 例如 setwd("C:/Users/Ling/Desktop/Lepus_analysis")

genes_508    <- readLines("KOBAS_508_shared.txt")
genes_eu     <- readLines("KOBAS_unique_europaeus.txt")
genes_ti     <- readLines("KOBAS_unique_timidus_clean.txt")

# 去除空行
genes_508 <- genes_508[genes_508 != ""]
genes_eu  <- genes_eu[genes_eu != ""]
genes_ti  <- genes_ti[genes_ti != ""]

cat("基因数量统计:\n")
cat("  508 shared DEGs:", length(genes_508), "\n")
cat("  Unique europaeus:", length(genes_eu), "\n")
cat("  Unique timidus (clean):", length(genes_ti), "\n")


# ─────────────────────────────────────────────
# SECTION 3: Gene Symbol → Entrez ID 转换
# （clusterProfiler的KEGG分析需要Entrez ID）
# ─────────────────────────────────────────────

convert_to_entrez <- function(gene_symbols, db = org.Mm.eg.db) {
  result <- bitr(gene_symbols,
                 fromType = "SYMBOL",
                 toType   = "ENTREZID",
                 OrgDb    = db)
  cat("  转换成功:", nrow(result), "/", length(gene_symbols), "个基因\n")
  return(result$ENTREZID)
}

cat("\n转换 508 shared genes...\n")
entrez_508 <- convert_to_entrez(genes_508)

cat("转换 unique europaeus genes...\n")
entrez_eu  <- convert_to_entrez(genes_eu)

cat("转换 unique timidus genes...\n")
entrez_ti  <- convert_to_entrez(genes_ti)


# ─────────────────────────────────────────────
# SECTION 4: GO 富集分析函数
# ─────────────────────────────────────────────

run_GO <- function(entrez_ids, label, ontology = "ALL") {
  cat("\n正在运行 GO 富集分析:", label, "| Ontology:", ontology, "\n")
  
  ego <- enrichGO(
    gene          = entrez_ids,
    OrgDb         = org.Mm.eg.db,
    keyType       = "ENTREZID",
    ont           = ontology,    # "BP", "MF", "CC", 或 "ALL"
    pAdjustMethod = "BH",        # Benjamini-Hochberg多重检验校正
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    readable      = TRUE         # 将Entrez ID转回基因Symbol显示
  )
  
  cat("  显著GO term数量:", nrow(as.data.frame(ego)), "\n")
  return(ego)
}


# ─────────────────────────────────────────────
# SECTION 5: KEGG 通路富集分析函数
# ─────────────────────────────────────────────

run_KEGG <- function(entrez_ids, label) {
  cat("\n正在运行 KEGG 富集分析:", label, "\n")
  
  kk <- enrichKEGG(
    gene         = entrez_ids,
    organism     = "mmu",        # mmu = Mus musculus
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )
  
  # 将KEGG ID转换为可读基因名
  kk <- setReadable(kk, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
  
  cat("  显著KEGG通路数量:", nrow(as.data.frame(kk)), "\n")
  return(kk)
}


# ─────────────────────────────────────────────
# SECTION 6: 运行所有分析
# ─────────────────────────────────────────────

# GO 分析（三组分别运行）
go_508 <- run_GO(entrez_508, "508_Shared_DEGs")
go_eu  <- run_GO(entrez_eu,  "Unique_europaeus")
go_ti  <- run_GO(entrez_ti,  "Unique_timidus")

# KEGG 分析（三组分别运行）
kegg_508 <- run_KEGG(entrez_508, "508_Shared_DEGs")
kegg_eu  <- run_KEGG(entrez_eu,  "Unique_europaeus")
kegg_ti  <- run_KEGG(entrez_ti,  "Unique_timidus")


# ─────────────────────────────────────────────
# SECTION 7: 可视化 — GO 富集图
# ─────────────────────────────────────────────

dir.create("figures", showWarnings = FALSE)
dir.create("tables",  showWarnings = FALSE)

# 函数：为某一基因集画GO dotplot + barplot
plot_GO <- function(ego, label, top_n = 20) {
  
  if (nrow(as.data.frame(ego)) == 0) {
    cat("  [警告]", label, "无显著GO term，跳过作图\n")
    return(NULL)
  }
  
  # --- Dotplot（气泡图）---
  p_dot <- dotplot(ego, showCategory = top_n, split = "ONTOLOGY") +
    facet_grid(ONTOLOGY ~ ., scale = "free") +
    ggtitle(paste0("GO Enrichment - ", label)) +
    theme(
      plot.title   = element_text(face = "bold", size = 13, hjust = 0.5),
      axis.text.y  = element_text(size = 9),
      strip.text   = element_text(face = "bold", size = 10)
    )
  
  ggsave(paste0("figures/GO_dotplot_", label, ".pdf"),
         p_dot, width = 10, height = 12)
  ggsave(paste0("figures/GO_dotplot_", label, ".png"),
         p_dot, width = 10, height = 12, dpi = 300)
  
  # --- Barplot（条形图）---
  p_bar <- barplot(ego, showCategory = top_n) +
    ggtitle(paste0("GO Barplot - ", label)) +
    theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5))
  
  ggsave(paste0("figures/GO_barplot_", label, ".pdf"),
         p_bar, width = 10, height = 10)
  
  # --- 网络图（GO term间关系）---
  p_emap <- emapplot(pairwise_termsim(ego), showCategory = top_n) +
    ggtitle(paste0("GO Network - ", label))
  
  ggsave(paste0("figures/GO_network_", label, ".pdf"),
         p_emap, width = 12, height = 10)
  
  cat("  图表已保存:", label, "\n")
}

plot_GO(go_508, "508_Shared_DEGs")
plot_GO(go_eu,  "Unique_europaeus")
plot_GO(go_ti,  "Unique_timidus")


# ─────────────────────────────────────────────
# SECTION 8: 可视化 — KEGG 通路图
# ─────────────────────────────────────────────

plot_KEGG <- function(kk, label, top_n = 20) {
  
  if (nrow(as.data.frame(kk)) == 0) {
    cat("  [警告]", label, "无显著KEGG通路，跳过作图\n")
    return(NULL)
  }
  
  # --- Dotplot ---
  p_dot <- dotplot(kk, showCategory = top_n) +
    ggtitle(paste0("KEGG Pathway Enrichment - ", label)) +
    theme(
      plot.title  = element_text(face = "bold", size = 13, hjust = 0.5),
      axis.text.y = element_text(size = 9)
    )
  
  ggsave(paste0("figures/KEGG_dotplot_", label, ".pdf"),
         p_dot, width = 10, height = 9)
  ggsave(paste0("figures/KEGG_dotplot_", label, ".png"),
         p_dot, width = 10, height = 9, dpi = 300)
  
  # --- Barplot ---
  p_bar <- barplot(kk, showCategory = top_n) +
    ggtitle(paste0("KEGG Barplot - ", label))
  
  ggsave(paste0("figures/KEGG_barplot_", label, ".pdf"),
         p_bar, width = 10, height = 9)
  
  cat("  图表已保存:", label, "\n")
}

plot_KEGG(kegg_508, "508_Shared_DEGs")
plot_KEGG(kegg_eu,  "Unique_europaeus")
plot_KEGG(kegg_ti,  "Unique_timidus")


# ─────────────────────────────────────────────
# SECTION 9: 三组对比图（Comparison Plot）
# 展示三组基因集的富集通路差异
# ─────────────────────────────────────────────

# 合并三组基因为命名列表（用于compareCluster）
gene_list_GO <- list(
  "Shared_508"        = entrez_508,
  "Unique_europaeus"  = entrez_eu,
  "Unique_timidus"    = entrez_ti
)

# GO 比较分析
cat("\n运行三组GO比较分析...\n")
ck_GO <- compareCluster(
  geneClusters  = gene_list_GO,
  fun           = "enrichGO",
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",           # 先做Biological Process
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)

p_compare_GO <- dotplot(ck_GO, showCategory = 15) +
  ggtitle("GO-BP Enrichment Comparison (3 Gene Sets)") +
  theme(
    plot.title   = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.text.x  = element_text(angle = 30, hjust = 1, size = 10),
    axis.text.y  = element_text(size = 9)
  )

ggsave("figures/GO_comparison_3groups.pdf",  p_compare_GO, width = 14, height = 12)
ggsave("figures/GO_comparison_3groups.png",  p_compare_GO, width = 14, height = 12, dpi = 300)

# KEGG 比较分析
cat("运行三组KEGG比较分析...\n")
ck_KEGG <- compareCluster(
  geneClusters  = gene_list_GO,
  fun           = "enrichKEGG",
  organism      = "mmu",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05
)
ck_KEGG <- setReadable(ck_KEGG, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")

p_compare_KEGG <- dotplot(ck_KEGG, showCategory = 15) +
  ggtitle("KEGG Pathway Enrichment Comparison (3 Gene Sets)") +
  theme(
    plot.title   = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.text.x  = element_text(angle = 30, hjust = 1, size = 10),
    axis.text.y  = element_text(size = 9)
  )

ggsave("figures/KEGG_comparison_3groups.pdf",  p_compare_KEGG, width = 14, height = 10)
ggsave("figures/KEGG_comparison_3groups.png",  p_compare_KEGG, width = 14, height = 10, dpi = 300)

cat("比较图已保存！\n")


# ─────────────────────────────────────────────
# SECTION 10: 导出结果表格（Excel格式）
# ─────────────────────────────────────────────

cat("\n导出结果表格...\n")

wb <- createWorkbook()

save_enrichment_table <- function(wb, result_obj, sheet_name) {
  df <- as.data.frame(result_obj)
  if (nrow(df) == 0) {
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, data.frame(Note = "No significant results"))
  } else {
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, df)
    # 冻结首行
    freezePane(wb, sheet_name, firstRow = TRUE)
    # 格式化首行
    addStyle(wb, sheet_name,
             style = createStyle(fontColour = "#FFFFFF", fgFill = "#4472C4",
                                 halign = "center", textDecoration = "bold"),
             rows = 1, cols = 1:ncol(df), gridExpand = TRUE)
    # 自动列宽
    setColWidths(wb, sheet_name, cols = 1:ncol(df), widths = "auto")
  }
}

# GO结果
save_enrichment_table(wb, go_508,  "GO_508_Shared")
save_enrichment_table(wb, go_eu,   "GO_Unique_europaeus")
save_enrichment_table(wb, go_ti,   "GO_Unique_timidus")

# KEGG结果
save_enrichment_table(wb, kegg_508,  "KEGG_508_Shared")
save_enrichment_table(wb, kegg_eu,   "KEGG_Unique_europaeus")
save_enrichment_table(wb, kegg_ti,   "KEGG_Unique_timidus")

saveWorkbook(wb, "tables/Enrichment_Results_All.xlsx", overwrite = TRUE)
cat("结果表格已保存至: tables/Enrichment_Results_All.xlsx\n")


# ─────────────────────────────────────────────
# SECTION 11: Session 信息（用于论文方法部分）
# ─────────────────────────────────────────────

cat("\n========== Session Info ==========\n")
sessionInfo()

# =============================================================================
# 输出文件清单 (Output Files):
#
# figures/
#   GO_dotplot_508_Shared_DEGs.pdf/.png      <- 508共享基因GO气泡图
#   GO_dotplot_Unique_europaeus.pdf/.png     <- europaeus特有基因GO气泡图
#   GO_dotplot_Unique_timidus.pdf/.png       <- timidus特有基因GO气泡图
#   GO_barplot_*.pdf                         <- 各组GO条形图
#   GO_network_*.pdf                         <- 各组GO网络图
#   KEGG_dotplot_*.pdf/.png                  <- 各组KEGG气泡图
#   KEGG_barplot_*.pdf                       <- 各组KEGG条形图
#   GO_comparison_3groups.pdf/.png           <- 三组GO对比图 ★ 论文主图
#   KEGG_comparison_3groups.pdf/.png         <- 三组KEGG对比图 ★ 论文主图
#
# tables/
#   Enrichment_Results_All.xlsx              <- 所有富集结果表（6个sheet）
# =============================================================================
