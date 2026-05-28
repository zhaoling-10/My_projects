# 运行这段 R 代码，生成分组条形图 PDF

library(ggplot2)
library(dplyr)

# 数据
data <- data.frame(
  Category = c("Shared only", "LE only", "LT only", 
               "Shared+LE", "Shared+LT", "LE+LT", "All three"),
  Count = c(28, 9, 2, 1, 2, 4, 1)
)

# 设置因子顺序（从多到少，或从少到多）
data$Category <- factor(data$Category, 
                        levels = c("All three", "Shared+LT", "LE+LT", "Shared+LE", 
                                   "LT only", "LE only", "Shared only"))

# 颜色
colors <- c("Shared only" = "#1B9E77", 
            "LE only" = "#D95F02", 
            "LT only" = "#7570B3",
            "Shared+LE" = "#E7298A",
            "Shared+LT" = "#66A61E", 
            "LE+LT" = "#E6AB02",
            "All three" = "#A6761D")

# 绘图
p <- ggplot(data, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = Count), vjust = -0.5, size = 5) +
  scale_fill_manual(values = colors) +
  labs(
    title = "KEGG Pathway Overlap Across Gene Sets (mmu reference)",
    subtitle = "Number of pathways uniquely or jointly enriched in Shared / LE-unique / LT-unique",
    x = NULL,
    y = "Number of pathways",
    caption = "FDR ≤ 0.05"
  ) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, size = 11),
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "none"
  )

# 保存
ggsave("KEGG_pathway_overlap_barplot_presentation.pdf", p, width = 10, height = 6, dpi = 300)
ggsave("KEGG_pathway_overlap_barplot_presentation.png", p, width = 10, height = 6, dpi = 300)

print(p)