# ---- read + keep Geneid and sample columns ----

file <- "gene_counts_featureCounts_DFS_LE.txt"

# Option A (recommended): data.table::fread is very robust & fast
library(data.table)

dt <- fread(
  file,
  sep = "\t",
  header = TRUE,
  data.table = FALSE,
  quote = ""   # avoids issues if there are unexpected quotes
)

# Check columns read correctly
print(names(dt))

# Keep Geneid + all columns that look like your sample BAM paths
sample_cols <- grep("\\.bam$", names(dt), value = TRUE)   # columns ending with .bam
# If your file might not end with .bam, use: grep("^aln/", names(dt), value = TRUE)

counts <- dt[, c("Geneid", sample_cols), drop = FALSE]

# Make Geneid rownames (optional; many tools like DESeq2 accept rownames)
counts_mat <- as.data.frame(counts)
rownames(counts_mat) <- counts_mat$Geneid
counts_mat$Geneid <- NULL

# Convert to numeric matrix (important for DESeq2/edgeR)
counts_mat[] <- lapply(counts_mat, function(x) as.numeric(x))

# Optional: simplify sample names
colnames(counts_mat) <- sub("^aln/", "", colnames(counts_mat))
colnames(counts_mat) <- sub("\\.Aligned\\.sortedByCoord\\.out\\.bam$", "", colnames(counts_mat))

# Quick sanity checks
dim(counts_mat)
head(counts_mat[, 1:3, drop = FALSE])

# Save cleaned counts table if you want
write.table(
  cbind(Geneid = rownames(counts_mat), counts_mat),
  file = "counts_only.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
