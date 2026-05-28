#!/usr/bin/env bash
set -euo pipefail
source ./config.sh

module load subread 2>/dev/null || true

cd "$WORKDIR"

featureCounts \
  -T "$THREADS" \
  -p \
  -s 0 \
  -a ref/genes.gtf \
  -o counts/gene_counts_featureCounts.txt \
  aln/*.Aligned.sortedByCoord.out.bam \
  2>&1 | tee logs/09_featurecounts.log

echo "featureCounts done: counts/gene_counts_featureCounts.txt"
echo "NOTE: -s 0 assumes unstranded library. If stranded, use -s 1 or -s 2."
