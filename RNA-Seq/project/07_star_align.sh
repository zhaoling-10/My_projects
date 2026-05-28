#!/usr/bin/env bash
set -euo pipefail
source ./config.sh

module load star 2>/dev/null || true
module load samtools 2>/dev/null || true

cd "$WORKDIR"

while read -r srr; do
  echo "==> STAR align $srr"

  R1="trimmed/${srr}_1.fastq.gz"
  R2="trimmed/${srr}_2.fastq.gz"
  test -f "$R1" || { echo "Missing $R1"; exit 1; }
  test -f "$R2" || { echo "Missing $R2"; exit 1; }

  BAM="aln/${srr}.Aligned.sortedByCoord.out.bam"
  if [ ! -f "$BAM" ]; then
    STAR \
      --runThreadN "$THREADS" \
      --genomeDir "$WORKDIR/star_index" \
      --readFilesIn "$R1" "$R2" \
      --readFilesCommand zcat \
      --outFileNamePrefix "aln/${srr}." \
      --outSAMtype BAM SortedByCoordinate \
      2>&1 | tee -a "logs/07_star_align.log"
  else
    echo "  BAM exists, skipping: $srr"
  fi

  if [ ! -f "${BAM}.bai" ]; then
    samtools index "$BAM"
  fi
done < list.txt

echo "STAR alignment done."