#!/usr/bin/env bash
set -euo pipefail
source ./config.sh

module load star 2>/dev/null || true
module load samtools 2>/dev/null || true

test -f "$WORKDIR/ref/genome.fa" || { echo "Missing ref/genome.fa (run 04_prepare_reference.sh)"; exit 1; }
test -f "$WORKDIR/ref/genes.gtf" || { echo "Missing ref/genes.gtf (run 04_prepare_reference.sh)"; exit 1; }

cd "$WORKDIR"

# STAR index is created once
if [ -f "$WORKDIR/star_index/Genome" ]; then
  echo "STAR index already exists. Skipping."
  exit 0
fi

STAR \
  --runMode genomeGenerate \
  --runThreadN "$THREADS" \
  --genomeDir "$WORKDIR/star_index" \
  --genomeFastaFiles "$WORKDIR/ref/genome.fa" \
  --sjdbGTFfile "$WORKDIR/ref/genes.gtf" \
  --sjdbOverhang "$SJDB_OVERHANG" \
  2>&1 | tee "$WORKDIR/logs/06_star_index.log"

echo "STAR index done."