#!/usr/bin/env bash
set -euo pipefail
source ./config.sh

module load sratoolkit/3.0.0

cd "$WORKDIR"

while read -r acc; do
  echo "==> fasterq-dump $acc"

  # Skip if already converted
  if [[ -f "fastq/${acc}_1.fastq.gz" && -f "fastq/${acc}_2.fastq.gz" ]]; then
    echo "  FASTQ gz exists, skipping: $acc"
    continue
  fi

  fasterq-dump "$acc" --split-files -e "$THREADS" -O "$WORKDIR/fastq" 2>&1 | tee -a "$WORKDIR/logs/03_fasterq_dump.log"

  # Compress (pigz if available, else gzip)
  if command -v pigz >/dev/null 2>&1; then
    pigz -p "$THREADS" -f "$WORKDIR/fastq/${acc}_1.fastq" "$WORKDIR/fastq/${acc}_2.fastq"
  else
    gzip -f "$WORKDIR/fastq/${acc}_1.fastq" "$WORKDIR/fastq/${acc}_2.fastq"
  fi
done < "$WORKDIR/list.txt"

echo "FASTQ conversion done."