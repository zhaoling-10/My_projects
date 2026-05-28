#!/usr/bin/env bash
set -euo pipefail
source ./config.sh

module load fastp 2>/dev/null || true

cd "$WORKDIR"

while read -r srr; do
  echo "==> fastp $srr"

  IN1="fastq/${srr}_1.fastq.gz"
  IN2="fastq/${srr}_2.fastq.gz"
  OUT1="trimmed/${srr}_1.fastq.gz"
  OUT2="trimmed/${srr}_2.fastq.gz"

  test -f "$IN1" || { echo "Missing $IN1"; exit 1; }
  test -f "$IN2" || { echo "Missing $IN2"; exit 1; }

  # Skip if already trimmed
  if [[ -f "$OUT1" && -f "$OUT2" ]]; then
    echo "  Trimmed exists, skipping: $srr"
    continue
  fi

  fastp \
    --in1 "$IN1" \
    --in2 "$IN2" \
    --out1 "$OUT1" \
    --out2 "$OUT2" \
    --detect_adapter_for_pe \
    --length_required 50 \
    --thread "$THREADS" \
    --html "qc/${srr}.fastp.html" \
    --json "qc/${srr}.fastp.json" \
    2>&1 | tee -a "logs/05_fastp.log"
done < list.txt

echo "fastp done. Reports in qc/"