#!/usr/bin/env bash
set -euo pipefail
source ./config.sh

module load samtools 2>/dev/null || true

cd "$WORKDIR"

# flagstat
while read -r srr; do
  BAM="aln/${srr}.Aligned.sortedByCoord.out.bam"
  test -f "$BAM" || { echo "Missing $BAM"; exit 1; }
  samtools flagstat "$BAM" > "qc/${srr}.flagstat.txt"
done < list.txt

# STAR mapping summary
echo -e "sample\tuniquely_mapped_pct" > qc/star_mapping_summary.tsv
while read -r srr; do
  LOG="aln/${srr}.Log.final.out"
  test -f "$LOG" || { echo "Missing $LOG"; exit 1; }
  uniq=$(grep "Uniquely mapped reads %" "$LOG" | awk -F'|\t' '{print $2}' | xargs)
  echo -e "${srr}\t${uniq}" >> qc/star_mapping_summary.tsv
done < list.txt

echo "Post-alignment QC done."