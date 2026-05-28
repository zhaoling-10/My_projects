#!/usr/bin/env bash
set -euo pipefail
source ./config.sh

module load samtools 2>/dev/null || true

test -f "$WORKDIR/ref_src/$REF_FA_GZ" || { echo "Missing $WORKDIR/ref_src/$REF_FA_GZ"; exit 1; }
test -f "$WORKDIR/ref_src/$REF_GTF_GZ" || { echo "Missing $WORKDIR/ref_src/$REF_GTF_GZ"; exit 1; }

echo "Decompressing reference to $WORKDIR/ref/ ..."
gunzip -c "$WORKDIR/ref_src/$REF_FA_GZ" > "$WORKDIR/ref/genome.fa"
gunzip -c "$WORKDIR/ref_src/$REF_GTF_GZ" > "$WORKDIR/ref/genes.gtf"

echo "Indexing FASTA..."
samtools faidx "$WORKDIR/ref/genome.fa"

echo "Reference ready:"
ls -lh "$WORKDIR/ref/genome.fa" "$WORKDIR/ref/genes.gtf" "$WORKDIR/ref/genome.fa.fai"