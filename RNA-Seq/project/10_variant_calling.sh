#!/usr/bin/env bash
set -euo pipefail
source ./config.sh

module load bcftools 2>/dev/null || true

cd "$WORKDIR"

bcftools mpileup \
  -Ou \
  -f ref/genome.fa \
  -q 20 -Q 20 \
  -a INFO/DP,FORMAT/DP,FORMAT/AD \
  aln/*.Aligned.sortedByCoord.out.bam \
| bcftools call -mv -Ou \
| bcftools view -Oz -o variants/all.raw.vcf.gz

bcftools index -f variants/all.raw.vcf.gz

bcftools filter \
  -i 'QUAL>100 && MQ>20 && INFO/DP>10' \
  -Oz -o variants/all.filtered.vcf.gz \
  variants/all.raw.vcf.gz

bcftools index -f variants/all.filtered.vcf.gz

echo -n "Filtered variant count: "
bcftools view -H variants/all.filtered.vcf.gz | wc -l