#!/bin/bash
set -euo pipefail

cd /scratch/project_2002674/RNAseq_hares/scripts/RNA-Seq_trial/RNA-Seq_PRJNA826339
module load samtools

echo "========== START $(date) =========="

# 就这么简单：8个BAM文件一起处理
bcftools mpileup \
  -Ou \
  -f ref/genome.fa \
  -q 20 -Q 20 \
  aln/SRR18740835.Aligned.sortedByCoord.out.bam \
  aln/SRR18740836.Aligned.sortedByCoord.out.bam \
  aln/SRR18740837.Aligned.sortedByCoord.out.bam \
  aln/SRR18740838.Aligned.sortedByCoord.out.bam \
  aln/SRR18740839.Aligned.sortedByCoord.out.bam \
  aln/SRR18740840.Aligned.sortedByCoord.out.bam \
  aln/SRR18740841.Aligned.sortedByCoord.out.bam \
  aln/SRR18740842.Aligned.sortedByCoord.out.bam | \
bcftools call -mv -Ou | \
bcftools view -Oz -o variants/all.raw.vcf.gz

echo "Indexing..."
bcftools index -f variants/all.raw.vcf.gz

echo "Filtering..."
bcftools filter \
  -i 'QUAL>100 && MQ>20 && INFO/DP>10' \
  -Oz -o variants/all.filtered.vcf.gz \
  variants/all.raw.vcf.gz

bcftools index -f variants/all.filtered.vcf.gz

echo "Raw: $(bcftools view -H variants/all.raw.vcf.gz | wc -l)"
echo "Filtered: $(bcftools view -H variants/all.filtered.vcf.gz | wc -l)"
echo "Files:"
ls -lh variants/*.vcf.gz*

echo "========== DONE $(date) =========="

