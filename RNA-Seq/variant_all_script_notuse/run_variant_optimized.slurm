#!/bin/bash
#SBATCH --job-name=vcf_optimized
#SBATCH --account=project_2002674
#SBATCH --partition=small
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=256G
#SBATCH --output=logs/vcf_optimized_%j.log

set -euo pipefail

echo "=========================================="
echo "Variant Calling - Optimized"
echo "Start: $(date)"
echo "=========================================="

cd /scratch/project_2002674/RNAseq_hares/scripts/RNA-Seq_trial/RNA-Seq_PRJNA826339

module load samtools
module load miniconda
conda activate rna-seq

mkdir -p variants logs

echo ""
echo "Running mpileup + call..."
echo "This may take 8-10 hours..."
echo ""

# 使用分染色体处理以加快速度
GENOME_FILE="ref/genome.fa"

# 获取所有染色体名
echo "Getting chromosome list..."
samtools view -H aln/SRR18740835.Aligned.sortedByCoord.out.bam | \
  grep "^@SQ" | \
  cut -f2 | sed 's/SN://' > chrom_list.txt

# 为每条染色体运行 mpileup
echo "Processing chromosomes..."
i=0
while read chr; do
  i=$((i+1))
  echo "[$i] Processing $chr..."
  
  bcftools mpileup \
    -Ou \
    -r $chr \
    -f $GENOME_FILE \
    -q 20 -Q 20 \
    -a INFO/DP,FORMAT/DP,FORMAT/AD \
    aln/*.Aligned.sortedByCoord.out.bam | \
  bcftools call -mv -Ou | \
  bcftools view -Oz -o variants/${chr}.raw.vcf.gz
  
  # 验证每条染色体的文件
  if ! gunzip -t variants/${chr}.raw.vcf.gz 2>/dev/null; then
    echo "ERROR: ${chr}.raw.vcf.gz is corrupted"
    exit 1
  fi
  echo "  ✓ ${chr} completed"
done < chrom_list.txt

echo ""
echo "Merging VCF files..."
# 合并所有染色体
bcftools concat -Oz -o variants/all.raw.vcf.gz variants/CM*.raw.vcf.gz variants/JBED*.raw.vcf.gz

# 验证合并的文件
if ! gunzip -t variants/all.raw.vcf.gz 2>/dev/null; then
  echo "ERROR: all.raw.vcf.gz is corrupted"
  exit 1
fi
echo "✓ Merged successfully"

echo ""
echo "Creating index..."
bcftools index -f variants/all.raw.vcf.gz
echo "✓ Index created"

echo ""
echo "Filtering variants..."
bcftools filter \
  -i 'QUAL>100 && MQ>20 && INFO/DP>10' \
  -Oz -o variants/all.filtered.vcf.gz \
  variants/all.raw.vcf.gz
echo "✓ Filtered"

echo ""
echo "Creating index for filtered VCF..."
bcftools index -f variants/all.filtered.vcf.gz
echo "✓ Index created"

echo ""
echo "=========================================="
echo "Summary"
echo "=========================================="
RAW=$(bcftools view -H variants/all.raw.vcf.gz | wc -l)
FILTERED=$(bcftools view -H variants/all.filtered.vcf.gz | wc -l)

echo "Raw variants: $RAW"
echo "Filtered variants (QUAL>100, MQ>20, DP>10): $FILTERED"
echo ""
echo "Final files:"
ls -lh variants/*.vcf.gz*

echo ""
echo "End: $(date)"
echo "=========================================="
echo "✓✓✓ Completed successfully!"

