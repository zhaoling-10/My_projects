#!/bin/bash
#SBATCH --job-name=variant_fixed
#SBATCH --account=project_2002674
#SBATCH --partition=small
#SBATCH --time=36:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=256G
#SBATCH --output=logs/variant_fixed_%j.log
#SBATCH --error=logs/variant_fixed_%j.err

set -euo pipefail

WORKDIR="/scratch/project_2002674/RNAseq_hares/scripts/RNA-Seq_trial/RNA-Seq_PRJNA826339"
cd "$WORKDIR"

echo "=========================================="
echo "Job started: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "=========================================="

module load samtools 2>&1

# 删除旧的失败文件
rm -f variants/*.raw.vcf.gz* variants/.*.pileup variants/.*.called

bash variant_calling_per_sample_FIXED.sh

echo ""
echo "=========================================="
echo "Job finished: $(date)"
echo "=========================================="

