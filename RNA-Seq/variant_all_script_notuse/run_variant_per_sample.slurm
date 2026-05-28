#!/bin/bash
#SBATCH --job-name=variant_per_sample
#SBATCH --account=project_2002674
#SBATCH --partition=small
#SBATCH --time=36:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=256G
#SBATCH --output=logs/variant_per_sample_%j.log
#SBATCH --error=logs/variant_per_sample_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lingzhao@example.com

set -euo pipefail

WORKDIR="/scratch/project_2002674/RNAseq_hares/scripts/RNA-Seq_trial/RNA-Seq_PRJNA826339"
cd "$WORKDIR"

echo "=========================================="
echo "SBATCH Job started at $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $HOSTNAME"
echo "=========================================="

# 加载模块
module load samtools 2>&1 || exit 1

# 运行脚本
bash variant_calling_per_sample.sh 2>&1

EXIT_CODE=$?

echo ""
echo "=========================================="
echo "SBATCH Job finished at $(date)"
echo "Exit code: $EXIT_CODE"
echo "=========================================="

exit $EXIT_CODE

