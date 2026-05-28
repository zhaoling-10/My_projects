#!/bin/bash
#SBATCH --job-name=variant
#SBATCH --account=project_2002674
#SBATCH --partition=small
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=256G
#SBATCH --output=logs/simple_%j.log

cd /scratch/project_2002674/RNAseq_hares/scripts/RNA-Seq_trial/RNA-Seq_PRJNA826339
module load samtools
bash simple_variant.sh

