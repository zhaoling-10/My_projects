#######*************************** bash

#!/bin/bash
#SBATCH --account=project_2002674
#SBATCH --job-name=Rtf01_trim_qc
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x.%j.out
#SBATCH --error=logs/%x.%j.err

set -euo pipefail

# ---- Conda inside batch jobs ----
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate reindeer_core

# ---- Load project variables (OMNI_R1 / OMNI_R2 etc.) ----
source 00_meta/config.sh

# ---- Output folders ----
TRIM_DIR="${QC_DIR}/trimmed"
QC_TRIM_DIR="${QC_DIR}/fastqc_trimmed"
mkdir -p "${TRIM_DIR}" "${QC_TRIM_DIR}"

# ---- Output filenames ----
R1_TRIM="${TRIM_DIR}/Rtf01_R1.trim.fq.gz"
R2_TRIM="${TRIM_DIR}/Rtf01_R2.trim.fq.gz"
REPORT="${TRIM_DIR}/cutadapt_report.txt"

echo "=== Input reads ==="
ls -lh "${OMNI_R1}" "${OMNI_R2}"

echo "=== Trimming with cutadapt (adapters + polyG/polyA tails) ==="
# Notes:
# - -q 20,20  : trim low quality from both ends (Q<20)
# - --minimum-length 50 : discard too-short reads
# - -a/-A    : standard Illumina adapter trimming (works well in practice)
# - --nextseq-trim 20 : helps with NovaSeq/2-color chemistry polyG artifacts
# - --trim-n : remove Ns at ends
cutadapt \
  -j ${SLURM_CPUS_PER_TASK} \
  -q 20,20 \
  --minimum-length 50 \
  --trim-n \
  --nextseq-trim 20 \
  -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
  -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
  -o "${R1_TRIM}" \
  -p "${R2_TRIM}" \
  "${OMNI_R1}" "${OMNI_R2}" \
  > "${REPORT}"

echo "=== Trimming finished. Outputs: ==="
ls -lh "${R1_TRIM}" "${R2_TRIM}" "${REPORT}"

echo "=== FastQC on trimmed reads ==="
fastqc -t ${SLURM_CPUS_PER_TASK} -o "${QC_TRIM_DIR}" "${R1_TRIM}" "${R2_TRIM}"

echo "=== MultiQC summary (trimmed reads) ==="
multiqc -o "${QC_TRIM_DIR}" "${QC_TRIM_DIR}"

echo "DONE."