#!/usr/bin/env bash
set -euo pipefail
source ./config.sh

module load multiqc 2>/dev/null || true
cd "$WORKDIR"

multiqc -o qc qc aln logs 2>&1 | tee logs/11_multiqc.log
echo "MultiQC report in qc/"