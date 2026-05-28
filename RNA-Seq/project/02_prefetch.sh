#!/usr/bin/env bash
set -euo pipefail
source ./config.sh

module load sratoolkit/3.0.0

cd "$WORKDIR/sra"
prefetch --option-file "$WORKDIR/list.txt" --output-directory . 2>&1 | tee "$WORKDIR/logs/02_prefetch.log"

echo "prefetch done. Check $WORKDIR/sra"