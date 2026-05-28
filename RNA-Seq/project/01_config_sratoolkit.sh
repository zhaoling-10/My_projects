#!/usr/bin/env bash
set -euo pipefail
source ./config.sh

module load sratoolkit/3.0.0

mkdir -p "$HOME/.ncbi"
mkdir -p "$HOME/ncbi/public/sra"

vdb-config --restore-defaults
vdb-config --set /repository/user/main/public/root="$HOME/ncbi/public"

echo "✅ SRA Toolkit configured."

# ---- Smoke test (recommended) ----
# Test ONLY the first SRR in list.txt (fast)
TEST_SRR=$(head -n 1 "$WORKDIR/list.txt")
echo "Testing download with: $TEST_SRR"
prefetch -v "$TEST_SRR" | tee "$WORKDIR/logs/01_prefetch_test.log"

echo "Test OK. Now run the full download with: bash 02_prefetch.sh"
