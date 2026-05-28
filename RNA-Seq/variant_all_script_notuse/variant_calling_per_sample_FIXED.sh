#!/usr/bin/env bash

set -euo pipefail

WORKDIR="/scratch/project_2002674/RNAseq_hares/scripts/RNA-Seq_trial/RNA-Seq_PRJNA826339"
LOG_DIR="${WORKDIR}/logs"
VARIANTS_DIR="${WORKDIR}/variants"
mkdir -p "$LOG_DIR" "$VARIANTS_DIR"

MAIN_LOG="${LOG_DIR}/variant_calling_fixed_$(date +%Y%m%d_%H%M%S).log"

log_msg() {
  local level=$1
  shift
  local msg="$@"
  local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
  echo "[$timestamp] [$level] $msg" | tee -a "$MAIN_LOG"
}

error_exit() {
  log_msg "ERROR" "$@"
  exit 1
}

log_msg "INFO" "========== START =========="
log_msg "INFO" "Time: $(date)"

module load samtools 2>&1 || error_exit "Failed to load samtools"
log_msg "INFO" "✓ Modules loaded"

cd "$WORKDIR"

# 处理每个样本
for BAM_FILE in aln/*.Aligned.sortedByCoord.out.bam; do
  SAMPLE_NAME=$(basename "$BAM_FILE" .Aligned.sortedByCoord.out.bam)
  SAMPLE_LOG="${LOG_DIR}/variant_calling_${SAMPLE_NAME}.log"
  SAMPLE_VCF="${VARIANTS_DIR}/${SAMPLE_NAME}.raw.vcf.gz"
  
  # 中间文件
  PILEUP_FILE="${VARIANTS_DIR}/.${SAMPLE_NAME}.pileup"
  CALLED_FILE="${VARIANTS_DIR}/.${SAMPLE_NAME}.called"
  
  # 检查是否已完成
  if [ -f "$SAMPLE_VCF" ] && gunzip -t "$SAMPLE_VCF" 2>/dev/null > /dev/null; then
    log_msg "INFO" "✓ Sample $SAMPLE_NAME already completed"
    continue
  fi
  
  rm -f "$SAMPLE_VCF" "${SAMPLE_VCF}.csi" "$PILEUP_FILE" "$CALLED_FILE"
  
  log_msg "INFO" "=========================================="
  log_msg "INFO" "Processing sample: $SAMPLE_NAME"
  log_msg "INFO" "BAM file: $BAM_FILE"
  log_msg "INFO" "=========================================="
  
  (
    set -euo pipefail
    
    echo "$(date '+%Y-%m-%d %H:%M:%S') [INFO] START: $SAMPLE_NAME"
    
    # 第1步：运行 mpileup（保存为中间文件）
    echo "$(date '+%Y-%m-%d %H:%M:%S') [INFO] Running mpileup for $SAMPLE_NAME..."
    if ! bcftools mpileup \
      -Ou \
      -f "${WORKDIR}/ref/genome.fa" \
      -q 20 -Q 20 \
      -a INFO/DP,FORMAT/DP,FORMAT/AD \
      "$BAM_FILE" \
      > "$PILEUP_FILE" 2>&1; then
      echo "$(date '+%Y-%m-%d %H:%M:%S') [ERROR] mpileup failed for $SAMPLE_NAME"
      echo "$(date '+%Y-%m-%d %H:%M:%S') [ERROR] BAM file: $BAM_FILE"
      exit 1
    fi
    
    if [ ! -f "$PILEUP_FILE" ] || [ ! -s "$PILEUP_FILE" ]; then
      echo "$(date '+%Y-%m-%d %H:%M:%S') [ERROR] mpileup output is empty"
      exit 1
    fi
    echo "$(date '+%Y-%m-%d %H:%M:%S') [INFO] ✓ mpileup completed, pileup size: $(du -h $PILEUP_FILE | cut -f1)"
    
    # 第2步：运行 call（从中间文件读取）
    echo "$(date '+%Y-%m-%d %H:%M:%S') [INFO] Running call for $SAMPLE_NAME..."
    if ! bcftools call \
      -mv \
      -Ou \
      "$PILEUP_FILE" \
      > "$CALLED_FILE" 2>&1; then
      echo "$(date '+%Y-%m-%d %H:%M:%S') [ERROR] call failed for $SAMPLE_NAME"
      exit 1
    fi
    
    if [ ! -f "$CALLED_FILE" ] || [ ! -s "$CALLED_FILE" ]; then
      echo "$(date '+%Y-%m-%d %H:%M:%S') [ERROR] call output is empty"
      exit 1
    fi
    echo "$(date '+%Y-%m-%d %H:%M:%S') [INFO] ✓ call completed"
    
    # 第3步：转换为压缩 VCF（从中间文件读取）
    echo "$(date '+%Y-%m-%d %H:%M:%S') [INFO] Converting to VCF for $SAMPLE_NAME..."
    if ! bcftools view \
      -Oz \
      "$CALLED_FILE" \
      > "$SAMPLE_VCF" 2>&1; then
      echo "$(date '+%Y-%m-%d %H:%M:%S') [ERROR] view failed for $SAMPLE_NAME"
      exit 1
    fi
    
    if [ ! -f "$SAMPLE_VCF" ] || [ ! -s "$SAMPLE_VCF" ]; then
      echo "$(date '+%Y-%m-%d %H:%M:%S') [ERROR] VCF file is empty"
      exit 1
    fi
    echo "$(date '+%Y-%m-%d %H:%M:%S') [INFO] ✓ VCF conversion completed"
    
    # 验证
    echo "$(date '+%Y-%m-%d %H:%M:%S') [INFO] Verifying VCF..."
    if ! gunzip -t "$SAMPLE_VCF" 2>&1; then
      echo "$(date '+%Y-%m-%d %H:%M:%S') [ERROR] VCF is corrupted"
      exit 1
    fi
    echo "$(date '+%Y-%m-%d %H:%M:%S') [INFO] ✓ VCF verified"
    
    # 索引
    echo "$(date '+%Y-%m-%d %H:%M:%S') [INFO] Creating index..."
    bcftools index -f "$SAMPLE_VCF"
    echo "$(date '+%Y-%m-%d %H:%M:%S') [INFO] ✓ Index created"
    
    # 统计
    COUNT=$(bcftools view -H "$SAMPLE_VCF" 2>/dev/null | wc -l)
    echo "$(date '+%Y-%m-%d %H:%M:%S') [INFO] Variant count: $COUNT"
    
    # 清理中间文件
    echo "$(date '+%Y-%m-%d %H:%M:%S') [INFO] Cleaning up intermediate files..."
    rm -f "$PILEUP_FILE" "$CALLED_FILE"
    
    echo "$(date '+%Y-%m-%d %H:%M:%S') [INFO] DONE: $SAMPLE_NAME (variants: $COUNT)"
    
  ) >> "$SAMPLE_LOG" 2>&1
  
  if [ $? -eq 0 ]; then
    log_msg "INFO" "✓ $SAMPLE_NAME completed successfully"
  else
    log_msg "ERROR" "✗ $SAMPLE_NAME failed - check $SAMPLE_LOG"
    log_msg "ERROR" "Sample log content:"
    cat "$SAMPLE_LOG" >> "$MAIN_LOG"
  fi
  
done

# 合并 VCF
log_msg "INFO" "=========================================="
log_msg "INFO" "Merging VCF files..."

COMPLETED_COUNT=$(ls "${VARIANTS_DIR}"/*.raw.vcf.gz 2>/dev/null | wc -l)
log_msg "INFO" "Completed samples: $COMPLETED_COUNT"

if [ "$COMPLETED_COUNT" -eq 0 ]; then
  error_exit "No completed samples found!"
fi

RAW_VCF="${VARIANTS_DIR}/all.raw.vcf.gz"
if [ ! -f "$RAW_VCF" ]; then
  (
    echo "$(date '+%Y-%m-%d %H:%M:%S') [INFO] Starting VCF merge..."
    bcftools concat -Oz -o "$RAW_VCF" $(ls "${VARIANTS_DIR}"/*.raw.vcf.gz | grep -v all.raw | sort)
    echo "$(date '+%Y-%m-%d %H:%M:%S') [INFO] Verifying..."
    gunzip -t "$RAW_VCF" || exit 1
    echo "$(date '+%Y-%m-%d %H:%M:%S') [INFO] Creating index..."
    bcftools index -f "$RAW_VCF"
    echo "$(date '+%Y-%m-%d %H:%M:%S') [INFO] ✓ Merge completed"
  ) >> "${LOG_DIR}/vcf_merge.log" 2>&1
  
  if [ $? -eq 0 ]; then
    log_msg "INFO" "✓ VCF merge completed"
  else
    error_exit "VCF merge failed. Check ${LOG_DIR}/vcf_merge.log"
  fi
fi

# 过滤
log_msg "INFO" "=========================================="
log_msg "INFO" "Filtering variants..."

FILTERED_VCF="${VARIANTS_DIR}/all.filtered.vcf.gz"
if [ ! -f "$FILTERED_VCF" ]; then
  (
    echo "$(date '+%Y-%m-%d %H:%M:%S') [INFO] Starting filter..."
    bcftools filter \
      -i 'QUAL>100 && MQ>20 && INFO/DP>10' \
      -Oz -o "$FILTERED_VCF" \
      "$RAW_VCF"
    echo "$(date '+%Y-%m-%d %H:%M:%S') [INFO] Creating index..."
    bcftools index -f "$FILTERED_VCF"
    echo "$(date '+%Y-%m-%d %H:%M:%S') [INFO] ✓ Filter completed"
  ) >> "${LOG_DIR}/vcf_filter.log" 2>&1
  
  if [ $? -eq 0 ]; then
    log_msg "INFO" "✓ Filtering completed"
  else
    error_exit "Filtering failed"
  fi
fi

# 最终统计
log_msg "INFO" "=========================================="
log_msg "INFO" "FINAL RESULTS"
log_msg "INFO" "=========================================="

RAW_COUNT=$(bcftools view -H "$RAW_VCF" 2>/dev/null | wc -l)
FILTERED_COUNT=$(bcftools view -H "$FILTERED_VCF" 2>/dev/null | wc -l)

log_msg "INFO" "Raw variants: $RAW_COUNT"
log_msg "INFO" "Filtered variants: $FILTERED_COUNT"
log_msg "INFO" ""
log_msg "INFO" "Final files:"
ls -lh "${VARIANTS_DIR}"/*.vcf.gz* >> "$MAIN_LOG" 2>&1

log_msg "INFO" ""
log_msg "INFO" "=========================================="
log_msg "INFO" "✓✓✓ ALL STEPS COMPLETED SUCCESSFULLY"
log_msg "INFO" "=========================================="

