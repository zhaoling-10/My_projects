#!/usr/bin/env bash

# 详细的日志记录脚本
# 每个样本独立处理，有完整的日志和断点恢复

set -euo pipefail

# ========== 配置 ==========
WORKDIR="/scratch/project_2002674/RNAseq_hares/scripts/RNA-Seq_trial/RNA-Seq_PRJNA826339"
LOG_DIR="${WORKDIR}/logs"
VARIANTS_DIR="${WORKDIR}/variants"
mkdir -p "$LOG_DIR" "$VARIANTS_DIR"

# 主日志文件
MAIN_LOG="${LOG_DIR}/variant_calling_per_sample_$(date +%Y%m%d_%H%M%S).log"

# 记录日志的函数
log_msg() {
  local level=$1
  shift
  local msg="$@"
  local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
  echo "[$timestamp] [$level] $msg" | tee -a "$MAIN_LOG"
}

# 错误处理函数
error_exit() {
  log_msg "ERROR" "$@"
  exit 1
}

# ========== 开始 ==========
log_msg "INFO" "=========================================="
log_msg "INFO" "Starting variant calling - per sample"
log_msg "INFO" "Timestamp: $(date '+%Y-%m-%d %H:%M:%S')"
log_msg "INFO" "=========================================="

# 加载模块
log_msg "INFO" "Loading modules..."
module load samtools 2>&1 | tee -a "$MAIN_LOG" || error_exit "Failed to load samtools"

# 验证前置文件
log_msg "INFO" "Verifying input files..."
test -f "${WORKDIR}/ref/genome.fa" || error_exit "Missing genome.fa"
log_msg "INFO" "✓ Reference genome found"

cd "$WORKDIR"
log_msg "INFO" "Changed to working directory: $WORKDIR"

# 获取所有 BAM 文件
log_msg "INFO" "Scanning BAM files..."
BAM_FILES=(aln/*.Aligned.sortedByCoord.out.bam)
log_msg "INFO" "Found ${#BAM_FILES[@]} BAM files"

# ========== 每个样本的处理 ==========
for BAM_FILE in "${BAM_FILES[@]}"; do
  # 提取样本名
  SAMPLE_NAME=$(basename "$BAM_FILE" .Aligned.sortedByCoord.out.bam)
  SAMPLE_LOG="${LOG_DIR}/variant_calling_${SAMPLE_NAME}.log"
  SAMPLE_VCF="${VARIANTS_DIR}/${SAMPLE_NAME}.raw.vcf.gz"
  
  # 检查是否已完成
  if [ -f "$SAMPLE_VCF" ]; then
    # 验证文件完整性
    if gunzip -t "$SAMPLE_VCF" 2>/dev/null > /dev/null; then
      log_msg "INFO" "Sample $SAMPLE_NAME already completed and valid"
      continue
    else
      log_msg "WARN" "Sample $SAMPLE_NAME VCF is corrupted, reprocessing..."
      rm -f "$SAMPLE_VCF" "${SAMPLE_VCF}.csi"
    fi
  fi
  
  log_msg "INFO" "=========================================="
  log_msg "INFO" "Processing sample: $SAMPLE_NAME"
  log_msg "INFO" "BAM file: $BAM_FILE"
  log_msg "INFO" "Output: $SAMPLE_VCF"
  log_msg "INFO" "=========================================="
  
  # 记录样本级别的详细日志
  (
    set -euo pipefail
    
    echo "$(date '+%Y-%m-%d %H:%M:%S') Starting mpileup + call for $SAMPLE_NAME"
    
    # 运行 mpileup + call
    bcftools mpileup \
      -Ou \
      -f "${WORKDIR}/ref/genome.fa" \
      -q 20 -Q 20 \
      -a INFO/DP,FORMAT/DP,FORMAT/AD \
      "$BAM_FILE" 2>&1 | \
    bcftools call -mv -Ou 2>&1 | \
    bcftools view -Oz -o "$SAMPLE_VCF" 2>&1
    
    echo "$(date '+%Y-%m-%d %H:%M:%S') mpileup + call completed for $SAMPLE_NAME"
    
    # 验证文件
    echo "$(date '+%Y-%m-%d %H:%M:%S') Verifying $SAMPLE_VCF..."
    if ! gunzip -t "$SAMPLE_VCF" 2>&1; then
      echo "ERROR: $SAMPLE_VCF is corrupted!"
      exit 1
    fi
    echo "$(date '+%Y-%m-%d %H:%M:%S') ✓ File verified successfully"
    
    # 建立索引
    echo "$(date '+%Y-%m-%d %H:%M:%S') Creating index for $SAMPLE_NAME..."
    bcftools index -f "$SAMPLE_VCF" 2>&1
    echo "$(date '+%Y-%m-%d %H:%M:%S') ✓ Index created"
    
    # 统计变异
    VARIANT_COUNT=$(bcftools view -H "$SAMPLE_VCF" 2>/dev/null | wc -l)
    echo "$(date '+%Y-%m-%d %H:%M:%S') Variant count for $SAMPLE_NAME: $VARIANT_COUNT"
    
  ) >> "$SAMPLE_LOG" 2>&1
  
  if [ $? -eq 0 ]; then
    log_msg "INFO" "✓ Sample $SAMPLE_NAME completed successfully"
  else
    log_msg "ERROR" "✗ Sample $SAMPLE_NAME failed"
    log_msg "ERROR" "Check log: $SAMPLE_LOG"
    # 不中断，继续处理其他样本
  fi
  
done

# ========== 合并所有 VCF ==========
log_msg "INFO" "=========================================="
log_msg "INFO" "Merging individual sample VCFs..."
log_msg "INFO" "=========================================="

# 检查是否所有样本都完成
COMPLETED_COUNT=$(ls "${VARIANTS_DIR}"/*.raw.vcf.gz 2>/dev/null | wc -l)
log_msg "INFO" "Completed samples: $COMPLETED_COUNT / ${#BAM_FILES[@]}"

if [ "$COMPLETED_COUNT" -eq 0 ]; then
  error_exit "No completed samples found!"
fi

# 创建样本列表文件
SAMPLE_LIST="${LOG_DIR}/sample_list.txt"
ls "${VARIANTS_DIR}"/*.raw.vcf.gz | sort > "$SAMPLE_LIST"

log_msg "INFO" "Samples to merge:"
cat "$SAMPLE_LIST" | tee -a "$MAIN_LOG"

# 合并 VCF
RAW_VCF="${VARIANTS_DIR}/all.raw.vcf.gz"
if [ ! -f "$RAW_VCF" ]; then
  (
    echo "$(date '+%Y-%m-%d %H:%M:%S') Starting VCF merge..."
    bcftools concat -Oz -o "$RAW_VCF" $(cat "$SAMPLE_LIST") 2>&1
    echo "$(date '+%Y-%m-%d %H:%M:%S') ✓ VCF merge completed"
    
    echo "$(date '+%Y-%m-%d %H:%M:%S') Verifying merged VCF..."
    if ! gunzip -t "$RAW_VCF" 2>&1; then
      echo "ERROR: Merged VCF is corrupted!"
      exit 1
    fi
    echo "$(date '+%Y-%m-%d %H:%M:%S') ✓ Merged VCF verified"
    
    echo "$(date '+%Y-%m-%d %H:%M:%S') Creating index..."
    bcftools index -f "$RAW_VCF" 2>&1
    echo "$(date '+%Y-%m-%d %H:%M:%S') ✓ Index created"
    
  ) >> "${LOG_DIR}/vcf_merge.log" 2>&1
  
  if [ $? -eq 0 ]; then
    log_msg "INFO" "✓ VCF merge completed successfully"
  else
    error_exit "VCF merge failed. Check ${LOG_DIR}/vcf_merge.log"
  fi
else
  log_msg "INFO" "Merged VCF already exists, skipping merge"
fi

# ========== 过滤 ==========
log_msg "INFO" "=========================================="
log_msg "INFO" "Filtering variants..."
log_msg "INFO" "=========================================="

FILTERED_VCF="${VARIANTS_DIR}/all.filtered.vcf.gz"
if [ ! -f "$FILTERED_VCF" ]; then
  (
    echo "$(date '+%Y-%m-%d %H:%M:%S') Starting filtering..."
    bcftools filter \
      -i 'QUAL>100 && MQ>20 && INFO/DP>10' \
      -Oz -o "$FILTERED_VCF" \
      "$RAW_VCF" 2>&1
    echo "$(date '+%Y-%m-%d %H:%M:%S') ✓ Filtering completed"
    
    echo "$(date '+%Y-%m-%d %H:%M:%S') Creating index..."
    bcftools index -f "$FILTERED_VCF" 2>&1
    echo "$(date '+%Y-%m-%d %H:%M:%S') ✓ Index created"
    
  ) >> "${LOG_DIR}/vcf_filter.log" 2>&1
  
  if [ $? -eq 0 ]; then
    log_msg "INFO" "✓ Filtering completed successfully"
  else
    error_exit "Filtering failed. Check ${LOG_DIR}/vcf_filter.log"
  fi
else
  log_msg "INFO" "Filtered VCF already exists, skipping filtering"
fi

# ========== 最终统计 ==========
log_msg "INFO" "=========================================="
log_msg "INFO" "FINAL SUMMARY"
log_msg "INFO" "=========================================="

RAW_COUNT=$(bcftools view -H "$RAW_VCF" 2>/dev/null | wc -l)
FILTERED_COUNT=$(bcftools view -H "$FILTERED_VCF" 2>/dev/null | wc -l)

log_msg "INFO" "Raw variants: $RAW_COUNT"
log_msg "INFO" "Filtered variants (QUAL>100, MQ>20, DP>10): $FILTERED_COUNT"
log_msg "INFO" ""
log_msg "INFO" "Final files:"
log_msg "INFO" "$(ls -lh ${VARIANTS_DIR}/*.vcf.gz* | tail -4)"

log_msg "INFO" ""
log_msg "INFO" "=========================================="
log_msg "INFO" "✓✓✓ ALL STEPS COMPLETED SUCCESSFULLY"
log_msg "INFO" "Completion time: $(date '+%Y-%m-%d %H:%M:%S')"
log_msg "INFO" "=========================================="

log_msg "INFO" ""
log_msg "INFO" "Log files:"
log_msg "INFO" "  Main log: $MAIN_LOG"
log_msg "INFO" "  Sample logs: ${LOG_DIR}/variant_calling_*.log"
log_msg "INFO" "  Merge log: ${LOG_DIR}/vcf_merge.log"
log_msg "INFO" "  Filter log: ${LOG_DIR}/vcf_filter.log"

