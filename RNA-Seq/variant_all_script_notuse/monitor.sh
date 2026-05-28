#!/bin/bash
while true; do
  clear
  echo "=== $(date) ==="
  echo "【文件】"
  ls -lh /scratch/project_2002674/RNAseq_hares/scripts/RNA-Seq_trial/RNA-Seq_PRJNA826339/variants/ 2>/dev/null | tail -1 || echo "No files"
  echo ""
  echo "【进程】"
  ps aux | grep bcftools | grep -v grep | wc -l
  echo "bcftools processes"
  echo ""
  echo "【内存】"
  free -h | grep Mem
  echo ""
  echo "【任务】"
  squeue -u lingzhao | grep bash || echo "srun task not found"
  echo ""
  sleep 60
done
