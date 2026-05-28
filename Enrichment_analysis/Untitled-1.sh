#!/usr/bin/env bash
# =============================================================================
# KOBAS Output Folder Comparison Script
# Compares your files vs colleague's files for both GO and KEGG results
# =============================================================================

set -euo pipefail

# ---- Colour codes ----
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
BOLD='\033[1m'
NC='\033[0m' # No Colour

# ---- Your folders ----
MY_GO="/home/zhaoling/Documents/Enrichment_analysis/KOBAS_output_mmu_GO_1"
MY_KEGG="/home/zhaoling/Documents/Enrichment_analysis/KOBAS_output_mmu_KEGG_1"

# ---- Colleague's folders ----
COL_GO="/home/zhaoling/Documents/Data_analysi_for_Ling/9.Enrichment analysis/9.2.Enrichment_analysis/2. GO enrichment: Mus musculus/1.int_ref_Genes508_adj"
COL_KEGG="/home/zhaoling/Documents/Data_analysi_for_Ling/9.Enrichment analysis/9.2.Enrichment_analysis/1. KEGG enrichment: Mus musculus/1.int_ref_Genes508_adj"

# ---- Report output ----
REPORT="kobas_comparison_report_$(date +%Y%m%d_%H%M%S).txt"

# ---- Counters ----
TOTAL_COMPARED=0
TOTAL_IDENTICAL=0
TOTAL_DIFFERENT=0
TOTAL_MISSING=0

# =============================================================================
# Helper functions
# =============================================================================

log()  { echo -e "$*" | tee -a "$REPORT"; }
hr()   { log "${BLUE}$(printf '=%.0s' {1..70})${NC}"; }
hr2()  { log "${BLUE}$(printf -- '-%.0s' {1..70})${NC}"; }

check_folder_exists() {
    local folder="$1"
    local label="$2"
    if [[ ! -d "$folder" ]]; then
        log "${RED}ERROR: Folder not found: ${label}${NC}"
        log "       Path: $folder"
        return 1
    fi
    log "${GREEN}OK${NC}  Found: $label"
    log "     Path: $folder"
    return 0
}

compare_two_files() {
    local file1="$1"   # your file
    local file2="$2"   # colleague's file
    local label="$3"   # display name

    TOTAL_COMPARED=$((TOTAL_COMPARED + 1))

    # --- Check existence ---
    local missing1=0 missing2=0
    [[ ! -f "$file1" ]] && missing1=1
    [[ ! -f "$file2" ]] && missing2=1

    if (( missing1 && missing2 )); then
        log "  ${YELLOW}MISSING BOTH${NC}  $label"
        TOTAL_MISSING=$((TOTAL_MISSING + 1))
        return
    elif (( missing1 )); then
        log "  ${YELLOW}MISSING YOURS ${NC} $label"
        log "              Your path:      $file1"
        TOTAL_MISSING=$((TOTAL_MISSING + 1))
        return
    elif (( missing2 )); then
        log "  ${YELLOW}MISSING THEIRS${NC} $label"
        log "              Colleague path: $file2"
        TOTAL_MISSING=$((TOTAL_MISSING + 1))
        return
    fi

    # --- File sizes ---
    local size1 size2
    size1=$(wc -c < "$file1")
    size2=$(wc -c < "$file2")

    # --- MD5 hash comparison (fast, definitive) ---
    local md5_1 md5_2
    md5_1=$(md5sum "$file1" 2>/dev/null | awk '{print $1}' || \
            md5      "$file1" 2>/dev/null | awk '{print $NF}')
    md5_2=$(md5sum "$file2" 2>/dev/null | awk '{print $1}' || \
            md5      "$file2" 2>/dev/null | awk '{print $NF}')

    if [[ "$md5_1" == "$md5_2" ]]; then
        log "  ${GREEN}IDENTICAL${NC}     $label  (${size1} bytes, MD5: ${md5_1:0:12}...)"
        TOTAL_IDENTICAL=$((TOTAL_IDENTICAL + 1))
    else
        TOTAL_DIFFERENT=$((TOTAL_DIFFERENT + 1))
        log "  ${RED}DIFFERENT${NC}     $label"
        log "              Your size:       ${size1} bytes  | MD5: $md5_1"
        log "              Colleague size:  ${size2} bytes  | MD5: $md5_2"

        # --- Line count diff ---
        local lines1 lines2
        lines1=$(wc -l < "$file1")
        lines2=$(wc -l < "$file2")
        log "              Your lines:      $lines1"
        log "              Colleague lines: $lines2"

        # --- Show first differing line (skip comment lines for clarity) ---
        local diff_line
        diff_line=$(diff --unified=0 \
                        <(grep -v '^#' "$file1" 2>/dev/null || cat "$file1") \
                        <(grep -v '^#' "$file2" 2>/dev/null || cat "$file2") \
                    | head -20 2>/dev/null || true)
        if [[ -n "$diff_line" ]]; then
            log "              First differences (- yours / + colleague, comments stripped):"
            while IFS= read -r dline; do
                log "                $dline"
            done <<< "$diff_line"
        fi
    fi
}

compare_folder_pair() {
    local my_folder="$1"
    local col_folder="$2"
    local analysis_type="$3"   # "GO" or "KEGG"
    local set_label="$4"       # e.g. "Shared_508"

    hr
    log "${BOLD}Comparing: ${analysis_type} — ${set_label}${NC}"
    log "  Your folder:      $my_folder"
    log "  Colleague folder: $col_folder"
    hr2

    # Determine prefix from analysis type
    # Your files use KOBAS_output_mmu_{TYPE}_1 prefix
    # Colleague files use KOBAS_output_mmu_{TYPE}_1 prefix too — we detect dynamically

    # List ALL files in your folder and compare each
    local my_files=()
    while IFS= read -r -d '' f; do
        my_files+=("$f")
    done < <(find "$my_folder" -maxdepth 1 -type f -print0 2>/dev/null | sort -z)

    if (( ${#my_files[@]} == 0 )); then
        log "  ${YELLOW}WARNING: Your folder is empty or not accessible${NC}"
        return
    fi

    # For each of your files, find the matching colleague file by basename
    for my_file in "${my_files[@]}"; do
        local basename
        basename=$(basename "$my_file")
        local col_file="${col_folder}/${basename}"
        compare_two_files "$my_file" "$col_file" "$basename"
    done

    # Also check if colleague has extra files you don't have
    log ""
    log "  Checking for extra files in colleague's folder..."
    local col_files=()
    while IFS= read -r -d '' f; do
        col_files+=("$f")
    done < <(find "$col_folder" -maxdepth 1 -type f -print0 2>/dev/null | sort -z)

    local extra_count=0
    for col_file in "${col_files[@]}"; do
        local bn
        bn=$(basename "$col_file")
        if [[ ! -f "${my_folder}/${bn}" ]]; then
            log "  ${YELLOW}EXTRA IN COLLEAGUE${NC}  $bn (not in your folder)"
            extra_count=$((extra_count + 1))
        fi
    done
    if (( extra_count == 0 )); then
        log "  ${GREEN}No extra files in colleague's folder${NC}"
    fi
}

# =============================================================================
# Content-level checks for key KOBAS files
# =============================================================================

content_check_identify() {
    local file1="$1"
    local file2="$2"
    local label="$3"

    if [[ ! -f "$file1" || ! -f "$file2" ]]; then return; fi

    log ""
    log "  ${BOLD}Content check: $label${NC}"

    # Count significant terms (non-comment, non-empty lines)
    local terms1 terms2
    terms1=$(grep -v '^#' "$file1" | grep -v '^$' | wc -l || echo 0)
    terms2=$(grep -v '^#' "$file2" | grep -v '^$' | wc -l || echo 0)
    log "    Enriched terms — Yours: $terms1 | Colleague: $terms2"

    if [[ "$terms1" != "$terms2" ]]; then
        log "    ${YELLOW}Term count DIFFERS${NC}"

        # Show terms in yours but not colleague's
        local only_yours only_col
        only_yours=$(comm -23 \
            <(grep -v '^#' "$file1" | awk -F'\t' '{print $1}' | sort) \
            <(grep -v '^#' "$file2" | awk -F'\t' '{print $1}' | sort) \
            2>/dev/null | head -10)
        only_col=$(comm -13 \
            <(grep -v '^#' "$file1" | awk -F'\t' '{print $1}' | sort) \
            <(grep -v '^#' "$file2" | awk -F'\t' '{print $1}' | sort) \
            2>/dev/null | head -10)

        if [[ -n "$only_yours" ]]; then
            log "    Terms ONLY in yours:"
            while IFS= read -r t; do log "      $t"; done <<< "$only_yours"
        fi
        if [[ -n "$only_col" ]]; then
            log "    Terms ONLY in colleague's:"
            while IFS= read -r t; do log "      $t"; done <<< "$only_col"
        fi
    else
        log "    ${GREEN}Term counts match${NC}"
    fi

    # Check top enriched term
    local top1 top2
    top1=$(grep -v '^#' "$file1" | grep -v '^$' | head -1 | awk -F'\t' '{print $1}' || true)
    top2=$(grep -v '^#' "$file2" | grep -v '^$' | head -1 | awk -F'\t' '{print $1}' || true)
    log "    Top term — Yours: '$top1'"
    log "               Colleague: '$top2'"
    if [[ "$top1" == "$top2" ]]; then
        log "    ${GREEN}Top term matches${NC}"
    else
        log "    ${YELLOW}Top term DIFFERS${NC}"
    fi
}

# =============================================================================
# MAIN
# =============================================================================

{
log ""
hr
log "${BOLD}  KOBAS OUTPUT COMPARISON REPORT${NC}"
log "  Generated: $(date)"
log "  Script:    compare_kobas_folders.sh"
hr

# --- Pre-flight: check all four folders exist ---
log ""
log "${BOLD}[0] FOLDER EXISTENCE CHECK${NC}"
log ""
all_exist=1
check_folder_exists "$MY_GO"   "Your GO folder"       || all_exist=0
check_folder_exists "$MY_KEGG" "Your KEGG folder"     || all_exist=0
check_folder_exists "$COL_GO"  "Colleague GO folder"  || all_exist=0
check_folder_exists "$COL_KEGG" "Colleague KEGG folder" || all_exist=0

if (( ! all_exist )); then
    log ""
    log "${RED}One or more folders not found. Please check paths and re-run.${NC}"
    log "Tip: Folder names with spaces need no escaping in this script —"
    log "     paths are stored as variables and passed correctly."
fi

# --- File-by-file comparison ---
log ""
log "${BOLD}[1] FILE-BY-FILE COMPARISON (GO analysis)${NC}"
compare_folder_pair "$MY_GO" "$COL_GO" "GO" "Shared_508 (mmu)"

log ""
log "${BOLD}[2] FILE-BY-FILE COMPARISON (KEGG analysis)${NC}"
compare_folder_pair "$MY_KEGG" "$COL_KEGG" "KEGG" "Shared_508 (mmu)"

# --- Deep content check on _identify.txt files ---
log ""
hr
log "${BOLD}[3] DEEP CONTENT CHECK (_identify.txt enrichment results)${NC}"

# Find identify files dynamically
GO_ID_MINE=$(  find "$MY_GO"   -name "*identify*" -type f 2>/dev/null | head -1)
GO_ID_COL=$(   find "$COL_GO"  -name "*identify*" -type f 2>/dev/null | head -1)
KEGG_ID_MINE=$(find "$MY_KEGG" -name "*identify*" -type f 2>/dev/null | head -1)
KEGG_ID_COL=$( find "$COL_KEGG" -name "*identify*" -type f 2>/dev/null | head -1)

content_check_identify "$GO_ID_MINE"   "$GO_ID_COL"   "GO _identify.txt"
content_check_identify "$KEGG_ID_MINE" "$KEGG_ID_COL" "KEGG _identify.txt"

# --- Summary ---
log ""
hr
log "${BOLD}[4] SUMMARY${NC}"
log ""
log "  Files compared:  $TOTAL_COMPARED"
log "  ${GREEN}Identical:       $TOTAL_IDENTICAL${NC}"
log "  ${RED}Different:       $TOTAL_DIFFERENT${NC}"
log "  ${YELLOW}Missing:         $TOTAL_MISSING${NC}"
log ""

if (( TOTAL_DIFFERENT == 0 && TOTAL_MISSING == 0 && TOTAL_COMPARED > 0 )); then
    log "  ${GREEN}${BOLD}CONCLUSION: All compared files are IDENTICAL.${NC}"
    log "  Your results and your colleague's results match perfectly."
elif (( TOTAL_DIFFERENT == 0 && TOTAL_MISSING > 0 )); then
    log "  ${YELLOW}${BOLD}CONCLUSION: Files present in both folders are identical,${NC}"
    log "  ${YELLOW}but some files are missing in one folder (see details above).${NC}"
elif (( TOTAL_DIFFERENT > 0 )); then
    log "  ${RED}${BOLD}CONCLUSION: Differences detected.${NC}"
    log "  Review the DIFFERENT entries above. Common causes:"
    log "    - Different input gene lists used"
    log "    - Different KOBAS run dates (database version may differ)"
    log "    - Different KOBAS parameters (e-value threshold, etc.)"
    log "    - File was modified after download"
else
    log "  ${YELLOW}No files were compared — check that folders exist and contain files.${NC}"
fi

hr
log "  Full report saved to: $REPORT"
hr

} 2>&1 | tee "$REPORT"

echo ""
echo "Done. Report saved to: $REPORT"
