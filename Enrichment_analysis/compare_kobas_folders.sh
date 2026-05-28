#!/usr/bin/env bash
# =============================================================================
# KOBAS Output Folder Comparison Script
# Compares your KOBAS results vs colleague's results across all 12 folder pairs
# (mmu/ocu × GO/KEGG × gene sets 508 / 158 / 133)
# =============================================================================

set -euo pipefail

# ---- Colour codes ----
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
BOLD='\033[1m'
NC='\033[0m' # No Colour

# ---- Base paths ----
MY_BASE="/home/zhaoling/Documents/Enrichment_analysis"
COL_BASE="/home/zhaoling/Documents/Data_analysi_for_Ling/9.Enrichment analysis/9.2.Enrichment_analysis"

# ---- All 12 folder pairs: my_folder|colleague_folder|analysis_type|label ----
# Colleague subfolders use 1./2./3.int_ref_* (not all 1.int_ref_*).
# GO set 3 maps to Genes133_adj (colleague folder name on disk).
FOLDER_PAIRS=(
    "${MY_BASE}/KOBAS_output_mmu_GO_1|${COL_BASE}/2. GO enrichment: Mus musculus/1.int_ref_Genes508_adj|GO|mmu GO — 508 genes"
    "${MY_BASE}/KOBAS_output_mmu_GO_2|${COL_BASE}/2. GO enrichment: Mus musculus/2.int_ref_Genes158_adj|GO|mmu GO — 158 genes"
    "${MY_BASE}/KOBAS_output_mmu_GO_3|${COL_BASE}/2. GO enrichment: Mus musculus/3.int_ref_Genes133_adj|GO|mmu GO — 133 genes"
    "${MY_BASE}/KOBAS_output_mmu_KEGG_1|${COL_BASE}/1. KEGG enrichment: Mus musculus/1.int_ref_Genes508_adj|KEGG|mmu KEGG — 508 genes"
    "${MY_BASE}/KOBAS_output_mmu_KEGG_2|${COL_BASE}/1. KEGG enrichment: Mus musculus/2.int_ref_Genes158_adj|KEGG|mmu KEGG — 158 genes"
    "${MY_BASE}/KOBAS_output_mmu_KEGG_3|${COL_BASE}/1. KEGG enrichment: Mus musculus/3.int_ref_Genes133_adj|KEGG|mmu KEGG — 133 genes"
    "${MY_BASE}/KOBAS_output_ocu_GO_1|${COL_BASE}/4.GO Oryctolagus cuniculus_(rabbit)/1.int_ref_Genes508_adj|GO|ocu GO — 508 genes"
    "${MY_BASE}/KOBAS_output_ocu_GO_2|${COL_BASE}/4.GO Oryctolagus cuniculus_(rabbit)/2.int_ref_Genes158_adj|GO|ocu GO — 158 genes"
    "${MY_BASE}/KOBAS_output_ocu_GO_3|${COL_BASE}/4.GO Oryctolagus cuniculus_(rabbit)/3.int_ref_Genes133_adj|GO|ocu GO — 133 genes"
    "${MY_BASE}/KOBAS_output_ocu_KEGG_1|${COL_BASE}/3.KEGG Oryctolagus cuniculus_(rabbit)/1.int_ref_Genes508_adj|KEGG|ocu KEGG — 508 genes"
    "${MY_BASE}/KOBAS_output_ocu_KEGG_2|${COL_BASE}/3.KEGG Oryctolagus cuniculus_(rabbit)/2.int_ref_Genes158_adj|KEGG|ocu KEGG — 158 genes"
    "${MY_BASE}/KOBAS_output_ocu_KEGG_3|${COL_BASE}/3.KEGG Oryctolagus cuniculus_(rabbit)/3.int_ref_Genes133_adj|KEGG|ocu KEGG — 133 genes"
)

# ---- Report output (written under Enrichment_analysis) ----
REPORT="${MY_BASE}/kobas_comparison_report_$(date +%Y%m%d_%H%M%S).txt"

# ---- Counters ----
TOTAL_COMPARED=0
TOTAL_IDENTICAL=0
TOTAL_DIFFERENT=0
TOTAL_MISSING=0

# Per-pair summary lines (filled during run)
declare -a PAIR_SUMMARIES=()

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

# Map your filename to colleague's filename (prefixes differ between folders)
map_to_colleague_name() {
    local my_basename="$1"
    local analysis_type="$2"   # GO or KEGG
    local species="$3"         # mmu or ocu

    case "$my_basename" in
        mmu.tsv|ocu.tsv)
            printf '%s\n' "$my_basename"
            return 0
            ;;
    esac

    if [[ "$analysis_type" == "GO" ]]; then
        case "$my_basename" in
            *_identify.txt)
                printf 'GO_output_%s_identify.txt\n' "$species"
                ;;
            *_KOBAS_acc_pathways.tsv)
                printf 'GO_output_%s_KOBAS_acc_pathways.tsv\n' "$species"
                ;;
            *_KOBAS_pathways_acc.tsv)
                printf 'GO_output_%s_KOBAS_pathways_acc.tsv\n' "$species"
                ;;
            KOBAS_output_"${species}"_GO_*.tsv)
                printf 'GO_output_%s\n' "$species"
                ;;
            *)
                return 1
                ;;
        esac
    else
        case "$my_basename" in
            *_identify.txt)
                printf 'KOBAS_output_%s_identify.txt\n' "$species"
                ;;
            *_KOBAS_acc_pathways.tsv)
                printf 'KOBAS_output_%s_KOBAS_acc_pathways.tsv\n' "$species"
                ;;
            *_KOBAS_pathways_acc.tsv)
                printf 'KOBAS_output_%s_KOBAS_pathways_acc.tsv\n' "$species"
                ;;
            KOBAS_output_"${species}"_KEGG_*.tsv)
                printf 'KOBAS_output_%s\n' "$species"
                ;;
            *)
                return 1
                ;;
        esac
    fi
}

# Infer species (mmu/ocu) from folder basename
infer_species() {
    local folder="$1"
    local base
    base=$(basename "$folder")
    if [[ "$base" == *"_mmu_"* || "$base" == *"_mmu" ]]; then
        echo mmu
    elif [[ "$base" == *"_ocu_"* || "$base" == *"_ocu" ]]; then
        echo ocu
    else
        echo ""
    fi
}

compare_two_files() {
    local file1="$1"   # your file
    local file2="$2"   # colleague's file
    local label="$3"   # display name

    TOTAL_COMPARED=$((TOTAL_COMPARED + 1))

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

    local size1 size2
    size1=$(wc -c < "$file1")
    size2=$(wc -c < "$file2")

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

        local lines1 lines2
        lines1=$(wc -l < "$file1")
        lines2=$(wc -l < "$file2")
        log "              Your lines:      $lines1"
        log "              Colleague lines: $lines2"

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
    local analysis_type="$3"
    local set_label="$4"

    local species
    species=$(infer_species "$my_folder")
    if [[ -z "$species" ]]; then
        log "  ${RED}ERROR: Cannot infer species (mmu/ocu) from folder: $my_folder${NC}"
        return 1
    fi

    local before_compared=$TOTAL_COMPARED
    local before_identical=$TOTAL_IDENTICAL
    local before_different=$TOTAL_DIFFERENT
    local before_missing=$TOTAL_MISSING

    hr
    log "${BOLD}Comparing: ${set_label}${NC}"
    log "  Your folder:      $my_folder"
    log "  Colleague folder: $col_folder"
    log "  Type: ${analysis_type} | Species: ${species}"
    hr2

    local my_files=()
    while IFS= read -r -d '' f; do
        my_files+=("$f")
    done < <(find "$my_folder" -maxdepth 1 -type f -print0 2>/dev/null | sort -z)

    if (( ${#my_files[@]} == 0 )); then
        log "  ${YELLOW}WARNING: Your folder is empty or not accessible${NC}"
        PAIR_SUMMARIES+=("${set_label}|EMPTY|0|0|0|0")
        return
    fi

    for my_file in "${my_files[@]}"; do
        local basename col_name col_file
        basename=$(basename "$my_file")

        if ! col_name=$(map_to_colleague_name "$basename" "$analysis_type" "$species"); then
            log "  ${YELLOW}SKIP (no mapping)${NC}  $basename"
            continue
        fi

        col_file="${col_folder}/${col_name}"
        log "  Pair: ${basename}  <->  ${col_name}"
        compare_two_files "$my_file" "$col_file" "${basename}  <->  ${col_name}"
    done

    local p_compared p_identical p_different p_missing
    p_compared=$((TOTAL_COMPARED - before_compared))
    p_identical=$((TOTAL_IDENTICAL - before_identical))
    p_different=$((TOTAL_DIFFERENT - before_different))
    p_missing=$((TOTAL_MISSING - before_missing))

    log ""
    log "  Checking for unmapped files in colleague's folder..."
    local col_files=()
    while IFS= read -r -d '' f; do
        col_files+=("$f")
    done < <(find "$col_folder" -maxdepth 1 -type f -print0 2>/dev/null | sort -z)

    local extra_count=0
    for col_file in "${col_files[@]}"; do
        local bn mapped=0
        bn=$(basename "$col_file")
        for my_file in "${my_files[@]}"; do
            local my_bn expected_col
            my_bn=$(basename "$my_file")
            expected_col=$(map_to_colleague_name "$my_bn" "$analysis_type" "$species" 2>/dev/null) || continue
            if [[ "$bn" == "$expected_col" ]]; then
                mapped=1
                break
            fi
        done
        if (( ! mapped )); then
            log "  ${YELLOW}EXTRA IN COLLEAGUE${NC}  $bn (no mapped file in your folder)"
            extra_count=$((extra_count + 1))
        fi
    done
    if (( extra_count == 0 )); then
        log "  ${GREEN}No unmapped extra files in colleague's folder${NC}"
    fi

    local pair_status="OK"
    if (( p_different > 0 )); then
        pair_status="DIFFER"
    elif (( p_missing > 0 )); then
        pair_status="MISSING"
    elif (( p_compared == 0 )); then
        pair_status="NO_FILES"
    fi
    PAIR_SUMMARIES+=("${set_label}|${pair_status}|${p_compared}|${p_identical}|${p_different}|${p_missing}")

    # Deep content check on _identify files for this pair
    local id_mine id_col
    id_mine=$(find "$my_folder" -name "*identify*" -type f 2>/dev/null | head -1)
    if [[ -n "$id_mine" ]]; then
        local id_bn id_col_name
        id_bn=$(basename "$id_mine")
        id_col_name=$(map_to_colleague_name "$id_bn" "$analysis_type" "$species")
        id_col="${col_folder}/${id_col_name}"
        content_check_identify "$id_mine" "$id_col" "${set_label} — identify"
    fi
}

content_check_identify() {
    local file1="$1"
    local file2="$2"
    local label="$3"

    if [[ ! -f "$file1" || ! -f "$file2" ]]; then return; fi

    log ""
    log "  ${BOLD}Content check: $label${NC}"

    local terms1 terms2
    terms1=$(grep -v '^#' "$file1" | grep -v '^$' | wc -l || echo 0)
    terms2=$(grep -v '^#' "$file2" | grep -v '^$' | wc -l || echo 0)
    log "    Enriched terms — Yours: $terms1 | Colleague: $terms2"

    if [[ "$terms1" != "$terms2" ]]; then
        log "    ${YELLOW}Term count DIFFERS${NC}"

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
log "${BOLD}  KOBAS OUTPUT COMPARISON REPORT (12 folder pairs)${NC}"
log "  Generated: $(date)"
log "  Script:    compare_kobas_folders.sh"
hr

log ""
log "${BOLD}[0] FOLDER EXISTENCE CHECK (all 12 pairs)${NC}"
log ""
all_exist=1
for entry in "${FOLDER_PAIRS[@]}"; do
    IFS='|' read -r my_folder col_folder analysis_type set_label <<< "$entry"
    check_folder_exists "$my_folder" "Yours: ${set_label}"       || all_exist=0
    check_folder_exists "$col_folder" "Colleague: ${set_label}"  || all_exist=0
    log ""
done

if (( ! all_exist )); then
    log "${YELLOW}One or more folders not found — comparisons for missing paths will show errors.${NC}"
    log "Note: Colleague gene-set folders use 1./2./3.int_ref_* prefixes on disk"
    log "      (e.g. set 2 -> 2.int_ref_Genes158_adj, not 1.int_ref_Genes158_adj)."
    log ""
fi

log "${BOLD}[1] FILE-BY-FILE COMPARISON (all pairs)${NC}"
PAIR_NUM=0
for entry in "${FOLDER_PAIRS[@]}"; do
    IFS='|' read -r my_folder col_folder analysis_type set_label <<< "$entry"
    PAIR_NUM=$((PAIR_NUM + 1))
    log ""
    log "${BOLD}--- Pair ${PAIR_NUM}/12: ${set_label} ---${NC}"
    if [[ -d "$my_folder" && -d "$col_folder" ]]; then
        compare_folder_pair "$my_folder" "$col_folder" "$analysis_type" "$set_label"
    else
        log "  ${RED}Skipped (folder missing)${NC}"
        PAIR_SUMMARIES+=("${set_label}|SKIPPED|0|0|0|0")
    fi
done

log ""
hr
log "${BOLD}[2] PER-PAIR SUMMARY${NC}"
log ""
printf "  %-32s  %-10s  compared  identical  different  missing\n" "Pair" "Status"
hr2
for summary in "${PAIR_SUMMARIES[@]}"; do
    IFS='|' read -r lbl status pc pi pd pm <<< "$summary"
    case "$status" in
        OK)       status_color="${GREEN}${status}${NC}" ;;
        DIFFER)   status_color="${RED}${status}${NC}" ;;
        MISSING)  status_color="${YELLOW}${status}${NC}" ;;
        *)        status_color="${YELLOW}${status}${NC}" ;;
    esac
    log "$(printf '  %-32s  ' "$lbl")${status_color}  $(printf '%8d %10d %10d %8d' "$pc" "$pi" "$pd" "$pm")"
done

log ""
hr
log "${BOLD}[3] OVERALL SUMMARY${NC}"
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
