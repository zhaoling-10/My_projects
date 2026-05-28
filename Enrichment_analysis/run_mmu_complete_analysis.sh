#!/usr/bin/env bash
# =============================================================================
# Complete Mus musculus–focused analysis pipeline
#   1. Compare your mmu KOBAS folders vs colleague's (6 pairs)
#   2. Classify up/down regulation from DESeq2 and link to mmu enrichment
#   3. Dot/bubble/bar plots of mmu KOBAS enrichment
#   4. Ranked candidate tables + seasonal module (supervisor-ready)
#   5. Expression direction for all KOBAS enriched genes
#   6. Upregulated-only KOBAS plots (+ comparison vs all genes)
#   7. Multi-panel regulation plots (volcano, heatmap, network, ...)
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

echo "=============================================="
echo " Step 1: KOBAS folder comparison (mmu only)"
echo "=============================================="
bash "$SCRIPT_DIR/compare_kobas_folders_mmu.sh"

echo ""
echo "=============================================="
echo " Step 2: Expression + mmu enrichment integration"
echo "=============================================="
if ! command -v Rscript >/dev/null 2>&1; then
    echo "ERROR: Rscript not found. Install R to run mmu_enrichment_expression_analysis.R"
    exit 1
fi
Rscript "$SCRIPT_DIR/mmu_enrichment_expression_analysis.R"

echo ""
echo "=============================================="
echo " Step 3: KOBAS enrichment plots (dot / bubble / bar)"
echo "=============================================="
Rscript "$SCRIPT_DIR/visualize_mmu_kobas_results.R"

echo ""
echo "=============================================="
echo " Step 4: Ranked candidates + seasonal module tables"
echo "=============================================="
Rscript "$SCRIPT_DIR/rank_candidates_and_module_tables.R"

echo ""
echo "=============================================="
echo " Step 5: Expression direction (all KOBAS enriched genes)"
echo "=============================================="
Rscript "$SCRIPT_DIR/kobas_enrichment_expression_direction.R"

echo ""
echo "=============================================="
echo " Step 6: Upregulated-only plots (comparison figures)"
echo "=============================================="
Rscript "$SCRIPT_DIR/visualize_kobas_upregulated_only.R"

echo ""
echo "=============================================="
echo " Step 7: Multi-panel LT vs LE regulation plots"
echo "=============================================="
Rscript "$SCRIPT_DIR/expression_regulation_multipanel_plots.R"

echo ""
echo "=============================================="
echo " All done."
echo "  - Comparison report: kobas_comparison_report_mmu_*.txt"
echo "  - Expression tables:  mmu_analysis/"
echo "  - KOBAS figures:      mmu_analysis/figures/"
echo "  - Meeting tables:     mmu_analysis/tables/"
echo "    (start with seasonal_module_summary.txt)"
echo "  - All pathway genes + direction:"
echo "    kobas_all_pathway_genes_with_direction.csv"
echo "  - Upregulated-only figures: mmu_analysis/figures/upregulated_only/"
echo "  - Regulation multipanel:    mmu_analysis/figures/regulation_multipanel/"
echo "=============================================="
