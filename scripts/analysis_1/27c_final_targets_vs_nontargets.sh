#!/bin/bash
#SBATCH --job-name=a1_27c_final
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --time=1:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --output=logs/27c_final_targets_vs_nontargets.out
#SBATCH --error=logs/27c_final_targets_vs_nontargets.err

# =============================================================================
# FINAL PROPER COMPARISON: TARGETS vs NON-TARGETS (Hypermethylated)
# =============================================================================
#
# This script compares MUTUALLY EXCLUSIVE groups:
# - TARGETS: 3,052 hypermethylated genes WITH TES/TEAD1 binding (within gene+10kb)
# - NON-TARGETS: 3,156 hypermethylated genes WITHOUT any binding
#
# Previous comparisons showed similar binding because groups were overlapping.
# This corrects that by using properly separated populations.
#
# =============================================================================

echo "=========================================="
echo "FINAL TARGETS vs NON-TARGETS COMPARISON"
echo "=========================================="
echo "Started: $(date)"
echo ""

cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_1

OUTDIR="output/27_non_targets_extended_flanks"

# -----------------------------------------------------------------------------
# FILE PATHS
# -----------------------------------------------------------------------------

# BED files - properly separated groups
TARGETS_BED="${OUTDIR}/TARGETS_hyper_with_binding.bed"
NONTARGETS_BED="${OUTDIR}/strict_non_targets_hyper.bed"

# BigWig files
TES_METH="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP/results/05_bigwig/TES_combined_RPKM.bw"
GFP_METH="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP/results/05_bigwig/GFP_combined_RPKM.bw"
TES_BIND="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_CUTandTAG/results/06_bigwig/TES_comb.bw"
TEAD1_BIND="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_CUTandTAG/results/06_bigwig/TEAD1_comb.bw"

# Verify files exist
echo "Verifying input files..."
for f in "$TARGETS_BED" "$NONTARGETS_BED" "$TES_METH" "$GFP_METH" "$TES_BIND" "$TEAD1_BIND"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: File not found: $f"
        exit 1
    fi
done
echo "  All files present."
echo ""

# Count genes
N_TARGETS=$(wc -l < ${TARGETS_BED})
N_NONTARGETS=$(wc -l < ${NONTARGETS_BED})

echo "Gene counts (mutually exclusive):"
echo "  TARGETS (hypermeth + binding): ${N_TARGETS}"
echo "  NON-TARGETS (hypermeth + NO binding): ${N_NONTARGETS}"
echo "  Total: $((N_TARGETS + N_NONTARGETS))"
echo ""

# -----------------------------------------------------------------------------
# ACTIVATE ENVIRONMENT
# -----------------------------------------------------------------------------

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate tg

# -----------------------------------------------------------------------------
# MAIN COMPARISON PLOT
# -----------------------------------------------------------------------------

echo "=== Creating main comparison plot ==="

computeMatrix scale-regions \
    -S $TES_METH $GFP_METH $TES_BIND $TEAD1_BIND \
    -R ${TARGETS_BED} ${NONTARGETS_BED} \
    --beforeRegionStartLength 10000 \
    --afterRegionStartLength 10000 \
    --regionBodyLength 5000 \
    --binSize 100 \
    --skipZeros \
    --missingDataAsZero \
    -o ${OUTDIR}/FINAL_targets_nontargets_matrix.gz \
    -p 8 \
    2>&1 | grep -v "Skipping\|did not match"

plotProfile -m ${OUTDIR}/FINAL_targets_nontargets_matrix.gz \
    -out ${OUTDIR}/FINAL_targets_vs_nontargets_hyper.png \
    --perGroup \
    --colors "#7B3294" "#636363" "#E31A1C" "#377EB8" \
    --startLabel "TSS" \
    --endLabel "TES" \
    --samplesLabel "TES meth" "GFP meth" "TES bind" "TEAD1 bind" \
    --regionsLabel "TARGETS (n=${N_TARGETS})" "NON-TARGETS (n=${N_NONTARGETS})" \
    --plotTitle "Hypermethylated Genes: TARGETS (bound) vs NON-TARGETS (unbound)" \
    --plotHeight 12 \
    --plotWidth 20 \
    --legendLocation "upper-left" \
    --yMin 0 \
    --dpi 300

echo "  Done: FINAL_targets_vs_nontargets_hyper.png"

# -----------------------------------------------------------------------------
# BINDING-ONLY PLOT (clearer view)
# -----------------------------------------------------------------------------

echo ""
echo "=== Creating binding-only comparison ==="

computeMatrix scale-regions \
    -S $TES_BIND $TEAD1_BIND \
    -R ${TARGETS_BED} ${NONTARGETS_BED} \
    --beforeRegionStartLength 10000 \
    --afterRegionStartLength 10000 \
    --regionBodyLength 5000 \
    --binSize 100 \
    --skipZeros \
    --missingDataAsZero \
    -o ${OUTDIR}/FINAL_binding_only_matrix.gz \
    -p 8 \
    2>&1 | grep -v "Skipping\|did not match"

plotProfile -m ${OUTDIR}/FINAL_binding_only_matrix.gz \
    -out ${OUTDIR}/FINAL_binding_targets_vs_nontargets.png \
    --perGroup \
    --colors "#E31A1C" "#377EB8" \
    --startLabel "TSS" \
    --endLabel "TES" \
    --samplesLabel "TES bind" "TEAD1 bind" \
    --regionsLabel "TARGETS (n=${N_TARGETS})" "NON-TARGETS (n=${N_NONTARGETS})" \
    --plotTitle "TES/TEAD1 Binding: TARGETS vs NON-TARGETS (Hypermethylated)" \
    --plotHeight 10 \
    --plotWidth 18 \
    --legendLocation "upper-right" \
    --yMin 0 \
    --dpi 300

echo "  Done: FINAL_binding_targets_vs_nontargets.png"

# -----------------------------------------------------------------------------
# METHYLATION-ONLY PLOT
# -----------------------------------------------------------------------------

echo ""
echo "=== Creating methylation-only comparison ==="

computeMatrix scale-regions \
    -S $TES_METH $GFP_METH \
    -R ${TARGETS_BED} ${NONTARGETS_BED} \
    --beforeRegionStartLength 10000 \
    --afterRegionStartLength 10000 \
    --regionBodyLength 5000 \
    --binSize 100 \
    --skipZeros \
    --missingDataAsZero \
    -o ${OUTDIR}/FINAL_methylation_only_matrix.gz \
    -p 8 \
    2>&1 | grep -v "Skipping\|did not match"

plotProfile -m ${OUTDIR}/FINAL_methylation_only_matrix.gz \
    -out ${OUTDIR}/FINAL_methylation_targets_vs_nontargets.png \
    --perGroup \
    --colors "#7B3294" "#636363" \
    --startLabel "TSS" \
    --endLabel "TES" \
    --samplesLabel "TES meth" "GFP meth" \
    --regionsLabel "TARGETS (n=${N_TARGETS})" "NON-TARGETS (n=${N_NONTARGETS})" \
    --plotTitle "Methylation: TARGETS vs NON-TARGETS (Both Hypermethylated)" \
    --plotHeight 10 \
    --plotWidth 18 \
    --legendLocation "lower-right" \
    --yMin 0 \
    --dpi 300

echo "  Done: FINAL_methylation_targets_vs_nontargets.png"

# -----------------------------------------------------------------------------
# SUMMARY
# -----------------------------------------------------------------------------

echo ""
echo "=========================================="
echo "COMPLETE"
echo "=========================================="
echo "Finished: $(date)"
echo ""
echo "Output files:"
ls -lh ${OUTDIR}/FINAL_*.png 2>/dev/null
echo ""
echo "Key questions these plots should answer:"
echo "  1. Do TARGETS show elevated TES/TEAD1 binding vs NON-TARGETS? (expected: YES)"
echo "  2. Do both groups show similar hypermethylation patterns? (testing hypothesis)"
echo "  3. Is methylation independent of TES/TEAD1 binding at non-targets?"
echo ""
