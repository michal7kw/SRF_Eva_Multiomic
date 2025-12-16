#!/bin/bash
#SBATCH --job-name=a1_27_non_targets
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --time=3:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --output=logs/27_non_targets_extended_flanks.out
#SBATCH --error=logs/27_non_targets_extended_flanks.err

# =============================================================================
# NON-TARGET GENES WITH EXTENDED FLANKING REGIONS
# =============================================================================
#
# This script analyzes genes that are NOT bound by TES or TEAD1 to see:
# 1. Do methylation changes occur at non-target genes?
# 2. How far do methylation differences extend from gene body?
#
# We test multiple flanking region sizes to find where differences disappear.
#
# =============================================================================

echo "=========================================="
echo "NON-TARGET GENES - EXTENDED FLANKS"
echo "=========================================="
echo "Started: $(date)"
echo ""

cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_1

# -----------------------------------------------------------------------------
# STEP 1: RUN R SCRIPT TO DEFINE NON-TARGET GENE SETS
# -----------------------------------------------------------------------------
echo "=== Step 1: Defining non-target gene sets ==="

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate r_chipseq_env

Rscript 27_non_targets_extended_flanks.R

# Check if BED files were created
OUTDIR="output/27_non_targets_extended_flanks"
if [ ! -f "${OUTDIR}/non_targets_hypermethylated.bed" ]; then
    echo "ERROR: BED files not created. Exiting."
    exit 1
fi

# -----------------------------------------------------------------------------
# STEP 2: SWITCH TO DEEPTOOLS ENVIRONMENT
# -----------------------------------------------------------------------------
echo ""
echo "=== Step 2: Setting up visualization ==="

conda activate tg

# BigWig files
TES_METH="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP/results/05_bigwig/TES_combined_RPKM.bw"
GFP_METH="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP/results/05_bigwig/GFP_combined_RPKM.bw"
TES_BIND="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_CUTandTAG/results/06_bigwig/TES_comb.bw"
TEAD1_BIND="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_CUTandTAG/results/06_bigwig/TEAD1_comb.bw"

# Count genes
N_HYPER=$(wc -l < ${OUTDIR}/non_targets_hypermethylated.bed)
N_HYPO=$(wc -l < ${OUTDIR}/non_targets_hypomethylated.bed)
N_NODMR=$(wc -l < ${OUTDIR}/non_targets_no_dmr.bed)

echo "Non-target gene counts:"
echo "  Hypermethylated: ${N_HYPER}"
echo "  Hypomethylated: ${N_HYPO}"
echo "  No DMR: ${N_NODMR}"
echo ""

# -----------------------------------------------------------------------------
# STEP 3: GENERATE PROFILES WITH DIFFERENT FLANK SIZES
# -----------------------------------------------------------------------------

# Test multiple flank sizes to find where differences disappear
FLANK_SIZES="5000 10000 15000"

for FLANK in ${FLANK_SIZES}; do
    FLANK_KB=$((FLANK / 1000))
    echo ""
    echo "=== Processing with ${FLANK_KB}kb flanks ==="

    # Create comparison matrix (all 3 gene sets)
    echo "  Creating comparison matrix..."
    computeMatrix scale-regions \
        -S $TES_METH $GFP_METH $TES_BIND $TEAD1_BIND \
        -R ${OUTDIR}/non_targets_hypermethylated.bed \
           ${OUTDIR}/non_targets_hypomethylated.bed \
           ${OUTDIR}/non_targets_no_dmr.bed \
        --beforeRegionStartLength ${FLANK} \
        --afterRegionStartLength ${FLANK} \
        --regionBodyLength 5000 \
        --binSize 100 \
        --skipZeros \
        --missingDataAsZero \
        -o ${OUTDIR}/comparison_${FLANK_KB}kb_matrix.gz \
        -p 16 \
        2>&1 | grep -v "Skipping\|did not match"

    # Create profile plot
    echo "  Creating profile plot..."
    plotProfile -m ${OUTDIR}/comparison_${FLANK_KB}kb_matrix.gz \
        -out ${OUTDIR}/non_targets_comparison_${FLANK_KB}kb_flanks.png \
        --perGroup \
        --colors "#7B3294" "#636363" "#E31A1C" "#377EB8" \
        --startLabel "TSS" \
        --endLabel "TES" \
        --samplesLabel "TES meth" "GFP meth" "TES bind" "TEAD1 bind" \
        --regionsLabel "Hyper (n=${N_HYPER})" "Hypo (n=${N_HYPO})" "No DMR (n=${N_NODMR})" \
        --plotTitle "NON-TARGETS: Gene Body Profiles (${FLANK_KB}kb flanks)" \
        --plotHeight 12 \
        --plotWidth 20 \
        --legendLocation "upper-left" \
        --yMin 0 \
        --dpi 300

    echo "  Done: non_targets_comparison_${FLANK_KB}kb_flanks.png"
done

# -----------------------------------------------------------------------------
# STEP 4: HYPERMETHYLATED NON-TARGETS ONLY (DETAILED VIEW)
# -----------------------------------------------------------------------------
echo ""
echo "=== Hypermethylated non-targets detailed view ==="

# 10kb flanks for detailed view
computeMatrix scale-regions \
    -S $TES_METH $GFP_METH $TES_BIND $TEAD1_BIND \
    -R ${OUTDIR}/non_targets_hypermethylated.bed \
    --beforeRegionStartLength 10000 \
    --afterRegionStartLength 10000 \
    --regionBodyLength 5000 \
    --binSize 100 \
    --skipZeros \
    --missingDataAsZero \
    -o ${OUTDIR}/hyper_non_targets_10kb_matrix.gz \
    -p 16 \
    2>&1 | grep -v "Skipping\|did not match"

plotProfile -m ${OUTDIR}/hyper_non_targets_10kb_matrix.gz \
    -out ${OUTDIR}/hyper_non_targets_10kb_profile.png \
    --perGroup \
    --colors "#7B3294" "#636363" "#E31A1C" "#377EB8" \
    --startLabel "TSS" \
    --endLabel "TES" \
    --samplesLabel "TES meth" "GFP meth" "TES bind" "TEAD1 bind" \
    --plotTitle "NON-TARGETS with Hypermethylated DMRs (n=${N_HYPER})" \
    --plotHeight 10 \
    --plotWidth 18 \
    --legendLocation "lower-right" \
    --yMin 0 \
    --dpi 300

echo "  Done: hyper_non_targets_10kb_profile.png"

# -----------------------------------------------------------------------------
# STEP 5: COMPARE TARGETS VS NON-TARGETS (HYPERMETHYLATED)
# -----------------------------------------------------------------------------
echo ""
echo "=== Comparing Targets vs Non-Targets (Hypermethylated) ==="

# Check if target gene BED exists from script 25
TARGET_BED="output/25_genebody_methylation_DMR_genes/genes_with_hypermethylated_DMRs.bed"

if [ -f "${TARGET_BED}" ]; then
    N_TARGETS=$(wc -l < ${TARGET_BED})

    computeMatrix scale-regions \
        -S $TES_METH $GFP_METH $TES_BIND $TEAD1_BIND \
        -R ${TARGET_BED} \
           ${OUTDIR}/non_targets_hypermethylated.bed \
        --beforeRegionStartLength 10000 \
        --afterRegionStartLength 10000 \
        --regionBodyLength 5000 \
        --binSize 100 \
        --skipZeros \
        --missingDataAsZero \
        -o ${OUTDIR}/targets_vs_non_targets_matrix.gz \
        -p 16 \
        2>&1 | grep -v "Skipping\|did not match"

    plotProfile -m ${OUTDIR}/targets_vs_non_targets_matrix.gz \
        -out ${OUTDIR}/targets_vs_non_targets_hyper.png \
        --perGroup \
        --colors "#7B3294" "#636363" "#E31A1C" "#377EB8" \
        --startLabel "TSS" \
        --endLabel "TES" \
        --samplesLabel "TES meth" "GFP meth" "TES bind" "TEAD1 bind" \
        --regionsLabel "All Hyper (n=${N_TARGETS})" "Non-target Hyper (n=${N_HYPER})" \
        --plotTitle "Hypermethylated: All Genes vs Non-Targets Only" \
        --plotHeight 12 \
        --plotWidth 20 \
        --legendLocation "upper-left" \
        --yMin 0 \
        --dpi 300

    echo "  Done: targets_vs_non_targets_hyper.png"
else
    echo "  Skipping comparison - target BED file not found"
fi

# -----------------------------------------------------------------------------
# SUMMARY
# -----------------------------------------------------------------------------
echo ""
echo "=========================================="
echo "COMPLETE"
echo "=========================================="
echo "Finished: $(date)"
echo ""
echo "Output directory: ${OUTDIR}"
echo ""
echo "PNG files generated:"
ls -lh ${OUTDIR}/*.png
echo ""
echo "Key findings to look for:"
echo "  1. Do non-targets show methylation differences similar to targets?"
echo "  2. At what distance from gene body do differences disappear?"
echo "  3. Is there any residual TES/TEAD1 binding at non-targets?"
echo ""
