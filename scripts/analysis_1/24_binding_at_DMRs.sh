#!/bin/bash
#SBATCH --job-name=a1_24_binding_DMRs
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --time=4:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --output=logs/24_binding_at_DMRs.out
#SBATCH --error=logs/24_binding_at_DMRs.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it

echo "=========================================="
echo "TES/TEAD1 BINDING AT DMRs"
echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Started: $(date)"
echo ""
echo "Purpose: Show that TES binding correlates with methylation increase"
echo "         by visualizing binding signal CENTERED ON DMRs"
echo ""

# Set working directory
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_1

# Output directory
OUTDIR="output/24_binding_at_DMRs"
mkdir -p ${OUTDIR}/{matrices,profiles,heatmaps}
mkdir -p logs

# =============================================================================
# DATA PATHS
# =============================================================================

# DMR BED files from meDIP analysis
DMR_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP/results/17_dmr_heatmaps/beds"
HYPER_BED="${DMR_DIR}/hypermethylated_dmrs.bed"
HYPO_BED="${DMR_DIR}/hypomethylated_dmrs.bed"
ALL_DMR_BED="${DMR_DIR}/all_dmrs.bed"
STRINGENT_BED="${DMR_DIR}/stringent_dmrs.bed"

# Cut&Tag BigWig files (binding)
CUTANDTAG_BIGWIG="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_CUTandTAG/results/06_bigwig"
TES_BIND="${CUTANDTAG_BIGWIG}/TES_comb.bw"
TEAD1_BIND="${CUTANDTAG_BIGWIG}/TEAD1_comb.bw"

# meDIP BigWig files (methylation)
MEDIP_BIGWIG="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP/results/05_bigwig"
TES_METH="${MEDIP_BIGWIG}/TES_combined_RPKM.bw"
GFP_METH="${MEDIP_BIGWIG}/GFP_combined_RPKM.bw"

# =============================================================================
# VERIFY INPUT FILES
# =============================================================================

echo "=== Verifying Input Files ==="
echo ""

# Check DMR files
echo "DMR BED files:"
for f in "$HYPER_BED" "$HYPO_BED" "$ALL_DMR_BED" "$STRINGENT_BED"; do
    if [[ -f "$f" ]]; then
        count=$(wc -l < "$f")
        echo "  Found: $(basename $f) ($count regions)"
    else
        echo "  WARNING: Not found: $f"
    fi
done

echo ""
echo "BigWig files:"

# Check BigWig files
for bw in "$TES_BIND" "$TEAD1_BIND" "$TES_METH" "$GFP_METH"; do
    if [ -f "$bw" ]; then
        echo "  Found: $(basename $bw)"
    else
        echo "  ERROR: BigWig file not found: $bw"
        exit 1
    fi
done

# =============================================================================
# ACTIVATE ENVIRONMENT
# =============================================================================

echo ""
echo "=== Activating Environment ==="
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate tg  # deepTools environment

# Verify deepTools is available
if ! command -v computeMatrix &> /dev/null; then
    echo "ERROR: deepTools not found"
    exit 1
fi
echo "deepTools available"

# =============================================================================
# FUNCTION: Create binding profile at DMRs
# =============================================================================

create_dmr_profile() {
    local dmr_type=$1
    local bed_file=$2
    local title=$3

    if [[ ! -f "$bed_file" ]]; then
        echo "  Skipping $dmr_type - BED file not found"
        return
    fi

    local region_count=$(wc -l < "$bed_file")
    echo ""
    echo "=== Processing $dmr_type ($region_count regions) ==="

    # -------------------------------------------------------------------------
    # Matrix 1: BINDING ONLY at DMRs
    # -------------------------------------------------------------------------
    echo "  Computing binding matrix centered on DMRs..."
    computeMatrix reference-point \
        --referencePoint center \
        -b 5000 -a 5000 \
        -R "$bed_file" \
        -S "$TES_BIND" "$TEAD1_BIND" \
        --samplesLabel "TES binding" "TEAD1 binding" \
        --skipZeros \
        --binSize 50 \
        --numberOfProcessors 16 \
        -o "${OUTDIR}/matrices/${dmr_type}_binding_matrix.gz" \
        2>&1 | grep -v "Skipping"

    # Binding profile plot
    echo "  Creating binding profile..."
    plotProfile -m "${OUTDIR}/matrices/${dmr_type}_binding_matrix.gz" \
        -out "${OUTDIR}/profiles/${dmr_type}_binding_profile.pdf" \
        --perGroup \
        --colors "#E31A1C" "#377EB8" \
        --refPointLabel "DMR center" \
        --plotTitle "TES/TEAD1 Binding at $title" \
        --yAxisLabel "Mean CPM" \
        --plotHeight 8 \
        --plotWidth 12 \
        2>&1 | grep -v "Skipping"

    # Binding heatmap
    echo "  Creating binding heatmap..."
    plotHeatmap -m "${OUTDIR}/matrices/${dmr_type}_binding_matrix.gz" \
        -out "${OUTDIR}/heatmaps/${dmr_type}_binding_heatmap.pdf" \
        --colorMap Blues Reds \
        --sortRegions descend \
        --sortUsing mean \
        --heatmapHeight 15 \
        --heatmapWidth 4 \
        --refPointLabel "DMR" \
        --plotTitle "TES/TEAD1 Binding at $title" \
        --legendLocation upper-right \
        2>&1 | grep -v "Skipping"

    # -------------------------------------------------------------------------
    # Matrix 2: METHYLATION at DMRs (TES vs GFP)
    # -------------------------------------------------------------------------
    echo "  Computing methylation matrix..."
    computeMatrix reference-point \
        --referencePoint center \
        -b 5000 -a 5000 \
        -R "$bed_file" \
        -S "$TES_METH" "$GFP_METH" \
        --samplesLabel "TES meDIP" "GFP meDIP" \
        --skipZeros \
        --binSize 50 \
        --numberOfProcessors 16 \
        -o "${OUTDIR}/matrices/${dmr_type}_methylation_matrix.gz" \
        2>&1 | grep -v "did not match"

    # Methylation profile plot
    echo "  Creating methylation profile..."
    plotProfile -m "${OUTDIR}/matrices/${dmr_type}_methylation_matrix.gz" \
        -out "${OUTDIR}/profiles/${dmr_type}_methylation_profile.pdf" \
        --perGroup \
        --colors "#7B3294" "#636363" \
        --refPointLabel "DMR center" \
        --plotTitle "DNA Methylation at $title" \
        --yAxisLabel "Mean RPKM" \
        --plotHeight 8 \
        --plotWidth 12 \
        2>&1 | grep -v "did not match"

    # Methylation heatmap
    echo "  Creating methylation heatmap..."
    plotHeatmap -m "${OUTDIR}/matrices/${dmr_type}_methylation_matrix.gz" \
        -out "${OUTDIR}/heatmaps/${dmr_type}_methylation_heatmap.pdf" \
        --colorMap Purples Greys \
        --sortRegions descend \
        --sortUsing mean \
        --heatmapHeight 15 \
        --heatmapWidth 4 \
        --refPointLabel "DMR" \
        --plotTitle "Methylation at $title" \
        --legendLocation upper-right \
        2>&1 | grep -v "did not match"

    # -------------------------------------------------------------------------
    # Matrix 3: COMBINED (binding + methylation)
    # -------------------------------------------------------------------------
    echo "  Computing combined matrix (binding + methylation)..."
    computeMatrix reference-point \
        --referencePoint center \
        -b 5000 -a 5000 \
        -R "$bed_file" \
        -S "$TES_BIND" "$TEAD1_BIND" "$TES_METH" "$GFP_METH" \
        --samplesLabel "TES bind" "TEAD1 bind" "TES meth" "GFP meth" \
        --skipZeros \
        --binSize 50 \
        --numberOfProcessors 16 \
        -o "${OUTDIR}/matrices/${dmr_type}_combined_matrix.gz" \
        2>&1 | grep -v "Skipping\|did not match"

    # Combined profile (separate Y-axes for binding vs methylation)
    echo "  Creating combined profile (per-sample panels)..."
    plotProfile -m "${OUTDIR}/matrices/${dmr_type}_combined_matrix.gz" \
        -out "${OUTDIR}/profiles/${dmr_type}_combined_profile.pdf" \
        --perSample \
        --colors "#E31A1C" "#377EB8" "#7B3294" "#636363" \
        --refPointLabel "DMR center" \
        --plotTitle "Binding & Methylation at $title" \
        --plotHeight 8 \
        --plotWidth 14 \
        2>&1 | grep -v "Skipping\|did not match"

    echo "  Done with $dmr_type"
}

# =============================================================================
# PROCESS EACH DMR CATEGORY
# =============================================================================

echo ""
echo "=========================================="
echo "CREATING DMR-CENTERED PROFILES"
echo "=========================================="

# MOST IMPORTANT: Hypermethylated DMRs (TES > GFP)
# If TES causes methylation, we should see TES binding enriched here
create_dmr_profile "hypermethylated" "$HYPER_BED" "Hypermethylated DMRs (TES > GFP)"

# Hypomethylated DMRs (control comparison)
create_dmr_profile "hypomethylated" "$HYPO_BED" "Hypomethylated DMRs (TES < GFP)"

# Stringent DMRs (highest confidence)
create_dmr_profile "stringent" "$STRINGENT_BED" "Stringent DMRs (FDR<0.01, |FC|>4)"

# =============================================================================
# COMPARISON: HYPER vs HYPO DMRs
# =============================================================================

echo ""
echo "=========================================="
echo "CREATING COMPARISON: HYPER vs HYPO"
echo "=========================================="

if [[ -f "$HYPER_BED" ]] && [[ -f "$HYPO_BED" ]]; then

    # Binding comparison
    echo "  Computing binding comparison matrix..."
    computeMatrix reference-point \
        --referencePoint center \
        -b 5000 -a 5000 \
        -R "$HYPER_BED" "$HYPO_BED" \
        -S "$TES_BIND" "$TEAD1_BIND" \
        --samplesLabel "TES binding" "TEAD1 binding" \
        --skipZeros \
        --binSize 50 \
        --numberOfProcessors 16 \
        -o "${OUTDIR}/matrices/hyper_vs_hypo_binding_matrix.gz" \
        2>&1 | grep -v "Skipping"

    # Binding comparison profile
    echo "  Creating binding comparison profile..."
    plotProfile -m "${OUTDIR}/matrices/hyper_vs_hypo_binding_matrix.gz" \
        -out "${OUTDIR}/profiles/hyper_vs_hypo_binding_profile.pdf" \
        --perGroup \
        --colors "#E31A1C" "#377EB8" \
        --refPointLabel "DMR center" \
        --regionsLabel "Hypermethylated" "Hypomethylated" \
        --plotTitle "TES/TEAD1 Binding: Hypermethylated vs Hypomethylated DMRs" \
        --yAxisLabel "Mean CPM" \
        --plotHeight 8 \
        --plotWidth 14 \
        --legendLocation upper-right \
        2>&1 | grep -v "Skipping"

    # Binding comparison heatmap
    echo "  Creating binding comparison heatmap..."
    plotHeatmap -m "${OUTDIR}/matrices/hyper_vs_hypo_binding_matrix.gz" \
        -out "${OUTDIR}/heatmaps/hyper_vs_hypo_binding_heatmap.pdf" \
        --colorMap Blues Reds \
        --sortRegions descend \
        --sortUsing mean \
        --heatmapHeight 15 \
        --heatmapWidth 4 \
        --refPointLabel "DMR" \
        --regionsLabel "Hypermethylated" "Hypomethylated" \
        --plotTitle "Binding at Hyper vs Hypo DMRs" \
        2>&1 | grep -v "Skipping"

    # Methylation comparison
    echo "  Computing methylation comparison matrix..."
    computeMatrix reference-point \
        --referencePoint center \
        -b 5000 -a 5000 \
        -R "$HYPER_BED" "$HYPO_BED" \
        -S "$TES_METH" "$GFP_METH" \
        --samplesLabel "TES meDIP" "GFP meDIP" \
        --skipZeros \
        --binSize 50 \
        --numberOfProcessors 16 \
        -o "${OUTDIR}/matrices/hyper_vs_hypo_methylation_matrix.gz" \
        2>&1 | grep -v "did not match"

    # Methylation comparison profile
    echo "  Creating methylation comparison profile..."
    plotProfile -m "${OUTDIR}/matrices/hyper_vs_hypo_methylation_matrix.gz" \
        -out "${OUTDIR}/profiles/hyper_vs_hypo_methylation_profile.pdf" \
        --perGroup \
        --colors "#7B3294" "#636363" \
        --refPointLabel "DMR center" \
        --regionsLabel "Hypermethylated" "Hypomethylated" \
        --plotTitle "DNA Methylation: Hypermethylated vs Hypomethylated DMRs" \
        --yAxisLabel "Mean RPKM" \
        --plotHeight 8 \
        --plotWidth 14 \
        --legendLocation upper-right \
        2>&1 | grep -v "did not match"

    echo "  Done with hyper vs hypo comparison"
fi

# =============================================================================
# COMBINED 4-SIGNAL COMPARISON (most informative plot)
# =============================================================================

echo ""
echo "=========================================="
echo "CREATING COMBINED 4-SIGNAL COMPARISON"
echo "=========================================="

if [[ -f "$HYPER_BED" ]] && [[ -f "$HYPO_BED" ]]; then

    echo "  Computing combined 4-signal matrix..."
    computeMatrix reference-point \
        --referencePoint center \
        -b 5000 -a 5000 \
        -R "$HYPER_BED" "$HYPO_BED" \
        -S "$TES_BIND" "$TEAD1_BIND" "$TES_METH" "$GFP_METH" \
        --samplesLabel "TES bind" "TEAD1 bind" "TES meth" "GFP meth" \
        --skipZeros \
        --binSize 50 \
        --numberOfProcessors 16 \
        -o "${OUTDIR}/matrices/comparison_all_signals_matrix.gz" \
        2>&1 | grep -v "Skipping\|did not match"

    # Per-sample profile (each signal gets its own Y-axis)
    echo "  Creating per-sample profile..."
    plotProfile -m "${OUTDIR}/matrices/comparison_all_signals_matrix.gz" \
        -out "${OUTDIR}/profiles/comparison_all_signals_profile.pdf" \
        --perSample \
        --colors "#E31A1C" "#377EB8" "#7B3294" "#636363" \
        --refPointLabel "DMR center" \
        --regionsLabel "Hypermethylated" "Hypomethylated" \
        --plotTitle "Binding & Methylation: Hyper vs Hypo DMRs" \
        --plotHeight 10 \
        --plotWidth 16 \
        --legendLocation upper-right \
        2>&1 | grep -v "Skipping\|did not match"

    # Combined heatmap with all signals
    echo "  Creating combined heatmap..."
    plotHeatmap -m "${OUTDIR}/matrices/comparison_all_signals_matrix.gz" \
        -out "${OUTDIR}/heatmaps/comparison_all_signals_heatmap.pdf" \
        --colorList "white,#2166AC" "white,#B2182B" "white,#762A83" "white,#636363" \
        --zMin 0 0 0 0 \
        --zMax 2 2 150 150 \
        --sortRegions descend \
        --sortUsing mean \
        --sortUsingSamples 1 2 \
        --heatmapHeight 15 \
        --heatmapWidth 3 \
        --refPointLabel "DMR" \
        --regionsLabel "Hypermethylated" "Hypomethylated" \
        --plotTitle "All Signals at Hyper vs Hypo DMRs" \
        2>&1 | grep -v "Skipping\|did not match"

    echo "  Done with combined comparison"
fi

# =============================================================================
# COMPLETION
# =============================================================================

echo ""
echo "=========================================="
echo "TES/TEAD1 BINDING AT DMRs COMPLETE"
echo "=========================================="
echo "Completed: $(date)"
echo ""
echo "Output directory: ${OUTDIR}/"
echo ""
echo "Generated profiles:"
ls -lh ${OUTDIR}/profiles/*.pdf 2>/dev/null
echo ""
echo "Generated heatmaps:"
ls -lh ${OUTDIR}/heatmaps/*.pdf 2>/dev/null
echo ""
echo "=========================================="
echo "KEY INTERPRETATION:"
echo "=========================================="
echo ""
echo "If TES causes methylation, you should see:"
echo "  1. TES binding ENRICHED at hypermethylated DMRs"
echo "  2. Higher TES binding at hyper vs hypo DMRs"
echo "  3. Clear separation between TES and GFP methylation at hyper DMRs"
echo ""
echo "Key plots to examine:"
echo "  - hypermethylated_binding_profile.pdf: TES/TEAD1 at sites that gained methylation"
echo "  - hyper_vs_hypo_binding_profile.pdf: Compare binding at gain vs loss sites"
echo "  - comparison_all_signals_profile.pdf: All 4 signals showing full picture"
echo ""
