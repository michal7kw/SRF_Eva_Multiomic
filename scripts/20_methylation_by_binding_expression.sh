#!/bin/bash
#SBATCH --job-name=methylation_binding_expr
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --mem=48G
#SBATCH --cpus-per-task=8
#SBATCH --time=6:00:00
#SBATCH --output=logs/20_methylation_binding_expression.out
#SBATCH --error=logs/20_methylation_binding_expression.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it

# ============================================================================
# Methylation Analysis by Binding and Expression Status
# ============================================================================
# Purpose: Comprehensive methylation visualization showing:
#   - meDIP signals at PROMOTERS for genes stratified by:
#       - Binding (TES-bound, TEAD1-bound, both, neither)
#       - Expression (upregulated, downregulated)
#   - meDIP signals at GENE BODIES for same categories
#
# This allows assessment of whether Cut&Tag binding correlates with
# methylation changes and gene expression changes.
# ============================================================================

set -e
set -o pipefail

cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis

# Create directories
OUTDIR="results/20_methylation_binding_expression"
mkdir -p "$OUTDIR"/{matrices,promoter_heatmaps,genebody_heatmaps,profiles}
mkdir -p logs

echo "=============================================="
echo "Methylation by Binding/Expression Analysis"
echo "Started: $(date)"
echo "=============================================="

# ============================================================================
# Step 1: Run R script to classify genes and create BED files
# ============================================================================

echo ""
echo "=== Step 1: Gene classification and BED file creation ==="

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate r_chipseq_env

# Run R classification script
if [[ ! -d "$OUTDIR/beds" ]] || [[ $(ls "$OUTDIR/beds"/*.bed 2>/dev/null | wc -l) -lt 5 ]]; then
    echo "Running R classification script..."
    Rscript scripts/20_methylation_by_binding_expression.R
else
    echo "BED files already exist. Skipping R script."
fi

# Switch to deepTools environment
conda activate peak_calling_new

if ! command -v computeMatrix &> /dev/null; then
    echo "ERROR: deepTools not found"
    exit 1
fi

# ============================================================================
# Define BigWig files
# ============================================================================

# meDIP BigWig files
MEDIP_DIR="../meDIP/results/05_bigwig"
GFP_MEDIP_BW1="$MEDIP_DIR/GFP-1-IP_RPKM.bw"
GFP_MEDIP_BW2="$MEDIP_DIR/GFP-2-IP_RPKM.bw"
TES_MEDIP_BW1="$MEDIP_DIR/TES-1-IP_RPKM.bw"
TES_MEDIP_BW2="$MEDIP_DIR/TES-2-IP_RPKM.bw"

# Cut&Tag BigWig files (for binding signal comparison)
TES_BW="../SRF_Eva_CUTandTAG/results/06_bigwig/TES_comb.bw"
TEAD1_BW="../SRF_Eva_CUTandTAG/results/06_bigwig/TEAD1_comb.bw"

# Verify files
echo ""
echo "=== Checking BigWig files ==="
for f in "$GFP_MEDIP_BW1" "$GFP_MEDIP_BW2" "$TES_MEDIP_BW1" "$TES_MEDIP_BW2" "$TES_BW" "$TEAD1_BW"; do
    if [[ -f "$f" ]]; then
        echo "  Found: $f"
    else
        echo "  WARNING: Not found: $f"
    fi
done

# ============================================================================
# BED files
# ============================================================================

BED_DIR="$OUTDIR/beds"

echo ""
echo "=== Available BED files ==="
for bed in "$BED_DIR"/*_promoter.bed; do
    if [[ -f "$bed" ]]; then
        count=$(wc -l < "$bed")
        name=$(basename "$bed" _promoter.bed)
        echo "  $name: $count genes"
    fi
done

# ============================================================================
# Function to create methylation heatmap
# ============================================================================

create_methylation_heatmap() {
    local category=$1
    local prom_bed="$BED_DIR/${category}_promoter.bed"
    local body_bed="$BED_DIR/${category}_gene_body.bed"

    if [[ ! -f "$prom_bed" ]]; then
        echo "  Skipping $category - promoter BED not found"
        return
    fi

    local region_count=$(wc -l < "$prom_bed")
    if [[ $region_count -lt 10 ]]; then
        echo "  Skipping $category - too few regions ($region_count)"
        return
    fi

    echo ""
    echo "=== Processing $category ($region_count genes) ==="

    # -------------------------------------------------------------------------
    # PROMOTER analysis (reference-point mode)
    # -------------------------------------------------------------------------
    echo "  Creating promoter methylation matrix..."

    # meDIP signals at promoter
    computeMatrix reference-point \
        --referencePoint center \
        -b 5000 -a 5000 \
        -R "$prom_bed" \
        -S "$GFP_MEDIP_BW1" "$GFP_MEDIP_BW2" "$TES_MEDIP_BW1" "$TES_MEDIP_BW2" \
        --samplesLabel "GFP-1" "GFP-2" "TES-1" "TES-2" \
        --skipZeros \
        --binSize 50 \
        --numberOfProcessors 8 \
        -o "$OUTDIR/matrices/${category}_promoter_medip_matrix.gz" \
        2>/dev/null || {
            echo "  WARNING: Promoter matrix failed for $category"
            return
        }

    # Promoter heatmap
    echo "  Creating promoter heatmap..."
    plotHeatmap \
        -m "$OUTDIR/matrices/${category}_promoter_medip_matrix.gz" \
        -o "$OUTDIR/promoter_heatmaps/${category}_promoter_heatmap.pdf" \
        --colorMap YlGn YlGn Purples Purples \
        --sortRegions descend \
        --sortUsing mean \
        --heatmapHeight 15 \
        --heatmapWidth 2.5 \
        --xAxisLabel "Distance from TSS (bp)" \
        --refPointLabel "TSS" \
        --plotTitle "meDIP at Promoters: $category" \
        2>/dev/null

    # Promoter profile
    plotProfile \
        -m "$OUTDIR/matrices/${category}_promoter_medip_matrix.gz" \
        -o "$OUTDIR/profiles/${category}_promoter_profile.pdf" \
        --plotTitle "meDIP at Promoters: $category" \
        --perGroup \
        --colors lightgreen darkgreen plum purple \
        --legendLocation upper-right \
        2>/dev/null

    # -------------------------------------------------------------------------
    # GENE BODY analysis (scale-regions mode)
    # -------------------------------------------------------------------------
    if [[ -f "$body_bed" ]]; then
        echo "  Creating gene body methylation matrix..."

        computeMatrix scale-regions \
            -b 2000 -a 2000 \
            --regionBodyLength 5000 \
            -R "$body_bed" \
            -S "$GFP_MEDIP_BW1" "$GFP_MEDIP_BW2" "$TES_MEDIP_BW1" "$TES_MEDIP_BW2" \
            --samplesLabel "GFP-1" "GFP-2" "TES-1" "TES-2" \
            --skipZeros \
            --binSize 50 \
            --numberOfProcessors 8 \
            -o "$OUTDIR/matrices/${category}_genebody_medip_matrix.gz" \
            2>/dev/null || {
                echo "  WARNING: Gene body matrix failed for $category"
                return
            }

        # Gene body heatmap
        echo "  Creating gene body heatmap..."
        plotHeatmap \
            -m "$OUTDIR/matrices/${category}_genebody_medip_matrix.gz" \
            -o "$OUTDIR/genebody_heatmaps/${category}_genebody_heatmap.pdf" \
            --colorMap YlGn YlGn Purples Purples \
            --sortRegions descend \
            --sortUsing mean \
            --heatmapHeight 15 \
            --heatmapWidth 2.5 \
            --startLabel "TSS" \
            --endLabel "TES" \
            --plotTitle "meDIP at Gene Body: $category" \
            2>/dev/null

        # Gene body profile
        plotProfile \
            -m "$OUTDIR/matrices/${category}_genebody_medip_matrix.gz" \
            -o "$OUTDIR/profiles/${category}_genebody_profile.pdf" \
            --plotTitle "meDIP at Gene Body: $category" \
            --perGroup \
            --colors lightgreen darkgreen plum purple \
            --legendLocation upper-right \
            --startLabel "TSS" \
            --endLabel "TES" \
            2>/dev/null
    fi

    echo "  Done with $category"
}

# ============================================================================
# Process each category
# ============================================================================

echo ""
echo "=== Processing individual categories ==="

# Binding categories
for cat in "TES_only_bound" "TEAD1_only_bound" "TES_TEAD1_bound" "Neither_bound"; do
    create_methylation_heatmap "$cat"
done

# Expression categories
for cat in "upregulated" "downregulated" "unchanged"; do
    create_methylation_heatmap "$cat"
done

# Combined binding + expression categories
for cat in "TES_bound_upregulated" "TES_bound_downregulated" \
           "TEAD1_bound_upregulated" "TEAD1_bound_downregulated" \
           "Both_bound_upregulated" "Both_bound_downregulated" \
           "Indirect_upregulated" "Indirect_downregulated"; do
    create_methylation_heatmap "$cat"
done

# ============================================================================
# Comparative visualizations
# ============================================================================

echo ""
echo "=== Creating comparative visualizations ==="

# -------------------------------------------------------------------------
# Comparison 1: Binding status (TES vs TEAD1 vs Both vs Neither)
# -------------------------------------------------------------------------
echo ""
echo "Comparison: By binding status..."

BINDING_BEDS=""
BINDING_LABELS=""
for cat in "TES_only_bound" "TEAD1_only_bound" "TES_TEAD1_bound" "Neither_bound"; do
    if [[ -f "$BED_DIR/${cat}_promoter.bed" ]]; then
        count=$(wc -l < "$BED_DIR/${cat}_promoter.bed")
        if [[ $count -ge 10 ]]; then
            BINDING_BEDS="$BINDING_BEDS $BED_DIR/${cat}_promoter.bed"
            BINDING_LABELS="$BINDING_LABELS $cat"
        fi
    fi
done

if [[ -n "$BINDING_BEDS" ]]; then
    echo "  Creating binding comparison matrix..."
    computeMatrix reference-point \
        --referencePoint center \
        -b 5000 -a 5000 \
        -R $BINDING_BEDS \
        -S "$GFP_MEDIP_BW1" "$TES_MEDIP_BW1" \
        --samplesLabel "GFP_meDIP" "TES_meDIP" \
        --skipZeros \
        --binSize 50 \
        --numberOfProcessors 8 \
        -o "$OUTDIR/matrices/binding_comparison_promoter_matrix.gz" \
        2>/dev/null

    if [[ -f "$OUTDIR/matrices/binding_comparison_promoter_matrix.gz" ]]; then
        plotProfile \
            -m "$OUTDIR/matrices/binding_comparison_promoter_matrix.gz" \
            -o "$OUTDIR/profiles/binding_comparison_promoter_profile.pdf" \
            --plotTitle "Promoter Methylation by Binding Status" \
            --perGroup \
            --colors green purple \
            --regionsLabel "TES_only" "TEAD1_only" "Both" "Neither" \
            --legendLocation upper-right \
            2>/dev/null

        plotHeatmap \
            -m "$OUTDIR/matrices/binding_comparison_promoter_matrix.gz" \
            -o "$OUTDIR/promoter_heatmaps/binding_comparison_heatmap.pdf" \
            --colorMap YlGn Purples \
            --sortRegions descend \
            --sortUsing mean \
            --heatmapHeight 20 \
            --heatmapWidth 3 \
            --regionsLabel "TES_only" "TEAD1_only" "Both" "Neither" \
            --plotTitle "Promoter Methylation by Binding Status" \
            2>/dev/null
    fi
fi

# -------------------------------------------------------------------------
# Comparison 2: Expression status (Up vs Down)
# -------------------------------------------------------------------------
echo ""
echo "Comparison: By expression status..."

if [[ -f "$BED_DIR/upregulated_promoter.bed" ]] && [[ -f "$BED_DIR/downregulated_promoter.bed" ]]; then
    up_count=$(wc -l < "$BED_DIR/upregulated_promoter.bed")
    down_count=$(wc -l < "$BED_DIR/downregulated_promoter.bed")

    if [[ $up_count -ge 10 ]] && [[ $down_count -ge 10 ]]; then
        echo "  Creating expression comparison matrix..."
        computeMatrix reference-point \
            --referencePoint center \
            -b 5000 -a 5000 \
            -R "$BED_DIR/upregulated_promoter.bed" "$BED_DIR/downregulated_promoter.bed" \
            -S "$GFP_MEDIP_BW1" "$TES_MEDIP_BW1" \
            --samplesLabel "GFP_meDIP" "TES_meDIP" \
            --skipZeros \
            --binSize 50 \
            --numberOfProcessors 8 \
            -o "$OUTDIR/matrices/expression_comparison_promoter_matrix.gz" \
            2>/dev/null

        if [[ -f "$OUTDIR/matrices/expression_comparison_promoter_matrix.gz" ]]; then
            plotProfile \
                -m "$OUTDIR/matrices/expression_comparison_promoter_matrix.gz" \
                -o "$OUTDIR/profiles/expression_comparison_promoter_profile.pdf" \
                --plotTitle "Promoter Methylation: Upregulated vs Downregulated" \
                --perGroup \
                --colors green purple \
                --regionsLabel "Upregulated" "Downregulated" \
                --legendLocation upper-right \
                2>/dev/null

            plotHeatmap \
                -m "$OUTDIR/matrices/expression_comparison_promoter_matrix.gz" \
                -o "$OUTDIR/promoter_heatmaps/expression_comparison_heatmap.pdf" \
                --colorMap YlGn Purples \
                --sortRegions descend \
                --sortUsing mean \
                --heatmapHeight 15 \
                --heatmapWidth 3 \
                --regionsLabel "Upregulated" "Downregulated" \
                --plotTitle "Promoter Methylation: Up vs Down DEGs" \
                2>/dev/null
        fi
    fi
fi

# -------------------------------------------------------------------------
# Comparison 3: TES-bound up vs down
# -------------------------------------------------------------------------
echo ""
echo "Comparison: TES-bound up vs down..."

if [[ -f "$BED_DIR/TES_bound_upregulated_promoter.bed" ]] && \
   [[ -f "$BED_DIR/TES_bound_downregulated_promoter.bed" ]]; then

    up_count=$(wc -l < "$BED_DIR/TES_bound_upregulated_promoter.bed" 2>/dev/null || echo "0")
    down_count=$(wc -l < "$BED_DIR/TES_bound_downregulated_promoter.bed" 2>/dev/null || echo "0")

    if [[ $up_count -ge 10 ]] && [[ $down_count -ge 10 ]]; then
        echo "  Creating TES-bound comparison matrix..."
        computeMatrix reference-point \
            --referencePoint center \
            -b 5000 -a 5000 \
            -R "$BED_DIR/TES_bound_upregulated_promoter.bed" "$BED_DIR/TES_bound_downregulated_promoter.bed" \
            -S "$TES_BW" "$GFP_MEDIP_BW1" "$TES_MEDIP_BW1" \
            --samplesLabel "TES_binding" "GFP_meDIP" "TES_meDIP" \
            --skipZeros \
            --binSize 50 \
            --numberOfProcessors 8 \
            -o "$OUTDIR/matrices/TES_bound_updown_matrix.gz" \
            2>/dev/null

        if [[ -f "$OUTDIR/matrices/TES_bound_updown_matrix.gz" ]]; then
            plotProfile \
                -m "$OUTDIR/matrices/TES_bound_updown_matrix.gz" \
                -o "$OUTDIR/profiles/TES_bound_updown_profile.pdf" \
                --plotTitle "TES-Bound Genes: Binding & Methylation" \
                --perGroup \
                --colors darkblue green purple \
                --regionsLabel "TES_up" "TES_down" \
                --legendLocation upper-right \
                2>/dev/null

            plotHeatmap \
                -m "$OUTDIR/matrices/TES_bound_updown_matrix.gz" \
                -o "$OUTDIR/promoter_heatmaps/TES_bound_updown_heatmap.pdf" \
                --colorList "white,darkblue" "white,green" "white,purple" \
                --sortRegions descend \
                --sortUsing mean \
                --heatmapHeight 15 \
                --heatmapWidth 3 \
                --regionsLabel "TES_up" "TES_down" \
                --plotTitle "TES-Bound: Binding + Methylation" \
                2>/dev/null
        fi
    fi
fi

# ============================================================================
# Summary
# ============================================================================

echo ""
echo "=============================================="
echo "Methylation by Binding/Expression Complete"
echo "=============================================="
echo ""
echo "Output directory: $OUTDIR/"
echo ""
echo "Contents:"
echo "  - beds/: Gene set BED files"
echo "  - matrices/: deepTools matrices"
echo "  - promoter_heatmaps/: Methylation at promoters (PDF)"
echo "  - genebody_heatmaps/: Methylation at gene bodies (PDF)"
echo "  - profiles/: Metagene profiles (PDF)"
echo ""
echo "Key visualizations:"
echo "  - Individual categories: Promoter and gene body methylation"
echo "  - Binding comparison: TES vs TEAD1 vs Both vs Neither"
echo "  - Expression comparison: Upregulated vs Downregulated"
echo "  - TES-bound comparison: Up vs Down with binding overlay"
echo ""
echo "Interpretation guide:"
echo "  - Compare GFP vs TES meDIP signals to see methylation changes"
echo "  - Promoter methylation typically correlates with gene repression"
echo "  - Gene body methylation may correlate with active transcription"
echo "  - Compare bound vs unbound genes to assess binding-methylation correlation"
echo ""
echo "Finished: $(date)"
