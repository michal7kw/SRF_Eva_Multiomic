#!/bin/bash
#SBATCH --job-name=metagene_degs
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --mem=48G
#SBATCH --cpus-per-task=8
#SBATCH --time=4:00:00
#SBATCH --output=logs/19_metagene_degs.out
#SBATCH --error=logs/19_metagene_degs.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it

# ============================================================================
# Metagene Analysis for DEG Gene Sets
# ============================================================================
# Purpose: Create metagene profiles overlaying:
#   - TES Cut&Tag binding signal
#   - TEAD1 Cut&Tag binding signal
#   - meDIP GFP methylation signal
#   - meDIP TES methylation signal
#
# Gene sets:
#   - Downregulated DEGs (from RNA-seq)
#   - TES_degs.txt custom gene list
#   - Upregulated DEGs (from RNA-seq)
#
# Output:
#   - Metagene profiles (TSS-centered)
#   - Scale-region profiles (gene body)
#   - Heatmaps with signal overlay
# ============================================================================

set -e
set -o pipefail

cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis

# Create directories
OUTDIR="results/19_metagene_degs"
mkdir -p "$OUTDIR"/{matrices,profiles,heatmaps}
mkdir -p logs

echo "=============================================="
echo "Metagene DEG Analysis"
echo "Started: $(date)"
echo "=============================================="

# ============================================================================
# Step 1: Prepare gene set BED files
# ============================================================================
echo ""
echo "=== Step 1: Preparing gene set BED files ==="

# Activate R environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate r_chipseq_env

# Run R script to create BED files
if [[ ! -d "results/19_metagene_beds" ]] || [[ $(ls results/19_metagene_beds/*.bed 2>/dev/null | wc -l) -lt 2 ]]; then
    echo "Creating gene set BED files..."
    Rscript scripts/19_prepare_geneset_beds.R
else
    echo "Gene set BED files already exist. Skipping R script."
fi

# Switch to deepTools environment
conda activate peak_calling_new

# Verify deepTools is available
if ! command -v computeMatrix &> /dev/null; then
    echo "ERROR: deepTools not found"
    exit 1
fi

# ============================================================================
# Define BigWig files
# ============================================================================

# Cut&Tag BigWig files
TES_BW="../SRF_Eva_CUTandTAG/results/06_bigwig/TES_comb.bw"
TEAD1_BW="../SRF_Eva_CUTandTAG/results/06_bigwig/TEAD1_comb.bw"

# Alternative: use average BigWigs from combined replicates
TES_AVG_BW="../SRF_Eva_CUTandTAG/results/11_combined_replicates_narrow/bigwig/TES_average.bw"
TEAD1_AVG_BW="../SRF_Eva_CUTandTAG/results/11_combined_replicates_narrow/bigwig/TEAD1_average.bw"

# Select best available BigWig
if [[ -f "$TES_AVG_BW" ]]; then
    TES_BW_USE="$TES_AVG_BW"
else
    TES_BW_USE="$TES_BW"
fi

if [[ -f "$TEAD1_AVG_BW" ]]; then
    TEAD1_BW_USE="$TEAD1_AVG_BW"
else
    TEAD1_BW_USE="$TEAD1_BW"
fi

# meDIP BigWig files
MEDIP_DIR="../meDIP/results/05_bigwig"
GFP_MEDIP_BW1="$MEDIP_DIR/GFP-1-IP_RPKM.bw"
GFP_MEDIP_BW2="$MEDIP_DIR/GFP-2-IP_RPKM.bw"
TES_MEDIP_BW1="$MEDIP_DIR/TES-1-IP_RPKM.bw"
TES_MEDIP_BW2="$MEDIP_DIR/TES-2-IP_RPKM.bw"

# Verify files
echo ""
echo "=== Checking BigWig files ==="
for f in "$TES_BW_USE" "$TEAD1_BW_USE" "$GFP_MEDIP_BW1" "$TES_MEDIP_BW1"; do
    if [[ -f "$f" ]]; then
        echo "  Found: $f"
    else
        echo "  WARNING: Not found: $f"
    fi
done

# ============================================================================
# Define BED files
# ============================================================================

BED_DIR="results/19_metagene_beds"

# Gene body BEDs
TES_DEGS_BODY="$BED_DIR/TES_degs_gene_body.bed"
DOWNREG_BODY="$BED_DIR/downregulated_degs_gene_body.bed"
UPREG_BODY="$BED_DIR/upregulated_degs_gene_body.bed"

# Promoter BEDs
TES_DEGS_PROM="$BED_DIR/TES_degs_promoter.bed"
DOWNREG_PROM="$BED_DIR/downregulated_degs_promoter.bed"
UPREG_PROM="$BED_DIR/upregulated_degs_promoter.bed"

# Check BED files
echo ""
echo "=== Checking BED files ==="
for f in "$TES_DEGS_BODY" "$TES_DEGS_PROM" "$DOWNREG_BODY" "$DOWNREG_PROM"; do
    if [[ -f "$f" ]]; then
        count=$(wc -l < "$f")
        echo "  Found: $f ($count regions)"
    else
        echo "  Not found: $f"
    fi
done

# ============================================================================
# Function to create metagene profiles
# ============================================================================

create_metagene_analysis() {
    local name=$1
    local body_bed=$2
    local prom_bed=$3
    local title=$4

    if [[ ! -f "$body_bed" ]] && [[ ! -f "$prom_bed" ]]; then
        echo "  Skipping $name - no BED files found"
        return
    fi

    echo ""
    echo "=== Processing: $name ==="

    # -------------------------------------------------------------------------
    # Reference-point mode: TSS-centered (promoter)
    # -------------------------------------------------------------------------
    if [[ -f "$prom_bed" ]]; then
        local prom_count=$(wc -l < "$prom_bed")
        if [[ $prom_count -lt 10 ]]; then
            echo "  Skipping promoter analysis - too few regions ($prom_count)"
        else
            echo "  Creating TSS-centered matrix ($prom_count genes)..."

            # Matrix with all 4 signals
            computeMatrix reference-point \
                --referencePoint center \
                -b 5000 -a 5000 \
                -R "$prom_bed" \
                -S "$TES_BW_USE" "$TEAD1_BW_USE" "$GFP_MEDIP_BW1" "$TES_MEDIP_BW1" \
                --samplesLabel "TES_binding" "TEAD1_binding" "GFP_meDIP" "TES_meDIP" \
                --skipZeros \
                --binSize 50 \
                --numberOfProcessors 8 \
                -o "$OUTDIR/matrices/${name}_TSS_matrix.gz" \
                2>/dev/null || echo "  WARNING: TSS matrix computation failed"

            if [[ -f "$OUTDIR/matrices/${name}_TSS_matrix.gz" ]]; then
                # Profile plot
                echo "  Creating TSS profile..."
                plotProfile \
                    -m "$OUTDIR/matrices/${name}_TSS_matrix.gz" \
                    -o "$OUTDIR/profiles/${name}_TSS_profile.pdf" \
                    --plotTitle "$title (TSS ± 5kb)" \
                    --perGroup \
                    --colors darkblue darkred green purple \
                    --legendLocation upper-right \
                    --yAxisLabel "Signal" \
                    2>/dev/null

                # Heatmap
                echo "  Creating TSS heatmap..."
                plotHeatmap \
                    -m "$OUTDIR/matrices/${name}_TSS_matrix.gz" \
                    -o "$OUTDIR/heatmaps/${name}_TSS_heatmap.pdf" \
                    --colorList "white,darkblue" "white,darkred" "white,green" "white,purple" \
                    --sortRegions descend \
                    --sortUsing mean \
                    --sortUsingSamples 1 2 \
                    --heatmapHeight 15 \
                    --heatmapWidth 3 \
                    --xAxisLabel "Distance from TSS (bp)" \
                    --refPointLabel "TSS" \
                    --plotTitle "$title (TSS ± 5kb)" \
                    2>/dev/null
            fi
        fi
    fi

    # -------------------------------------------------------------------------
    # Scale-regions mode: Gene body
    # -------------------------------------------------------------------------
    if [[ -f "$body_bed" ]]; then
        local body_count=$(wc -l < "$body_bed")
        if [[ $body_count -lt 10 ]]; then
            echo "  Skipping gene body analysis - too few regions ($body_count)"
        else
            echo "  Creating scaled gene body matrix ($body_count genes)..."

            # Matrix with all 4 signals
            computeMatrix scale-regions \
                -b 2000 -a 2000 \
                --regionBodyLength 5000 \
                -R "$body_bed" \
                -S "$TES_BW_USE" "$TEAD1_BW_USE" "$GFP_MEDIP_BW1" "$TES_MEDIP_BW1" \
                --samplesLabel "TES_binding" "TEAD1_binding" "GFP_meDIP" "TES_meDIP" \
                --skipZeros \
                --binSize 50 \
                --numberOfProcessors 8 \
                -o "$OUTDIR/matrices/${name}_genebody_matrix.gz" \
                2>/dev/null || echo "  WARNING: Gene body matrix computation failed"

            if [[ -f "$OUTDIR/matrices/${name}_genebody_matrix.gz" ]]; then
                # Profile plot
                echo "  Creating gene body profile..."
                plotProfile \
                    -m "$OUTDIR/matrices/${name}_genebody_matrix.gz" \
                    -o "$OUTDIR/profiles/${name}_genebody_profile.pdf" \
                    --plotTitle "$title (Gene Body)" \
                    --perGroup \
                    --colors darkblue darkred green purple \
                    --legendLocation upper-right \
                    --yAxisLabel "Signal" \
                    --startLabel "TSS" \
                    --endLabel "TES" \
                    2>/dev/null

                # Heatmap
                echo "  Creating gene body heatmap..."
                plotHeatmap \
                    -m "$OUTDIR/matrices/${name}_genebody_matrix.gz" \
                    -o "$OUTDIR/heatmaps/${name}_genebody_heatmap.pdf" \
                    --colorList "white,darkblue" "white,darkred" "white,green" "white,purple" \
                    --sortRegions descend \
                    --sortUsing mean \
                    --sortUsingSamples 1 2 \
                    --heatmapHeight 15 \
                    --heatmapWidth 3 \
                    --startLabel "TSS" \
                    --endLabel "TES" \
                    --plotTitle "$title (Gene Body)" \
                    2>/dev/null
            fi
        fi
    fi

    echo "  Done with $name"
}

# ============================================================================
# Process each gene set
# ============================================================================

echo ""
echo "=== Creating metagene analyses ==="

# TES_degs gene list
create_metagene_analysis "TES_degs" "$TES_DEGS_BODY" "$TES_DEGS_PROM" "TES Target Genes"

# Downregulated DEGs
create_metagene_analysis "downregulated_degs" "$DOWNREG_BODY" "$DOWNREG_PROM" "Downregulated DEGs"

# Upregulated DEGs
create_metagene_analysis "upregulated_degs" "$UPREG_BODY" "$UPREG_PROM" "Upregulated DEGs"

# ============================================================================
# Comparative analysis: Downreg vs Upreg
# ============================================================================

echo ""
echo "=== Creating comparative analysis ==="

if [[ -f "$DOWNREG_PROM" ]] && [[ -f "$UPREG_PROM" ]]; then

    downreg_count=$(wc -l < "$DOWNREG_PROM")
    upreg_count=$(wc -l < "$UPREG_PROM")

    if [[ $downreg_count -ge 10 ]] && [[ $upreg_count -ge 10 ]]; then

        echo "  Creating comparative TSS matrix..."
        computeMatrix reference-point \
            --referencePoint center \
            -b 5000 -a 5000 \
            -R "$DOWNREG_PROM" "$UPREG_PROM" \
            -S "$TES_BW_USE" "$TEAD1_BW_USE" "$GFP_MEDIP_BW1" "$TES_MEDIP_BW1" \
            --samplesLabel "TES_binding" "TEAD1_binding" "GFP_meDIP" "TES_meDIP" \
            --skipZeros \
            --binSize 50 \
            --numberOfProcessors 8 \
            -o "$OUTDIR/matrices/downreg_vs_upreg_TSS_matrix.gz" \
            2>/dev/null

        if [[ -f "$OUTDIR/matrices/downreg_vs_upreg_TSS_matrix.gz" ]]; then
            echo "  Creating comparative profile..."
            plotProfile \
                -m "$OUTDIR/matrices/downreg_vs_upreg_TSS_matrix.gz" \
                -o "$OUTDIR/profiles/downreg_vs_upreg_TSS_profile.pdf" \
                --plotTitle "Downregulated vs Upregulated DEGs (TSS ± 5kb)" \
                --perGroup \
                --colors darkblue darkred green purple \
                --regionsLabel "Downregulated" "Upregulated" \
                --legendLocation upper-right \
                2>/dev/null

            echo "  Creating comparative heatmap..."
            plotHeatmap \
                -m "$OUTDIR/matrices/downreg_vs_upreg_TSS_matrix.gz" \
                -o "$OUTDIR/heatmaps/downreg_vs_upreg_TSS_heatmap.pdf" \
                --colorList "white,darkblue" "white,darkred" "white,green" "white,purple" \
                --sortRegions descend \
                --sortUsing mean \
                --heatmapHeight 15 \
                --heatmapWidth 3 \
                --regionsLabel "Downregulated" "Upregulated" \
                --plotTitle "Downregulated vs Upregulated DEGs" \
                2>/dev/null
        fi
    else
        echo "  Skipping comparative analysis - insufficient regions"
    fi
fi

# ============================================================================
# Binding-only comparative (TES vs TEAD1)
# ============================================================================

echo ""
echo "=== Creating binding comparison ==="

if [[ -f "$DOWNREG_PROM" ]]; then
    echo "  Creating TES vs TEAD1 binding comparison at downregulated genes..."

    computeMatrix reference-point \
        --referencePoint center \
        -b 5000 -a 5000 \
        -R "$DOWNREG_PROM" \
        -S "$TES_BW_USE" "$TEAD1_BW_USE" \
        --samplesLabel "TES" "TEAD1" \
        --skipZeros \
        --binSize 50 \
        --numberOfProcessors 8 \
        -o "$OUTDIR/matrices/binding_comparison_downreg_matrix.gz" \
        2>/dev/null

    if [[ -f "$OUTDIR/matrices/binding_comparison_downreg_matrix.gz" ]]; then
        plotProfile \
            -m "$OUTDIR/matrices/binding_comparison_downreg_matrix.gz" \
            -o "$OUTDIR/profiles/binding_comparison_downreg_profile.pdf" \
            --plotTitle "TES vs TEAD1 Binding at Downregulated DEGs" \
            --perGroup \
            --colors darkblue darkred \
            --legendLocation upper-right \
            2>/dev/null
    fi
fi

# ============================================================================
# Summary
# ============================================================================

echo ""
echo "=============================================="
echo "Metagene Analysis Complete"
echo "=============================================="
echo ""
echo "Output directory: $OUTDIR/"
echo ""
echo "Contents:"
echo "  - matrices/: deepTools matrices"
echo "  - profiles/: Metagene profile plots (PDF)"
echo "  - heatmaps/: Metagene heatmaps (PDF)"
echo ""
echo "Gene sets analyzed:"

for name in "TES_degs" "downregulated_degs" "upregulated_degs"; do
    if [[ -f "$BED_DIR/${name}_gene_body.bed" ]]; then
        count=$(wc -l < "$BED_DIR/${name}_gene_body.bed")
        echo "  - $name: $count genes"
    fi
done

echo ""
echo "Interpretation:"
echo "  - TSS profiles: Signal around transcription start site"
echo "  - Gene body profiles: Signal across entire gene (scaled)"
echo "  - Compare TES vs TEAD1 binding at DEGs"
echo "  - Compare meDIP signals between conditions (GFP vs TES)"
echo ""
echo "Finished: $(date)"
