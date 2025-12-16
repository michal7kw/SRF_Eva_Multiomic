#!/bin/bash
#SBATCH --job-name=a1_25_genebody_meth
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --time=4:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --output=logs/25_genebody_methylation_DMR_genes.out
#SBATCH --error=logs/25_genebody_methylation_DMR_genes.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it

echo "=========================================="
echo "GENE BODY METHYLATION AT DMR-ASSOCIATED GENES"
echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Started: $(date)"
echo ""
echo "Purpose: Show gene body methylation profiles stratified by DMR status"
echo "         to demonstrate that TES causes methylation at specific genes"
echo ""

# Set working directory
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_1

# Output directory
OUTDIR="output/25_genebody_methylation_DMR_genes"
mkdir -p ${OUTDIR}
mkdir -p logs

# =============================================================================
# DATA PATHS
# =============================================================================

# Cut&Tag BigWig files (binding)
CUTANDTAG_BIGWIG="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_CUTandTAG/results/06_bigwig"
TES_BIND="${CUTANDTAG_BIGWIG}/TES_comb.bw"
TEAD1_BIND="${CUTANDTAG_BIGWIG}/TEAD1_comb.bw"

# meDIP BigWig files (methylation)
MEDIP_BIGWIG="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP/results/05_bigwig"
TES_METH="${MEDIP_BIGWIG}/TES_combined_RPKM.bw"
GFP_METH="${MEDIP_BIGWIG}/GFP_combined_RPKM.bw"

# =============================================================================
# STEP 1: PREPARE GENE LISTS (R script)
# =============================================================================

echo ""
echo "=== STEP 1: Preparing Gene Lists by DMR Status ==="
echo ""

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate r_chipseq_env

Rscript 25_genebody_methylation_DMR_genes.R

if [ $? -ne 0 ]; then
    echo "ERROR: Gene list preparation failed!"
    exit 1
fi

# Verify BED files were created
if [ ! -f "${OUTDIR}/genes_with_hypermethylated_DMRs.bed" ]; then
    echo "ERROR: BED files not created!"
    exit 1
fi

echo ""
echo "BED files created successfully"

# =============================================================================
# STEP 2: COMPUTE MATRICES WITH DEEPTOOLS (scale-regions mode)
# =============================================================================

echo ""
echo "=== STEP 2: Computing Gene Body Coverage Matrices ==="
echo ""

# Switch to deepTools environment
conda activate tg

# Check BigWig files exist
echo "Checking BigWig files..."
for bw in $TES_BIND $TEAD1_BIND $TES_METH $GFP_METH; do
    if [ ! -f "$bw" ]; then
        echo "ERROR: BigWig file not found: $bw"
        exit 1
    fi
done
echo "All BigWig files found"

# -----------------------------------------------------------------------------
# METHYLATION at genes with HYPERMETHYLATED DMRs
# This is the KEY plot showing TES increases methylation
# -----------------------------------------------------------------------------

echo ""
echo "Computing METHYLATION at genes with hypermethylated DMRs..."
computeMatrix scale-regions \
    -S $TES_METH $GFP_METH \
    -R ${OUTDIR}/genes_with_hypermethylated_DMRs.bed \
    --beforeRegionStartLength 3000 \
    --afterRegionStartLength 3000 \
    --regionBodyLength 5000 \
    --binSize 50 \
    --skipZeros \
    --missingDataAsZero \
    -o ${OUTDIR}/methylation_hypermethylated_genes_matrix.gz \
    -p 16 \
    2>&1 | grep -v "did not match"

echo "Creating methylation profile plot..."
plotProfile -m ${OUTDIR}/methylation_hypermethylated_genes_matrix.gz \
    -out ${OUTDIR}/methylation_hypermethylated_genes_profile.pdf \
    --perGroup \
    --colors "#7B3294" "#636363" \
    --startLabel "TSS" \
    --endLabel "TES" \
    --samplesLabel "TES meDIP" "GFP meDIP" \
    --plotTitle "DNA Methylation at Genes with Hypermethylated DMRs" \
    --yAxisLabel "Mean RPKM" \
    --plotHeight 8 \
    --plotWidth 12

echo "Creating methylation heatmap..."
plotHeatmap -m ${OUTDIR}/methylation_hypermethylated_genes_matrix.gz \
    -out ${OUTDIR}/methylation_hypermethylated_genes_heatmap.pdf \
    --colorMap Purples Greys \
    --sortRegions descend \
    --sortUsing mean \
    --sortUsingSamples 1 \
    --heatmapHeight 15 \
    --heatmapWidth 4 \
    --startLabel "TSS" \
    --endLabel "TES" \
    --plotTitle "Methylation at Hypermethylated Genes"

# -----------------------------------------------------------------------------
# METHYLATION at genes with HYPOMETHYLATED DMRs (comparison)
# -----------------------------------------------------------------------------

echo ""
echo "Computing METHYLATION at genes with hypomethylated DMRs..."
computeMatrix scale-regions \
    -S $TES_METH $GFP_METH \
    -R ${OUTDIR}/genes_with_hypomethylated_DMRs.bed \
    --beforeRegionStartLength 3000 \
    --afterRegionStartLength 3000 \
    --regionBodyLength 5000 \
    --binSize 50 \
    --skipZeros \
    --missingDataAsZero \
    -o ${OUTDIR}/methylation_hypomethylated_genes_matrix.gz \
    -p 16 \
    2>&1 | grep -v "did not match"

plotProfile -m ${OUTDIR}/methylation_hypomethylated_genes_matrix.gz \
    -out ${OUTDIR}/methylation_hypomethylated_genes_profile.pdf \
    --perGroup \
    --colors "#7B3294" "#636363" \
    --startLabel "TSS" \
    --endLabel "TES" \
    --samplesLabel "TES meDIP" "GFP meDIP" \
    --plotTitle "DNA Methylation at Genes with Hypomethylated DMRs" \
    --yAxisLabel "Mean RPKM" \
    --plotHeight 8 \
    --plotWidth 12

# -----------------------------------------------------------------------------
# METHYLATION at genes WITHOUT DMRs (control)
# -----------------------------------------------------------------------------

echo ""
echo "Computing METHYLATION at genes without DMRs (control)..."
computeMatrix scale-regions \
    -S $TES_METH $GFP_METH \
    -R ${OUTDIR}/genes_without_DMRs.bed \
    --beforeRegionStartLength 3000 \
    --afterRegionStartLength 3000 \
    --regionBodyLength 5000 \
    --binSize 50 \
    --skipZeros \
    --missingDataAsZero \
    -o ${OUTDIR}/methylation_no_DMR_genes_matrix.gz \
    -p 16 \
    2>&1 | grep -v "did not match"

plotProfile -m ${OUTDIR}/methylation_no_DMR_genes_matrix.gz \
    -out ${OUTDIR}/methylation_no_DMR_genes_profile.pdf \
    --perGroup \
    --colors "#7B3294" "#636363" \
    --startLabel "TSS" \
    --endLabel "TES" \
    --samplesLabel "TES meDIP" "GFP meDIP" \
    --plotTitle "DNA Methylation at Genes without DMRs (Control)" \
    --yAxisLabel "Mean RPKM" \
    --plotHeight 8 \
    --plotWidth 12

# -----------------------------------------------------------------------------
# COMPARISON: Hyper vs Hypo vs No-DMR genes (same plot)
# -----------------------------------------------------------------------------

echo ""
echo "Creating COMPARISON methylation matrix..."
computeMatrix scale-regions \
    -S $TES_METH $GFP_METH \
    -R ${OUTDIR}/genes_with_hypermethylated_DMRs.bed \
       ${OUTDIR}/genes_with_hypomethylated_DMRs.bed \
       ${OUTDIR}/genes_without_DMRs.bed \
    --beforeRegionStartLength 3000 \
    --afterRegionStartLength 3000 \
    --regionBodyLength 5000 \
    --binSize 50 \
    --skipZeros \
    --missingDataAsZero \
    -o ${OUTDIR}/methylation_comparison_matrix.gz \
    -p 16 \
    2>&1 | grep -v "did not match"

echo "Creating comparison profile..."
plotProfile -m ${OUTDIR}/methylation_comparison_matrix.gz \
    -out ${OUTDIR}/methylation_comparison_profile.pdf \
    --perGroup \
    --colors "#7B3294" "#636363" \
    --startLabel "TSS" \
    --endLabel "TES" \
    --samplesLabel "TES meDIP" "GFP meDIP" \
    --regionsLabel "Hypermethylated" "Hypomethylated" "No DMR" \
    --plotTitle "Gene Body Methylation by DMR Status" \
    --yAxisLabel "Mean RPKM" \
    --plotHeight 8 \
    --plotWidth 14 \
    --legendLocation upper-right

# -----------------------------------------------------------------------------
# BINDING at genes with HYPERMETHYLATED DMRs
# Shows TES/TEAD1 binding at genes that gained methylation
# -----------------------------------------------------------------------------

echo ""
echo "Computing BINDING at genes with hypermethylated DMRs..."
computeMatrix scale-regions \
    -S $TES_BIND $TEAD1_BIND \
    -R ${OUTDIR}/genes_with_hypermethylated_DMRs.bed \
    --beforeRegionStartLength 3000 \
    --afterRegionStartLength 3000 \
    --regionBodyLength 5000 \
    --binSize 50 \
    --skipZeros \
    --missingDataAsZero \
    -o ${OUTDIR}/binding_hypermethylated_genes_matrix.gz \
    -p 16 \
    2>&1 | grep -v "Skipping"

plotProfile -m ${OUTDIR}/binding_hypermethylated_genes_matrix.gz \
    -out ${OUTDIR}/binding_hypermethylated_genes_profile.pdf \
    --perGroup \
    --colors "#E31A1C" "#377EB8" \
    --startLabel "TSS" \
    --endLabel "TES" \
    --samplesLabel "TES binding" "TEAD1 binding" \
    --plotTitle "TES/TEAD1 Binding at Genes with Hypermethylated DMRs" \
    --yAxisLabel "Mean CPM" \
    --plotHeight 8 \
    --plotWidth 12

# -----------------------------------------------------------------------------
# BINDING COMPARISON: Hyper vs Hypo vs No-DMR genes
# -----------------------------------------------------------------------------

echo ""
echo "Creating BINDING comparison matrix..."
computeMatrix scale-regions \
    -S $TES_BIND $TEAD1_BIND \
    -R ${OUTDIR}/genes_with_hypermethylated_DMRs.bed \
       ${OUTDIR}/genes_with_hypomethylated_DMRs.bed \
       ${OUTDIR}/genes_without_DMRs.bed \
    --beforeRegionStartLength 3000 \
    --afterRegionStartLength 3000 \
    --regionBodyLength 5000 \
    --binSize 50 \
    --skipZeros \
    --missingDataAsZero \
    -o ${OUTDIR}/binding_comparison_matrix.gz \
    -p 16 \
    2>&1 | grep -v "Skipping"

echo "Creating binding comparison profile..."
plotProfile -m ${OUTDIR}/binding_comparison_matrix.gz \
    -out ${OUTDIR}/binding_comparison_profile.pdf \
    --perGroup \
    --colors "#E31A1C" "#377EB8" \
    --startLabel "TSS" \
    --endLabel "TES" \
    --samplesLabel "TES binding" "TEAD1 binding" \
    --regionsLabel "Hypermethylated" "Hypomethylated" "No DMR" \
    --plotTitle "TES/TEAD1 Binding by Gene DMR Status" \
    --yAxisLabel "Mean CPM" \
    --plotHeight 8 \
    --plotWidth 14 \
    --legendLocation upper-right

# -----------------------------------------------------------------------------
# DEGs DOWN with hypermethylated DMRs (if file exists)
# -----------------------------------------------------------------------------

if [ -f "${OUTDIR}/DEGs_DOWN_hypermethylated.bed" ]; then
    echo ""
    echo "Computing profiles for DEGs DOWN with hypermethylated DMRs..."

    # Methylation
    computeMatrix scale-regions \
        -S $TES_METH $GFP_METH \
        -R ${OUTDIR}/DEGs_DOWN_hypermethylated.bed \
        --beforeRegionStartLength 3000 \
        --afterRegionStartLength 3000 \
        --regionBodyLength 5000 \
        --binSize 50 \
        --skipZeros \
        --missingDataAsZero \
        -o ${OUTDIR}/methylation_DEGs_DOWN_hyper_matrix.gz \
        -p 16 \
        2>&1 | grep -v "did not match"

    plotProfile -m ${OUTDIR}/methylation_DEGs_DOWN_hyper_matrix.gz \
        -out ${OUTDIR}/methylation_DEGs_DOWN_hyper_profile.pdf \
        --perGroup \
        --colors "#7B3294" "#636363" \
        --startLabel "TSS" \
        --endLabel "TES" \
        --samplesLabel "TES meDIP" "GFP meDIP" \
        --plotTitle "Methylation at Downregulated DEGs with Hypermethylated DMRs" \
        --yAxisLabel "Mean RPKM" \
        --plotHeight 8 \
        --plotWidth 12

    # Binding
    computeMatrix scale-regions \
        -S $TES_BIND $TEAD1_BIND \
        -R ${OUTDIR}/DEGs_DOWN_hypermethylated.bed \
        --beforeRegionStartLength 3000 \
        --afterRegionStartLength 3000 \
        --regionBodyLength 5000 \
        --binSize 50 \
        --skipZeros \
        --missingDataAsZero \
        -o ${OUTDIR}/binding_DEGs_DOWN_hyper_matrix.gz \
        -p 16 \
        2>&1 | grep -v "Skipping"

    plotProfile -m ${OUTDIR}/binding_DEGs_DOWN_hyper_matrix.gz \
        -out ${OUTDIR}/binding_DEGs_DOWN_hyper_profile.pdf \
        --perGroup \
        --colors "#E31A1C" "#377EB8" \
        --startLabel "TSS" \
        --endLabel "TES" \
        --samplesLabel "TES binding" "TEAD1 binding" \
        --plotTitle "Binding at Downregulated DEGs with Hypermethylated DMRs" \
        --yAxisLabel "Mean CPM" \
        --plotHeight 8 \
        --plotWidth 12

    # Combined binding + methylation for DEGs DOWN hyper
    computeMatrix scale-regions \
        -S $TES_BIND $TEAD1_BIND $TES_METH $GFP_METH \
        -R ${OUTDIR}/DEGs_DOWN_hypermethylated.bed \
        --beforeRegionStartLength 3000 \
        --afterRegionStartLength 3000 \
        --regionBodyLength 5000 \
        --binSize 50 \
        --skipZeros \
        --missingDataAsZero \
        -o ${OUTDIR}/combined_DEGs_DOWN_hyper_matrix.gz \
        -p 16 \
        2>&1 | grep -v "Skipping\|did not match"

    plotHeatmap -m ${OUTDIR}/combined_DEGs_DOWN_hyper_matrix.gz \
        -out ${OUTDIR}/combined_DEGs_DOWN_hyper_heatmap.pdf \
        --colorList "white,#2166AC" "white,#B2182B" "white,#762A83" "white,#636363" \
        --zMin 0 0 0 0 \
        --zMax 3 3 200 200 \
        --sortRegions descend \
        --sortUsing mean \
        --sortUsingSamples 1 2 \
        --heatmapHeight 15 \
        --heatmapWidth 3 \
        --startLabel "TSS" \
        --endLabel "TES" \
        --plotTitle "DEGs DOWN with Hypermethylation: Binding + Methylation"
fi

# =============================================================================
# COMPLETION
# =============================================================================

echo ""
echo "=========================================="
echo "GENE BODY METHYLATION ANALYSIS COMPLETE"
echo "=========================================="
echo "Completed: $(date)"
echo ""
echo "Output directory: ${OUTDIR}/"
echo ""
echo "Generated files:"
ls -lh ${OUTDIR}/*.pdf 2>/dev/null
echo ""
echo "=========================================="
echo "KEY INTERPRETATION:"
echo "=========================================="
echo ""
echo "These plots show gene body methylation (TSS → gene body → TES)"
echo "stratified by whether genes contain DMRs."
echo ""
echo "Key plots:"
echo "  1. methylation_hypermethylated_genes_profile.pdf"
echo "     - Should show TES meDIP HIGHER than GFP meDIP across gene body"
echo "     - This demonstrates TES increases methylation at these genes"
echo ""
echo "  2. methylation_comparison_profile.pdf"
echo "     - Compares: Hypermethylated vs Hypomethylated vs No-DMR genes"
echo "     - TES-GFP difference should be largest for hypermethylated genes"
echo ""
echo "  3. binding_hypermethylated_genes_profile.pdf"
echo "     - Shows TES/TEAD1 binding at genes that gained methylation"
echo ""
echo "  4. DEGs_DOWN_hyper plots (if generated)"
echo "     - Most direct evidence: genes downregulated by TES AND hypermethylated"
echo ""
