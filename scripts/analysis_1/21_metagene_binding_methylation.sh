#!/bin/bash
#SBATCH --job-name=a1_21_metagene_bind_meth
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --time=4:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --output=logs/21_metagene_binding_methylation.out
#SBATCH --error=logs/21_metagene_binding_methylation.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it

echo "=========================================="
echo "METAGENE PROFILES: BINDING + METHYLATION"
echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Started: $(date)"
echo ""
echo "Task: Generate metagene profiles showing TES/TEAD1 binding"
echo "      and DNA methylation (TES vs GFP) at gene promoters"
echo ""

# Set working directory
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_1

# Output directory
OUTDIR="output/21_metagene_binding_methylation"
mkdir -p ${OUTDIR}
mkdir -p logs

# =============================================================================
# DATA PATHS
# =============================================================================

# Cut&Tag BigWig files (TES/TEAD1 binding)
CUTANDTAG_BIGWIG="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_CUTandTAG/results/06_bigwig"

# meDIP BigWig files (DNA methylation)
MEDIP_BIGWIG="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP/results/05_bigwig"

# =============================================================================
# STEP 1: PREPARE GENE LISTS (R script)
# =============================================================================

echo ""
echo "=== STEP 1: Preparing Gene Lists ==="
echo ""

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate r_chipseq_env

Rscript 21_metagene_binding_methylation.R

if [ $? -ne 0 ]; then
    echo "ERROR: Gene list preparation failed!"
    exit 1
fi

# Verify BED files were created
if [ ! -f "${OUTDIR}/all_expressed_genes_TSS.bed" ] || [ ! -f "${OUTDIR}/DEGs_DOWN_TSS.bed" ]; then
    echo "ERROR: BED files not created!"
    exit 1
fi

echo ""
echo "BED files created successfully"
echo "  - All expressed genes: $(wc -l < ${OUTDIR}/all_expressed_genes_TSS.bed) TSS regions"
echo "  - DEGs DOWN: $(wc -l < ${OUTDIR}/DEGs_DOWN_TSS.bed) TSS regions"

# =============================================================================
# STEP 2: COMPUTE MATRICES WITH DEEPTOOLS
# =============================================================================

echo ""
echo "=== STEP 2: Computing Coverage Matrices ==="
echo ""

# Switch to deepTools environment
conda activate tg

# Define BigWig files
# Cut&Tag (binding) - individual replicates
TES_BINDING="${CUTANDTAG_BIGWIG}/TES-1_CPM.bw ${CUTANDTAG_BIGWIG}/TES-2_CPM.bw ${CUTANDTAG_BIGWIG}/TES-3_CPM.bw"
TEAD1_BINDING="${CUTANDTAG_BIGWIG}/TEAD1-1_CPM.bw ${CUTANDTAG_BIGWIG}/TEAD1-2_CPM.bw ${CUTANDTAG_BIGWIG}/TEAD1-3_CPM.bw"

# meDIP (methylation) - individual replicates
TES_METH="${MEDIP_BIGWIG}/TES-1-IP_RPKM.bw ${MEDIP_BIGWIG}/TES-2-IP_RPKM.bw"
GFP_METH="${MEDIP_BIGWIG}/GFP-1-IP_RPKM.bw ${MEDIP_BIGWIG}/GFP-2-IP_RPKM.bw"

# Verify files exist
echo "Checking BigWig files..."
for bw in $TES_BINDING $TEAD1_BINDING $TES_METH $GFP_METH; do
    if [ ! -f "$bw" ]; then
        echo "ERROR: BigWig file not found: $bw"
        exit 1
    fi
done
echo "✓ All BigWig files found"

# -----------------------------------------------------------------------------
# MATRIX 1: BINDING at ALL EXPRESSED GENES
# -----------------------------------------------------------------------------
echo ""
echo "Computing binding matrix at all expressed genes..."
computeMatrix reference-point \
    --referencePoint TSS \
    -b 10000 -a 10000 \
    --binSize 50 \
    -R ${OUTDIR}/all_expressed_genes_TSS.bed \
    -S $TES_BINDING $TEAD1_BINDING \
    --skipZeros \
    --missingDataAsZero \
    -o ${OUTDIR}/binding_all_genes_matrix.gz \
    -p 16 \
    2>&1 | tee ${OUTDIR}/binding_all_genes_matrix.log

# -----------------------------------------------------------------------------
# MATRIX 2: BINDING at DEGs DOWN
# -----------------------------------------------------------------------------
echo ""
echo "Computing binding matrix at DEGs DOWN..."
computeMatrix reference-point \
    --referencePoint TSS \
    -b 10000 -a 10000 \
    --binSize 50 \
    -R ${OUTDIR}/DEGs_DOWN_TSS.bed \
    -S $TES_BINDING $TEAD1_BINDING \
    --skipZeros \
    --missingDataAsZero \
    -o ${OUTDIR}/binding_DEGs_DOWN_matrix.gz \
    -p 16 \
    2>&1 | tee ${OUTDIR}/binding_DEGs_DOWN_matrix.log

# -----------------------------------------------------------------------------
# MATRIX 3: METHYLATION at ALL EXPRESSED GENES
# -----------------------------------------------------------------------------
echo ""
echo "Computing methylation matrix at all expressed genes..."
computeMatrix reference-point \
    --referencePoint TSS \
    -b 10000 -a 10000 \
    --binSize 50 \
    -R ${OUTDIR}/all_expressed_genes_TSS.bed \
    -S $TES_METH $GFP_METH \
    --skipZeros \
    --missingDataAsZero \
    -o ${OUTDIR}/methylation_all_genes_matrix.gz \
    -p 16 \
    2>&1 | tee ${OUTDIR}/methylation_all_genes_matrix.log

# -----------------------------------------------------------------------------
# MATRIX 4: METHYLATION at DEGs DOWN
# -----------------------------------------------------------------------------
echo ""
echo "Computing methylation matrix at DEGs DOWN..."
computeMatrix reference-point \
    --referencePoint TSS \
    -b 10000 -a 10000 \
    --binSize 50 \
    -R ${OUTDIR}/DEGs_DOWN_TSS.bed \
    -S $TES_METH $GFP_METH \
    --skipZeros \
    --missingDataAsZero \
    -o ${OUTDIR}/methylation_DEGs_DOWN_matrix.gz \
    -p 16 \
    2>&1 | tee ${OUTDIR}/methylation_DEGs_DOWN_matrix.log

# =============================================================================
# STEP 3: GENERATE PROFILE PLOTS
# =============================================================================

echo ""
echo "=== STEP 3: Generating Profile Plots ==="
echo ""

# Color schemes
# Binding: TES (reds/oranges), TEAD1 (blues/greens)
BINDING_COLORS="#E31A1C #FB9A99 #FDBF6F #1F78B4 #A6CEE3 #B2DF8A"

# Methylation: TES (purples), GFP (grays)
METH_COLORS="#7B3294 #C2A5CF #636363 #BDBDBD"

# -----------------------------------------------------------------------------
# PLOT 1: BINDING PROFILE - ALL GENES
# -----------------------------------------------------------------------------
echo "Creating binding profile at all genes..."
plotProfile -m ${OUTDIR}/binding_all_genes_matrix.gz \
    -out ${OUTDIR}/binding_profile_all_genes.pdf \
    --perGroup \
    --colors $BINDING_COLORS \
    --refPointLabel "TSS" \
    --plotTitle "TES/TEAD1 Binding Profile at All Expressed Genes" \
    --samplesLabel "TES-1" "TES-2" "TES-3" "TEAD1-1" "TEAD1-2" "TEAD1-3" \
    --yAxisLabel "Mean CPM" \
    --plotHeight 6 \
    --plotWidth 8 \
    2>&1 | tee ${OUTDIR}/binding_profile_all_genes.log

# -----------------------------------------------------------------------------
# PLOT 2: BINDING PROFILE - DEGs DOWN
# -----------------------------------------------------------------------------
echo "Creating binding profile at DEGs DOWN..."
plotProfile -m ${OUTDIR}/binding_DEGs_DOWN_matrix.gz \
    -out ${OUTDIR}/binding_profile_DEGs_DOWN.pdf \
    --perGroup \
    --colors $BINDING_COLORS \
    --refPointLabel "TSS" \
    --plotTitle "TES/TEAD1 Binding Profile at Downregulated DEGs" \
    --samplesLabel "TES-1" "TES-2" "TES-3" "TEAD1-1" "TEAD1-2" "TEAD1-3" \
    --yAxisLabel "Mean CPM" \
    --plotHeight 6 \
    --plotWidth 8 \
    2>&1 | tee ${OUTDIR}/binding_profile_DEGs_DOWN.log

# -----------------------------------------------------------------------------
# PLOT 3: METHYLATION PROFILE - ALL GENES
# -----------------------------------------------------------------------------
echo "Creating methylation profile at all genes..."
plotProfile -m ${OUTDIR}/methylation_all_genes_matrix.gz \
    -out ${OUTDIR}/methylation_profile_all_genes.pdf \
    --perGroup \
    --colors $METH_COLORS \
    --refPointLabel "TSS" \
    --plotTitle "DNA Methylation Profile at All Expressed Genes (TES vs GFP)" \
    --samplesLabel "TES-1" "TES-2" "GFP-1" "GFP-2" \
    --yAxisLabel "Mean RPKM" \
    --plotHeight 6 \
    --plotWidth 8 \
    2>&1 | tee ${OUTDIR}/methylation_profile_all_genes.log

# -----------------------------------------------------------------------------
# PLOT 4: METHYLATION PROFILE - DEGs DOWN
# -----------------------------------------------------------------------------
echo "Creating methylation profile at DEGs DOWN..."
plotProfile -m ${OUTDIR}/methylation_DEGs_DOWN_matrix.gz \
    -out ${OUTDIR}/methylation_profile_DEGs_DOWN.pdf \
    --perGroup \
    --colors $METH_COLORS \
    --refPointLabel "TSS" \
    --plotTitle "DNA Methylation Profile at Downregulated DEGs (TES vs GFP)" \
    --samplesLabel "TES-1" "TES-2" "GFP-1" "GFP-2" \
    --yAxisLabel "Mean RPKM" \
    --plotHeight 6 \
    --plotWidth 8 \
    2>&1 | tee ${OUTDIR}/methylation_profile_DEGs_DOWN.log

# =============================================================================
# STEP 4: GENERATE COMBINED PANEL PLOTS
# =============================================================================

echo ""
echo "=== STEP 4: Generating Combined Panel Plots ==="
echo ""

# Create a combined matrix for ALL GENES (binding + methylation)
echo "Creating combined plot for all expressed genes..."
computeMatrix reference-point \
    --referencePoint TSS \
    -b 10000 -a 10000 \
    --binSize 50 \
    -R ${OUTDIR}/all_expressed_genes_TSS.bed \
    -S $TES_BINDING $TEAD1_BINDING $TES_METH $GFP_METH \
    --skipZeros \
    --missingDataAsZero \
    -o ${OUTDIR}/combined_all_genes_matrix.gz \
    -p 16

plotProfile -m ${OUTDIR}/combined_all_genes_matrix.gz \
    -out ${OUTDIR}/combined_profile_all_genes.pdf \
    --perGroup \
    --colors "#E31A1C" "#FB9A99" "#FDBF6F" "#1F78B4" "#A6CEE3" "#B2DF8A" "#7B3294" "#C2A5CF" "#636363" "#BDBDBD" \
    --refPointLabel "TSS" \
    --plotTitle "Binding + Methylation at All Expressed Genes" \
    --samplesLabel "TES-1" "TES-2" "TES-3" "TEAD1-1" "TEAD1-2" "TEAD1-3" "mTES-1" "mTES-2" "mGFP-1" "mGFP-2" \
    --plotHeight 8 \
    --plotWidth 10 \
    2>&1 | tee ${OUTDIR}/combined_profile_all_genes.log

# Create combined matrix for DEGs DOWN
echo "Creating combined plot for DEGs DOWN..."
computeMatrix reference-point \
    --referencePoint TSS \
    -b 10000 -a 10000 \
    --binSize 50 \
    -R ${OUTDIR}/DEGs_DOWN_TSS.bed \
    -S $TES_BINDING $TEAD1_BINDING $TES_METH $GFP_METH \
    --skipZeros \
    --missingDataAsZero \
    -o ${OUTDIR}/combined_DEGs_DOWN_matrix.gz \
    -p 16

plotProfile -m ${OUTDIR}/combined_DEGs_DOWN_matrix.gz \
    -out ${OUTDIR}/combined_profile_DEGs_DOWN.pdf \
    --perGroup \
    --colors "#E31A1C" "#FB9A99" "#FDBF6F" "#1F78B4" "#A6CEE3" "#B2DF8A" "#7B3294" "#C2A5CF" "#636363" "#BDBDBD" \
    --refPointLabel "TSS" \
    --plotTitle "Binding + Methylation at Downregulated DEGs" \
    --samplesLabel "TES-1" "TES-2" "TES-3" "TEAD1-1" "TEAD1-2" "TEAD1-3" "mTES-1" "mTES-2" "mGFP-1" "mGFP-2" \
    --plotHeight 8 \
    --plotWidth 10 \
    2>&1 | tee ${OUTDIR}/combined_profile_DEGs_DOWN.log

# =============================================================================
# STEP 5: GENERATE HEATMAPS
# =============================================================================

echo ""
echo "=== STEP 5: Generating Heatmaps ==="
echo ""

# Binding heatmap at DEGs DOWN
echo "Creating binding heatmap at DEGs DOWN..."
plotHeatmap -m ${OUTDIR}/binding_DEGs_DOWN_matrix.gz \
    -out ${OUTDIR}/binding_heatmap_DEGs_DOWN.pdf \
    --colorMap RdYlBu_r \
    --whatToShow 'heatmap and colorbar' \
    --zMin 0 --zMax 5 \
    --heatmapHeight 15 \
    --heatmapWidth 4 \
    --refPointLabel "TSS" \
    --samplesLabel "TES-1" "TES-2" "TES-3" "TEAD1-1" "TEAD1-2" "TEAD1-3" \
    --plotTitle "TES/TEAD1 Binding at Downregulated DEGs" \
    2>&1 | tee ${OUTDIR}/binding_heatmap_DEGs_DOWN.log

# Methylation heatmap at DEGs DOWN
echo "Creating methylation heatmap at DEGs DOWN..."
plotHeatmap -m ${OUTDIR}/methylation_DEGs_DOWN_matrix.gz \
    -out ${OUTDIR}/methylation_heatmap_DEGs_DOWN.pdf \
    --colorMap PuOr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 15 \
    --heatmapWidth 4 \
    --refPointLabel "TSS" \
    --samplesLabel "TES-1" "TES-2" "GFP-1" "GFP-2" \
    --plotTitle "DNA Methylation at Downregulated DEGs" \
    2>&1 | tee ${OUTDIR}/methylation_heatmap_DEGs_DOWN.log

# =============================================================================
# COMPLETION
# =============================================================================

echo ""
echo "=========================================="
echo "METAGENE PROFILES COMPLETE"
echo "=========================================="
echo "Completed: $(date)"
echo ""
echo "Output directory: ${OUTDIR}/"
echo ""
echo "Generated files:"
echo ""
echo "PROFILE PLOTS:"
ls -lh ${OUTDIR}/*profile*.pdf 2>/dev/null
echo ""
echo "HEATMAPS:"
ls -lh ${OUTDIR}/*heatmap*.pdf 2>/dev/null
echo ""
echo "Key outputs:"
echo "  - binding_profile_all_genes.pdf    : TES/TEAD1 binding at all genes"
echo "  - binding_profile_DEGs_DOWN.pdf    : TES/TEAD1 binding at downregulated DEGs"
echo "  - methylation_profile_all_genes.pdf: Methylation at all genes (TES vs GFP)"
echo "  - methylation_profile_DEGs_DOWN.pdf: Methylation at downregulated DEGs"
echo "  - combined_profile_*.pdf           : Binding + Methylation combined"
echo "  - *_heatmap_*.pdf                  : Heatmap visualizations"
