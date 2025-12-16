#!/bin/bash
#SBATCH --job-name=a1_23_genebody
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --time=4:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --output=logs/23_genebody_metagene.out
#SBATCH --error=logs/23_genebody_metagene.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it

echo "=========================================="
echo "GENE BODY METAGENE PROFILES"
echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Started: $(date)"
echo ""
echo "Task: Create metagene profiles with gene body structure"
echo "      (Upstream -2kb -> Gene Body (scaled) -> Downstream +1kb)"
echo ""

# Set working directory
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_1

# Output directory
OUTDIR="output/23_genebody_metagene"
mkdir -p ${OUTDIR}
mkdir -p logs

# =============================================================================
# DATA PATHS
# =============================================================================

# Cut&Tag BigWig files (binding)
CUTANDTAG_BIGWIG="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_CUTandTAG/results/06_bigwig"

# meDIP BigWig files (methylation)
MEDIP_BIGWIG="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP/results/05_bigwig"

# =============================================================================
# STEP 1: PREPARE GENE LISTS (R script)
# =============================================================================

echo ""
echo "=== STEP 1: Preparing Gene Lists ==="
echo ""

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate r_chipseq_env

Rscript 23_genebody_metagene.R

if [ $? -ne 0 ]; then
    echo "ERROR: Gene list preparation failed!"
    exit 1
fi

# Verify BED files were created
if [ ! -f "${OUTDIR}/all_expressed_genes.bed" ] || [ ! -f "${OUTDIR}/DEGs_DOWN_genes.bed" ]; then
    echo "ERROR: BED files not created!"
    exit 1
fi

echo ""
echo "BED files created successfully"

# =============================================================================
# STEP 2: COMPUTE MATRICES WITH DEEPTOOLS (scale-regions mode)
# =============================================================================

echo ""
echo "=== STEP 2: Computing Coverage Matrices (scale-regions mode) ==="
echo ""

# Switch to deepTools environment
conda activate tg

# Define BigWig files
# Binding
TES_BIND="${CUTANDTAG_BIGWIG}/TES_comb.bw"
TEAD1_BIND="${CUTANDTAG_BIGWIG}/TEAD1_comb.bw"

# Methylation
TES_METH="${MEDIP_BIGWIG}/TES_combined_RPKM.bw"
GFP_METH="${MEDIP_BIGWIG}/GFP_combined_RPKM.bw"

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
# BINDING MATRICES
# -----------------------------------------------------------------------------

echo ""
echo "Computing BINDING matrices..."

# All genes - Binding
echo "  Binding at all genes..."
computeMatrix scale-regions \
    -S $TES_BIND $TEAD1_BIND \
    -R ${OUTDIR}/all_expressed_genes.bed \
    --beforeRegionStartLength 2000 \
    --afterRegionStartLength 1000 \
    --regionBodyLength 5000 \
    --binSize 50 \
    --skipZeros \
    --missingDataAsZero \
    -o ${OUTDIR}/binding_all_genes_matrix.gz \
    -p 16 \
    2>&1 | tee ${OUTDIR}/binding_all_genes_matrix.log

# DEGs DOWN - Binding
echo "  Binding at DEGs DOWN..."
computeMatrix scale-regions \
    -S $TES_BIND $TEAD1_BIND \
    -R ${OUTDIR}/DEGs_DOWN_genes.bed \
    --beforeRegionStartLength 2000 \
    --afterRegionStartLength 1000 \
    --regionBodyLength 5000 \
    --binSize 50 \
    --skipZeros \
    --missingDataAsZero \
    -o ${OUTDIR}/binding_DEGs_DOWN_matrix.gz \
    -p 16 \
    2>&1 | tee ${OUTDIR}/binding_DEGs_DOWN_matrix.log

# -----------------------------------------------------------------------------
# METHYLATION MATRICES
# -----------------------------------------------------------------------------

echo ""
echo "Computing METHYLATION matrices..."

# All genes - Methylation
echo "  Methylation at all genes..."
computeMatrix scale-regions \
    -S $TES_METH $GFP_METH \
    -R ${OUTDIR}/all_expressed_genes.bed \
    --beforeRegionStartLength 2000 \
    --afterRegionStartLength 1000 \
    --regionBodyLength 5000 \
    --binSize 50 \
    --skipZeros \
    --missingDataAsZero \
    -o ${OUTDIR}/methylation_all_genes_matrix.gz \
    -p 16 \
    2>&1 | tee ${OUTDIR}/methylation_all_genes_matrix.log

# DEGs DOWN - Methylation
echo "  Methylation at DEGs DOWN..."
computeMatrix scale-regions \
    -S $TES_METH $GFP_METH \
    -R ${OUTDIR}/DEGs_DOWN_genes.bed \
    --beforeRegionStartLength 2000 \
    --afterRegionStartLength 1000 \
    --regionBodyLength 5000 \
    --binSize 50 \
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
BINDING_COLORS="#E31A1C #377EB8"  # TES (red), TEAD1 (blue)
METH_COLORS="#7B3294 #636363"     # TES (purple), GFP (gray)

# -----------------------------------------------------------------------------
# BINDING PROFILES
# -----------------------------------------------------------------------------

echo "Creating binding profile plots..."

# All genes
plotProfile -m ${OUTDIR}/binding_all_genes_matrix.gz \
    -out ${OUTDIR}/binding_all_genes_profile.pdf \
    --perGroup \
    --colors $BINDING_COLORS \
    --startLabel "TSS" \
    --endLabel "TES" \
    --samplesLabel "TES binding" "TEAD1 binding" \
    --plotTitle "TES/TEAD1 Binding Profile - All Expressed Genes" \
    --yAxisLabel "Mean CPM" \
    --plotHeight 8 \
    --plotWidth 12 \
    2>&1 | tee ${OUTDIR}/binding_all_genes_profile.log

# DEGs DOWN
plotProfile -m ${OUTDIR}/binding_DEGs_DOWN_matrix.gz \
    -out ${OUTDIR}/binding_DEGs_DOWN_profile.pdf \
    --perGroup \
    --colors $BINDING_COLORS \
    --startLabel "TSS" \
    --endLabel "TES" \
    --samplesLabel "TES binding" "TEAD1 binding" \
    --plotTitle "TES/TEAD1 Binding Profile - Downregulated DEGs" \
    --yAxisLabel "Mean CPM" \
    --plotHeight 8 \
    --plotWidth 12 \
    2>&1 | tee ${OUTDIR}/binding_DEGs_DOWN_profile.log

echo "  Created binding profile plots"

# -----------------------------------------------------------------------------
# METHYLATION PROFILES
# -----------------------------------------------------------------------------

echo "Creating methylation profile plots..."

# All genes
plotProfile -m ${OUTDIR}/methylation_all_genes_matrix.gz \
    -out ${OUTDIR}/methylation_all_genes_profile.pdf \
    --perGroup \
    --colors $METH_COLORS \
    --startLabel "TSS" \
    --endLabel "TES" \
    --samplesLabel "TES meDIP" "GFP meDIP" \
    --plotTitle "DNA Methylation Profile - All Expressed Genes" \
    --yAxisLabel "Mean RPKM" \
    --plotHeight 8 \
    --plotWidth 12 \
    2>&1 | tee ${OUTDIR}/methylation_all_genes_profile.log

# DEGs DOWN
plotProfile -m ${OUTDIR}/methylation_DEGs_DOWN_matrix.gz \
    -out ${OUTDIR}/methylation_DEGs_DOWN_profile.pdf \
    --perGroup \
    --colors $METH_COLORS \
    --startLabel "TSS" \
    --endLabel "TES" \
    --samplesLabel "TES meDIP" "GFP meDIP" \
    --plotTitle "DNA Methylation Profile - Downregulated DEGs" \
    --yAxisLabel "Mean RPKM" \
    --plotHeight 8 \
    --plotWidth 12 \
    2>&1 | tee ${OUTDIR}/methylation_DEGs_DOWN_profile.log

echo "  Created methylation profile plots"

# =============================================================================
# STEP 4: GENERATE COMBINED COMPARISON PLOTS
# =============================================================================

echo ""
echo "=== STEP 4: Generating Combined Plots ==="
echo ""

# Combined matrix for comparison - Binding
echo "Creating combined binding comparison..."
computeMatrix scale-regions \
    -S $TES_BIND $TEAD1_BIND \
    -R ${OUTDIR}/all_expressed_genes.bed ${OUTDIR}/DEGs_DOWN_genes.bed \
    --beforeRegionStartLength 2000 \
    --afterRegionStartLength 1000 \
    --regionBodyLength 5000 \
    --binSize 50 \
    --skipZeros \
    --missingDataAsZero \
    -o ${OUTDIR}/binding_comparison_matrix.gz \
    -p 16

plotProfile -m ${OUTDIR}/binding_comparison_matrix.gz \
    -out ${OUTDIR}/binding_comparison_profile.pdf \
    --perGroup \
    --colors "#E31A1C" "#377EB8" \
    --startLabel "TSS" \
    --endLabel "TES" \
    --samplesLabel "TES binding" "TEAD1 binding" \
    --regionsLabel "All Genes" "DEGs DOWN" \
    --plotTitle "TES/TEAD1 Binding: All Genes vs DEGs DOWN" \
    --yAxisLabel "Mean CPM" \
    --plotHeight 8 \
    --plotWidth 14 \
    2>&1 | tee ${OUTDIR}/binding_comparison_profile.log

# Combined matrix for comparison - Methylation
echo "Creating combined methylation comparison..."
computeMatrix scale-regions \
    -S $TES_METH $GFP_METH \
    -R ${OUTDIR}/all_expressed_genes.bed ${OUTDIR}/DEGs_DOWN_genes.bed \
    --beforeRegionStartLength 2000 \
    --afterRegionStartLength 1000 \
    --regionBodyLength 5000 \
    --binSize 50 \
    --skipZeros \
    --missingDataAsZero \
    -o ${OUTDIR}/methylation_comparison_matrix.gz \
    -p 16

plotProfile -m ${OUTDIR}/methylation_comparison_matrix.gz \
    -out ${OUTDIR}/methylation_comparison_profile.pdf \
    --perGroup \
    --colors "#7B3294" "#636363" \
    --startLabel "TSS" \
    --endLabel "TES" \
    --samplesLabel "TES meDIP" "GFP meDIP" \
    --regionsLabel "All Genes" "DEGs DOWN" \
    --plotTitle "DNA Methylation: All Genes vs DEGs DOWN" \
    --yAxisLabel "Mean RPKM" \
    --plotHeight 8 \
    --plotWidth 14 \
    2>&1 | tee ${OUTDIR}/methylation_comparison_profile.log

echo "  Created comparison profile plots"

# =============================================================================
# STEP 5: ALL-IN-ONE COMBINED PLOT (Binding + Methylation)
# =============================================================================

echo ""
echo "=== STEP 5: Creating All-in-One Combined Plot ==="
echo ""

# Create matrix with all signals for DEGs DOWN
computeMatrix scale-regions \
    -S $TES_BIND $TEAD1_BIND $TES_METH $GFP_METH \
    -R ${OUTDIR}/DEGs_DOWN_genes.bed \
    --beforeRegionStartLength 2000 \
    --afterRegionStartLength 1000 \
    --regionBodyLength 5000 \
    --binSize 50 \
    --skipZeros \
    --missingDataAsZero \
    -o ${OUTDIR}/combined_DEGs_DOWN_matrix.gz \
    -p 16

plotProfile -m ${OUTDIR}/combined_DEGs_DOWN_matrix.gz \
    -out ${OUTDIR}/combined_DEGs_DOWN_profile.pdf \
    --perGroup \
    --colors "#E31A1C" "#377EB8" "#7B3294" "#636363" \
    --startLabel "TSS" \
    --endLabel "TES" \
    --samplesLabel "TES bind" "TEAD1 bind" "TES meth" "GFP meth" \
    --plotTitle "Binding + Methylation at Downregulated DEGs" \
    --plotHeight 8 \
    --plotWidth 14 \
    2>&1 | tee ${OUTDIR}/combined_DEGs_DOWN_profile.log

echo "  Created combined profile plot"

# =============================================================================
# COMPLETION
# =============================================================================

echo ""
echo "=========================================="
echo "GENE BODY METAGENE PROFILES COMPLETE"
echo "=========================================="
echo "Completed: $(date)"
echo ""
echo "Output directory: ${OUTDIR}/"
echo ""
echo "Generated files:"
echo ""
echo "BINDING PROFILES:"
ls -lh ${OUTDIR}/binding*profile*.pdf 2>/dev/null
echo ""
echo "METHYLATION PROFILES:"
ls -lh ${OUTDIR}/methylation*profile*.pdf 2>/dev/null
echo ""
echo "COMPARISON PROFILES:"
ls -lh ${OUTDIR}/*comparison*profile*.pdf 2>/dev/null
echo ""
echo "COMBINED PROFILE:"
ls -lh ${OUTDIR}/combined*profile*.pdf 2>/dev/null
echo ""
echo "Key outputs:"
echo "  - binding_all_genes_profile.pdf"
echo "  - binding_DEGs_DOWN_profile.pdf"
echo "  - methylation_all_genes_profile.pdf"
echo "  - methylation_DEGs_DOWN_profile.pdf"
echo "  - binding_comparison_profile.pdf"
echo "  - methylation_comparison_profile.pdf"
echo "  - combined_DEGs_DOWN_profile.pdf"
