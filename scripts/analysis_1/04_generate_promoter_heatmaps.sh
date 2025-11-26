#!/bin/bash
#SBATCH --job-name=a1_04_generate_promoter_heatmaps
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --time=4:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --output=logs/04_generate_promoter_heatmaps.out
#SBATCH --error=logs/04_generate_promoter_heatmaps.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it

echo "=========================================="
echo "PROMOTER HEATMAP GENERATION"
echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Started: $(date)"
echo ""
echo "Task: Visualize TES/TEAD1 binding at UP vs DOWN gene promoters"
echo ""

# Set working directory
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_1

# Output directory
OUTDIR="output/04_generate_promoter_heatmaps"
mkdir -p ${OUTDIR}
mkdir -p logs

# BigWig files location (absolute path)
CUTANDTAG_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_CUTandTAG"
BIGWIG_DIR="${CUTANDTAG_DIR}/results/06_bigwig"

# =============================================================================
# STEP 1: PREPARE GENE LISTS FROM INTEGRATIVE ANALYSIS
# =============================================================================

echo ""
echo "=== STEP 1: Preparing Gene Lists ==="
echo ""

# Run R script to create BED files
echo "Running R script for gene list preparation..."
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate r_chipseq_env

Rscript 04_generate_promoter_heatmaps.R

if [ $? -ne 0 ]; then
    echo "ERROR: Gene list preparation failed!"
    exit 1
fi

# =============================================================================
# STEP 2: COMPUTE COVERAGE MATRICES WITH DEEPTOOLS
# =============================================================================

echo ""
echo "=== STEP 2: Computing Coverage Matrices ==="
echo ""

# Switch to deepTools environment
echo "Switching to deepTools environment..."
conda activate tg

# Check if BigWig files exist
if [ ! -d "$BIGWIG_DIR" ]; then
    echo "ERROR: BigWig directory not found: $BIGWIG_DIR"
    exit 1
fi

# Find BigWig files (pattern: TES-1_CPM.bw, TES-2_CPM.bw, etc.)
TES_BIGWIGS=$(ls ${BIGWIG_DIR}/TES-[0-9]*_CPM.bw 2>/dev/null | tr '\n' ' ')
TEAD1_BIGWIGS=$(ls ${BIGWIG_DIR}/TEAD1-[0-9]*_CPM.bw 2>/dev/null | tr '\n' ' ')

if [ -z "$TES_BIGWIGS" ] || [ -z "$TEAD1_BIGWIGS" ]; then
    echo "ERROR: BigWig files not found!"
    echo "TES BigWigs: $TES_BIGWIGS"
    echo "TEAD1 BigWigs: $TEAD1_BIGWIGS"
    exit 1
fi

echo "TES BigWigs: $TES_BIGWIGS"
echo "TEAD1 BigWigs: $TEAD1_BIGWIGS"
echo ""

# Combine BED files for single analysis
cat ${OUTDIR}/TES_UP_promoters.bed > ${OUTDIR}/TES_promoters_combined.bed
cat ${OUTDIR}/TES_DOWN_promoters.bed >> ${OUTDIR}/TES_promoters_combined.bed

cat ${OUTDIR}/TEAD1_UP_promoters.bed > ${OUTDIR}/TEAD1_promoters_combined.bed
cat ${OUTDIR}/TEAD1_DOWN_promoters.bed >> ${OUTDIR}/TEAD1_promoters_combined.bed

# Compute matrix for TES targets
echo "Computing matrix for TES targets..."
computeMatrix reference-point \
    --referencePoint TSS \
    -b 10000 -a 10000 \
    --binSize 10 \
    -R ${OUTDIR}/TES_UP_promoters.bed ${OUTDIR}/TES_DOWN_promoters.bed \
    -S $TES_BIGWIGS $TEAD1_BIGWIGS \
    --skipZeros \
    --missingDataAsZero \
    -o ${OUTDIR}/TES_targets_matrix.gz \
    --outFileSortedRegions ${OUTDIR}/TES_targets_sorted.bed \
    -p 16 \
    2>&1 | tee ${OUTDIR}/computeMatrix_TES.log

if [ $? -ne 0 ]; then
    echo "WARNING: TES matrix computation had issues, continuing..."
fi

# Compute matrix for TEAD1 targets
echo "Computing matrix for TEAD1 targets..."
computeMatrix reference-point \
    --referencePoint TSS \
    -b 10000 -a 10000 \
    --binSize 10 \
    -R ${OUTDIR}/TEAD1_UP_promoters.bed ${OUTDIR}/TEAD1_DOWN_promoters.bed \
    -S $TES_BIGWIGS $TEAD1_BIGWIGS \
    --skipZeros \
    --missingDataAsZero \
    -o ${OUTDIR}/TEAD1_targets_matrix.gz \
    --outFileSortedRegions ${OUTDIR}/TEAD1_targets_sorted.bed \
    -p 16 \
    2>&1 | tee ${OUTDIR}/computeMatrix_TEAD1.log

if [ $? -ne 0 ]; then
    echo "WARNING: TEAD1 matrix computation had issues, continuing..."
fi

# =============================================================================
# STEP 3: GENERATE HEATMAPS
# =============================================================================

echo ""
echo "=== STEP 3: Generating Heatmaps ==="
echo ""

# Heatmap for TES targets
if [ -f "${OUTDIR}/TES_targets_matrix.gz" ]; then
    echo "Creating TES targets heatmap..."
    plotHeatmap -m ${OUTDIR}/TES_targets_matrix.gz \
        -out ${OUTDIR}/TES_targets_heatmap.pdf \
        --colorMap RdYlBu_r \
        --whatToShow 'heatmap and colorbar' \
        --zMin 0 --zMax 10 \
        --heatmapHeight 15 \
        --refPointLabel "TSS" \
        --regionsLabel "TES UP" "TES DOWN" \
        --samplesLabel "TES-1" "TES-2" "TES-3" "TEAD1-1" "TEAD1-2" "TEAD1-3" \
        --plotTitle "TES/TEAD1 Binding at TES Target Promoters (UP vs DOWN)" \
        2>&1 | tee ${OUTDIR}/plotHeatmap_TES.log

    # Profile plot for TES targets
    echo "Creating TES targets profile plot..."
    plotProfile -m ${OUTDIR}/TES_targets_matrix.gz \
        -out ${OUTDIR}/TES_targets_profile.pdf \
        --perGroup \
        --colors "#E31A1C" "#FB9A99" "#FDBF6F" "#1F78B4" "#A6CEE3" "#B2DF8A" \
        --refPointLabel "TSS" \
        --regionsLabel "TES UP" "TES DOWN" \
        --plotTitle "Average TES/TEAD1 Binding Profile at TES Targets" \
        2>&1 | tee ${OUTDIR}/plotProfile_TES.log
fi

# Heatmap for TEAD1 targets
if [ -f "${OUTDIR}/TEAD1_targets_matrix.gz" ]; then
    echo "Creating TEAD1 targets heatmap..."
    plotHeatmap -m ${OUTDIR}/TEAD1_targets_matrix.gz \
        -out ${OUTDIR}/TEAD1_targets_heatmap.pdf \
        --colorMap RdYlBu_r \
        --whatToShow 'heatmap and colorbar' \
        --zMin 0 --zMax 10 \
        --heatmapHeight 15 \
        --refPointLabel "TSS" \
        --regionsLabel "TEAD1 UP" "TEAD1 DOWN" \
        --samplesLabel "TES-1" "TES-2" "TES-3" "TEAD1-1" "TEAD1-2" "TEAD1-3" \
        --plotTitle "TES/TEAD1 Binding at TEAD1 Target Promoters (UP vs DOWN)" \
        2>&1 | tee ${OUTDIR}/plotHeatmap_TEAD1.log

    # Profile plot for TEAD1 targets
    echo "Creating TEAD1 targets profile plot..."
    plotProfile -m ${OUTDIR}/TEAD1_targets_matrix.gz \
        -out ${OUTDIR}/TEAD1_targets_profile.pdf \
        --perGroup \
        --colors "#E31A1C" "#FB9A99" "#FDBF6F" "#1F78B4" "#A6CEE3" "#B2DF8A" \
        --refPointLabel "TSS" \
        --regionsLabel "TEAD1 UP" "TEAD1 DOWN" \
        --plotTitle "Average TES/TEAD1 Binding Profile at TEAD1 Targets" \
        2>&1 | tee ${OUTDIR}/plotProfile_TEAD1.log
fi

# =============================================================================
# COMPLETION
# =============================================================================

echo ""
echo "=========================================="
echo "PROMOTER HEATMAP GENERATION COMPLETE"
echo "=========================================="
echo "Completed: $(date)"
echo ""
echo "Output directory: ${OUTDIR}/"
echo ""
echo "Generated files:"
ls -lh ${OUTDIR}/*.pdf 2>/dev/null || echo "  (Check for PDF files manually)"
echo ""
echo "Key outputs:"
echo "  - TES_targets_heatmap.pdf"
echo "  - TES_targets_profile.pdf"
echo "  - TEAD1_targets_heatmap.pdf"
echo "  - TEAD1_targets_profile.pdf"
