#!/bin/bash
#SBATCH --job-name=a1_30_peaked_metagene
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --time=4:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --output=logs/30_peaked_genes_metagene.out
#SBATCH --error=logs/30_peaked_genes_metagene.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it

echo "=========================================="
echo "PEAKED GENES METAGENE PROFILES"
echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Started: $(date)"
echo ""
echo "Task: Create metagene profiles for genes with Cut&Tag peaks"
echo "      Gene body structure: TSS -> Gene Body (scaled) -> TTS"
echo "      Flanking: +/- 5kb"
echo ""
echo "Version 1: All genes with peaks (no DEG filter)"
echo "Version 2: DEGs (UP/DOWN) with peaks"
echo ""

# Set working directory
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_1

# Output directory
OUTDIR="output/30_peaked_genes_metagene"
mkdir -p ${OUTDIR}
mkdir -p logs

# =============================================================================
# DATA PATHS
# =============================================================================

# Cut&Tag BigWig files (binding) - Combined/mean signal
CUTANDTAG_BIGWIG="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_CUTandTAG/results/06_bigwig"
TES_BIND="${CUTANDTAG_BIGWIG}/TES_comb.bw"
TEAD1_BIND="${CUTANDTAG_BIGWIG}/TEAD1_comb.bw"

# =============================================================================
# STEP 1: PREPARE GENE LISTS (R script)
# =============================================================================

echo ""
echo "=== STEP 1: Preparing Gene Lists ==="
echo ""

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate r_chipseq_env

Rscript 30_peaked_genes_metagene.R

if [ $? -ne 0 ]; then
    echo "ERROR: Gene list preparation failed!"
    exit 1
fi

# Verify BED files were created
echo ""
echo "Checking BED files..."
for bed in TES_peaked_genes.bed TEAD1_peaked_genes.bed \
           TES_DEGs_UP_peaked.bed TES_DEGs_DOWN_peaked.bed \
           TEAD1_DEGs_UP_peaked.bed TEAD1_DEGs_DOWN_peaked.bed; do
    if [ -f "${OUTDIR}/${bed}" ]; then
        count=$(wc -l < "${OUTDIR}/${bed}")
        echo "  ${bed}: ${count} genes"
    else
        echo "  WARNING: ${bed} not found!"
    fi
done

# =============================================================================
# STEP 2: COMPUTE MATRICES WITH DEEPTOOLS (scale-regions mode)
# =============================================================================

echo ""
echo "=== STEP 2: Computing Coverage Matrices (scale-regions mode) ==="
echo ""

# Switch to deepTools environment
conda activate tg

# Check BigWig files exist
echo "Checking BigWig files..."
for bw in $TES_BIND $TEAD1_BIND; do
    if [ ! -f "$bw" ]; then
        echo "ERROR: BigWig file not found: $bw"
        exit 1
    fi
done
echo "All BigWig files found"

# -----------------------------------------------------------------------------
# VERSION 1: ALL PEAKED GENES
# -----------------------------------------------------------------------------

echo ""
echo "=== VERSION 1: All Peaked Genes ==="
echo ""

# TES-peaked genes with TES binding
echo "Computing TES binding at TES-peaked genes..."
computeMatrix scale-regions \
    -S ${TES_BIND} \
    -R ${OUTDIR}/TES_peaked_genes.bed \
    --beforeRegionStartLength 5000 \
    --afterRegionStartLength 5000 \
    --regionBodyLength 5000 \
    --binSize 50 \
    --skipZeros \
    --missingDataAsZero \
    -o ${OUTDIR}/TES_peaked_genes_matrix.gz \
    -p 16 \
    2>&1 | tee ${OUTDIR}/TES_peaked_genes_matrix.log

# TEAD1-peaked genes with TEAD1 binding
echo "Computing TEAD1 binding at TEAD1-peaked genes..."
computeMatrix scale-regions \
    -S ${TEAD1_BIND} \
    -R ${OUTDIR}/TEAD1_peaked_genes.bed \
    --beforeRegionStartLength 5000 \
    --afterRegionStartLength 5000 \
    --regionBodyLength 5000 \
    --binSize 50 \
    --skipZeros \
    --missingDataAsZero \
    -o ${OUTDIR}/TEAD1_peaked_genes_matrix.gz \
    -p 16 \
    2>&1 | tee ${OUTDIR}/TEAD1_peaked_genes_matrix.log

# -----------------------------------------------------------------------------
# VERSION 2: DEGs WITH PEAKS (UP/DOWN separated)
# -----------------------------------------------------------------------------

echo ""
echo "=== VERSION 2: DEGs with Peaks ==="
echo ""

# TES DEGs (UP and DOWN in one matrix for comparison)
echo "Computing TES binding at TES-peaked DEGs..."
computeMatrix scale-regions \
    -S ${TES_BIND} \
    -R ${OUTDIR}/TES_DEGs_UP_peaked.bed ${OUTDIR}/TES_DEGs_DOWN_peaked.bed \
    --beforeRegionStartLength 5000 \
    --afterRegionStartLength 5000 \
    --regionBodyLength 5000 \
    --binSize 50 \
    --skipZeros \
    --missingDataAsZero \
    -o ${OUTDIR}/TES_DEGs_peaked_matrix.gz \
    -p 16 \
    2>&1 | tee ${OUTDIR}/TES_DEGs_peaked_matrix.log

# TEAD1 DEGs (UP and DOWN in one matrix for comparison)
echo "Computing TEAD1 binding at TEAD1-peaked DEGs..."
computeMatrix scale-regions \
    -S ${TEAD1_BIND} \
    -R ${OUTDIR}/TEAD1_DEGs_UP_peaked.bed ${OUTDIR}/TEAD1_DEGs_DOWN_peaked.bed \
    --beforeRegionStartLength 5000 \
    --afterRegionStartLength 5000 \
    --regionBodyLength 5000 \
    --binSize 50 \
    --skipZeros \
    --missingDataAsZero \
    -o ${OUTDIR}/TEAD1_DEGs_peaked_matrix.gz \
    -p 16 \
    2>&1 | tee ${OUTDIR}/TEAD1_DEGs_peaked_matrix.log

# =============================================================================
# STEP 3: GENERATE PROFILE PLOTS
# =============================================================================

echo ""
echo "=== STEP 3: Generating Profile Plots ==="
echo ""

# Color schemes
TES_COLOR="#E31A1C"    # Red for TES
TEAD1_COLOR="#377EB8"  # Blue for TEAD1

# -----------------------------------------------------------------------------
# VERSION 1: ALL PEAKED GENES
# -----------------------------------------------------------------------------

echo "Creating VERSION 1 profile plots..."

# TES binding at TES-peaked genes
plotProfile -m ${OUTDIR}/TES_peaked_genes_matrix.gz \
    -out ${OUTDIR}/TES_peaked_genes_profile.pdf \
    --perGroup \
    --colors ${TES_COLOR} \
    --startLabel "TSS" \
    --endLabel "TTS" \
    --samplesLabel "TES binding" \
    --plotTitle "TES Binding at Genes with TES Peaks" \
    --yAxisLabel "Mean CPM" \
    --plotHeight 8 \
    --plotWidth 12 \
    2>&1 | tee ${OUTDIR}/TES_peaked_genes_profile.log

# TEAD1 binding at TEAD1-peaked genes
plotProfile -m ${OUTDIR}/TEAD1_peaked_genes_matrix.gz \
    -out ${OUTDIR}/TEAD1_peaked_genes_profile.pdf \
    --perGroup \
    --colors ${TEAD1_COLOR} \
    --startLabel "TSS" \
    --endLabel "TTS" \
    --samplesLabel "TEAD1 binding" \
    --plotTitle "TEAD1 Binding at Genes with TEAD1 Peaks" \
    --yAxisLabel "Mean CPM" \
    --plotHeight 8 \
    --plotWidth 12 \
    2>&1 | tee ${OUTDIR}/TEAD1_peaked_genes_profile.log

echo "  Created VERSION 1 profile plots"

# -----------------------------------------------------------------------------
# VERSION 2: DEGs WITH PEAKS
# -----------------------------------------------------------------------------

echo "Creating VERSION 2 profile plots..."

# TES binding at TES-peaked DEGs (UP vs DOWN)
# Note: Need 2 colors for 2 region groups (UP and DOWN)
plotProfile -m ${OUTDIR}/TES_DEGs_peaked_matrix.gz \
    -out ${OUTDIR}/TES_DEGs_peaked_profile.pdf \
    --perGroup \
    --colors "${TES_COLOR}" "${TES_COLOR}" \
    --startLabel "TSS" \
    --endLabel "TTS" \
    --samplesLabel "TES binding" \
    --regionsLabel "DEGs UP" "DEGs DOWN" \
    --plotTitle "TES Binding at DEGs with TES Peaks" \
    --yAxisLabel "Mean CPM" \
    --plotHeight 8 \
    --plotWidth 12 \
    2>&1 | tee ${OUTDIR}/TES_DEGs_peaked_profile.log

# TEAD1 binding at TEAD1-peaked DEGs (UP vs DOWN)
# Note: Need 2 colors for 2 region groups (UP and DOWN)
plotProfile -m ${OUTDIR}/TEAD1_DEGs_peaked_matrix.gz \
    -out ${OUTDIR}/TEAD1_DEGs_peaked_profile.pdf \
    --perGroup \
    --colors "${TEAD1_COLOR}" "${TEAD1_COLOR}" \
    --startLabel "TSS" \
    --endLabel "TTS" \
    --samplesLabel "TEAD1 binding" \
    --regionsLabel "DEGs UP" "DEGs DOWN" \
    --plotTitle "TEAD1 Binding at DEGs with TEAD1 Peaks" \
    --yAxisLabel "Mean CPM" \
    --plotHeight 8 \
    --plotWidth 12 \
    2>&1 | tee ${OUTDIR}/TEAD1_DEGs_peaked_profile.log

echo "  Created VERSION 2 profile plots"

# =============================================================================
# COMPLETION
# =============================================================================

echo ""
echo "=========================================="
echo "PEAKED GENES METAGENE PROFILES COMPLETE"
echo "=========================================="
echo "Completed: $(date)"
echo ""
echo "Output directory: ${OUTDIR}/"
echo ""
echo "Generated files:"
echo ""
echo "VERSION 1 - All Peaked Genes:"
ls -lh ${OUTDIR}/TES_peaked_genes_profile.pdf 2>/dev/null
ls -lh ${OUTDIR}/TEAD1_peaked_genes_profile.pdf 2>/dev/null
echo ""
echo "VERSION 2 - DEGs with Peaks:"
ls -lh ${OUTDIR}/TES_DEGs_peaked_profile.pdf 2>/dev/null
ls -lh ${OUTDIR}/TEAD1_DEGs_peaked_profile.pdf 2>/dev/null
echo ""
echo "Key outputs:"
echo "  VERSION 1:"
echo "    - TES_peaked_genes_profile.pdf (TES binding at all TES-peaked genes)"
echo "    - TEAD1_peaked_genes_profile.pdf (TEAD1 binding at all TEAD1-peaked genes)"
echo "  VERSION 2:"
echo "    - TES_DEGs_peaked_profile.pdf (TES binding at TES-peaked DEGs, UP vs DOWN)"
echo "    - TEAD1_DEGs_peaked_profile.pdf (TEAD1 binding at TEAD1-peaked DEGs, UP vs DOWN)"
echo ""
