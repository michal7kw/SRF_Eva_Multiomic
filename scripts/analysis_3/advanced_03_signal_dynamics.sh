#!/bin/bash
#SBATCH --job-name=signal_dynamics_a3_
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --cpus-per-task=16
#SBATCH --time=04:00:00
#SBATCH --partition=workq
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --output=logs/advanced_03_signal_dynamics.out
#SBATCH --error=logs/advanced_03_signal_dynamics.err

echo "=================================================="
echo "Phase 1.3: Binding Site Signal Dynamics"
echo "Start time: $(date)"
echo "=================================================="

# Change to project directory
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top

# Set paths
INPUT_DIR="SRF_Eva_integrated_analysis/scripts/analysis_3/results/01_binding_classification"
OUTPUT_DIR="SRF_Eva_integrated_analysis/scripts/analysis_3/results/03_signal_dynamics"
BIGWIG_DIR="SRF_Eva_CUTandTAG/results/06_bigwig"

# Create output directory
mkdir -p ${OUTPUT_DIR}

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate analysis3_env 

# Check if classification results exist
if [ ! -d "${INPUT_DIR}" ]; then
    echo "ERROR: Binding classification results not found!"
    echo "Please run Phase 1.1 first: sbatch scripts/advanced_01_binding_classification.sh"
    exit 1
fi

echo ""
echo "Step 1: Preparing BED files for each category..."
echo ""

# Categories to analyze - dynamically detect available BED files
echo "Detecting available categories..."
CATEGORIES=()
CATEGORY_LABELS=""
for CAT in TES_unique TEAD1_unique Shared_high Shared_TES_dominant Shared_TEAD1_dominant Shared_equivalent; do
    if [ -f "${INPUT_DIR}/${CAT}.bed" ]; then
        CATEGORIES+=("${CAT}")
        if [ -z "$CATEGORY_LABELS" ]; then
            CATEGORY_LABELS="\"${CAT}\""
        else
            CATEGORY_LABELS="${CATEGORY_LABELS} \"${CAT}\""
        fi
    fi
done
echo "Found categories: ${CATEGORIES[*]}"

# Check if BED files exist
for CATEGORY in "${CATEGORIES[@]}"; do
    BED_FILE="${INPUT_DIR}/${CATEGORY}.bed"
    if [ ! -f "${BED_FILE}" ]; then
        echo "WARNING: ${BED_FILE} not found"
    else
        N_PEAKS=$(wc -l < ${BED_FILE})
        echo "  ${CATEGORY}: ${N_PEAKS} peaks"
    fi
done

echo ""
echo "Step 2: Computing signal matrices with deepTools..."
echo ""

# Define BigWig files (note: files have _CPM suffix)
TES_BIGWIGS="${BIGWIG_DIR}/TES-1_CPM.bw ${BIGWIG_DIR}/TES-2_CPM.bw ${BIGWIG_DIR}/TES-3_CPM.bw"
TEAD1_BIGWIGS="${BIGWIG_DIR}/TEAD1-1_CPM.bw ${BIGWIG_DIR}/TEAD1-2_CPM.bw ${BIGWIG_DIR}/TEAD1-3_CPM.bw"
ALL_BIGWIGS="${TES_BIGWIGS} ${TEAD1_BIGWIGS}"

# Create region file list
REGION_FILES=""
for CATEGORY in "${CATEGORIES[@]}"; do
    BED_FILE="${INPUT_DIR}/${CATEGORY}.bed"
    if [ -f "${BED_FILE}" ]; then
        REGION_FILES="${REGION_FILES} ${BED_FILE}"
    fi
done

echo "Computing matrix for all categories..."
computeMatrix reference-point \
    --referencePoint center \
    -S ${ALL_BIGWIGS} \
    -R ${REGION_FILES} \
    -a 3000 -b 3000 \
    --binSize 10 \
    --missingDataAsZero \
    --numberOfProcessors 16 \
    -o ${OUTPUT_DIR}/signal_matrix_all_categories.mat.gz

echo ""
echo "Step 3: Creating heatmaps and profile plots..."
echo ""

# Heatmap with clustering
echo "  Generating clustered heatmap..."
plotHeatmap \
    -m ${OUTPUT_DIR}/signal_matrix_all_categories.mat.gz \
    -o ${OUTPUT_DIR}/binding_signal_heatmap_clustered.pdf \
    --colorMap RdYlBu_r \
    --whatToShow 'heatmap and colorbar' \
    --kmeans ${#CATEGORIES[@]} \
    --outFileSortedRegions ${OUTPUT_DIR}/clustered_regions.bed \
    --heatmapHeight 20 \
    --heatmapWidth 6 \
    --legendLocation upper-left \
    --refPointLabel "Peak Center" \
    --regionsLabel ${CATEGORY_LABELS}

# Profile plot
echo "  Generating profile plot..."
plotProfile \
    -m ${OUTPUT_DIR}/signal_matrix_all_categories.mat.gz \
    -o ${OUTPUT_DIR}/binding_signal_profiles.pdf \
    --perGroup \
    --plotHeight 15 \
    --plotWidth 20 \
    --refPointLabel "Peak Center" \
    --legendLocation upper-right \
    --regionsLabel ${CATEGORY_LABELS}

# Heatmap without clustering (preserve category order)
echo "  Generating non-clustered heatmap..."

# Create shorter labels to avoid y-axis overlap
SHORT_LABELS=""
for CAT in "${CATEGORIES[@]}"; do
    # Abbreviate long category names
    SHORT_CAT=$(echo "$CAT" | sed 's/Shared_TES_dominant/Shared_TES/g' | \
                              sed 's/Shared_TEAD1_dominant/Shared_TEAD1/g' | \
                              sed 's/Shared_equivalent/Shared_Eq/g' | \
                              sed 's/Shared_high/Shared_Hi/g')
    if [ -z "$SHORT_LABELS" ]; then
        SHORT_LABELS="\"${SHORT_CAT}\""
    else
        SHORT_LABELS="${SHORT_LABELS} \"${SHORT_CAT}\""
    fi
done

plotHeatmap \
    -m ${OUTPUT_DIR}/signal_matrix_all_categories.mat.gz \
    -o ${OUTPUT_DIR}/binding_signal_heatmap_by_category.pdf \
    --colorMap RdYlBu_r \
    --whatToShow 'heatmap and colorbar' \
    --sortRegions keep \
    --heatmapHeight 30 \
    --heatmapWidth 8 \
    --legendLocation upper-left \
    --refPointLabel "Peak Center" \
    --regionsLabel ${SHORT_LABELS}

echo ""
echo "Step 4: Separate analysis for TES-unique and TEAD1-unique peaks..."
echo ""

# TES-unique peaks
if [ -f "${INPUT_DIR}/TES_unique.bed" ]; then
    echo "  Analyzing TES-unique peaks..."
    computeMatrix reference-point \
        --referencePoint center \
        -S ${ALL_BIGWIGS} \
        -R ${INPUT_DIR}/TES_unique.bed \
        -a 3000 -b 3000 \
        --binSize 10 \
        --missingDataAsZero \
        --numberOfProcessors 16 \
        -o ${OUTPUT_DIR}/TES_unique_matrix.mat.gz

    plotHeatmap \
        -m ${OUTPUT_DIR}/TES_unique_matrix.mat.gz \
        -o ${OUTPUT_DIR}/TES_unique_heatmap.pdf \
        --colorMap Reds \
        --whatToShow 'heatmap and colorbar' \
        --sortUsing mean \
        --heatmapHeight 15 \
        --refPointLabel "Peak Center"
fi

# TEAD1-unique peaks
if [ -f "${INPUT_DIR}/TEAD1_unique.bed" ]; then
    echo "  Analyzing TEAD1-unique peaks..."
    computeMatrix reference-point \
        --referencePoint center \
        -S ${ALL_BIGWIGS} \
        -R ${INPUT_DIR}/TEAD1_unique.bed \
        -a 3000 -b 3000 \
        --binSize 10 \
        --missingDataAsZero \
        --numberOfProcessors 16 \
        -o ${OUTPUT_DIR}/TEAD1_unique_matrix.mat.gz

    plotHeatmap \
        -m ${OUTPUT_DIR}/TEAD1_unique_matrix.mat.gz \
        -o ${OUTPUT_DIR}/TEAD1_unique_heatmap.pdf \
        --colorMap Blues \
        --whatToShow 'heatmap and colorbar' \
        --sortUsing mean \
        --heatmapHeight 15 \
        --refPointLabel "Peak Center"
fi

echo ""
echo "Step 5: Quantitative signal analysis with R..."
echo ""

Rscript ./SRF_Eva_integrated_analysis/scripts/analysis_3/advanced_03_signal_dynamics.R

echo ""
echo "=================================================="
echo "Signal dynamics analysis complete: $(date)"
echo "=================================================="
echo ""
echo "Results location: ${OUTPUT_DIR}"
echo ""
echo "Key outputs:"
echo "  - binding_signal_heatmap_clustered.pdf"
echo "  - binding_signal_profiles.pdf"
echo "  - signal_ratio_distributions.pdf"
echo "  - signal_quantification.csv"
echo ""
echo "Next steps:"
echo "  1. Review heatmaps and profile plots"
echo "  2. Check signal quantification statistics"
echo "  3. Proceed to Phase 2 (Expression Analysis)"
