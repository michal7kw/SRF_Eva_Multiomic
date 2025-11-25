#!/bin/bash
#SBATCH --job-name=signal_dynamics
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
conda activate cutntag

# Check if classification results exist
if [ ! -d "${INPUT_DIR}" ]; then
    echo "ERROR: Binding classification results not found!"
    echo "Please run Phase 1.1 first: sbatch scripts/advanced_01_binding_classification.sh"
    exit 1
fi

echo ""
echo "Step 1: Preparing BED files for each category..."
echo ""

# Categories to analyze
CATEGORIES=(
    "TES_unique"
    "TEAD1_unique"
    "Shared_high"
    "Shared_TES_dominant"
    "Shared_TEAD1_dominant"
    "Shared_equivalent"
)

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

# Define BigWig files
TES_BIGWIGS="${BIGWIG_DIR}/TES-1.bw ${BIGWIG_DIR}/TES-2.bw ${BIGWIG_DIR}/TES-3.bw"
TEAD1_BIGWIGS="${BIGWIG_DIR}/TEAD1-1.bw ${BIGWIG_DIR}/TEAD1-2.bw ${BIGWIG_DIR}/TEAD1-3.bw"
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
    --kmeans 6 \
    --outFileSortedRegions ${OUTPUT_DIR}/clustered_regions.bed \
    --heatmapHeight 20 \
    --heatmapWidth 6 \
    --legendLocation upper-left \
    --refPointLabel "Peak Center" \
    --regionsLabel "TES_unique" "TEAD1_unique" "Shared_high" "Shared_TES_dominant" "Shared_TEAD1_dominant" "Shared_equivalent"

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
    --regionsLabel "TES_unique" "TEAD1_unique" "Shared_high" "Shared_TES_dominant" "Shared_TEAD1_dominant" "Shared_equivalent"

# Heatmap without clustering (preserve category order)
echo "  Generating non-clustered heatmap..."
plotHeatmap \
    -m ${OUTPUT_DIR}/signal_matrix_all_categories.mat.gz \
    -o ${OUTPUT_DIR}/binding_signal_heatmap_by_category.pdf \
    --colorMap RdYlBu_r \
    --whatToShow 'heatmap and colorbar' \
    --sortRegions keep \
    --heatmapHeight 20 \
    --heatmapWidth 6 \
    --legendLocation upper-left \
    --refPointLabel "Peak Center" \
    --regionsLabel "TES_unique" "TEAD1_unique" "Shared_high" "Shared_TES_dominant" "Shared_TEAD1_dominant" "Shared_equivalent"

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

conda activate r_chipseq_env

# Create R script for quantitative analysis
cat > ${OUTPUT_DIR}/quantify_signals.R << 'EOF'
#!/usr/bin/env Rscript

library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(ggplot2)
library(ggpubr)

setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top")

INPUT_DIR <- "SRF_Eva_integrated_analysis/scripts/analysis_3/results/01_binding_classification"
OUTPUT_DIR <- "SRF_Eva_integrated_analysis/scripts/analysis_3/results/03_signal_dynamics"

message("Loading binding classification data...")
load(file.path(INPUT_DIR, "binding_classification_data.RData"))

# Calculate additional signal metrics
message("Calculating signal metrics...")

peaks_df$log2_signal_ratio <- log2((peaks_df$tes_signal + 1) / (peaks_df$tead1_signal + 1))
peaks_df$signal_difference <- peaks_df$tes_signal - peaks_df$tead1_signal
peaks_df$signal_sum <- peaks_df$tes_signal + peaks_df$tead1_signal

# Statistical comparisons
message("Performing statistical comparisons...")

# Compare signal ratios between categories
stat_results <- compare_means(
  log2_signal_ratio ~ category,
  data = peaks_df,
  method = "wilcox.test",
  ref.group = "Shared_equivalent"
)

write.csv(stat_results,
          file.path(OUTPUT_DIR, "signal_ratio_statistics.csv"),
          row.names = FALSE)

# Visualizations
message("Creating visualizations...")

# 1. Signal ratio distribution
pdf(file.path(OUTPUT_DIR, "signal_ratio_distributions.pdf"), width = 12, height = 6)
ggplot(peaks_df, aes(x = category, y = log2_signal_ratio, fill = category)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.2, alpha = 0.5, outlier.shape = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_fill_brewer(palette = "Set2") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  labs(
    title = "Signal Ratio Distribution (TES/TEAD1) by Category",
    x = "Category",
    y = "log2(TES Signal / TEAD1 Signal)"
  )
dev.off()

# 2. Signal strength comparison
pdf(file.path(OUTPUT_DIR, "signal_strength_comparison.pdf"), width = 12, height = 8)
p1 <- ggplot(peaks_df, aes(x = category, y = log2(tes_signal + 1), fill = category)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.2, alpha = 0.5) +
  scale_fill_brewer(palette = "Set2") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  labs(title = "TES Signal", x = "", y = "log2(TES Signal + 1)")

p2 <- ggplot(peaks_df, aes(x = category, y = log2(tead1_signal + 1), fill = category)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.2, alpha = 0.5) +
  scale_fill_brewer(palette = "Set2") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  labs(title = "TEAD1 Signal", x = "Category", y = "log2(TEAD1 Signal + 1)")

gridExtra::grid.arrange(p1, p2, ncol = 1)
dev.off()

# 3. Peak width vs signal
pdf(file.path(OUTPUT_DIR, "peak_width_vs_signal.pdf"), width = 10, height = 8)
ggplot(peaks_df, aes(x = peak_width, y = signal_sum, color = category)) +
  geom_point(alpha = 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_brewer(palette = "Set2") +
  theme_classic() +
  labs(
    title = "Peak Width vs Total Signal",
    x = "Peak Width (bp, log scale)",
    y = "Total Signal (TES + TEAD1, log scale)",
    color = "Category"
  )
dev.off()

# Export quantified signals
write.csv(peaks_df,
          file.path(OUTPUT_DIR, "signal_quantification.csv"),
          row.names = FALSE)

message("Signal dynamics analysis complete!")
EOF

Rscript ${OUTPUT_DIR}/quantify_signals.R

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
