#!/bin/bash
#SBATCH --job-name=a1_29_degs_binding_meth
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --time=2:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --output=logs/29_degs_down_binding_methylation.out
#SBATCH --error=logs/29_degs_down_binding_methylation.err

# =============================================================================
# DEGs DOWN STRATIFIED BY TES/TEAD1 BINDING - METAGENE PROFILES
# =============================================================================
#
# Purpose: Test whether TES-induced methylation requires direct binding
#
# Three gene sets:
#   1. DEGs DOWN WITH TES/TEAD1 binding at promoter (±2kb)
#      → Expected: DIFFERENTIAL MeDIP (TES > GFP)
#
#   2. DEGs DOWN WITHOUT TES/TEAD1 binding at promoter
#      → Expected: SIMILAR MeDIP if methylation requires binding
#
#   3. Random expressed genes (same N) WITHOUT binding
#      → Expected: SIMILAR MeDIP (baseline control)
#
# =============================================================================

echo "=========================================="
echo "DEGs DOWN BINDING-STRATIFIED METAGENE ANALYSIS"
echo "=========================================="
echo "Started: $(date)"
echo ""

cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_1

OUTDIR="output/29_degs_down_binding_methylation"
mkdir -p ${OUTDIR}
mkdir -p logs

# =============================================================================
# STEP 1: CREATE GENE SET BED FILES
# =============================================================================

echo "=== STEP 1: Creating gene set BED files ==="
echo ""

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate r_chipseq_env

Rscript 29_degs_down_binding_methylation.R

# Check if BED files were created
if [ ! -f "${OUTDIR}/DEGs_DOWN_with_binding.bed" ]; then
    echo "ERROR: BED files not created. Check R script output."
    exit 1
fi

echo ""
echo "BED files created successfully."
echo ""

# =============================================================================
# STEP 2: SETUP FOR DEEPTOOLS
# =============================================================================

echo "=== STEP 2: Setting up deepTools ==="

conda activate tg

# BigWig files - Binding (Cut&Tag)
TES_BIND="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_CUTandTAG/results/06_bigwig/TES_comb.bw"
TEAD1_BIND="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_CUTandTAG/results/06_bigwig/TEAD1_comb.bw"

# BigWig files - Methylation (meDIP)
TES_METH="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP/results/05_bigwig/TES_combined_RPKM.bw"
GFP_METH="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP/results/05_bigwig/GFP_combined_RPKM.bw"

# BED files
BED_WITH_BINDING="${OUTDIR}/DEGs_DOWN_with_binding.bed"
BED_WITHOUT_BINDING="${OUTDIR}/DEGs_DOWN_without_binding.bed"
BED_RANDOM_CONTROL="${OUTDIR}/random_control_no_binding.bed"

# Verify all files exist
echo "Verifying input files..."
for f in "$TES_BIND" "$TEAD1_BIND" "$TES_METH" "$GFP_METH" \
         "$BED_WITH_BINDING" "$BED_WITHOUT_BINDING" "$BED_RANDOM_CONTROL"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: File not found: $f"
        exit 1
    fi
done
echo "  All files present."

# Get gene counts
N_WITH=$(wc -l < ${BED_WITH_BINDING})
N_WITHOUT=$(wc -l < ${BED_WITHOUT_BINDING})
N_RANDOM=$(wc -l < ${BED_RANDOM_CONTROL})

echo ""
echo "Gene counts:"
echo "  GROUP 1 - DEGs DOWN with binding: ${N_WITH}"
echo "  GROUP 2 - DEGs DOWN without binding: ${N_WITHOUT}"
echo "  GROUP 3 - Random control (no binding): ${N_RANDOM}"
echo ""

# =============================================================================
# STEP 3: MAIN COMPARISON - ALL 4 SIGNALS
# =============================================================================

echo "=== STEP 3: Creating main comparison (all signals) ==="

computeMatrix scale-regions \
    -S $TES_METH $GFP_METH $TES_BIND $TEAD1_BIND \
    -R ${BED_WITH_BINDING} ${BED_WITHOUT_BINDING} ${BED_RANDOM_CONTROL} \
    --beforeRegionStartLength 10000 \
    --afterRegionStartLength 10000 \
    --regionBodyLength 5000 \
    --binSize 100 \
    --skipZeros \
    --missingDataAsZero \
    -o ${OUTDIR}/main_comparison_matrix.gz \
    -p 16 \
    2>&1 | grep -v "Skipping\|did not match"

plotProfile -m ${OUTDIR}/main_comparison_matrix.gz \
    -out ${OUTDIR}/MAIN_degs_down_binding_comparison.png \
    --perGroup \
    --colors "#7B3294" "#636363" "#E31A1C" "#377EB8" \
    --startLabel "TSS" \
    --endLabel "TES" \
    --samplesLabel "TES meth" "GFP meth" "TES bind" "TEAD1 bind" \
    --regionsLabel "DEGs DOWN + Binding (n=${N_WITH})" \
                   "DEGs DOWN - NO Binding (n=${N_WITHOUT})" \
                   "Random Control (n=${N_RANDOM})" \
    --plotTitle "DEGs DOWN: Binding vs Methylation Profile" \
    --plotHeight 14 \
    --plotWidth 22 \
    --legendLocation "upper-left" \
    --yMin 0 \
    --dpi 300

echo "  Created: MAIN_degs_down_binding_comparison.png"

# Also create heatmap
plotHeatmap -m ${OUTDIR}/main_comparison_matrix.gz \
    -out ${OUTDIR}/MAIN_degs_down_binding_heatmap.png \
    --colorList "white,#762A83" "white,#636363" "white,#B2182B" "white,#2166AC" \
    --zMin 0 0 0 0 \
    --zMax 200 200 5 5 \
    --sortRegions descend \
    --sortUsing mean \
    --sortUsingSamples 3 \
    --whatToShow "heatmap and colorbar" \
    --heatmapHeight 20 \
    --heatmapWidth 3 \
    --startLabel "TSS" \
    --endLabel "TES" \
    --regionsLabel "DEGs DOWN + Binding" "DEGs DOWN - NO Binding" "Random Control" \
    --samplesLabel "TES meth" "GFP meth" "TES bind" "TEAD1 bind" \
    --dpi 300

echo "  Created: MAIN_degs_down_binding_heatmap.png"

# =============================================================================
# STEP 4: METHYLATION-ONLY COMPARISON
# =============================================================================

echo ""
echo "=== STEP 4: Creating methylation-only comparison ==="

computeMatrix scale-regions \
    -S $TES_METH $GFP_METH \
    -R ${BED_WITH_BINDING} ${BED_WITHOUT_BINDING} ${BED_RANDOM_CONTROL} \
    --beforeRegionStartLength 10000 \
    --afterRegionStartLength 10000 \
    --regionBodyLength 5000 \
    --binSize 100 \
    --skipZeros \
    --missingDataAsZero \
    -o ${OUTDIR}/methylation_only_matrix.gz \
    -p 16 \
    2>&1 | grep -v "Skipping\|did not match"

plotProfile -m ${OUTDIR}/methylation_only_matrix.gz \
    -out ${OUTDIR}/METHYLATION_degs_down_comparison.png \
    --perGroup \
    --colors "#7B3294" "#636363" \
    --startLabel "TSS" \
    --endLabel "TES" \
    --samplesLabel "TES meth" "GFP meth" \
    --regionsLabel "DEGs DOWN + Binding (n=${N_WITH})" \
                   "DEGs DOWN - NO Binding (n=${N_WITHOUT})" \
                   "Random Control (n=${N_RANDOM})" \
    --plotTitle "Methylation Profile: DEGs DOWN Stratified by Binding" \
    --plotHeight 12 \
    --plotWidth 20 \
    --legendLocation "lower-right" \
    --yMin 0 \
    --dpi 300

echo "  Created: METHYLATION_degs_down_comparison.png"

# =============================================================================
# STEP 5: BINDING-ONLY COMPARISON
# =============================================================================

echo ""
echo "=== STEP 5: Creating binding-only comparison ==="

computeMatrix scale-regions \
    -S $TES_BIND $TEAD1_BIND \
    -R ${BED_WITH_BINDING} ${BED_WITHOUT_BINDING} ${BED_RANDOM_CONTROL} \
    --beforeRegionStartLength 10000 \
    --afterRegionStartLength 10000 \
    --regionBodyLength 5000 \
    --binSize 100 \
    --skipZeros \
    --missingDataAsZero \
    -o ${OUTDIR}/binding_only_matrix.gz \
    -p 16 \
    2>&1 | grep -v "Skipping\|did not match"

plotProfile -m ${OUTDIR}/binding_only_matrix.gz \
    -out ${OUTDIR}/BINDING_degs_down_comparison.png \
    --perGroup \
    --colors "#E31A1C" "#377EB8" \
    --startLabel "TSS" \
    --endLabel "TES" \
    --samplesLabel "TES bind" "TEAD1 bind" \
    --regionsLabel "DEGs DOWN + Binding (n=${N_WITH})" \
                   "DEGs DOWN - NO Binding (n=${N_WITHOUT})" \
                   "Random Control (n=${N_RANDOM})" \
    --plotTitle "Binding Profile: DEGs DOWN Stratified by Binding Status" \
    --plotHeight 12 \
    --plotWidth 20 \
    --legendLocation "upper-right" \
    --yMin 0 \
    --dpi 300

echo "  Created: BINDING_degs_down_comparison.png"

# =============================================================================
# STEP 6: QUANTIFY METHYLATION DIFFERENCES
# =============================================================================

echo ""
echo "=== STEP 6: Quantifying methylation differences ==="

conda activate r_chipseq_env

Rscript - << 'RSCRIPT_QUANTIFY'
suppressPackageStartupMessages({
    library(data.table)
    library(jsonlite)
    library(dplyr)
})

setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_1")

OUTPUT_DIR <- "output/29_degs_down_binding_methylation"

cat("=== Quantifying Methylation Differences ===\n\n")

# Read the methylation-only matrix
matrix_file <- file.path(OUTPUT_DIR, "methylation_only_matrix.gz")

# Read header (JSON)
con <- gzfile(matrix_file, "rt")
header_line <- readLines(con, n = 1)
close(con)

header_json <- fromJSON(gsub("^@", "", header_line))

cat("Matrix metadata:\n")
cat(sprintf("  Samples: %s\n", paste(header_json$sample_labels, collapse = ", ")))
cat(sprintf("  Group boundaries: %s\n", paste(header_json$group_boundaries, collapse = ", ")))
cat("\n")

# Read data
mat <- fread(cmd = paste("zcat", matrix_file, "| tail -n +2"), header = FALSE)

# Calculate bins
n_bins <- (header_json$upstream + header_json$body + header_json$downstream) / header_json$`bin size`

# Gene body bins (avoid TSS dip) - using bins 110-150 out of 250 (gene body region)
body_bins <- 110:150

# Column indices
tes_meth_body <- 6 + body_bins
gfp_meth_body <- 6 + n_bins + body_bins

# Group boundaries
group_bounds <- header_json$group_boundaries
n_group1 <- group_bounds[2]
n_group2 <- group_bounds[3] - group_bounds[2]
n_group3 <- nrow(mat) - group_bounds[3]

group1_rows <- 1:n_group1
group2_rows <- (n_group1 + 1):(n_group1 + n_group2)
group3_rows <- (n_group1 + n_group2 + 1):nrow(mat)

cat(sprintf("Groups:\n"))
cat(sprintf("  Group 1 (DEGs DOWN + binding): rows 1-%d (%d genes)\n", n_group1, length(group1_rows)))
cat(sprintf("  Group 2 (DEGs DOWN - no binding): rows %d-%d (%d genes)\n",
            n_group1 + 1, n_group1 + n_group2, length(group2_rows)))
cat(sprintf("  Group 3 (Random control): rows %d-%d (%d genes)\n",
            n_group1 + n_group2 + 1, nrow(mat), length(group3_rows)))
cat("\n")

# Calculate methylation metrics
calculate_meth <- function(rows) {
    tes <- rowMeans(as.matrix(mat[rows, ..tes_meth_body]), na.rm = TRUE)
    gfp <- rowMeans(as.matrix(mat[rows, ..gfp_meth_body]), na.rm = TRUE)
    data.frame(tes_meth = tes, gfp_meth = gfp, diff = tes - gfp)
}

group1_meth <- calculate_meth(group1_rows)
group2_meth <- calculate_meth(group2_rows)
group3_meth <- calculate_meth(group3_rows)

# Summary statistics
cat("Gene Body Methylation (RPKM):\n")
cat("==============================\n\n")

summary_df <- data.frame(
    Group = c("DEGs_DOWN_with_binding", "DEGs_DOWN_without_binding", "Random_control"),
    N = c(length(group1_rows), length(group2_rows), length(group3_rows)),
    TES_meth_mean = c(mean(group1_meth$tes_meth, na.rm = TRUE),
                      mean(group2_meth$tes_meth, na.rm = TRUE),
                      mean(group3_meth$tes_meth, na.rm = TRUE)),
    TES_meth_median = c(median(group1_meth$tes_meth, na.rm = TRUE),
                        median(group2_meth$tes_meth, na.rm = TRUE),
                        median(group3_meth$tes_meth, na.rm = TRUE)),
    GFP_meth_mean = c(mean(group1_meth$gfp_meth, na.rm = TRUE),
                      mean(group2_meth$gfp_meth, na.rm = TRUE),
                      mean(group3_meth$gfp_meth, na.rm = TRUE)),
    GFP_meth_median = c(median(group1_meth$gfp_meth, na.rm = TRUE),
                        median(group2_meth$gfp_meth, na.rm = TRUE),
                        median(group3_meth$gfp_meth, na.rm = TRUE)),
    Diff_mean = c(mean(group1_meth$diff, na.rm = TRUE),
                  mean(group2_meth$diff, na.rm = TRUE),
                  mean(group3_meth$diff, na.rm = TRUE)),
    Diff_median = c(median(group1_meth$diff, na.rm = TRUE),
                    median(group2_meth$diff, na.rm = TRUE),
                    median(group3_meth$diff, na.rm = TRUE))
)

print(summary_df, row.names = FALSE)
cat("\n")

# Statistical tests
cat("Statistical Comparisons (Wilcoxon rank-sum test):\n")
cat("================================================\n\n")

# Compare Group 1 vs Group 2
wilcox_1v2 <- wilcox.test(group1_meth$diff, group2_meth$diff)
cat(sprintf("Group 1 vs Group 2 (DEGs with vs without binding):\n"))
cat(sprintf("  p-value: %.4e\n", wilcox_1v2$p.value))
cat(sprintf("  Interpretation: %s\n\n",
            ifelse(wilcox_1v2$p.value < 0.05, "SIGNIFICANT DIFFERENCE", "No significant difference")))

# Compare Group 1 vs Group 3
wilcox_1v3 <- wilcox.test(group1_meth$diff, group3_meth$diff)
cat(sprintf("Group 1 vs Group 3 (DEGs with binding vs Random control):\n"))
cat(sprintf("  p-value: %.4e\n", wilcox_1v3$p.value))
cat(sprintf("  Interpretation: %s\n\n",
            ifelse(wilcox_1v3$p.value < 0.05, "SIGNIFICANT DIFFERENCE", "No significant difference")))

# Compare Group 2 vs Group 3
wilcox_2v3 <- wilcox.test(group2_meth$diff, group3_meth$diff)
cat(sprintf("Group 2 vs Group 3 (DEGs without binding vs Random control):\n"))
cat(sprintf("  p-value: %.4e\n", wilcox_2v3$p.value))
cat(sprintf("  Interpretation: %s\n\n",
            ifelse(wilcox_2v3$p.value < 0.05, "SIGNIFICANT DIFFERENCE", "No significant difference")))

# Save results
summary_df$wilcox_pvalue_vs_group1 <- c(NA, wilcox_1v2$p.value, wilcox_1v3$p.value)
summary_df$wilcox_pvalue_vs_random <- c(wilcox_1v3$p.value, wilcox_2v3$p.value, NA)

write.csv(summary_df, file.path(OUTPUT_DIR, "methylation_quantification.csv"), row.names = FALSE)
cat("Saved: methylation_quantification.csv\n")

# Save per-gene data
all_meth <- rbind(
    data.frame(group1_meth, group = "DEGs_DOWN_with_binding"),
    data.frame(group2_meth, group = "DEGs_DOWN_without_binding"),
    data.frame(group3_meth, group = "Random_control")
)
write.csv(all_meth, file.path(OUTPUT_DIR, "methylation_per_gene.csv"), row.names = FALSE)
cat("Saved: methylation_per_gene.csv\n")

cat("\n=== Quantification Complete ===\n")
RSCRIPT_QUANTIFY

# =============================================================================
# SUMMARY
# =============================================================================

echo ""
echo "=========================================="
echo "ANALYSIS COMPLETE"
echo "=========================================="
echo "Finished: $(date)"
echo ""
echo "Output directory: ${OUTDIR}"
echo ""
echo "Generated plots:"
ls -lh ${OUTDIR}/*.png 2>/dev/null
echo ""
echo "Generated matrices (for re-plotting):"
ls -lh ${OUTDIR}/*.gz 2>/dev/null
echo ""
echo "Generated data files:"
ls -lh ${OUTDIR}/*.csv 2>/dev/null
echo ""
echo "Key results:"
echo "  - MAIN_degs_down_binding_comparison.png: All 4 signals comparison"
echo "  - METHYLATION_degs_down_comparison.png: TES vs GFP methylation"
echo "  - BINDING_degs_down_comparison.png: TES vs TEAD1 binding"
echo "  - methylation_quantification.csv: Statistical analysis"
echo ""
echo "Expected findings:"
echo "  - Group 1 (DEGs DOWN + binding): Differential MeDIP expected"
echo "  - Group 2 (DEGs DOWN - no binding): Similar MeDIP if binding required"
echo "  - Group 3 (Random control): Baseline MeDIP (no change expected)"
echo ""
