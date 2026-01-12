#!/bin/bash
#SBATCH --job-name=a1_31_enhancer_analysis
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --time=2:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --output=logs/31_enhancer_analysis.out
#SBATCH --error=logs/31_enhancer_analysis.err

# =============================================================================
# ENHANCER LEVEL ANALYSIS: BINDING & METHYLATION
# =============================================================================
#
# Purpose: Analyze Binding and Methylation profiles at Enhancer regions (Distal Peaks).
#          Stratified by target gene expression (DEGs DOWN vs Control).
#
# =============================================================================

echo "=========================================="
echo "ENHANCER LEVEL ANALYSIS"
echo "=========================================="
echo "Started: $(date)"
echo ""

cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_1

OUTDIR="output/31_enhancer_analysis"
mkdir -p ${OUTDIR}
mkdir -p logs

# =============================================================================
# STEP 1: CREATE BED FILES (R SCRIPT)
# =============================================================================

echo "=== STEP 1: Creating Enhancer BED files ==="
echo ""

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate r_chipseq_env

Rscript 31_enhancer_analysis.R

# Check if BED files were created
BED_DEGS="${OUTDIR}/Enhancers_DEGs_DOWN.bed"
BED_CONTROL="${OUTDIR}/Enhancers_Random_Control.bed"

if [[ ! -f "$BED_DEGS" || ! -f "$BED_CONTROL" ]]; then
    echo "ERROR: BED files not created. Check R script output."
    exit 1
fi

echo ""
echo "BED files created successfully."
echo ""

# Get counts
N_DEGS=$(wc -l < ${BED_DEGS})
N_CONTROL=$(wc -l < ${BED_CONTROL})

echo "Enhancer Counts:"
echo "  DEGs DOWN Enhancers: ${N_DEGS}"
echo "  Control Enhancers:   ${N_CONTROL}"
echo ""

# =============================================================================
# STEP 2: SETUP DEEPTOOLS
# =============================================================================

echo "=== STEP 2: Setting up deepTools ==="

conda activate tg

# BigWig files
TES_BIND="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_CUTandTAG/results/06_bigwig/TES_comb.bw"
TEAD1_BIND="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_CUTandTAG/results/06_bigwig/TEAD1_comb.bw"
TES_METH="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP/results/05_bigwig/TES_combined_RPKM.bw"
GFP_METH="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP/results/05_bigwig/GFP_combined_RPKM.bw"

# Verify inputs
for f in "$TES_BIND" "$TEAD1_BIND" "$TES_METH" "$GFP_METH"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: File not found: $f"
        exit 1
    fi
done

# =============================================================================
# STEP 3: COMPUTE MATRIX (Centered on Peak)
# =============================================================================

echo "=== STEP 3: Computing Matrix ==="

# Note: reference-point center (peaks are the reference)
computeMatrix reference-point \
    --referencePoint center \
    -S $TES_METH $GFP_METH $TES_BIND $TEAD1_BIND \
    -R ${BED_DEGS} ${BED_CONTROL} \
    --beforeRegionStartLength 5000 \
    --afterRegionStartLength 5000 \
    --binSize 50 \
    --skipZeros \
    --missingDataAsZero \
    -o ${OUTDIR}/enhancer_matrix.gz \
    -p 16 \
    2>&1 | grep -v "Skipping\|did not match"

# =============================================================================
# STEP 4: PLOT PROFILES
# =============================================================================

echo "=== STEP 4: Plotting Profiles ==="

# 1. Main Comparison (All Signals)
plotProfile -m ${OUTDIR}/enhancer_matrix.gz \
    -out ${OUTDIR}/MAIN_Enhancer_Comparison.png \
    --perGroup \
    --colors "#7B3294" "#636363" "#E31A1C" "#377EB8" \
    --refPointLabel "Peak Center" \
    --samplesLabel "TES meth" "GFP meth" "TES bind" "TEAD1 bind" \
    --regionsLabel "Enhancers (DEGs DOWN) (n=${N_DEGS})" \
                   "Enhancers (Control) (n=${N_CONTROL})" \
    --plotTitle "Enhancer Profile: Binding & Methylation (DEGs DOWN vs Control)" \
    --plotHeight 14 \
    --plotWidth 16 \
    --legendLocation "upper-left" \
    --yMin 0 \
    --dpi 300

echo "  Created: MAIN_Enhancer_Comparison.png"

# 2. Methylation Only
plotProfile -m ${OUTDIR}/enhancer_matrix.gz \
    -out ${OUTDIR}/METHYLATION_Enhancer_Comparison.png \
    --perGroup \
    --colors "#7B3294" "#636363" \
    --samplesLabel "TES meth" "GFP meth" \
    --regionsLabel "Enhancers (DEGs DOWN)" "Enhancers (Control)" \
    --plotTitle "Methylation at Enhancers" \
    --yMin 0 \
    --dpi 300

echo "  Created: METHYLATION_Enhancer_Comparison.png"

# =============================================================================
# STEP 5: QUANTIFY METHYLATION DIFFERENCES
# =============================================================================

echo "=== STEP 5: Quantifying Methylation via R ==="

conda activate r_chipseq_env

Rscript - << 'RSCRIPT_QUANTIFY'
suppressPackageStartupMessages({
    library(data.table)
    library(jsonlite)
    library(dplyr)
})

setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_1")
OUTPUT_DIR <- "output/31_enhancer_analysis"

# Read Matrix
matrix_file <- file.path(OUTPUT_DIR, "enhancer_matrix.gz")
con <- gzfile(matrix_file, "rt")
header_line <- readLines(con, n = 1)
close(con)
header_json <- fromJSON(gsub("^@", "", header_line))

mat <- fread(cmd = paste("zcat", matrix_file, "| tail -n +2"), header = FALSE)

# Groups
group_bounds <- header_json$group_boundaries
n_g1 <- group_bounds[2]
n_g2 <- group_bounds[3] - group_bounds[2]

g1_rows <- 1:n_g1
g2_rows <- (n_g1 + 1):(n_g1 + n_g2)

# Bins: Center is at 5000bp. Total 10000bp. 50bp bins -> 200 bins.
# Center bins: roughly 90-110 (±500bp from center)
center_bins <- 90:110

# Columns (first 2 samples are meth)
n_bins_per_sample <- (header_json$upstream + header_json$downstream) / header_json$`bin size`
tes_meth_cols <- 6 + center_bins
gfp_meth_cols <- 6 + n_bins_per_sample + center_bins

# Calc Function
calc_diff <- function(rows) {
    tes <- rowMeans(as.matrix(mat[rows, ..tes_meth_cols]), na.rm=TRUE)
    gfp <- rowMeans(as.matrix(mat[rows, ..gfp_meth_cols]), na.rm=TRUE)
    return(tes - gfp)
}

g1_diff <- calc_diff(g1_rows)
g2_diff <- calc_diff(g2_rows)

# Test
res <- wilcox.test(g1_diff, g2_diff)

cat(sprintf("\nMethylation Difference (TES - GFP) at Enhancers (Center ±500bp):\n"))
cat(sprintf("  Enhancers DEGs DOWN (mean diff): %.4f\n", mean(g1_diff, na.rm=TRUE)))
cat(sprintf("  Enhancers Control (mean diff):   %.4f\n", mean(g2_diff, na.rm=TRUE)))
cat(sprintf("  Wilcoxon p-value: %.4e\n", res$p.value))

sink(file.path(OUTPUT_DIR, "statistical_result.txt"))
cat("Methylation Difference (TES - GFP) at Enhancers:\n")
print(res)
sink()

RSCRIPT_QUANTIFY

echo ""
echo "=========================================="
echo "ANALYSIS COMPLETE"
echo "=========================================="
