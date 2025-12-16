#!/bin/bash
#SBATCH --job-name=a1_28_meth_quant
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --time=2:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --output=logs/28_methylation_quantification.out
#SBATCH --error=logs/28_methylation_quantification.err

# =============================================================================
# METHYLATION QUANTIFICATION AND DIRECT TARGETS ANALYSIS
# =============================================================================
#
# Task 1: Quantify methylation difference between TSS-bound vs NOT TSS-bound
# Task 2: Analyze direct TES targets (bound + DE) methylation pattern
#
# =============================================================================

echo "=========================================="
echo "METHYLATION QUANTIFICATION ANALYSIS"
echo "=========================================="
echo "Started: $(date)"
echo ""

cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_1

OUTDIR="output/27_non_targets_extended_flanks"
mkdir -p ${OUTDIR}

# =============================================================================
# TASK 1: QUANTIFY METHYLATION FROM EXISTING MATRIX
# =============================================================================

echo "=== TASK 1: Quantifying methylation from matrix ==="
echo ""

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate r_chipseq_env

Rscript - << 'RSCRIPT_TASK1'
suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(jsonlite)
})

setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_1")

OUTPUT_DIR <- "output/27_non_targets_extended_flanks"

cat("=== Task 1: Quantifying Methylation Difference ===\n\n")

# -----------------------------------------------------------------------------
# 1. Parse the deepTools matrix file
# -----------------------------------------------------------------------------

matrix_file <- file.path(OUTPUT_DIR, "TSS_proximal_matrix.gz")

# Read header (JSON)
con <- gzfile(matrix_file, "rt")
header_line <- readLines(con, n = 1)
close(con)

# Parse JSON header
header_json <- fromJSON(gsub("^@", "", header_line))
cat("Matrix metadata:\n")
cat(sprintf("  Samples: %s\n", paste(header_json$sample_labels, collapse = ", ")))
cat(sprintf("  Upstream: %d bp\n", header_json$upstream))
cat(sprintf("  Downstream: %d bp\n", header_json$downstream))
cat(sprintf("  Body: %d bp\n", header_json$body))
cat(sprintf("  Bin size: %d bp\n", header_json$`bin size`))
cat(sprintf("  Group boundaries: %s\n", paste(header_json$group_boundaries, collapse = ", ")))
cat("\n")

# Calculate number of bins per sample
bins_per_sample <- (header_json$upstream + header_json$body + header_json$downstream) / header_json$`bin size`
cat(sprintf("  Bins per sample: %.0f\n", bins_per_sample))

# Read data (skip header)
mat <- fread(cmd = paste("zcat", matrix_file, "| tail -n +2"), header = FALSE)
cat(sprintf("  Total genes in matrix: %d\n", nrow(mat)))
cat(sprintf("  Total columns: %d\n", ncol(mat)))
cat("\n")

# -----------------------------------------------------------------------------
# 2. Extract methylation columns
# -----------------------------------------------------------------------------

# Data columns start at column 7
# Sample order: TES_combined_RPKM, GFP_combined_RPKM, TES_comb, TEAD1_comb
# Each sample has 250 bins (100 upstream + 50 body + 100 downstream for 10kb+5kb+10kb)

# Recalculate based on actual matrix structure
n_bins <- as.integer(bins_per_sample)  # bins per sample

# Column indices for each sample (1-indexed, data starts at col 7)
tes_meth_cols <- 7:(7 + n_bins - 1)
gfp_meth_cols <- (7 + n_bins):(7 + 2*n_bins - 1)
tes_bind_cols <- (7 + 2*n_bins):(7 + 3*n_bins - 1)
tead1_bind_cols <- (7 + 3*n_bins):(7 + 4*n_bins - 1)

cat(sprintf("Column ranges (1-indexed):\n"))
cat(sprintf("  TES meth: cols %d-%d\n", min(tes_meth_cols), max(tes_meth_cols)))
cat(sprintf("  GFP meth: cols %d-%d\n", min(gfp_meth_cols), max(gfp_meth_cols)))
cat(sprintf("  TES bind: cols %d-%d\n", min(tes_bind_cols), max(tes_bind_cols)))
cat(sprintf("  TEAD1 bind: cols %d-%d\n", min(tead1_bind_cols), max(tead1_bind_cols)))
cat("\n")

# -----------------------------------------------------------------------------
# 3. Define gene body region (avoid TSS dip)
# -----------------------------------------------------------------------------

# With 10kb upstream, 5kb body, 10kb downstream at 100bp bins:
# Total bins per sample = (10000 + 5000 + 10000) / 100 = 250
# Upstream: bins 1-100 (positions -10kb to TSS)
# Body: bins 101-150 (positions TSS to TES, scaled to 5kb)
# Downstream: bins 151-250 (positions TES to +10kb)

# TSS is at bin 100 (end of upstream)
# TES is at bin 150 (end of body)
# Gene body = bins 101-150

# To avoid TSS dip, use bins 110-150 (middle to end of gene body)
body_bins_relative <- 110:150  # Relative to each sample's start

cat(sprintf("Gene body bins (avoiding TSS dip): %d-%d\n", min(body_bins_relative), max(body_bins_relative)))
cat("\n")

# Adjust to actual column indices
tes_meth_body <- tes_meth_cols[body_bins_relative]
gfp_meth_body <- gfp_meth_cols[body_bins_relative]

# -----------------------------------------------------------------------------
# 4. Split by groups and calculate metrics
# -----------------------------------------------------------------------------

# Group boundaries from header: [0, 700, 2800]
# Group 1: rows 1-700 (TSS-Bound)
# Group 2: rows 701-2800 (NOT TSS-Bound)

n_bound <- header_json$group_boundaries[2]  # 700
n_total <- nrow(mat)

tss_bound_rows <- 1:n_bound
not_bound_rows <- (n_bound + 1):n_total

cat(sprintf("Group sizes:\n"))
cat(sprintf("  TSS-Bound: rows 1-%d (%d genes)\n", n_bound, length(tss_bound_rows)))
cat(sprintf("  NOT TSS-Bound: rows %d-%d (%d genes)\n", n_bound + 1, n_total, length(not_bound_rows)))
cat("\n")

# Calculate per-gene gene body methylation
calculate_gene_body_meth <- function(row_indices, tes_cols, gfp_cols) {
    tes_values <- as.matrix(mat[row_indices, ..tes_cols])
    gfp_values <- as.matrix(mat[row_indices, ..gfp_cols])

    # Mean across gene body bins for each gene
    tes_mean <- rowMeans(tes_values, na.rm = TRUE)
    gfp_mean <- rowMeans(gfp_values, na.rm = TRUE)

    data.frame(
        tes_meth = tes_mean,
        gfp_meth = gfp_mean,
        difference = tes_mean - gfp_mean
    )
}

tss_bound_meth <- calculate_gene_body_meth(tss_bound_rows, tes_meth_body, gfp_meth_body)
tss_bound_meth$group <- "TSS-Bound"
tss_bound_meth$gene_name <- mat[tss_bound_rows, V4]

not_bound_meth <- calculate_gene_body_meth(not_bound_rows, tes_meth_body, gfp_meth_body)
not_bound_meth$group <- "NOT TSS-Bound"
not_bound_meth$gene_name <- mat[not_bound_rows, V4]

# -----------------------------------------------------------------------------
# 5. Statistical comparison
# -----------------------------------------------------------------------------

cat("=== RESULTS ===\n\n")

# Summary statistics
summary_stats <- rbind(
    data.frame(
        Group = "TSS-Bound",
        N = length(tss_bound_rows),
        TES_meth_mean = mean(tss_bound_meth$tes_meth, na.rm = TRUE),
        TES_meth_median = median(tss_bound_meth$tes_meth, na.rm = TRUE),
        GFP_meth_mean = mean(tss_bound_meth$gfp_meth, na.rm = TRUE),
        GFP_meth_median = median(tss_bound_meth$gfp_meth, na.rm = TRUE),
        Diff_mean = mean(tss_bound_meth$difference, na.rm = TRUE),
        Diff_median = median(tss_bound_meth$difference, na.rm = TRUE),
        Diff_SD = sd(tss_bound_meth$difference, na.rm = TRUE)
    ),
    data.frame(
        Group = "NOT TSS-Bound",
        N = length(not_bound_rows),
        TES_meth_mean = mean(not_bound_meth$tes_meth, na.rm = TRUE),
        TES_meth_median = median(not_bound_meth$tes_meth, na.rm = TRUE),
        GFP_meth_mean = mean(not_bound_meth$gfp_meth, na.rm = TRUE),
        GFP_meth_median = median(not_bound_meth$gfp_meth, na.rm = TRUE),
        Diff_mean = mean(not_bound_meth$difference, na.rm = TRUE),
        Diff_median = median(not_bound_meth$difference, na.rm = TRUE),
        Diff_SD = sd(not_bound_meth$difference, na.rm = TRUE)
    )
)

cat("Gene Body Methylation (RPKM):\n")
cat("------------------------------\n")
print(summary_stats, row.names = FALSE)
cat("\n")

# Wilcoxon test comparing TES-GFP differences between groups
wilcox_result <- wilcox.test(tss_bound_meth$difference, not_bound_meth$difference)
cat(sprintf("Wilcoxon rank-sum test (comparing TES-GFP difference between groups):\n"))
cat(sprintf("  W statistic: %.0f\n", wilcox_result$statistic))
cat(sprintf("  p-value: %.2e\n", wilcox_result$p.value))
cat("\n")

# T-test for comparison
ttest_result <- t.test(tss_bound_meth$difference, not_bound_meth$difference)
cat(sprintf("T-test (comparing TES-GFP difference between groups):\n"))
cat(sprintf("  t statistic: %.3f\n", ttest_result$statistic))
cat(sprintf("  p-value: %.2e\n", ttest_result$p.value))
cat(sprintf("  95%% CI: [%.4f, %.4f]\n", ttest_result$conf.int[1], ttest_result$conf.int[2]))
cat("\n")

# Effect size
effect_size <- (mean(not_bound_meth$difference, na.rm = TRUE) - mean(tss_bound_meth$difference, na.rm = TRUE)) /
               sd(c(tss_bound_meth$difference, not_bound_meth$difference), na.rm = TRUE)
cat(sprintf("Effect size (Cohen's d): %.3f\n", effect_size))
cat("\n")

# -----------------------------------------------------------------------------
# 6. Save results
# -----------------------------------------------------------------------------

# Combined per-gene data
all_meth <- rbind(tss_bound_meth, not_bound_meth)
write.csv(all_meth, file.path(OUTPUT_DIR, "methylation_per_gene.csv"), row.names = FALSE)
cat(sprintf("Saved per-gene data: methylation_per_gene.csv\n"))

# Summary statistics
summary_stats$Wilcoxon_pvalue <- wilcox_result$p.value
summary_stats$Ttest_pvalue <- ttest_result$p.value
summary_stats$Effect_size <- effect_size
write.csv(summary_stats, file.path(OUTPUT_DIR, "methylation_quantification.csv"), row.names = FALSE)
cat(sprintf("Saved summary: methylation_quantification.csv\n"))

cat("\n=== Task 1 Complete ===\n\n")
RSCRIPT_TASK1

# =============================================================================
# TASK 2: ANALYZE DIRECT TES TARGETS
# =============================================================================

echo ""
echo "=== TASK 2: Analyzing Direct TES Targets ==="
echo ""

Rscript - << 'RSCRIPT_TASK2'
suppressPackageStartupMessages({
    library(GenomicRanges)
    library(rtracklayer)
    library(dplyr)
})

setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_1")

OUTPUT_DIR <- "output/27_non_targets_extended_flanks"

cat("=== Task 2: Analyzing Direct TES Targets ===\n\n")

# -----------------------------------------------------------------------------
# 1. Load direct targets
# -----------------------------------------------------------------------------

direct_targets_file <- "output/10_final_integrative_analysis/direct_targets/TES_direct_targets.csv"

if (!file.exists(direct_targets_file)) {
    # Try alternative location
    direct_targets_file <- "output/only_degs/results/10_direct_targets/TES_direct_targets.csv"
}

if (!file.exists(direct_targets_file)) {
    cat("ERROR: Direct targets file not found. Checking available files...\n")
    system("find output -name '*direct*' -type f 2>/dev/null | head -20")
    stop("Cannot find direct targets file")
}

direct_targets <- read.csv(direct_targets_file)
cat(sprintf("Loaded direct targets: %d genes\n", nrow(direct_targets)))
cat(sprintf("Columns: %s\n", paste(colnames(direct_targets), collapse = ", ")))
cat("\n")

# -----------------------------------------------------------------------------
# 2. Load gene annotations (GTF)
# -----------------------------------------------------------------------------

GTF_FILE <- "/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/annotation/gencode.v44.annotation.gtf"
gtf <- rtracklayer::import(GTF_FILE)
genes_gtf <- gtf[gtf$type == "gene"]
genes_df <- as.data.frame(genes_gtf)

genes_df <- genes_df %>%
    filter(gene_type == "protein_coding") %>%
    filter(seqnames %in% paste0("chr", c(1:22, "X", "Y")))
genes_df$gene_id_clean <- gsub("\\..*", "", genes_df$gene_id)

cat(sprintf("Protein-coding genes in GTF: %d\n", nrow(genes_df)))

# -----------------------------------------------------------------------------
# 3. Match direct targets to gene coordinates
# -----------------------------------------------------------------------------

# Clean direct target gene IDs
if ("gene_id" %in% colnames(direct_targets)) {
    direct_targets$gene_id_clean <- gsub("\\..*", "", direct_targets$gene_id)
} else if ("ensembl_id" %in% colnames(direct_targets)) {
    direct_targets$gene_id_clean <- direct_targets$ensembl_id
} else {
    stop("Cannot find gene ID column in direct targets")
}

# Match to GTF
matched <- genes_df[genes_df$gene_id_clean %in% direct_targets$gene_id_clean, ]
cat(sprintf("Direct targets with GTF match: %d\n", nrow(matched)))

# Create BED file for direct targets
direct_bed <- data.frame(
    chr = matched$seqnames,
    start = matched$start - 1,
    end = matched$end,
    name = matched$gene_name,
    score = 0,
    strand = matched$strand
)
direct_bed <- direct_bed[!duplicated(paste(direct_bed$chr, direct_bed$start, direct_bed$end)), ]
direct_bed <- direct_bed[order(direct_bed$chr, direct_bed$start), ]

write.table(direct_bed, file.path(OUTPUT_DIR, "direct_targets.bed"),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
cat(sprintf("Created: direct_targets.bed (%d genes)\n", nrow(direct_bed)))

# -----------------------------------------------------------------------------
# 4. Create indirect targets BED (DE but NOT bound)
# -----------------------------------------------------------------------------

# Load DESeq2 results
DESEQ2_FILE <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_RNA/results/05_deseq2/deseq2_results_TES_vs_GFP.txt"
deseq2 <- read.delim(DESEQ2_FILE)
deseq2$gene_id_clean <- gsub("\\..*", "", deseq2$gene_id)

# DE genes (FDR < 0.05)
de_genes <- deseq2 %>% filter(padj < 0.05) %>% pull(gene_id_clean)
cat(sprintf("Differentially expressed genes: %d\n", length(de_genes)))

# Indirect = DE but NOT in direct targets
indirect_ids <- setdiff(de_genes, direct_targets$gene_id_clean)
cat(sprintf("Indirect targets (DE but not bound): %d\n", length(indirect_ids)))

# Match to GTF
indirect_matched <- genes_df[genes_df$gene_id_clean %in% indirect_ids, ]

indirect_bed <- data.frame(
    chr = indirect_matched$seqnames,
    start = indirect_matched$start - 1,
    end = indirect_matched$end,
    name = indirect_matched$gene_name,
    score = 0,
    strand = indirect_matched$strand
)
indirect_bed <- indirect_bed[!duplicated(paste(indirect_bed$chr, indirect_bed$start, indirect_bed$end)), ]
indirect_bed <- indirect_bed[order(indirect_bed$chr, indirect_bed$start), ]

write.table(indirect_bed, file.path(OUTPUT_DIR, "indirect_targets.bed"),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
cat(sprintf("Created: indirect_targets.bed (%d genes)\n", nrow(indirect_bed)))
cat("\n")

# Save gene counts for downstream analysis
counts <- data.frame(
    category = c("direct_targets", "indirect_targets"),
    n_genes = c(nrow(direct_bed), nrow(indirect_bed))
)
write.csv(counts, file.path(OUTPUT_DIR, "target_counts.csv"), row.names = FALSE)

cat("=== Task 2 Part 1 Complete (BED files created) ===\n\n")
RSCRIPT_TASK2

# -----------------------------------------------------------------------------
# Generate methylation profiles for direct vs indirect targets
# -----------------------------------------------------------------------------

echo "=== Generating methylation profiles ==="

conda activate tg

TES_METH="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP/results/05_bigwig/TES_combined_RPKM.bw"
GFP_METH="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP/results/05_bigwig/GFP_combined_RPKM.bw"
TES_BIND="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_CUTandTAG/results/06_bigwig/TES_comb.bw"
TEAD1_BIND="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_CUTandTAG/results/06_bigwig/TEAD1_comb.bw"

N_DIRECT=$(wc -l < ${OUTDIR}/direct_targets.bed)
N_INDIRECT=$(wc -l < ${OUTDIR}/indirect_targets.bed)

echo "  Direct targets: ${N_DIRECT}"
echo "  Indirect targets: ${N_INDIRECT}"

# Generate matrix
computeMatrix scale-regions \
    -S $TES_METH $GFP_METH $TES_BIND $TEAD1_BIND \
    -R ${OUTDIR}/direct_targets.bed ${OUTDIR}/indirect_targets.bed \
    --beforeRegionStartLength 10000 \
    --afterRegionStartLength 10000 \
    --regionBodyLength 5000 \
    --binSize 100 \
    --skipZeros \
    --missingDataAsZero \
    -o ${OUTDIR}/direct_indirect_matrix.gz \
    -p 8 \
    2>&1 | grep -v "Skipping\|did not match"

# Generate profile plot
plotProfile -m ${OUTDIR}/direct_indirect_matrix.gz \
    -out ${OUTDIR}/direct_vs_indirect_methylation.png \
    --perGroup \
    --colors "#7B3294" "#636363" "#E31A1C" "#377EB8" \
    --startLabel "TSS" \
    --endLabel "TES" \
    --samplesLabel "TES meth" "GFP meth" "TES bind" "TEAD1 bind" \
    --regionsLabel "Direct (n=${N_DIRECT})" "Indirect (n=${N_INDIRECT})" \
    --plotTitle "Direct vs Indirect TES Targets: Methylation Profile" \
    --plotHeight 12 \
    --plotWidth 20 \
    --legendLocation "upper-left" \
    --yMin 0 \
    --dpi 300

echo "  Created: direct_vs_indirect_methylation.png"

# Methylation-only plot
computeMatrix scale-regions \
    -S $TES_METH $GFP_METH \
    -R ${OUTDIR}/direct_targets.bed ${OUTDIR}/indirect_targets.bed \
    --beforeRegionStartLength 10000 \
    --afterRegionStartLength 10000 \
    --regionBodyLength 5000 \
    --binSize 100 \
    --skipZeros \
    --missingDataAsZero \
    -o ${OUTDIR}/direct_indirect_meth_only_matrix.gz \
    -p 8 \
    2>&1 | grep -v "Skipping\|did not match"

plotProfile -m ${OUTDIR}/direct_indirect_meth_only_matrix.gz \
    -out ${OUTDIR}/direct_vs_indirect_methylation_only.png \
    --perGroup \
    --colors "#7B3294" "#636363" \
    --startLabel "TSS" \
    --endLabel "TES" \
    --samplesLabel "TES meth" "GFP meth" \
    --regionsLabel "Direct (n=${N_DIRECT})" "Indirect (n=${N_INDIRECT})" \
    --plotTitle "Direct vs Indirect TES Targets: Methylation Only" \
    --plotHeight 10 \
    --plotWidth 18 \
    --legendLocation "lower-right" \
    --yMin 0 \
    --dpi 300

echo "  Created: direct_vs_indirect_methylation_only.png"

# =============================================================================
# QUANTIFY DIRECT vs INDIRECT METHYLATION
# =============================================================================

echo ""
echo "=== Quantifying direct vs indirect methylation ==="

conda activate r_chipseq_env

Rscript - << 'RSCRIPT_QUANTIFY'
suppressPackageStartupMessages({
    library(data.table)
    library(jsonlite)
})

setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_1")

OUTPUT_DIR <- "output/27_non_targets_extended_flanks"

cat("=== Quantifying Direct vs Indirect Methylation ===\n\n")

matrix_file <- file.path(OUTPUT_DIR, "direct_indirect_meth_only_matrix.gz")

# Read header
con <- gzfile(matrix_file, "rt")
header_line <- readLines(con, n = 1)
close(con)

header_json <- fromJSON(gsub("^@", "", header_line))

# Read data
mat <- fread(cmd = paste("zcat", matrix_file, "| tail -n +2"), header = FALSE)

# Get group boundaries
n_direct <- header_json$group_boundaries[2]
n_total <- nrow(mat)

cat(sprintf("Direct targets: %d\n", n_direct))
cat(sprintf("Indirect targets: %d\n", n_total - n_direct))
cat("\n")

# Calculate bins per sample
n_bins <- (header_json$upstream + header_json$body + header_json$downstream) / header_json$`bin size`

# Gene body bins (avoiding TSS dip)
body_bins <- 110:150

# Column indices
tes_meth_body <- 6 + body_bins
gfp_meth_body <- 6 + n_bins + body_bins

# Calculate per-gene metrics
calculate_meth <- function(rows) {
    tes <- rowMeans(as.matrix(mat[rows, ..tes_meth_body]), na.rm = TRUE)
    gfp <- rowMeans(as.matrix(mat[rows, ..gfp_meth_body]), na.rm = TRUE)
    data.frame(tes_meth = tes, gfp_meth = gfp, diff = tes - gfp)
}

direct_meth <- calculate_meth(1:n_direct)
indirect_meth <- calculate_meth((n_direct + 1):n_total)

# Summary
cat("Gene Body Methylation (RPKM):\n")
cat("------------------------------\n")
cat(sprintf("DIRECT targets (n=%d):\n", nrow(direct_meth)))
cat(sprintf("  TES meth: %.3f (mean), %.3f (median)\n",
            mean(direct_meth$tes_meth, na.rm = TRUE), median(direct_meth$tes_meth, na.rm = TRUE)))
cat(sprintf("  GFP meth: %.3f (mean), %.3f (median)\n",
            mean(direct_meth$gfp_meth, na.rm = TRUE), median(direct_meth$gfp_meth, na.rm = TRUE)))
cat(sprintf("  TES-GFP diff: %.3f (mean), %.3f (median)\n",
            mean(direct_meth$diff, na.rm = TRUE), median(direct_meth$diff, na.rm = TRUE)))
cat("\n")

cat(sprintf("INDIRECT targets (n=%d):\n", nrow(indirect_meth)))
cat(sprintf("  TES meth: %.3f (mean), %.3f (median)\n",
            mean(indirect_meth$tes_meth, na.rm = TRUE), median(indirect_meth$tes_meth, na.rm = TRUE)))
cat(sprintf("  GFP meth: %.3f (mean), %.3f (median)\n",
            mean(indirect_meth$gfp_meth, na.rm = TRUE), median(indirect_meth$gfp_meth, na.rm = TRUE)))
cat(sprintf("  TES-GFP diff: %.3f (mean), %.3f (median)\n",
            mean(indirect_meth$diff, na.rm = TRUE), median(indirect_meth$diff, na.rm = TRUE)))
cat("\n")

# Statistical test
wilcox <- wilcox.test(direct_meth$diff, indirect_meth$diff)
cat(sprintf("Wilcoxon test (comparing TES-GFP difference):\n"))
cat(sprintf("  p-value: %.2e\n", wilcox$p.value))
cat("\n")

# Save results
summary_df <- data.frame(
    Group = c("Direct", "Indirect"),
    N = c(nrow(direct_meth), nrow(indirect_meth)),
    TES_meth_mean = c(mean(direct_meth$tes_meth, na.rm = TRUE), mean(indirect_meth$tes_meth, na.rm = TRUE)),
    GFP_meth_mean = c(mean(direct_meth$gfp_meth, na.rm = TRUE), mean(indirect_meth$gfp_meth, na.rm = TRUE)),
    Diff_mean = c(mean(direct_meth$diff, na.rm = TRUE), mean(indirect_meth$diff, na.rm = TRUE)),
    Diff_median = c(median(direct_meth$diff, na.rm = TRUE), median(indirect_meth$diff, na.rm = TRUE))
)
summary_df$Wilcoxon_pvalue <- wilcox$p.value

write.csv(summary_df, file.path(OUTPUT_DIR, "direct_targets_methylation_stats.csv"), row.names = FALSE)
cat("Saved: direct_targets_methylation_stats.csv\n")

cat("\n=== Task 2 Complete ===\n")
RSCRIPT_QUANTIFY

# =============================================================================
# SUMMARY
# =============================================================================

echo ""
echo "=========================================="
echo "COMPLETE"
echo "=========================================="
echo "Finished: $(date)"
echo ""
echo "Output files:"
ls -lh ${OUTDIR}/methylation_*.csv ${OUTDIR}/direct_*.csv ${OUTDIR}/direct_*.png ${OUTDIR}/target_counts.csv 2>/dev/null
echo ""
echo "Key results:"
echo "  - methylation_quantification.csv: TSS-bound vs NOT TSS-bound comparison"
echo "  - direct_targets_methylation_stats.csv: Direct vs indirect targets comparison"
echo "  - direct_vs_indirect_methylation.png: Profile plot"
echo ""
