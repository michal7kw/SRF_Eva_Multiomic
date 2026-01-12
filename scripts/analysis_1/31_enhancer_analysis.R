#!/usr/bin/env Rscript
#
# ENHANCER ANALYSIS: DEGs DOWN vs CONTROL
# =============================================================================
#
# Purpose: Define "Enhancer" regions as TES/TEAD1 peaks that are NOT in promoters.
#          Stratify these enhancers by the expression status of their target gene.
#
# Three sets of regions:
#   1. Enhancers of DEGs DOWN (TES vs GFP)
#      → Peaks annotated to genes that are downregulated.
#   2. Enhancers of Random Control Genes
#      → Peaks annotated to expressed genes that are NOT differentially expressed.
#      → Subsampled to match the number of peaks in Group 1 (optional, but good for vis).
#
# =============================================================================

suppressPackageStartupMessages({
    library(GenomicRanges)
    library(rtracklayer)
    library(ChIPseeker)
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    library(org.Hs.eg.db)
    library(dplyr)
})

setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_1")

cat("==========================================================\n")
cat("ENHANCER ANALYSIS: DEGs DOWN vs CONTROL\n")
cat("==========================================================\n")
cat("Analysis started:", as.character(Sys.time()), "\n\n")

# =============================================================================
# PATH CONFIGURATION
# =============================================================================

# Peak files
TES_PEAKS <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_CUTandTAG/results/05_peaks_narrow/TES_peaks.narrowPeak"
TEAD1_PEAKS <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_CUTandTAG/results/05_peaks_narrow/TEAD1_peaks.narrowPeak"

# DESeq2 results
DESEQ2_FILE <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_RNA/results/05_deseq2/deseq2_results_TES_vs_GFP.txt"

# Output directory
OUTPUT_BASE <- "output/31_enhancer_analysis"
dir.create(OUTPUT_BASE, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# PHASE 1: LOAD AND ANNOTATE PEAKS
# =============================================================================

cat("=== PHASE 1: Loading and Annotating Peaks ===\n")

# Load peaks function
load_peaks <- function(peak_file) {
    peaks <- read.table(peak_file,
        header = FALSE, stringsAsFactors = FALSE,
        col.names = c(
            "chr", "start", "end", "name", "score",
            "strand", "signalValue", "pValue", "qValue", "peak"
        )
    )
    # Fix chr names
    if (!grepl("^chr", peaks$chr[1])) peaks$chr <- paste0("chr", peaks$chr)

    GRanges(seqnames = peaks$chr, ranges = IRanges(peaks$start + 1, peaks$end))
}

tes_gr <- load_peaks(TES_PEAKS)
tead1_gr <- load_peaks(TEAD1_PEAKS)

# Combine all binding peaks
all_peaks <- unique(c(tes_gr, tead1_gr))
cat(sprintf("  Total unique TES/TEAD1 peaks: %d\n", length(all_peaks)))

# Annotate peaks
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
cat("  Annotating peaks to nearest genes...\n")
peak_anno <- annotatePeak(all_peaks, TxDb = txdb, annoDb = "org.Hs.eg.db", verbose = FALSE)
peak_anno_df <- as.data.frame(peak_anno)

# Filter for Enhancers (Exclude Promoters)
# ChIPseeker "Promoter" is usually TSS -3kb to +3kb by default in annotatePeak if not specified,
# or explicit "Promoter (<=1kb)", "Promoter (1-2kb)", etc. in the annotation column.
# We will exclude anything starting with "Promoter".

# Check annotation categories
cat("  Annotation categories found:\n")
print(table(gsub(" \\(.*", "", peak_anno_df$annotation)))

enhancer_idx <- !grepl("^Promoter", peak_anno_df$annotation)
enhancer_peaks <- peak_anno_df[enhancer_idx, ]

cat(sprintf("\n  Defined 'Enhancers' (Non-Promoter Peaks): %d\n", nrow(enhancer_peaks)))
cat(sprintf("  Removed 'Promoter' Peaks: %d\n", nrow(peak_anno_df) - nrow(enhancer_peaks)))

# =============================================================================
# PHASE 2: INTEGRATE RNA-SEQ (DESEQ2)
# =============================================================================

cat("\n=== PHASE 2: Integrating RNA-seq Data ===\n")

deseq2 <- read.delim(DESEQ2_FILE, stringsAsFactors = FALSE)
cat(sprintf("  Total genes in DESeq2: %d\n", nrow(deseq2)))

# Define gene sets based on SYMBOLS
# Note: DESeq2 has 'gene_symbol', ChIPseeker has 'SYMBOL'.

# 1. DEGs DOWN
degs_down_genes <- deseq2 %>%
    filter(padj < 0.05, log2FoldChange < 0) %>%
    pull(gene_symbol)

# 2. Control Genes (Expressed but NOT DE)
# Expressed: padj is not NA (meaning it passed independent filtering usually implies expression) or baseMean > 10
control_genes <- deseq2 %>%
    filter(!is.na(padj), padj >= 0.05) %>%
    pull(gene_symbol)

cat(sprintf("  DEGs DOWN genes: %d\n", length(degs_down_genes)))
cat(sprintf("  Control genes (Expressed, Non-DE): %d\n", length(control_genes)))

# =============================================================================
# PHASE 3: ASSIGN ENHANCERS TO GROUPS
# =============================================================================

cat("\n=== PHASE 3: Assigning Enhancers to Groups ===\n")

# Enhancers of DEGs DOWN
enhancers_degs_down <- enhancer_peaks %>%
    filter(SYMBOL %in% degs_down_genes)

cat(sprintf("  Enhancers assigned to DEGs DOWN: %d\n", nrow(enhancers_degs_down)))

# Enhancers of Control Genes
enhancers_control_pool <- enhancer_peaks %>%
    filter(SYMBOL %in% control_genes)

cat(sprintf("  Enhancers assigned to Control Genes (Pool): %d\n", nrow(enhancers_control_pool)))

# Match sample size
n_target <- nrow(enhancers_degs_down)
set.seed(42)

if (nrow(enhancers_control_pool) > n_target) {
    enhancers_control <- enhancers_control_pool[sample(nrow(enhancers_control_pool), n_target), ]
    cat(sprintf("  Subsampled Control Enhancers to match N: %d\n", nrow(enhancers_control)))
} else {
    enhancers_control <- enhancers_control_pool
    cat("  WARNING: Control pool smaller than target. Using all available.\n")
}

# =============================================================================
# PHASE 4: EXPORT BED FILES
# =============================================================================

cat("\n=== PHASE 4: Exporting BED Files ===\n")

write_bed <- function(df, filename) {
    # BED format: chr, start, end
    # Note: ChIPseeker 'start'/'end' are 1-based (from GRanges). BED is 0-based.
    # makeGRangesFromDataFrame to standard bed export

    bed_df <- data.frame(
        chr = df$seqnames,
        start = df$start - 1,
        end = df$end,
        name = df$SYMBOL, # Use gene symbol as name
        score = 0,
        strand = df$strand
    )

    # Sort
    bed_df <- bed_df[order(bed_df$chr, bed_df$start), ]

    write.table(bed_df, filename, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
    cat(sprintf("  Created: %s\n", basename(filename)))
}

write_bed(enhancers_degs_down, file.path(OUTPUT_BASE, "Enhancers_DEGs_DOWN.bed"))
write_bed(enhancers_control, file.path(OUTPUT_BASE, "Enhancers_Random_Control.bed"))

# Save summary stats
summary_df <- data.frame(
    Group = c("Enhancers_DEGs_DOWN", "Enhancers_Random_Control"),
    N_Peaks = c(nrow(enhancers_degs_down), nrow(enhancers_control)),
    Description = c("Non-promoter peaks linked to DEGs DOWN", "Non-promoter peaks linked to Non-DE genes")
)
write.csv(summary_df, file.path(OUTPUT_BASE, "enhancer_counts.csv"), row.names = FALSE)

cat("\nAnalysis Complete.\n")
