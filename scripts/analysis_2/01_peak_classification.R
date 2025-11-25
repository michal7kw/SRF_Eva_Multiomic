#!/usr/bin/env Rscript

# 01_peak_classification.R
# Module 1: Refined Peak Classification for TES and TEAD1
# Part of the SRF_Eva_integrated_analysis pipeline
#
# FIXED: Path inconsistencies, removed unused ggVennDiagram

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(dplyr)
})

# ===================== Configuration =====================

# Base directories - use absolute paths throughout
BASE_DIR <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top"
SCRIPT_DIR <- file.path(BASE_DIR, "SRF_Eva_integrated_analysis/scripts/analysis_2")

# Input files (Consensus peaks from Cut&Tag pipeline)
TES_PEAKS_FILE <- file.path(BASE_DIR, "SRF_Eva_CUTandTAG/results/11_combined_replicates_narrow/peaks/TES_consensus_peaks.bed")
TEAD1_PEAKS_FILE <- file.path(BASE_DIR, "SRF_Eva_CUTandTAG/results/11_combined_replicates_narrow/peaks/TEAD1_consensus_peaks.bed")

# Output directory - use absolute path
OUTPUT_DIR <- file.path(SCRIPT_DIR, "results/01_peak_classification")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Annotation Database
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# ===================== Functions =====================

# Function to load peaks
load_peaks <- function(bed_file) {
  if (!file.exists(bed_file)) {
    stop(paste("File not found:", bed_file))
  }
  peaks <- import(bed_file, format = "BED")
  # Ensure standard chromosome style
  seqlevelsStyle(peaks) <- "UCSC"
  return(peaks)
}

# ===================== Main Analysis =====================

cat("=== Module 1: Peak Classification ===\n\n")

cat("Loading peaks...\n")
tes_peaks <- load_peaks(TES_PEAKS_FILE)
tead1_peaks <- load_peaks(TEAD1_PEAKS_FILE)

cat(sprintf("  TES peaks: %d\n", length(tes_peaks)))
cat(sprintf("  TEAD1 peaks: %d\n", length(tead1_peaks)))

# 1. Peak Intersection
cat("\nClassifying peaks...\n")

# Find overlaps
# Shared: TES peaks that overlap with TEAD1 (and vice versa)
# Strategy:
# - Shared: Intersection of TES and TEAD1
# - TES_Unique: TES peaks - TEAD1 peaks
# - TEAD1_Unique: TEAD1 peaks - TES peaks

# Use GenomicRanges set operations
# intersect() returns the genomic ranges that are common
shared_ranges <- GenomicRanges::intersect(tes_peaks, tead1_peaks)
# reduce() merges adjacent ranges to form a clean consensus
shared_ranges <- GenomicRanges::reduce(shared_ranges)

# Identify unique peaks
# subsetByOverlaps(A, B, invert=TRUE) returns ranges in A that do NOT overlap B
tes_unique <- subsetByOverlaps(tes_peaks, tead1_peaks, invert = TRUE)
tead1_unique <- subsetByOverlaps(tead1_peaks, tes_peaks, invert = TRUE)

cat(sprintf("  Shared regions: %d\n", length(shared_ranges)))
cat(sprintf("  TES-Unique peaks: %d\n", length(tes_unique)))
cat(sprintf("  TEAD1-Unique peaks: %d\n", length(tead1_unique)))

# 2. Export BED files
cat("\nExporting BED files...\n")
export(shared_ranges, file.path(OUTPUT_DIR, "Shared_TES_TEAD1.bed"))
export(tes_unique, file.path(OUTPUT_DIR, "TES_Unique.bed"))
export(tead1_unique, file.path(OUTPUT_DIR, "TEAD1_Unique.bed"))

# Create a combined Master Peak List with metadata
shared_ranges$category <- "Shared"
tes_unique$category <- "TES_Unique"
tead1_unique$category <- "TEAD1_Unique"

all_peaks <- c(shared_ranges, tes_unique, tead1_unique)
export(all_peaks, file.path(OUTPUT_DIR, "Master_Peak_List.bed"))

# 3. Annotation
cat("\nAnnotating peaks...\n")

annotate_peak_set <- function(peaks, name) {
  if (length(peaks) == 0) {
    warning(paste("No peaks to annotate for", name))
    return(NULL)
  }
  anno <- annotatePeak(peaks,
                       TxDb = txdb,
                       tssRegion = c(-2000, 2000),
                       verbose = FALSE)
  return(anno)
}

anno_shared <- annotate_peak_set(shared_ranges, "Shared")
anno_tes_unique <- annotate_peak_set(tes_unique, "TES_Unique")
anno_tead1_unique <- annotate_peak_set(tead1_unique, "TEAD1_Unique")

# Convert to data frames
df_shared <- as.data.frame(anno_shared)
df_shared$category <- "Shared"

df_tes_unique <- as.data.frame(anno_tes_unique)
df_tes_unique$category <- "TES_Unique"

df_tead1_unique <- as.data.frame(anno_tead1_unique)
df_tead1_unique$category <- "TEAD1_Unique"

# Combine annotation results
# Ensure columns match
common_cols <- intersect(names(df_shared), intersect(names(df_tes_unique), names(df_tead1_unique)))
combined_anno <- rbind(
  df_shared[, common_cols],
  df_tes_unique[, common_cols],
  df_tead1_unique[, common_cols]
)

# Add Gene Symbol
cat("  Adding gene symbols...\n")
if (!"SYMBOL" %in% colnames(combined_anno)) {
    symbols <- mapIds(org.Hs.eg.db,
                      keys = as.character(combined_anno$geneId),
                      column = "SYMBOL",
                      keytype = "ENTREZID",
                      multiVals = "first")
    combined_anno$SYMBOL <- symbols
}

# Save detailed annotation
write.csv(combined_anno, file.path(OUTPUT_DIR, "Master_Peak_Annotation.csv"), row.names = FALSE)

cat(sprintf("  Total annotated peaks: %d\n", nrow(combined_anno)))

# 4. Visualization
cat("\nGenerating visualizations...\n")

# Peak Category Counts (Bar Chart)
counts_df <- data.frame(
  Category = factor(c("TES_Unique", "Shared", "TEAD1_Unique"),
                    levels = c("TES_Unique", "Shared", "TEAD1_Unique")),
  Count = c(length(tes_unique), length(shared_ranges), length(tead1_unique))
)

pdf(file.path(OUTPUT_DIR, "Peak_Overlap_Venn.pdf"), width = 8, height = 6)
p1 <- ggplot(counts_df, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Count), vjust = -0.5, size = 5) +
  scale_fill_manual(values = c("TES_Unique" = "#E41A1C", "Shared" = "#984EA3", "TEAD1_Unique" = "#377EB8")) +
  theme_minimal(base_size = 14) +
  labs(title = "Peak Classification Counts",
       subtitle = "TES and TEAD1 binding site overlap",
       y = "Number of Peaks", x = "") +
  theme(legend.position = "none")
print(p1)
dev.off()

# Annotation Distribution Plot
pdf(file.path(OUTPUT_DIR, "Peak_Annotation_Distribution.pdf"), width = 10, height = 6)
anno_list <- list(Shared = anno_shared, TES_Unique = anno_tes_unique, TEAD1_Unique = anno_tead1_unique)
anno_list <- anno_list[!sapply(anno_list, is.null)]
if (length(anno_list) > 0) {
  p2 <- plotAnnoBar(anno_list)
  print(p2)
}
dev.off()

# Distance to TSS
pdf(file.path(OUTPUT_DIR, "Distance_to_TSS.pdf"), width = 10, height = 6)
if (length(anno_list) > 0) {
  p3 <- plotDistToTSS(anno_list)
  print(p3)
}
dev.off()

# Summary Statistics
cat("\n=== Summary Statistics ===\n")
cat(sprintf("Total TES peaks: %d\n", length(tes_peaks)))
cat(sprintf("Total TEAD1 peaks: %d\n", length(tead1_peaks)))
cat(sprintf("Shared (overlap): %d (%.1f%% of TES, %.1f%% of TEAD1)\n",
            length(shared_ranges),
            100 * length(shared_ranges) / length(tes_peaks),
            100 * length(shared_ranges) / length(tead1_peaks)))
cat(sprintf("TES-Unique: %d (%.1f%% of TES)\n",
            length(tes_unique),
            100 * length(tes_unique) / length(tes_peaks)))
cat(sprintf("TEAD1-Unique: %d (%.1f%% of TEAD1)\n",
            length(tead1_unique),
            100 * length(tead1_unique) / length(tead1_peaks)))

cat("\nAnalysis complete. Results saved to", OUTPUT_DIR, "\n")
