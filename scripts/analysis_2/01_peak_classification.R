#!/usr/bin/env Rscript

# 01_peak_classification.R
# Module 1: Refined Peak Classification for TES and TEAD1
# Part of the SRF_Eva_integrated_analysis pipeline
#
# Uses narrowPeak files from MACS2 (05_peaks_narrow) instead of
# problematic consensus peaks from failed IDR analysis

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(dplyr)
  library(VennDiagram)
  library(grid)
})

# ===================== Configuration =====================

# Base directories - use absolute paths throughout
BASE_DIR <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top"
SCRIPT_DIR <- file.path(BASE_DIR, "SRF_Eva_integrated_analysis/scripts/analysis_2")

# Input files - Use replicate narrowPeak files from MACS2
PEAK_DIR <- file.path(BASE_DIR, "SRF_Eva_CUTandTAG/results/05_peaks_narrow")

# Replicate peak files
TES_REPLICATES <- c(
  file.path(PEAK_DIR, "TES-1_peaks.narrowPeak"),
  file.path(PEAK_DIR, "TES-2_peaks.narrowPeak"),
  file.path(PEAK_DIR, "TES-3_peaks.narrowPeak")
)

TEAD1_REPLICATES <- c(
  file.path(PEAK_DIR, "TEAD1-1_peaks.narrowPeak"),
  file.path(PEAK_DIR, "TEAD1-2_peaks.narrowPeak"),
  file.path(PEAK_DIR, "TEAD1-3_peaks.narrowPeak")
)

# Minimum replicate overlap for consensus (2 out of 3 = reproducible)
MIN_REPLICATE_OVERLAP <- 2

# Output directory - use absolute path
OUTPUT_DIR <- file.path(SCRIPT_DIR, "results/01_peak_classification")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Annotation Database
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# ===================== Functions =====================

# Function to load narrowPeak files
load_narrowpeak <- function(peak_file, name = "peaks") {
  if (!file.exists(peak_file)) {
    stop(paste("File not found:", peak_file))
  }

  cat(sprintf("  Loading %s from %s\n", name, basename(peak_file)))

  # Read narrowPeak format (10 columns)
  # chr, start, end, name, score, strand, signalValue, pValue, qValue, peak
  peak_data <- read.table(peak_file, header = FALSE, sep = "\t",
                          stringsAsFactors = FALSE)

  # Handle different column counts (some files may have fewer columns)
  if (ncol(peak_data) >= 10) {
    colnames(peak_data)[1:10] <- c("chr", "start", "end", "name", "score",
                                    "strand", "signalValue", "pValue", "qValue", "peak")
  } else if (ncol(peak_data) >= 6) {
    colnames(peak_data)[1:6] <- c("chr", "start", "end", "name", "score", "strand")
  } else {
    colnames(peak_data)[1:3] <- c("chr", "start", "end")
  }

  # Create GRanges object
  peaks <- GRanges(
    seqnames = peak_data$chr,
    ranges = IRanges(start = peak_data$start + 1, end = peak_data$end),  # Convert 0-based to 1-based
    score = if ("score" %in% names(peak_data)) peak_data$score else NA,
    signalValue = if ("signalValue" %in% names(peak_data)) peak_data$signalValue else NA,
    pValue = if ("pValue" %in% names(peak_data)) peak_data$pValue else NA,
    qValue = if ("qValue" %in% names(peak_data)) peak_data$qValue else NA
  )

  # Add chr prefix if missing
  if (!any(grepl("^chr", seqnames(peaks)))) {
    seqlevels(peaks) <- paste0("chr", seqlevels(peaks))
  }

  # Standardize chromosome style
  seqlevelsStyle(peaks) <- "UCSC"

  # Filter to standard chromosomes only
  standard_chr <- paste0("chr", c(1:22, "X", "Y"))
  peaks <- keepSeqlevels(peaks, standard_chr[standard_chr %in% seqlevels(peaks)],
                         pruning.mode = "coarse")

  return(peaks)
}

# Function to create consensus peaks from replicates
# Requires peaks to be present in at least min_overlap replicates
create_consensus_peaks <- function(replicate_files, min_overlap = 2, name = "TF") {
  cat(sprintf("\n  Creating consensus peaks for %s (min %d/%d replicates)...\n",
              name, min_overlap, length(replicate_files)))

  # Load all replicate peaks
  rep_peaks <- list()
  for (i in seq_along(replicate_files)) {
    rep_name <- sprintf("%s_rep%d", name, i)
    rep_peaks[[rep_name]] <- load_narrowpeak(replicate_files[i], rep_name)
    cat(sprintf("    Replicate %d: %d peaks\n", i, length(rep_peaks[[rep_name]])))
  }

  # Merge all peaks into a union set
  # Use GRangesList and unlist to properly combine GRanges objects
  all_peaks <- unlist(GRangesList(rep_peaks))
  merged_peaks <- GenomicRanges::reduce(all_peaks)

  cat(sprintf("    Union of all peaks: %d regions\n", length(merged_peaks)))

  # Count how many replicates each merged region overlaps
  overlap_counts <- sapply(seq_along(merged_peaks), function(i) {
    region <- merged_peaks[i]
    sum(sapply(rep_peaks, function(rp) any(overlapsAny(region, rp))))
  })

  # Keep only regions present in >= min_overlap replicates
  consensus <- merged_peaks[overlap_counts >= min_overlap]
  consensus$replicate_count <- overlap_counts[overlap_counts >= min_overlap]

  cat(sprintf("    Consensus peaks (>=%d replicates): %d\n", min_overlap, length(consensus)))

  # Calculate reproducibility stats
  rep_stats <- table(overlap_counts)
  cat("    Replicate overlap distribution:\n")
  for (n in names(rep_stats)) {
    cat(sprintf("      %s replicates: %d peaks\n", n, rep_stats[n]))
  }

  return(consensus)
}

# Function to export GRanges to BED format
export_bed <- function(gr, filepath) {
  if (length(gr) == 0) {
    cat(sprintf("  Skipping %s (empty)\n", basename(filepath)))
    return(FALSE)
  }

  # Create data frame for BED export
  bed_df <- data.frame(
    chr = as.character(seqnames(gr)),
    start = start(gr) - 1,  # Convert back to 0-based for BED
    end = end(gr),
    name = if (!is.null(gr$name)) gr$name else paste0("peak_", seq_along(gr)),
    score = if (!is.null(gr$score)) gr$score else 0,
    strand = as.character(strand(gr))
  )
  bed_df$strand[bed_df$strand == "*"] <- "."

  write.table(bed_df, filepath, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
  cat(sprintf("  Exported %s (%d peaks)\n", basename(filepath), nrow(bed_df)))
  return(TRUE)
}

# ===================== Main Analysis =====================

cat("=== Module 1: Peak Classification ===\n\n")
cat(sprintf("Working directory: %s\n", getwd()))
cat(sprintf("Output directory: %s\n", OUTPUT_DIR))
cat(sprintf("Date: %s\n\n", Sys.time()))

# Create consensus peaks from replicates
cat("Creating consensus peaks from replicates...\n")

# Check that all replicate files exist
for (f in c(TES_REPLICATES, TEAD1_REPLICATES)) {
  if (!file.exists(f)) {
    stop(paste("File not found:", f))
  }
}

# Create consensus peaks (present in >= MIN_REPLICATE_OVERLAP replicates)
tes_peaks <- create_consensus_peaks(TES_REPLICATES, MIN_REPLICATE_OVERLAP, "TES")
tead1_peaks <- create_consensus_peaks(TEAD1_REPLICATES, MIN_REPLICATE_OVERLAP, "TEAD1")

cat(sprintf("\n  TES consensus peaks: %d\n", length(tes_peaks)))
cat(sprintf("  TEAD1 consensus peaks: %d\n", length(tead1_peaks)))

# Export individual consensus peak files for downstream use
export_bed(tes_peaks, file.path(OUTPUT_DIR, "TES_consensus_peaks.bed"))
export_bed(tead1_peaks, file.path(OUTPUT_DIR, "TEAD1_consensus_peaks.bed"))

# 1. Peak Classification
cat("\nClassifying peaks by overlap...\n")

# Find overlapping peaks
# Shared: TES peaks that overlap with TEAD1
tes_overlaps_tead1 <- overlapsAny(tes_peaks, tead1_peaks)
tead1_overlaps_tes <- overlapsAny(tead1_peaks, tes_peaks)

# Create classification
tes_shared <- tes_peaks[tes_overlaps_tead1]
tes_unique <- tes_peaks[!tes_overlaps_tead1]
tead1_unique <- tead1_peaks[!tead1_overlaps_tes]

# For shared regions, merge overlapping peaks from both TFs
shared_ranges <- GenomicRanges::reduce(c(tes_shared, tead1_peaks[tead1_overlaps_tes]))

cat(sprintf("  Shared regions (merged): %d\n", length(shared_ranges)))
cat(sprintf("  TES-Unique peaks: %d\n", length(tes_unique)))
cat(sprintf("  TEAD1-Unique peaks: %d\n", length(tead1_unique)))

# Sanity check
if (length(tes_unique) == 0 && length(tead1_unique) == 0 && length(shared_ranges) == 0) {
  stop("ERROR: All peak categories are empty. Check input files.")
}

# 2. Export BED files
cat("\nExporting BED files...\n")

# Add category metadata and export
shared_ranges$category <- "Shared"
tes_unique$category <- "TES_Unique"
tead1_unique$category <- "TEAD1_Unique"

export_bed(shared_ranges, file.path(OUTPUT_DIR, "Shared_TES_TEAD1.bed"))
export_bed(tes_unique, file.path(OUTPUT_DIR, "TES_Unique.bed"))
export_bed(tead1_unique, file.path(OUTPUT_DIR, "TEAD1_Unique.bed"))

# Create Master Peak List
all_peaks <- c(shared_ranges, tes_unique, tead1_unique)
export_bed(all_peaks, file.path(OUTPUT_DIR, "Master_Peak_List.bed"))

# 3. Annotation with ChIPseeker
cat("\nAnnotating peaks with ChIPseeker...\n")

annotate_peak_set <- function(peaks, name) {
  if (length(peaks) == 0) {
    cat(sprintf("  Skipping annotation for %s (empty)\n", name))
    return(NULL)
  }
  cat(sprintf("  Annotating %s (%d peaks)...\n", name, length(peaks)))
  anno <- annotatePeak(peaks,
                       TxDb = txdb,
                       tssRegion = c(-2000, 2000),
                       verbose = FALSE)
  return(anno)
}

anno_shared <- annotate_peak_set(shared_ranges, "Shared")
anno_tes_unique <- annotate_peak_set(tes_unique, "TES_Unique")
anno_tead1_unique <- annotate_peak_set(tead1_unique, "TEAD1_Unique")

# Convert annotations to data frames
df_list <- list()

if (!is.null(anno_shared)) {
  df_shared <- as.data.frame(anno_shared)
  df_shared$category <- "Shared"
  df_list[["shared"]] <- df_shared
}

if (!is.null(anno_tes_unique)) {
  df_tes_unique <- as.data.frame(anno_tes_unique)
  df_tes_unique$category <- "TES_Unique"
  df_list[["tes"]] <- df_tes_unique
}

if (!is.null(anno_tead1_unique)) {
  df_tead1_unique <- as.data.frame(anno_tead1_unique)
  df_tead1_unique$category <- "TEAD1_Unique"
  df_list[["tead1"]] <- df_tead1_unique
}

# Combine annotations
if (length(df_list) == 0) {
  stop("ERROR: No peaks were annotated. Cannot continue.")
}

common_cols <- Reduce(intersect, lapply(df_list, names))
combined_anno <- do.call(rbind, lapply(df_list, function(df) df[, common_cols]))

# Add Gene Symbols
cat("  Adding gene symbols...\n")
if (!"SYMBOL" %in% colnames(combined_anno) && "geneId" %in% colnames(combined_anno)) {
  symbols <- mapIds(org.Hs.eg.db,
                    keys = as.character(combined_anno$geneId),
                    column = "SYMBOL",
                    keytype = "ENTREZID",
                    multiVals = "first")
  combined_anno$SYMBOL <- symbols
}

# Save annotation
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

pdf(file.path(OUTPUT_DIR, "category_counts_barplot.pdf"), width = 8, height = 6)
p1 <- ggplot(counts_df, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Count), vjust = -0.5, size = 5) +
  scale_fill_manual(values = c("TES_Unique" = "#E41A1C",
                                "Shared" = "#984EA3",
                                "TEAD1_Unique" = "#377EB8")) +
  theme_minimal(base_size = 14) +
  labs(title = "Peak Classification Counts",
       subtitle = sprintf("TES consensus: %d | TEAD1 consensus: %d (>=%d/3 replicates)",
                          length(tes_peaks), length(tead1_peaks), MIN_REPLICATE_OVERLAP),
       y = "Number of Peaks", x = "") +
  theme(legend.position = "none") +
  ylim(0, max(counts_df$Count) * 1.15)
print(p1)
dev.off()
cat("  Saved category_counts_barplot.pdf\n")

# Annotation Distribution Plot
anno_list <- list(Shared = anno_shared,
                  TES_Unique = anno_tes_unique,
                  TEAD1_Unique = anno_tead1_unique)
anno_list <- anno_list[!sapply(anno_list, is.null)]

if (length(anno_list) > 0) {
  pdf(file.path(OUTPUT_DIR, "genomic_context_enrichment.pdf"), width = 10, height = 6)
  p2 <- plotAnnoBar(anno_list)
  print(p2)
  dev.off()
  cat("  Saved genomic_context_enrichment.pdf\n")

  pdf(file.path(OUTPUT_DIR, "distance_to_TSS.pdf"), width = 10, height = 6)
  p3 <- plotDistToTSS(anno_list)
  print(p3)
  dev.off()
  cat("  Saved distance_to_TSS.pdf\n")
}

# Venn diagram showing TES and TEAD1 peak overlap
pdf(file.path(OUTPUT_DIR, "binding_category_venn.pdf"), width = 8, height = 8)

# Suppress VennDiagram log file creation
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

venn.plot <- draw.pairwise.venn(
  area1 = length(tes_peaks),
  area2 = length(tead1_peaks),
  cross.area = sum(tes_overlaps_tead1),
  category = c("TES", "TEAD1"),
  fill = c("#E41A1C", "#377EB8"),
  alpha = 0.5,
  col = c("#E41A1C", "#377EB8"),
  lwd = 2,
  fontfamily = "sans",
  cat.fontfamily = "sans",
  cat.fontface = "bold",
  cat.cex = 1.5,
  cex = 1.3,
  cat.pos = c(-20, 20),
  cat.dist = 0.05,
  margin = 0.1,
  scaled = TRUE
)

# Add title
grid.text("TES and TEAD1 Binding Site Overlap",
          y = 0.95, gp = gpar(fontsize = 16, fontface = "bold"))
grid.text(sprintf("Consensus peaks (≥%d/3 replicates)", MIN_REPLICATE_OVERLAP),
          y = 0.90, gp = gpar(fontsize = 11))

dev.off()
cat("  Saved binding_category_venn.pdf\n")

# 5. Summary Statistics
cat("\n=== Summary Statistics ===\n")
cat(sprintf("Input: %d TES replicates, %d TEAD1 replicates\n",
            length(TES_REPLICATES), length(TEAD1_REPLICATES)))
cat(sprintf("Consensus threshold: >= %d replicates\n", MIN_REPLICATE_OVERLAP))
cat(sprintf("\nConsensus peaks:\n"))
cat(sprintf("  TES: %d\n", length(tes_peaks)))
cat(sprintf("  TEAD1: %d\n", length(tead1_peaks)))
cat(sprintf("\nClassification:\n"))
cat(sprintf("  Shared (both TES+TEAD1): %d (%.1f%% of TES, %.1f%% of TEAD1)\n",
            length(shared_ranges),
            100 * sum(tes_overlaps_tead1) / length(tes_peaks),
            100 * sum(tead1_overlaps_tes) / length(tead1_peaks)))
cat(sprintf("  TES-Unique: %d (%.1f%% of TES)\n",
            length(tes_unique),
            100 * length(tes_unique) / length(tes_peaks)))
cat(sprintf("  TEAD1-Unique: %d (%.1f%% of TEAD1)\n",
            length(tead1_unique),
            100 * length(tead1_unique) / length(tead1_peaks)))

# Save summary
summary_text <- sprintf("
=== Peak Classification Summary ===
Date: %s

Replicate-Based Consensus Peak Analysis
=======================================

Input Files:
  TES replicates:
    - %s
    - %s
    - %s
  TEAD1 replicates:
    - %s
    - %s
    - %s

Consensus Parameters:
  Minimum replicate overlap: %d out of 3

Consensus Peak Counts:
  TES consensus peaks: %d
  TEAD1 consensus peaks: %d

Classification Results:
  Shared (TES+TEAD1): %d merged regions
    - TES peaks with TEAD1 overlap: %d (%.1f%%)
    - TEAD1 peaks with TES overlap: %d (%.1f%%)
  TES-Unique: %d (%.1f%% of TES consensus)
  TEAD1-Unique: %d (%.1f%% of TEAD1 consensus)

Output Files:
  Consensus peaks:
    - TES_consensus_peaks.bed
    - TEAD1_consensus_peaks.bed
  Classification:
    - Shared_TES_TEAD1.bed
    - TES_Unique.bed
    - TEAD1_Unique.bed
    - Master_Peak_List.bed
    - Master_Peak_Annotation.csv
  Visualizations:
    - category_counts_barplot.pdf
    - genomic_context_enrichment.pdf
    - distance_to_TSS.pdf
    - binding_category_venn.pdf
",
Sys.time(),
basename(TES_REPLICATES[1]), basename(TES_REPLICATES[2]), basename(TES_REPLICATES[3]),
basename(TEAD1_REPLICATES[1]), basename(TEAD1_REPLICATES[2]), basename(TEAD1_REPLICATES[3]),
MIN_REPLICATE_OVERLAP,
length(tes_peaks),
length(tead1_peaks),
length(shared_ranges),
sum(tes_overlaps_tead1), 100 * sum(tes_overlaps_tead1) / length(tes_peaks),
sum(tead1_overlaps_tes), 100 * sum(tead1_overlaps_tes) / length(tead1_peaks),
length(tes_unique), 100 * length(tes_unique) / length(tes_peaks),
length(tead1_unique), 100 * length(tead1_unique) / length(tead1_peaks)
)

writeLines(summary_text, file.path(OUTPUT_DIR, "binding_category_summary.txt"))

cat("\n=== Module 1 Complete ===\n")
cat(sprintf("Results saved to: %s\n", OUTPUT_DIR))
