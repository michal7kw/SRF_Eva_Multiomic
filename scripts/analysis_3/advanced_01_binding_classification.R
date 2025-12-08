#!/usr/bin/env Rscript

################################################################################
# Advanced Binding Site Classification: TES vs TEAD1
# Phase 1.1 - Binding Site Characterization
#
# Purpose: Classify all binding sites into mutually exclusive categories
#          comparing TES and TEAD1 (excluding TESmut)
#
# Author: Advanced Multi-Omics Analysis Plan
# Date: 2025-01-24
################################################################################

message("=== Phase 1.1: Advanced Binding Site Classification ===")
message("Start time: ", Sys.time())

# Load required libraries
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(dplyr)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(VennDiagram)
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
})

# Set working directory and paths
setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top")

# Input files - using combined replicate peaks from MACS2
TES_PEAKS <- "SRF_Eva_CUTandTAG/results/05_peaks_narrow/TES_peaks.narrowPeak"
TEAD1_PEAKS <- "SRF_Eva_CUTandTAG/results/05_peaks_narrow/TEAD1_peaks.narrowPeak"

# BigWig files for signal quantification
TES_BIGWIGS <- c(
  "SRF_Eva_CUTandTAG/results/06_bigwig/TES-1_CPM.bw",
  "SRF_Eva_CUTandTAG/results/06_bigwig/TES-2_CPM.bw",
  "SRF_Eva_CUTandTAG/results/06_bigwig/TES-3_CPM.bw"
)

TEAD1_BIGWIGS <- c(
  "SRF_Eva_CUTandTAG/results/06_bigwig/TEAD1-1_CPM.bw",
  "SRF_Eva_CUTandTAG/results/06_bigwig/TEAD1-2_CPM.bw",
  "SRF_Eva_CUTandTAG/results/06_bigwig/TEAD1-3_CPM.bw"
)

# Output directory
OUTPUT_DIR <- "SRF_Eva_integrated_analysis/scripts/analysis_3/results/01_binding_classification"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Parameters
OVERLAP_THRESHOLD <- 500  # bp for considering peaks as overlapping
SIGNAL_RATIO_EQUIV <- 0.2 # log2 ratio threshold for "equivalent" signal
SIGNAL_PERCENTILE_HIGH <- 0.75 # Top 25% for "high" signal

################################################################################
# Step 1: Load and prepare peak data
################################################################################

message("\n[Step 1] Loading peak files...")

# Load TES peaks
if (!file.exists(TES_PEAKS)) {
  stop("TES peaks file not found: ", TES_PEAKS)
}
tes_peaks <- import(TES_PEAKS, format = "narrowPeak")
# NOTE: Do NOT convert to UCSC style here - BigWig files use Ensembl style (1, 2, 3...)
# We will convert to UCSC style AFTER signal extraction
message("  Loaded ", length(tes_peaks), " TES peaks")
message("  Chromosome style: ", paste(head(seqlevels(tes_peaks), 5), collapse = ", "), "...")

# Load TEAD1 peaks
if (!file.exists(TEAD1_PEAKS)) {
  stop("TEAD1 peaks file not found: ", TEAD1_PEAKS)
}
tead1_peaks <- import(TEAD1_PEAKS, format = "narrowPeak")
# NOTE: Do NOT convert to UCSC style here - keep Ensembl style for BigWig compatibility
message("  Loaded ", length(tead1_peaks), " TEAD1 peaks")
message("  Chromosome style: ", paste(head(seqlevels(tead1_peaks), 5), collapse = ", "), "...")

# Add unique IDs
tes_peaks$peak_id <- paste0("TES_", seq_along(tes_peaks))
tead1_peaks$peak_id <- paste0("TEAD1_", seq_along(tead1_peaks))

################################################################################
# Step 2: Quantify binding signal at all peaks
################################################################################

message("\n[Step 2] Quantifying binding signals from BigWig files...")

# Function to extract mean signal from BigWig
# Handles chromosome naming mismatches (chr1 vs 1)
extract_signal <- function(peaks, bigwig_files) {
  signals <- matrix(NA, nrow = length(peaks), ncol = length(bigwig_files))

  # Get peak chromosome style
  peak_chroms <- unique(as.character(seqnames(peaks)))
  peak_has_chr <- any(grepl("^chr", peak_chroms))
  message("    Peak chromosome style: ",
          ifelse(peak_has_chr, "UCSC (chr1)", "Ensembl (1)"))
  message("    Sample peak chromosomes: ",
          paste(head(peak_chroms, 5), collapse = ", "))

  for (i in seq_along(bigwig_files)) {
    message("    Processing: ", basename(bigwig_files[i]))

    # Check if file exists
    if (!file.exists(bigwig_files[i])) {
      message("      WARNING: File not found: ", bigwig_files[i])
      next
    }

    tryCatch({
      # Import BigWig as RleList - direct chromosome access
      bw <- import(bigwig_files[i], format = "BigWig", as = "RleList")
      bw_chroms <- names(bw)

      # Determine chromosome naming style in BigWig
      bw_has_chr <- any(grepl("^chr", bw_chroms))
      message("      BigWig chromosome style: ",
              ifelse(bw_has_chr, "UCSC (chr1)", "Ensembl (1)"))
      message("      Sample BigWig chromosomes: ",
              paste(head(bw_chroms, 5), collapse = ", "))

      # Build chromosome name mapping
      chr_map <- function(chr) {
        if (bw_has_chr && !grepl("^chr", chr)) {
          # BigWig uses chr1, peak uses 1 -> add chr prefix
          return(paste0("chr", chr))
        } else if (!bw_has_chr && grepl("^chr", chr)) {
          # BigWig uses 1, peak uses chr1 -> remove chr prefix
          return(gsub("^chr", "", chr))
        } else {
          # Same style, return as-is
          return(chr)
        }
      }

      # Test the mapping on first peak
      test_chr <- as.character(seqnames(peaks[1]))
      test_mapped <- chr_map(test_chr)
      message("      Chromosome mapping test: '", test_chr, "' -> '", test_mapped, "'")
      message("      Mapped chr in BigWig: ", test_mapped %in% bw_chroms)

      # Extract signal for each peak with progress
      n_success <- 0
      n_chr_miss <- 0
      n_error <- 0

      for (j in seq_along(peaks)) {
        chr <- as.character(seqnames(peaks[j]))
        target_chr <- chr_map(chr)

        # Check if chromosome exists in BigWig
        if (!target_chr %in% bw_chroms) {
          n_chr_miss <- n_chr_miss + 1
          signals[j, i] <- NA
          next
        }

        result <- tryCatch({
          st <- start(peaks[j])
          en <- end(peaks[j])

          # Get chromosome length from BigWig
          chr_len <- length(bw[[target_chr]])

          # Ensure indices are within bounds
          st <- max(1, st)
          en <- min(chr_len, en)

          if (st <= en && en > 0 && chr_len > 0) {
            region_signal <- as.numeric(bw[[target_chr]][st:en])
            val <- mean(region_signal, na.rm = TRUE)
            if (is.nan(val) || is.infinite(val)) NA else val
          } else {
            NA
          }
        }, error = function(e) {
          n_error <<- n_error + 1
          NA
        })

        if (!is.na(result)) n_success <- n_success + 1
        signals[j, i] <- result
      }

      # Report extraction results
      message("      Results: ", n_success, " successful, ",
              n_chr_miss, " chr mismatches, ", n_error, " errors")

    }, error = function(e) {
      message("      WARNING: Error reading BigWig: ", e$message)
    })
  }

  # Return mean signal across replicates
  result <- rowMeans(signals, na.rm = TRUE)
  # Convert NaN to NA
  result[is.nan(result)] <- NA

  # Final report
  n_valid <- sum(!is.na(result))
  message("  Final: ", n_valid, "/", length(result),
          " peaks have valid signal (", round(100*n_valid/length(result), 1), "%)")

  result
}

# Extract TES signal at TES peaks
message("  Extracting TES signal at TES peaks...")
tes_peaks$tes_signal <- extract_signal(tes_peaks, TES_BIGWIGS)

# Extract TEAD1 signal at TES peaks
message("  Extracting TEAD1 signal at TES peaks...")
tes_peaks$tead1_signal <- extract_signal(tes_peaks, TEAD1_BIGWIGS)

# Extract TES signal at TEAD1 peaks
message("  Extracting TES signal at TEAD1 peaks...")
tead1_peaks$tes_signal <- extract_signal(tead1_peaks, TES_BIGWIGS)

# Extract TEAD1 signal at TEAD1 peaks
message("  Extracting TEAD1 signal at TEAD1 peaks...")
tead1_peaks$tead1_signal <- extract_signal(tead1_peaks, TEAD1_BIGWIGS)

################################################################################
# Step 3: Find overlapping peaks
################################################################################

message("\n[Step 3] Identifying overlapping peaks...")

# Find overlaps with maximum gap
overlaps <- findOverlaps(
  resize(tes_peaks, width = width(tes_peaks) + 2*OVERLAP_THRESHOLD, fix = "center"),
  resize(tead1_peaks, width = width(tead1_peaks) + 2*OVERLAP_THRESHOLD, fix = "center")
)

message("  Found ", length(overlaps), " TES-TEAD1 peak overlaps")

# Get indices of overlapping peaks
tes_overlap_idx <- queryHits(overlaps)
tead1_overlap_idx <- subjectHits(overlaps)

################################################################################
# Step 4: Classify peaks into categories
################################################################################

message("\n[Step 4] Classifying peaks into binding categories...")

# Initialize classification
tes_peaks$category <- "TES_unique"
tead1_peaks$category <- "TEAD1_unique"

# For overlapping peaks, determine subcategory based on signal
if (length(overlaps) > 0) {
  # Get signal ratios for overlapping peaks
  tes_overlap_signal <- tes_peaks$tes_signal[tes_overlap_idx]
  tead1_overlap_signal <- tead1_peaks$tead1_signal[tead1_overlap_idx]

  # Calculate log2 ratio (TES/TEAD1)
  # Add pseudocount to avoid division by zero
  log2_ratio <- log2((tes_overlap_signal + 1) / (tead1_overlap_signal + 1))

  # Classify shared peaks
  for (i in seq_along(overlaps)) {
    tes_idx <- tes_overlap_idx[i]
    tead1_idx <- tead1_overlap_idx[i]
    ratio <- log2_ratio[i]

    # Handle NA values in ratio - default to equivalent if signal data missing
    if (is.na(ratio)) {
      category <- "Shared_equivalent"
    } else if (abs(ratio) <= SIGNAL_RATIO_EQUIV) {
      category <- "Shared_equivalent"
    } else if (ratio > SIGNAL_RATIO_EQUIV) {
      category <- "Shared_TES_dominant"
    } else {
      category <- "Shared_TEAD1_dominant"
    }

    # Check if both signals are high (top 25%) - only if signals are not NA
    tes_sig <- tes_peaks$tes_signal[tes_idx]
    tead1_sig <- tead1_peaks$tead1_signal[tead1_idx]

    if (!is.na(tes_sig) && !is.na(tead1_sig)) {
      tes_percentile <- ecdf(tes_peaks$tes_signal[!is.na(tes_peaks$tes_signal)])(tes_sig)
      tead1_percentile <- ecdf(tead1_peaks$tead1_signal[!is.na(tead1_peaks$tead1_signal)])(tead1_sig)

      if (tes_percentile >= SIGNAL_PERCENTILE_HIGH && tead1_percentile >= SIGNAL_PERCENTILE_HIGH) {
        category <- "Shared_high"
      }
    }

    tes_peaks$category[tes_idx] <- category
    tead1_peaks$category[tead1_idx] <- category
  }
}

# Create unified peak set
all_peaks <- c(tes_peaks, tead1_peaks)

# Remove duplicate shared peaks (keep TES version)
shared_tes_idx <- which(tes_peaks$category != "TES_unique")
shared_tead1_idx <- which(tead1_peaks$category != "TEAD1_unique")

if (length(shared_tead1_idx) > 0) {
  # Only keep TES peaks for shared regions
  all_peaks_unique <- c(
    tes_peaks,
    tead1_peaks[-shared_tead1_idx]
  )
} else {
  all_peaks_unique <- all_peaks
}

message("  Classification complete:")
message("    TES_unique: ", sum(all_peaks_unique$category == "TES_unique"))
message("    TEAD1_unique: ", sum(all_peaks_unique$category == "TEAD1_unique"))
message("    Shared_high: ", sum(all_peaks_unique$category == "Shared_high"))
message("    Shared_TES_dominant: ", sum(all_peaks_unique$category == "Shared_TES_dominant"))
message("    Shared_TEAD1_dominant: ", sum(all_peaks_unique$category == "Shared_TEAD1_dominant"))
message("    Shared_equivalent: ", sum(all_peaks_unique$category == "Shared_equivalent"))

################################################################################
# Step 5: Genomic annotation
################################################################################

message("\n[Step 5] Annotating genomic context...")

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Ensure chromosome naming style matches (UCSC style: chr1, chr2, ...)
suppressWarnings({
  seqlevelsStyle(all_peaks_unique) <- "UCSC"
})

# Annotate all peaks
peak_anno <- annotatePeak(
  all_peaks_unique,
  tssRegion = c(-3000, 3000),
  TxDb = txdb,
  annoDb = "org.Hs.eg.db"
)

# Extract annotation data
anno_df <- as.data.frame(peak_anno)

# Add to peaks object
all_peaks_unique$annotation <- anno_df$annotation
all_peaks_unique$gene_name <- anno_df$SYMBOL
all_peaks_unique$distance_to_tss <- anno_df$distanceToTSS

################################################################################
# Step 6: Calculate additional genomic features
################################################################################

message("\n[Step 6] Calculating genomic features...")

# Peak width
all_peaks_unique$peak_width <- width(all_peaks_unique)

# CpG content (simplified - actual CpG would need sequence data)
# For now, use CpG island annotation proximity as proxy
# This would need BSgenome package for actual CpG counting

# Conservation scores (would need phastCons/phyloP data)
# Placeholder for now
all_peaks_unique$conservation <- NA

################################################################################
# Step 7: Export results
################################################################################

message("\n[Step 7] Exporting results...")

# Convert to data frame for export
peaks_df <- as.data.frame(all_peaks_unique)

# Add simplified 3-category classification
peaks_df <- peaks_df %>%
  mutate(
    category_simple = case_when(
      category == "TES_unique" ~ "TES_Unique",
      category == "TEAD1_unique" ~ "TEAD1_Unique",
      TRUE ~ "Shared"  # All Shared_* categories become "Shared"
    )
  ) %>%
  dplyr::select(seqnames, start, end, width, peak_id,
         category, category_simple,
         tes_signal, tead1_signal,
         annotation, gene_name, distance_to_tss,
         peak_width)

# Write main classification file (includes both detailed and simple categories)
output_file <- file.path(OUTPUT_DIR, "binding_site_classification.csv")
write.csv(peaks_df, output_file, row.names = FALSE)
message("  Saved: ", output_file)

# Write separate BED files for each detailed category (6 categories)
categories <- unique(peaks_df$category)
for (cat in categories) {
  cat_peaks <- all_peaks_unique[all_peaks_unique$category == cat]
  bed_file <- file.path(OUTPUT_DIR, paste0(cat, ".bed"))
  export(cat_peaks, bed_file, format = "BED")
  message("  Saved: ", bed_file)
}

# Write separate BED files for simplified categories (3 categories)
message("  Creating simplified category BED files...")
for (cat_simple in c("TES_Unique", "Shared", "TEAD1_Unique")) {
  cat_peaks <- all_peaks_unique[peaks_df$category_simple == cat_simple]
  bed_file <- file.path(OUTPUT_DIR, paste0(cat_simple, "_simple.bed"))
  export(cat_peaks, bed_file, format = "BED")
  message("  Saved: ", bed_file)
}

################################################################################
# Step 8: Generate summary statistics
################################################################################

message("\n[Step 8] Generating summary statistics...")

# Category counts - detailed (6 categories)
category_counts <- table(peaks_df$category)

# Category counts - simplified (3 categories)
category_counts_simple <- table(peaks_df$category_simple)

# Genomic distribution by category
genomic_dist <- peaks_df %>%
  group_by(category, annotation) %>%
  summarise(count = n(), .groups = "drop")

# Signal statistics by detailed category (6 categories)
signal_stats <- peaks_df %>%
  group_by(category) %>%
  summarise(
    n_peaks = n(),
    mean_tes_signal = mean(tes_signal, na.rm = TRUE),
    median_tes_signal = median(tes_signal, na.rm = TRUE),
    mean_tead1_signal = mean(tead1_signal, na.rm = TRUE),
    median_tead1_signal = median(tead1_signal, na.rm = TRUE),
    mean_peak_width = mean(peak_width, na.rm = TRUE),
    .groups = "drop"
  )

# Signal statistics by simplified category (3 categories)
signal_stats_simple <- peaks_df %>%
  group_by(category_simple) %>%
  summarise(
    n_peaks = n(),
    mean_tes_signal = mean(tes_signal, na.rm = TRUE),
    median_tes_signal = median(tes_signal, na.rm = TRUE),
    mean_tead1_signal = mean(tead1_signal, na.rm = TRUE),
    median_tead1_signal = median(tead1_signal, na.rm = TRUE),
    mean_peak_width = mean(peak_width, na.rm = TRUE),
    .groups = "drop"
  )

# Write summary
summary_file <- file.path(OUTPUT_DIR, "binding_category_summary.txt")
sink(summary_file)
cat("=== Binding Site Classification Summary ===\n")
cat("Date:", as.character(Sys.time()), "\n\n")
cat("Input Files:\n")
cat("  TES peaks:", TES_PEAKS, "\n")
cat("  TEAD1 peaks:", TEAD1_PEAKS, "\n\n")
cat("Parameters:\n")
cat("  Overlap threshold:", OVERLAP_THRESHOLD, "bp\n")
cat("  Signal ratio threshold:", SIGNAL_RATIO_EQUIV, "(log2)\n")
cat("  High signal percentile:", SIGNAL_PERCENTILE_HIGH, "\n\n")

cat("=====================================\n")
cat("DETAILED CLASSIFICATION (6 categories)\n")
cat("=====================================\n")
cat("Category Counts:\n")
print(category_counts)
cat("\nSignal Statistics by Category:\n")
print(signal_stats)

cat("\n\n=====================================\n")
cat("SIMPLIFIED CLASSIFICATION (3 categories)\n")
cat("=====================================\n")
cat("Category Counts:\n")
print(category_counts_simple)
cat("\nSignal Statistics by Category:\n")
print(signal_stats_simple)

cat("\n\n=====================================\n")
cat("SHARED PEAKS BREAKDOWN\n")
cat("=====================================\n")
shared_only <- peaks_df %>% filter(category_simple == "Shared")
cat("Total Shared peaks:", nrow(shared_only), "\n")
cat("  - Shared_high:", sum(peaks_df$category == "Shared_high"), "\n")
cat("  - Shared_TES_dominant:", sum(peaks_df$category == "Shared_TES_dominant"), "\n")
cat("  - Shared_TEAD1_dominant:", sum(peaks_df$category == "Shared_TEAD1_dominant"), "\n")
cat("  - Shared_equivalent:", sum(peaks_df$category == "Shared_equivalent"), "\n")

cat("\n\nGenomic Distribution (detailed):\n")
print(genomic_dist)
sink()
message("  Saved: ", summary_file)

################################################################################
# Step 9: Generate visualizations
################################################################################

message("\n[Step 9] Generating visualizations...")

# 9.1: Category counts bar plot
pdf(file.path(OUTPUT_DIR, "category_counts_barplot.pdf"), width = 10, height = 6)
ggplot(peaks_df, aes(x = category, fill = category)) +
  geom_bar() +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5) +
  scale_fill_brewer(palette = "Set2") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  labs(
    title = "Binding Site Category Distribution",
    x = "Category",
    y = "Number of Peaks"
  )
dev.off()

# 9.2: Signal comparison scatter plot
pdf(file.path(OUTPUT_DIR, "signal_comparison_scatterplot.pdf"), width = 10, height = 8)
ggplot(peaks_df, aes(x = log2(tes_signal + 1), y = log2(tead1_signal + 1), color = category)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_brewer(palette = "Set2") +
  theme_classic() +
  theme(legend.position = "right") +
  labs(
    title = "TES vs TEAD1 Signal Comparison",
    x = "log2(TES Signal + 1)",
    y = "log2(TEAD1 Signal + 1)",
    color = "Category"
  )
dev.off()

# 9.3: Venn diagram with signal thresholds
pdf(file.path(OUTPUT_DIR, "binding_category_venn.pdf"), width = 8, height = 8)
venn.plot <- draw.pairwise.venn(
  area1 = length(tes_peaks),
  area2 = length(tead1_peaks),
  cross.area = length(overlaps),
  category = c("TES", "TEAD1"),
  fill = c("skyblue", "pink"),
  lty = "blank",
  cex = 2,
  cat.cex = 2,
  cat.pos = c(330, 30),
  cat.dist = 0.05,
  scaled = TRUE
)
grid.draw(venn.plot)
dev.off()

# 9.4: Genomic context enrichment
genomic_summary <- peaks_df %>%
  mutate(
    genomic_region = case_when(
      grepl("Promoter", annotation) ~ "Promoter",
      grepl("Exon|5' UTR|3' UTR", annotation) ~ "Gene Body",
      grepl("Intron", annotation) ~ "Intron",
      grepl("Downstream", annotation) ~ "Downstream",
      TRUE ~ "Intergenic"
    )
  ) %>%
  group_by(category, genomic_region) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(category) %>%
  mutate(percentage = count / sum(count) * 100)

pdf(file.path(OUTPUT_DIR, "genomic_context_enrichment.pdf"), width = 12, height = 6)
ggplot(genomic_summary, aes(x = category, y = percentage, fill = genomic_region)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_brewer(palette = "Set3") +
  scale_y_continuous(labels = scales::percent) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Genomic Context Distribution by Binding Category",
    x = "Category",
    y = "Percentage",
    fill = "Genomic Region"
  )
dev.off()

# 9.5: Peak width distribution
pdf(file.path(OUTPUT_DIR, "peak_width_distribution.pdf"), width = 10, height = 6)
ggplot(peaks_df, aes(x = category, y = peak_width, fill = category)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.2, alpha = 0.5, outlier.shape = NA) +
  scale_fill_brewer(palette = "Set2") +
  scale_y_log10() +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  labs(
    title = "Peak Width Distribution by Category",
    x = "Category",
    y = "Peak Width (bp, log scale)"
  )
dev.off()

# 9.6: Create simplified 3-category classification for comparison
message("  Creating simplified 3-category plots...")

peaks_df <- peaks_df %>%
  mutate(
    category_simple = case_when(
      category == "TES_unique" ~ "TES_Unique",
      category == "TEAD1_unique" ~ "TEAD1_Unique",
      TRUE ~ "Shared"  # All Shared_* categories become "Shared"
    )
  )

# 9.6a: Simplified 3-category bar plot (like analysis_2)
simple_counts <- peaks_df %>%
  group_by(category_simple) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(category_simple = factor(category_simple,
         levels = c("TES_Unique", "Shared", "TEAD1_Unique")))

pdf(file.path(OUTPUT_DIR, "category_counts_barplot_simple.pdf"), width = 8, height = 6)
ggplot(simple_counts, aes(x = category_simple, y = count, fill = category_simple)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = count), vjust = -0.5, size = 5) +
  scale_fill_manual(values = c("TES_Unique" = "#E41A1C",
                                "Shared" = "#984EA3",
                                "TEAD1_Unique" = "#377EB8")) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none") +
  labs(
    title = "Peak Classification (Simplified 3-Category)",
    subtitle = sprintf("TES peaks: %d | TEAD1 peaks: %d",
                       length(tes_peaks), length(tead1_peaks)),
    x = "", y = "Number of Peaks"
  ) +
  ylim(0, max(simple_counts$count) * 1.15)
dev.off()

# 9.6b: Side-by-side comparison of both classification schemes
detailed_counts <- peaks_df %>%
  group_by(category) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(scheme = "Detailed (6 categories)")

simple_counts_long <- simple_counts %>%
  rename(category = category_simple) %>%
  mutate(scheme = "Simplified (3 categories)")

comparison_df <- bind_rows(detailed_counts, simple_counts_long)

pdf(file.path(OUTPUT_DIR, "category_comparison_both_schemes.pdf"), width = 14, height = 6)
# Define colors for all 9 unique categories (6 detailed + 3 simplified)
# Using consistent color scheme: red tones for TES, blue tones for TEAD1, purple for shared
category_colors <- c(
  # Detailed categories (6)
  "TES_unique" = "#E41A1C",
  "TEAD1_unique" = "#377EB8",
  "Shared_high" = "#984EA3",
  "Shared_TES_dominant" = "#FF7F00",
  "Shared_TEAD1_dominant" = "#4DAF4A",
  "Shared_equivalent" = "#A65628",
  # Simplified categories (3)
  "TES_Unique" = "#E41A1C",
  "TEAD1_Unique" = "#377EB8",
  "Shared" = "#984EA3"
)
p_comparison <- ggplot(comparison_df, aes(x = category, y = count, fill = category)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = count), vjust = -0.5, size = 4) +
  facet_wrap(~scheme, scales = "free_x") +
  scale_fill_manual(values = category_colors) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    strip.text = element_text(size = 14, face = "bold")
  ) +
  labs(
    title = "Binding Site Classification: Detailed vs Simplified",
    x = "", y = "Number of Peaks"
  )
print(p_comparison)
dev.off()

# 9.6c: Stacked bar showing how Shared breaks down into subcategories
shared_breakdown <- peaks_df %>%
  filter(grepl("Shared", category)) %>%
  group_by(category) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(
    percentage = round(100 * count / sum(count), 1),
    label = paste0(category, "\n(", count, ", ", percentage, "%)")
  )

pdf(file.path(OUTPUT_DIR, "shared_category_breakdown.pdf"), width = 10, height = 8)
# Pie chart for Shared breakdown
ggplot(shared_breakdown, aes(x = "", y = count, fill = category)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  scale_fill_brewer(palette = "Pastel1") +
  theme_void(base_size = 14) +
  theme(legend.position = "right") +
  geom_text(aes(label = paste0(count, "\n(", percentage, "%)")),
            position = position_stack(vjust = 0.5), size = 4) +
  labs(
    title = "Breakdown of Shared Binding Sites by Signal Dominance",
    subtitle = sprintf("Total Shared peaks: %d", sum(shared_breakdown$count)),
    fill = "Category"
  )
dev.off()

# 9.6d: Signal scatter plot with simplified coloring
pdf(file.path(OUTPUT_DIR, "signal_comparison_simple.pdf"), width = 10, height = 8)
ggplot(peaks_df, aes(x = log2(tes_signal + 1), y = log2(tead1_signal + 1),
                     color = category_simple)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = c("TES_Unique" = "#E41A1C",
                                 "Shared" = "#984EA3",
                                 "TEAD1_Unique" = "#377EB8")) +
  theme_classic(base_size = 12) +
  theme(legend.position = "right") +
  labs(
    title = "TES vs TEAD1 Signal (Simplified Categories)",
    x = "log2(TES Signal + 1)",
    y = "log2(TEAD1 Signal + 1)",
    color = "Category"
  )
dev.off()

# Add simplified category to the main data frame for export
all_peaks_unique$category_simple <- peaks_df$category_simple

message("  Saved additional comparison plots")

################################################################################
# Step 10: Save R objects for downstream analysis
################################################################################

message("\n[Step 10] Saving R objects...")

save(
  all_peaks_unique,
  tes_peaks,
  tead1_peaks,
  peaks_df,
  signal_stats,
  category_counts,
  file = file.path(OUTPUT_DIR, "binding_classification_data.RData")
)

message("\n=== Analysis Complete ===")
message("Output directory: ", OUTPUT_DIR)
message("End time: ", Sys.time())
message("\nNext steps:")
message("  1. Review binding_site_classification.csv")
message("  2. Inspect visualization PDFs")
message("  3. Proceed to Phase 1.2 (Motif Analysis)")
