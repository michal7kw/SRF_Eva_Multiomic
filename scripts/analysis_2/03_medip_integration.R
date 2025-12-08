#!/usr/bin/env Rscript

# 03_medip_integration.R
# Module 3: Integration with DNA Methylation (meDIP)
# Part of the SRF_Eva_integrated_analysis pipeline
#
# FIXED: Signal extraction, paths, variable naming (TES vs GFP comparison)

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(ggplot2)
  library(dplyr)
})

# ===================== Configuration =====================

# Base directories - use absolute paths throughout
BASE_DIR <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top"
SCRIPT_DIR <- file.path(BASE_DIR, "SRF_Eva_integrated_analysis/scripts/analysis_2")

# Input: Peak Categories from Module 1 (using correct path)
PEAK_DIR <- file.path(SCRIPT_DIR, "results/01_peak_classification")
SHARED_PEAKS <- file.path(PEAK_DIR, "Shared_TES_TEAD1.bed")
TES_UNIQUE <- file.path(PEAK_DIR, "TES_Unique.bed")
TEAD1_UNIQUE <- file.path(PEAK_DIR, "TEAD1_Unique.bed")

# Input: meDIP Data
MEDIP_DIR <- file.path(BASE_DIR, "meDIP/results")
# TES condition meDIP (this is TES sample methylation, not an input control)
TES_MEDIP_BW <- file.path(MEDIP_DIR, "05_bigwig", "TES-1-IP_RPKM.bw")
# GFP control condition meDIP (the control comparison)
GFP_MEDIP_BW <- file.path(MEDIP_DIR, "05_bigwig", "GFP-1-IP_RPKM.bw")
# DMR file
DMR_FILE <- file.path(MEDIP_DIR, "07_differential_MEDIPS", "TES_vs_GFP_DMRs_FDR05.csv")

# Output
OUTPUT_DIR <- file.path(SCRIPT_DIR, "results/03_medip_integration")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ===================== Functions =====================

load_peaks <- function(bed_file) {
  if (!file.exists(bed_file)) {
    stop(paste("File not found:", bed_file))
  }
  import(bed_file, format = "BED")
}

# Function to compute average signal over regions - FIXED version
compute_signal <- function(regions, bw_file, name) {
  cat(sprintf("  Computing %s signal...\n", name))
  if (!file.exists(bw_file)) {
    warning(paste("BigWig not found:", bw_file))
    return(rep(NA, length(regions)))
  }

  tryCatch({
    # Get BigWig seqinfo
    bw_obj <- BigWigFile(bw_file)

    # Adjust chromosome style to match BigWig
    regions_adj <- regions

    # Check if styles match
    if (length(seqlevels(regions_adj)) > 0 && length(seqlevels(bw_obj)) > 0) {
      bw_style <- seqlevelsStyle(bw_obj)
      if (length(bw_style) > 0) {
        seqlevelsStyle(regions_adj) <- bw_style[1]
      }
    }

    # Filter regions to only those chromosomes present in BigWig
    valid_seqlevels <- intersect(seqlevels(regions_adj), seqlevels(bw_obj))
    if (length(valid_seqlevels) == 0) {
      warning("No overlapping chromosomes between regions and BigWig")
      return(rep(NA, length(regions)))
    }

    regions_adj <- keepSeqlevels(regions_adj, valid_seqlevels,
                                 pruning.mode = "coarse")

    # Import the signal for just these regions
    signal <- import(bw_file, format = "BigWig", which = regions_adj)

    # For each region, calculate mean signal
    scores <- sapply(seq_along(regions_adj), function(i) {
      region <- regions_adj[i]
      overlapping <- subsetByOverlaps(signal, region)
      if (length(overlapping) == 0) {
        return(NA)
      }

      # Weighted mean by overlap width
      ol <- pintersect(overlapping, region)
      overlap_widths <- width(ol)
      if (sum(overlap_widths) == 0) return(NA)

      weights <- overlap_widths / sum(overlap_widths)
      weighted_mean <- sum(overlapping$score * weights, na.rm = TRUE)
      return(weighted_mean)
    })

    # Handle case where some regions were filtered
    if (length(scores) < length(regions)) {
      full_scores <- rep(NA, length(regions))
      full_scores[1:length(scores)] <- scores
      return(full_scores)
    }

    return(scores)

  }, error = function(e) {
    warning(paste("Error computing signal:", e$message))
    return(rep(NA, length(regions)))
  })
}

# ===================== Main Analysis =====================

cat("=== Module 3: meDIP Integration ===\n\n")

# 1. Load Peak Categories
cat("Loading peak categories...\n")

# Check if peak files exist
if (!file.exists(SHARED_PEAKS)) {
  stop(paste("Peak file not found:", SHARED_PEAKS,
             "\nPlease run Module 1 (01_peak_classification.R) first."))
}

peaks_shared <- load_peaks(SHARED_PEAKS)
peaks_tes <- load_peaks(TES_UNIQUE)
peaks_tead1 <- load_peaks(TEAD1_UNIQUE)

cat(sprintf("  Shared peaks: %d\n", length(peaks_shared)))
cat(sprintf("  TES_Unique peaks: %d\n", length(peaks_tes)))
cat(sprintf("  TEAD1_Unique peaks: %d\n", length(peaks_tead1)))

# Combine into one list for easier processing
all_peaks <- list(
  Shared = peaks_shared,
  TES_Unique = peaks_tes,
  TEAD1_Unique = peaks_tead1
)

# 2. Compute Methylation Levels
cat("\nComputing methylation levels at binding sites...\n")

results_list <- list()

for (cat_name in names(all_peaks)) {
  p <- all_peaks[[cat_name]]
  cat(sprintf("\nProcessing %s (%d peaks)...\n", cat_name, length(p)))

  if (length(p) == 0) {
    cat("  Skipping (no peaks)\n")
    next
  }

  # Compute signal for TES and GFP conditions
  tes_meth <- compute_signal(p, TES_MEDIP_BW, "TES meDIP")
  gfp_meth <- compute_signal(p, GFP_MEDIP_BW, "GFP meDIP")

  # Create data frame
  df <- data.frame(
    chr = as.character(seqnames(p)),
    start = start(p),
    end = end(p),
    Category = cat_name,
    TES_Methylation = tes_meth,
    GFP_Methylation = gfp_meth,
    stringsAsFactors = FALSE
  )

  # Calculate difference (TES - GFP)
  df$Methylation_Diff <- df$TES_Methylation - df$GFP_Methylation

  results_list[[cat_name]] <- df
}

# Combine all results
results_df <- do.call(rbind, results_list)
rownames(results_df) <- NULL

# Save raw data
write.csv(results_df, file.path(OUTPUT_DIR, "Methylation_at_Binding_Sites.csv"),
          row.names = FALSE)

cat(sprintf("\nTotal binding sites analyzed: %d\n", nrow(results_df)))
cat(sprintf("Sites with valid TES methylation: %d\n",
            sum(!is.na(results_df$TES_Methylation))))
cat(sprintf("Sites with valid GFP methylation: %d\n",
            sum(!is.na(results_df$GFP_Methylation))))

# 3. Visualization: Boxplot of Methylation
cat("\nGenerating methylation boxplots...\n")

# Filter for valid values
plot_data <- results_df[!is.na(results_df$TES_Methylation), ]

if (nrow(plot_data) > 0) {
  pdf(file.path(OUTPUT_DIR, "Methylation_Boxplot.pdf"), width = 10, height = 6)

  # Plot 1: TES methylation by category
  p1 <- ggplot(plot_data, aes(x = Category, y = TES_Methylation, fill = Category)) +
    geom_boxplot(outlier.shape = NA) +
    coord_cartesian(ylim = quantile(plot_data$TES_Methylation,
                                    c(0.05, 0.95), na.rm = TRUE)) +
    theme_minimal(base_size = 14) +
    labs(
      title = "DNA Methylation at Binding Sites (TES Condition)",
      y = "Mean Methylation Signal (RPKM)",
      x = ""
    ) +
    scale_fill_manual(values = c("Shared" = "#984EA3",
                                 "TES_Unique" = "#E41A1C",
                                 "TEAD1_Unique" = "#377EB8")) +
    theme(legend.position = "none")
  print(p1)

  # Plot 2: Methylation difference (TES - GFP)
  plot_data_diff <- plot_data[!is.na(plot_data$Methylation_Diff), ]
  if (nrow(plot_data_diff) > 0) {
    p2 <- ggplot(plot_data_diff, aes(x = Category, y = Methylation_Diff,
                                     fill = Category)) +
      geom_boxplot(outlier.shape = NA) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      coord_cartesian(ylim = quantile(plot_data_diff$Methylation_Diff,
                                      c(0.05, 0.95), na.rm = TRUE)) +
      theme_minimal(base_size = 14) +
      labs(
        title = "Methylation Change at Binding Sites (TES vs GFP)",
        y = "Methylation Difference (TES - GFP)",
        x = ""
      ) +
      scale_fill_manual(values = c("Shared" = "#984EA3",
                                   "TES_Unique" = "#E41A1C",
                                   "TEAD1_Unique" = "#377EB8")) +
      theme(legend.position = "none")
    print(p2)
  }

  dev.off()
} else {
  warning("No valid methylation data to plot.")
}

# 4. DMR Overlap Analysis
cat("\nAnalyzing DMR overlap...\n")

if (file.exists(DMR_FILE)) {
  # Load DMRs from CSV
  dmr_df <- read.csv(DMR_FILE)

  # Check for required columns (the file uses 'stop' not 'end')
  required_cols <- c("chr", "start", "stop")
  if (all(required_cols %in% colnames(dmr_df))) {

    # Create GRanges - handle chromosome naming
    dmr_chr <- dmr_df$chr
    if (!grepl("^chr", dmr_chr[1])) {
      dmr_chr <- paste0("chr", dmr_chr)
    }

    dmrs <- GRanges(
      seqnames = dmr_chr,
      ranges = IRanges(start = dmr_df$start, end = dmr_df$stop)
    )

    # Add logFC if available
    if ("logFC" %in% colnames(dmr_df)) {
      dmrs$logFC <- dmr_df$logFC
    }

    cat(sprintf("  Loaded %d DMRs\n", length(dmrs)))

    # Calculate overlap statistics
    overlap_stats <- data.frame()

    for (cat_name in names(all_peaks)) {
      p <- all_peaks[[cat_name]]

      # Ensure same chromosome style
      seqlevelsStyle(p) <- "UCSC"

      hits <- findOverlaps(p, dmrs)
      n_overlapping <- length(unique(queryHits(hits)))
      pct <- n_overlapping / length(p) * 100

      overlap_stats <- rbind(overlap_stats, data.frame(
        Category = cat_name,
        Total_Peaks = length(p),
        Peaks_in_DMR = n_overlapping,
        DMR_Overlap_Pct = pct
      ))

      cat(sprintf("  %s: %d/%d peaks overlap DMRs (%.1f%%)\n",
                  cat_name, n_overlapping, length(p), pct))
    }

    write.csv(overlap_stats, file.path(OUTPUT_DIR, "DMR_Overlap_Stats.csv"),
              row.names = FALSE)

    # Plot
    pdf(file.path(OUTPUT_DIR, "DMR_Overlap_Barplot.pdf"), width = 8, height = 6)
    p_dmr <- ggplot(overlap_stats, aes(x = Category, y = DMR_Overlap_Pct,
                                       fill = Category)) +
      geom_bar(stat = "identity") +
      geom_text(aes(label = sprintf("%.1f%%", DMR_Overlap_Pct)), vjust = -0.5) +
      theme_minimal(base_size = 14) +
      labs(
        title = "Binding Sites Overlapping DMRs",
        subtitle = "TES vs GFP DMRs (FDR < 0.05)",
        y = "% of Peaks in DMRs",
        x = ""
      ) +
      scale_fill_manual(values = c("Shared" = "#984EA3",
                                   "TES_Unique" = "#E41A1C",
                                   "TEAD1_Unique" = "#377EB8")) +
      theme(legend.position = "none")
    print(p_dmr)
    dev.off()

  } else {
    warning(paste("DMR CSV missing required columns:",
                  paste(setdiff(required_cols, colnames(dmr_df)), collapse = ", ")))
  }
} else {
  warning("DMR file not found. Skipping DMR overlap analysis.")
  cat("  Expected file:", DMR_FILE, "\n")
}

# 5. Summary Statistics
cat("\n=== Summary Statistics ===\n")
for (cat_name in unique(results_df$Category)) {
  cat_data <- results_df[results_df$Category == cat_name, ]
  cat(sprintf("\n%s:\n", cat_name))
  cat(sprintf("  N peaks: %d\n", nrow(cat_data)))
  cat(sprintf("  Mean TES methylation: %.2f\n",
              mean(cat_data$TES_Methylation, na.rm = TRUE)))
  cat(sprintf("  Mean GFP methylation: %.2f\n",
              mean(cat_data$GFP_Methylation, na.rm = TRUE)))
  cat(sprintf("  Mean difference: %.2f\n",
              mean(cat_data$Methylation_Diff, na.rm = TRUE)))
}

cat("\n\nAnalysis complete. Results saved to", OUTPUT_DIR, "\n")
