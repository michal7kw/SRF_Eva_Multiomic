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

# Get available categories with sufficient data (at least 3 observations)
category_counts <- table(peaks_df$category)
valid_categories <- names(category_counts[category_counts >= 3])
message("  Categories with sufficient data: ", paste(valid_categories, collapse = ", "))

# Only run comparison if we have at least 2 valid categories
if (length(valid_categories) >= 2) {
  # Filter data to valid categories
  peaks_valid <- peaks_df[peaks_df$category %in% valid_categories, ]

  # Choose reference group: prefer Shared_equivalent if available, otherwise use first category
  ref_group <- if ("Shared_equivalent" %in% valid_categories) {
    "Shared_equivalent"
  } else {
    valid_categories[1]
  }
  message("  Reference group: ", ref_group)

  # Run comparison with error handling
  stat_results <- tryCatch({
    compare_means(
      log2_signal_ratio ~ category,
      data = peaks_valid,
      method = "wilcox.test",
      ref.group = ref_group
    )
  }, error = function(e) {
    message("  Warning: Statistical comparison failed - ", e$message)
    NULL
  })

  if (!is.null(stat_results)) {
    write.csv(stat_results,
      file.path(OUTPUT_DIR, "signal_ratio_statistics.csv"),
      row.names = FALSE
    )
    message("  Statistical results saved.")
  }
} else {
  message("  Skipping statistical comparison - fewer than 2 valid categories")
}

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
  row.names = FALSE
)

message("Signal dynamics analysis complete!")
