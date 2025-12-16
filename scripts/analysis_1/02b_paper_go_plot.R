#!/usr/bin/env Rscript
#
# PUBLICATION-READY GO ENRICHMENT PLOT
# Reuses precomputed results from 02_directional_go_enrichment.R
# Generates a cleaner, more readable plot for paper figures

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(scales)
})

setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_1")

cat("=== PUBLICATION-READY GO ENRICHMENT PLOT ===\n")
cat("Generating cleaner combined comparison plot for paper\n")
cat("Started:", as.character(Sys.time()), "\n\n")

# =============================================================================
# CONFIGURATION
# =============================================================================

# Number of top terms to show per direction (reduce for readability)
N_TERMS <- 10

# Maximum character length for term descriptions (truncate longer ones)
MAX_CHARS <- 50

# Output directory
output_dir <- "output/02_directional_go_enrichment"

# =============================================================================
# LOAD PRECOMPUTED RESULTS
# =============================================================================

cat("Loading precomputed GO enrichment results...\n")

up_go_bp <- read.csv(file.path(output_dir, "upregulated_GO_BP.csv"), stringsAsFactors = FALSE)
down_go_bp <- read.csv(file.path(output_dir, "downregulated_GO_BP.csv"), stringsAsFactors = FALSE)

cat(sprintf("  Upregulated GO BP: %d terms\n", nrow(up_go_bp)))
cat(sprintf("  Downregulated GO BP: %d terms\n\n", nrow(down_go_bp)))

# =============================================================================
# PREPARE DATA FOR PLOTTING
# =============================================================================

cat("Preparing data for publication plot...\n")

# Function to prepare data for one direction
prepare_go_data <- function(go_df, direction, n_terms = N_TERMS, max_chars = MAX_CHARS) {
  go_df %>%
    arrange(p.adjust) %>%
    head(n_terms) %>%
    mutate(
      direction = direction,
      # Calculate numeric gene ratio
      GeneRatio_numeric = sapply(
        strsplit(GeneRatio, "/"),
        function(x) as.numeric(x[1]) / as.numeric(x[2])
      ),
      is_significant = p.adjust < 0.05,
      # Truncate long descriptions with word-aware truncation
      Description_short = sapply(Description, function(desc) {
        if (nchar(desc) > max_chars) {
          # Try to break at word boundary
          truncated <- substr(desc, 1, max_chars - 3)
          # Find last space
          last_space <- max(gregexpr(" ", truncated)[[1]])
          if (last_space > max_chars * 0.6) {
            truncated <- substr(truncated, 1, last_space - 1)
          }
          paste0(truncated, "...")
        } else {
          desc
        }
      })
    )
}

up_top <- prepare_go_data(up_go_bp, "Upregulated")
down_top <- prepare_go_data(down_go_bp, "Downregulated")

combined <- rbind(up_top, down_top)

# Set factor levels for proper ordering
combined$direction <- factor(combined$direction, levels = c("Upregulated", "Downregulated"))

cat(sprintf("  Selected %d terms per direction\n\n", N_TERMS))

# =============================================================================
# CREATE PUBLICATION-READY PLOT
# =============================================================================

cat("Creating publication-ready plot...\n")

# Plot dimensions optimized for paper (wider to avoid shrinking)
plot_width <- 16  # inches (wider)
plot_height <- 6  # inches

p <- ggplot(combined, aes(
  x = GeneRatio_numeric,
  y = reorder(Description_short, GeneRatio_numeric),
  color = -log10(p.adjust),
  size = Count,
  shape = is_significant
)) +
  geom_point(stroke = 1.2) +
  facet_wrap(~direction, ncol = 2, scales = "free") +
  scale_color_gradient(
    low = "blue",
    high = "red",
    name = expression(-log[10](FDR))
  ) +
  scale_size_continuous(
    name = "Gene\ncount",
    range = c(3, 8),
    breaks = c(20, 40, 60, 80)
  ) +
  scale_shape_manual(
    values = c("TRUE" = 16, "FALSE" = 1),
    labels = c("TRUE" = "FDR < 0.05", "FALSE" = "FDR >= 0.05"),
    name = "Significance"
  ) +
  scale_x_continuous(
    n.breaks = 4,
    expand = expansion(mult = c(0.05, 0.15))  # Add space on right for bubbles
  ) +
  labs(
    x = "Gene ratio",
    y = NULL
  ) +
  theme_bw(base_size = 14) +
  theme(
    # Panel and strip styling - ensure full text visibility
    strip.background = element_rect(fill = "grey90", color = "grey50"),
    strip.text = element_text(face = "bold", size = 13, margin = margin(b = 5, t = 5)),
    strip.clip = "off",

    # Axis styling
    axis.text.y = element_text(size = 11, color = "black"),
    axis.text.x = element_text(size = 11, color = "black"),
    axis.title.x = element_text(size = 13, face = "bold", margin = margin(t = 10)),

    # Legend styling
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.8, "cm"),
    legend.position = "right",
    legend.box = "vertical",

    # Panel spacing
    panel.spacing = unit(2, "lines"),

    # Grid
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),

    # Plot margins
    plot.margin = margin(t = 15, r = 10, b = 10, l = 10)
  ) +
  guides(
    shape = guide_legend(order = 1),
    color = guide_colorbar(order = 2),
    size = guide_legend(order = 3)
  )

# Save as PDF
pdf(file.path(output_dir, "04_combined_comparison_paper.pdf"),
    width = plot_width, height = plot_height)
print(p)
dev.off()

cat(sprintf("  Saved: %s/04_combined_comparison_paper.pdf\n", output_dir))

# Also save as high-resolution PNG for easy preview
png(file.path(output_dir, "04_combined_comparison_paper.png"),
    width = plot_width, height = plot_height, units = "in", res = 300)
print(p)
dev.off()

cat(sprintf("  Saved: %s/04_combined_comparison_paper.png\n\n", output_dir))

# =============================================================================
# ALTERNATIVE: EVEN MORE COMPACT VERSION
# =============================================================================

cat("Creating compact version (8 terms per direction)...\n")

up_compact <- prepare_go_data(up_go_bp, "Upregulated", n_terms = 8, max_chars = 45)
down_compact <- prepare_go_data(down_go_bp, "Downregulated", n_terms = 8, max_chars = 45)

combined_compact <- rbind(up_compact, down_compact)
combined_compact$direction <- factor(combined_compact$direction, levels = c("Upregulated", "Downregulated"))

p_compact <- ggplot(combined_compact, aes(
  x = GeneRatio_numeric,
  y = reorder(Description_short, GeneRatio_numeric),
  color = -log10(p.adjust),
  size = Count,
  shape = is_significant
)) +
  geom_point(stroke = 1.1) +
  facet_wrap(~direction, ncol = 2, scales = "free") +
  scale_color_gradient(
    low = "blue",
    high = "red",
    name = expression(-log[10](FDR))
  ) +
  scale_size_continuous(
    name = "Gene\ncount",
    range = c(3, 7),
    breaks = c(20, 40, 60)
  ) +
  scale_shape_manual(
    values = c("TRUE" = 16, "FALSE" = 1),
    labels = c("TRUE" = "FDR < 0.05", "FALSE" = "FDR >= 0.05"),
    name = "Significance"
  ) +
  scale_x_continuous(
    n.breaks = 4,
    expand = expansion(mult = c(0.05, 0.15))  # Add space on right for bubbles
  ) +
  labs(
    x = "Gene ratio",
    y = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    strip.background = element_rect(fill = "grey90", color = "grey50"),
    strip.text = element_text(face = "bold", size = 12, margin = margin(b = 4, t = 4)),
    strip.clip = "off",
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 8)),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    legend.key.size = unit(0.6, "cm"),
    legend.position = "right",
    panel.spacing = unit(1.8, "lines"),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 12, r = 8, b = 8, l = 8)
  ) +
  guides(
    shape = guide_legend(order = 1),
    color = guide_colorbar(order = 2),
    size = guide_legend(order = 3)
  )

# Save compact version (wider)
pdf(file.path(output_dir, "04_combined_comparison_compact.pdf"),
    width = 14, height = 5)
print(p_compact)
dev.off()

png(file.path(output_dir, "04_combined_comparison_compact.png"),
    width = 14, height = 5, units = "in", res = 300)
print(p_compact)
dev.off()

cat(sprintf("  Saved: %s/04_combined_comparison_compact.pdf\n", output_dir))
cat(sprintf("  Saved: %s/04_combined_comparison_compact.png\n\n", output_dir))

# =============================================================================
# SUMMARY
# =============================================================================

cat("========================================\n")
cat("PUBLICATION PLOT GENERATION COMPLETE\n")
cat("========================================\n")
cat("Output files:\n")
cat("  1. 04_combined_comparison_paper.pdf  (16x6 inches, 10 terms/direction)\n")
cat("  2. 04_combined_comparison_paper.png  (300 DPI)\n")
cat("  3. 04_combined_comparison_compact.pdf (14x5 inches, 8 terms/direction)\n")
cat("  4. 04_combined_comparison_compact.png (300 DPI)\n\n")
cat("Note: Filled circles = significant (FDR < 0.05)\n")
cat("      Open circles = not significant (FDR >= 0.05)\n\n")
cat("Completed:", as.character(Sys.time()), "\n")
