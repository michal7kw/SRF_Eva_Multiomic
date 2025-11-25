#!/usr/bin/env Rscript
#
# ENHANCED GSEA VISUALIZATIONS
# Creates publication-quality plots with advanced visualization packages
#
# New packages used:
# - ggrepel: Better text label positioning
# - patchwork: Multi-panel figure composition
# - ggridges: Ridge plots for pathway distributions
# - ggsci: Scientific journal color palettes
# - ggpubr: Publication-ready themes and statistical comparisons
# - ComplexHeatmap: Advanced heatmaps with annotations
# - circlize: Color mapping functions
# - cowplot: Enhanced ggplot2 themes

suppressPackageStartupMessages({
  library(fgsea)
  library(dplyr)
  library(ggplot2)
  library(ggrepel) # Better label placement
  library(patchwork) # Combine plots elegantly
  library(ggridges) # Ridge plots
  library(ggsci) # Science journal color palettes
  library(ggpubr) # Publication themes
  library(ComplexHeatmap)
  library(circlize)
  library(cowplot)
  library(scales)
  library(tidyr)
})

setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_1")

cat("=== ENHANCED GSEA VISUALIZATIONS ===\n")
cat("Creating publication-quality figures\n\n")

# =============================================================================
# PATH CONFIGURATION
# =============================================================================
# Input paths from GSEA analysis (script 01)
INPUT_BASE <- "output/01_true_gsea_analysis"

# Output paths
OUTPUT_BASE <- file.path(INPUT_BASE, "16_enhanced_plots")
dir.create(OUTPUT_BASE, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# LOAD PREVIOUS GSEA RESULTS
# =============================================================================

cat("Loading GSEA results...\n")
cat(sprintf("Input directory: %s\n", INPUT_BASE))
fgsea_results <- read.csv(file.path(INPUT_BASE, "fgsea_all_pathways.csv"))
fgsea_sig <- read.csv(file.path(INPUT_BASE, "fgsea_significant_pathways.csv"))
fgsea_cancer <- read.csv(file.path(INPUT_BASE, "fgsea_cancer_pathways.csv"))

cat(sprintf(
  "Loaded %d total pathways, %d significant\n",
  nrow(fgsea_results), nrow(fgsea_sig)
))

# =============================================================================
# PUBLICATION-QUALITY THEME
# =============================================================================

theme_publication <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      # Text
      plot.title = element_text(face = "bold", size = rel(1.3), hjust = 0.5),
      plot.subtitle = element_text(size = rel(1.1), hjust = 0.5, margin = margin(b = 10)),
      axis.title = element_text(face = "bold", size = rel(1.1)),
      axis.text = element_text(size = rel(0.9)),
      legend.title = element_text(face = "bold", size = rel(1.0)),
      legend.text = element_text(size = rel(0.9)),

      # Panel
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),

      # Legend
      legend.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white", color = NA),

      # Strip (facet labels)
      strip.background = element_rect(fill = "grey95", color = "black", linewidth = 0.8),
      strip.text = element_text(face = "bold", size = rel(1.0))
    )
}

# =============================================================================
# PLOT 1: VOLCANO PLOT WITH LABELS (REPLACES BUBBLE PLOT)
# =============================================================================

cat("\n1. Creating enhanced volcano plot...\n")

if (nrow(fgsea_sig) > 0) {
  # Add significance labels
  fgsea_sig$significance <- ifelse(abs(fgsea_sig$NES) > 2 & fgsea_sig$padj < 0.001, "High",
    ifelse(abs(fgsea_sig$NES) > 1.5 & fgsea_sig$padj < 0.01, "Medium", "Low")
  )

  # Select top pathways to label
  top_pathways <- fgsea_sig %>%
    arrange(padj) %>%
    head(20)

  fgsea_sig$label <- ifelse(fgsea_sig$pathway %in% top_pathways$pathway,
    fgsea_sig$pathway, ""
  )

  p1 <- ggplot(fgsea_sig, aes(x = NES, y = -log10(padj))) +
    geom_point(aes(color = NES, size = size), alpha = 0.6) +
    geom_point(
      data = subset(fgsea_sig, label != ""),
      aes(color = NES, size = size), alpha = 0.9, shape = 21,
      stroke = 1.5, fill = NA
    ) +
    geom_text_repel(aes(label = label),
      size = 3, max.overlaps = 20,
      box.padding = 0.5, segment.color = "grey50"
    ) +
    geom_vline(xintercept = c(-1.5, 1.5), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "grey40") +
    scale_color_gradient2(
      low = "#2166AC", mid = "grey90", high = "#B2182B",
      midpoint = 0, name = "NES"
    ) +
    scale_size_continuous(range = c(1, 8), name = "Gene Set Size") +
    labs(
      title = "GSEA Enrichment Landscape",
      subtitle = "Top 20 pathways labeled | Dashed lines: |NES| = 1.5, FDR = 0.01",
      x = "Normalized Enrichment Score (NES)",
      y = "-log10(Adjusted p-value)"
    ) +
    theme_publication()

  ggsave(file.path(OUTPUT_BASE, "01_volcano_plot.pdf"), p1,
    width = 14, height = 10, device = cairo_pdf
  )
  ggsave(file.path(OUTPUT_BASE, "01_volcano_plot.png"), p1,
    width = 14, height = 10, dpi = 300
  )
}

# =============================================================================
# PLOT 2: LOLLIPOP PLOT FOR TOP PATHWAYS (BETTER THAN BAR PLOT)
# =============================================================================

cat("2. Creating lollipop plot for top pathways...\n")

if (nrow(fgsea_cancer) > 0) {
  top_up <- fgsea_cancer %>%
    filter(NES > 0) %>%
    arrange(desc(NES)) %>%
    head(15) %>%
    mutate(direction = "Upregulated")

  top_down <- fgsea_cancer %>%
    filter(NES < 0) %>%
    arrange(NES) %>%
    head(15) %>%
    mutate(direction = "Downregulated")

  top_combined <- rbind(top_up, top_down) %>%
    mutate(pathway_short = substr(pathway, 1, 60)) # Trim long names

  p2 <- ggplot(top_combined, aes(x = NES, y = reorder(pathway_short, NES))) +
    geom_segment(aes(x = 0, xend = NES, yend = pathway_short, color = direction),
      linewidth = 1.5, alpha = 0.8
    ) +
    geom_point(aes(size = -log10(padj), color = direction), alpha = 0.9) +
    scale_color_manual(
      values = c("Upregulated" = "#D6604D", "Downregulated" = "#4393C3"),
      name = "Direction"
    ) +
    scale_size_continuous(range = c(3, 10), name = "-log10(FDR)") +
    geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.8) +
    labs(
      title = "Top Cancer-Relevant Pathways",
      subtitle = "Top 15 enriched and depleted pathways",
      x = "Normalized Enrichment Score",
      y = NULL
    ) +
    theme_publication() +
    theme(axis.text.y = element_text(size = 9))

  ggsave(file.path(OUTPUT_BASE, "02_lollipop_plot.pdf"), p2,
    width = 14, height = 12, device = cairo_pdf
  )
  ggsave(file.path(OUTPUT_BASE, "02_lollipop_plot.png"), p2,
    width = 14, height = 12, dpi = 300
  )
}

# =============================================================================
# PLOT 3: HEATMAP OF PATHWAY ENRICHMENT BY CATEGORY
# =============================================================================

cat("3. Creating category enrichment heatmap...\n")

if (nrow(fgsea_cancer) > 0 && !all(is.na(fgsea_cancer$category))) {
  # Calculate category-level statistics
  category_stats <- fgsea_cancer %>%
    filter(!is.na(category)) %>%
    group_by(category) %>%
    summarise(
      mean_NES = mean(NES),
      median_NES = median(NES),
      n_pathways = n(),
      n_up = sum(NES > 0),
      n_down = sum(NES < 0),
      mean_pval = mean(-log10(padj)),
      .groups = "drop"
    ) %>%
    arrange(desc(mean_NES))

  # Create matrix for heatmap
  mat <- as.matrix(category_stats[, c("mean_NES", "median_NES", "mean_pval")])
  rownames(mat) <- category_stats$category
  colnames(mat) <- c("Mean NES", "Median NES", "Mean -log10(FDR)")

  # Column annotation
  ha_col <- HeatmapAnnotation(
    Metric = colnames(mat),
    col = list(Metric = c(
      "Mean NES" = "#E41A1C",
      "Median NES" = "#377EB8",
      "Mean -log10(FDR)" = "#4DAF4A"
    )),
    show_legend = FALSE
  )

  # Row annotation
  ha_row <- rowAnnotation(
    "N Pathways" = anno_barplot(category_stats$n_pathways,
      gp = gpar(fill = "#8DD3C7")
    ),
    "Up/Down" = anno_barplot(cbind(category_stats$n_up, category_stats$n_down),
      gp = gpar(fill = c("#FB8072", "#80B1D3"))
    ),
    width = unit(4, "cm")
  )

  pdf(file.path(OUTPUT_BASE, "03_category_heatmap.pdf"), width = 12, height = 8)
  ht <- Heatmap(mat,
    name = "Score",
    top_annotation = ha_col,
    right_annotation = ha_row,
    col = colorRamp2(
      c(min(mat), 0, max(mat)),
      c("#2166AC", "white", "#B2182B")
    ),
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    row_names_gp = gpar(fontsize = 11, fontface = "bold"),
    column_names_gp = gpar(fontsize = 11, fontface = "bold"),
    column_title = "Pathway Category Enrichment Statistics",
    column_title_gp = gpar(fontsize = 14, fontface = "bold"),
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(sprintf("%.2f", mat[i, j]), x, y,
        gp = gpar(fontsize = 9, fontface = "bold")
      )
    }
  )
  draw(ht)
  dev.off()

  png(file.path(OUTPUT_BASE, "03_category_heatmap.png"),
    width = 12, height = 8, units = "in", res = 300
  )
  draw(ht)
  dev.off()
}

# =============================================================================
# PLOT 4: RIDGE PLOT FOR NES DISTRIBUTION BY CATEGORY
# =============================================================================

cat("4. Creating ridge plot for NES distribution...\n")

if (nrow(fgsea_cancer) > 0 && !all(is.na(fgsea_cancer$category))) {
  p4 <- ggplot(
    fgsea_cancer %>% filter(!is.na(category)),
    aes(x = NES, y = reorder(category, NES, median), fill = stat(x))
  ) +
    geom_density_ridges_gradient(scale = 2, rel_min_height = 0.01, alpha = 0.8) +
    scale_fill_gradient2(
      low = "#2166AC", mid = "grey90", high = "#B2182B",
      midpoint = 0, name = "NES"
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
    labs(
      title = "Pathway Enrichment Distribution by Category",
      subtitle = "Ridge plot showing NES distribution across functional categories",
      x = "Normalized Enrichment Score (NES)",
      y = "Functional Category"
    ) +
    theme_publication() +
    theme(legend.position = "right")

  ggsave(file.path(OUTPUT_BASE, "04_ridge_plot.pdf"), p4,
    width = 12, height = 10, device = cairo_pdf
  )
  ggsave(file.path(OUTPUT_BASE, "04_ridge_plot.png"), p4,
    width = 12, height = 10, dpi = 300
  )
}

# =============================================================================
# PLOT 5: MULTI-PANEL SUMMARY FIGURE (PATCHWORK)
# =============================================================================

cat("5. Creating comprehensive multi-panel figure...\n")

if (nrow(fgsea_sig) > 0 && nrow(fgsea_cancer) > 0) {
  # Panel A: Summary statistics
  summary_data <- data.frame(
    Category = c(
      "Total Tested", "Significant (FDR<0.05)", "Cancer-Relevant",
      "Upregulated (NES>0)", "Downregulated (NES<0)"
    ),
    Count = c(
      nrow(fgsea_results), nrow(fgsea_sig), nrow(fgsea_cancer),
      sum(fgsea_sig$NES > 0), sum(fgsea_sig$NES < 0)
    )
  )

  pA <- ggplot(summary_data, aes(x = reorder(Category, Count), y = Count, fill = Category)) +
    geom_bar(stat = "identity", color = "black", linewidth = 0.5, alpha = 0.8) +
    geom_text(aes(label = Count), hjust = -0.2, size = 4, fontface = "bold") +
    coord_flip() +
    scale_fill_brewer(palette = "Set2") +
    labs(title = "A) GSEA Overview", x = NULL, y = "Number of Pathways") +
    theme_publication() +
    theme(legend.position = "none") +
    ylim(0, max(summary_data$Count) * 1.1)

  # Panel B: Top pathways
  top10_up <- fgsea_sig %>%
    filter(NES > 0) %>%
    arrange(desc(NES)) %>%
    head(10)
  top10_down <- fgsea_sig %>%
    filter(NES < 0) %>%
    arrange(NES) %>%
    head(10)
  top10 <- rbind(top10_up, top10_down) %>%
    mutate(
      pathway_label = substr(pathway, 1, 40),
      direction = ifelse(NES > 0, "UP", "DOWN")
    )

  pB <- ggplot(top10, aes(x = NES, y = reorder(pathway_label, NES), fill = direction)) +
    geom_col(color = "black", linewidth = 0.3) +
    scale_fill_manual(
      values = c("UP" = "#E41A1C", "DOWN" = "#377EB8"),
      name = "Direction"
    ) +
    labs(title = "B) Top 20 Enriched Pathways", x = "NES", y = NULL) +
    theme_publication() +
    theme(axis.text.y = element_text(size = 8))

  # Panel C: Significance distribution
  pC <- ggplot(fgsea_sig, aes(x = -log10(padj))) +
    geom_histogram(bins = 30, fill = "#377EB8", color = "black", alpha = 0.7) +
    geom_vline(
      xintercept = -log10(0.05), linetype = "dashed",
      color = "red", linewidth = 1
    ) +
    annotate("text",
      x = -log10(0.05), y = Inf, label = "FDR = 0.05",
      vjust = 2, hjust = -0.1, color = "red", fontface = "bold"
    ) +
    labs(
      title = "C) Significance Distribution",
      x = "-log10(Adjusted p-value)", y = "Count"
    ) +
    theme_publication()

  # Panel D: NES vs Size
  pD <- ggplot(fgsea_sig, aes(x = size, y = NES, color = -log10(padj))) +
    geom_point(alpha = 0.6, size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    scale_color_viridis_c(option = "plasma", name = "-log10(FDR)") +
    labs(
      title = "D) Gene Set Size vs Enrichment",
      x = "Gene Set Size", y = "NES"
    ) +
    theme_publication()

  # Combine with patchwork
  combined <- (pA | pB) / (pC | pD) +
    plot_annotation(
      title = "GSEA Analysis Summary - TES vs GFP",
      subtitle = sprintf(
        "Total: %d pathways tested | %d significant (FDR < 0.05)",
        nrow(fgsea_results), nrow(fgsea_sig)
      ),
      theme = theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5)
      )
    )

  ggsave(file.path(OUTPUT_BASE, "05_multipanel_summary.pdf"), combined,
    width = 18, height = 12, device = cairo_pdf
  )
  ggsave(file.path(OUTPUT_BASE, "05_multipanel_summary.png"), combined,
    width = 18, height = 12, dpi = 300
  )
}

# =============================================================================
# COMPLETION
# =============================================================================

cat("\n========================================\n")
cat("ENHANCED VISUALIZATIONS COMPLETE\n")
cat("========================================\n")
cat(sprintf("Input directory: %s\n", INPUT_BASE))
cat(sprintf("Output directory: %s\n", OUTPUT_BASE))
cat("\nGenerated files:\n")
cat("  1. 01_volcano_plot.pdf/png - Labeled volcano plot\n")
cat("  2. 02_lollipop_plot.pdf/png - Top pathways lollipop\n")
cat("  3. 03_category_heatmap.pdf/png - Category enrichment heatmap\n")
cat("  4. 04_ridge_plot.pdf/png - NES distribution ridge plot\n")
cat("  5. 05_multipanel_summary.pdf/png - Comprehensive 4-panel figure\n")
cat("\nAll plots use publication-quality themes and colors\n")
