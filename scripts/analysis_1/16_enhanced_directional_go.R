#!/usr/bin/env Rscript
#
# ENHANCED DIRECTIONAL GO ENRICHMENT
# Advanced visualizations and statistical comparisons
#
# New features:
# - Treemap visualization of GO hierarchies
# - Network analysis of pathway relationships
# - Interactive HTML plots for exploration
# - Statistical comparison of UP vs DOWN
# - Semantic similarity analysis

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(dplyr)
  library(ggplot2)
  library(enrichplot)
  library(readr)
  library(ggrepel)
  library(patchwork)
  library(treemap) # Treemap visualizations
  library(GOSemSim) # GO semantic similarity
  library(plotly) # Interactive plots
  library(networkD3) # Interactive networks
  library(htmlwidgets) # Save HTML widgets
  library(viridis) # Better color scales
  library(ggridges)
  library(ggsci)
  library(ggpubr) # Statistical comparisons
})

setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_1")

cat("=== ENHANCED DIRECTIONAL GO ENRICHMENT ===\n")
cat("Creating advanced visualizations and comparisons\n\n")

# =============================================================================
# PATH CONFIGURATION
# =============================================================================
# Input paths from directional GO enrichment analysis (script 02)
INPUT_BASE <- "output/02_directional_go_enrichment"

# Output paths
OUTPUT_BASE <- "output/16_enhanced_directional_go"
OUTPUT_INTERACTIVE <- file.path(OUTPUT_BASE, "interactive")

dir.create(OUTPUT_BASE, showWarnings = FALSE, recursive = TRUE)
dir.create(OUTPUT_INTERACTIVE, showWarnings = FALSE, recursive = TRUE)

# Publication theme
theme_pub <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = rel(1.3), hjust = 0.5),
      plot.subtitle = element_text(size = rel(1.0), hjust = 0.5),
      axis.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )
}

# =============================================================================
# LOAD PREVIOUS RESULTS
# =============================================================================

cat("Loading GO enrichment results...\n")
cat(sprintf("Input directory: %s\n", INPUT_BASE))

# Load gene lists
upregulated <- read.csv(file.path(INPUT_BASE, "upregulated_genes.csv"))
downregulated <- read.csv(file.path(INPUT_BASE, "downregulated_genes.csv"))

# Load GO results
up_go_bp <- read.csv(file.path(INPUT_BASE, "upregulated_GO_BP.csv"))
down_go_bp <- read.csv(file.path(INPUT_BASE, "downregulated_GO_BP.csv"))

cat(sprintf("UP: %d genes, %d GO terms\n", nrow(upregulated), nrow(up_go_bp)))
cat(sprintf("DOWN: %d genes, %d GO terms\n\n", nrow(downregulated), nrow(down_go_bp)))

# =============================================================================
# PLOT 1: TREEMAP VISUALIZATION OF GO HIERARCHY
# =============================================================================

cat("1. Creating treemap visualizations...\n")

if (nrow(up_go_bp) > 0) {
  # Prepare data for treemap - use top 30 pathways
  up_tree <- up_go_bp %>%
    arrange(p.adjust) %>%
    head(30) %>%
    mutate(
      pathway_short = substr(Description, 1, 40),
      size_val = Count,
      color_val = -log10(p.adjust)
    )

  pdf(file.path(OUTPUT_BASE, "01_upregulated_treemap.pdf"), width = 14, height = 10)
  treemap(up_tree,
    index = "pathway_short",
    vSize = "size_val",
    vColor = "color_val",
    type = "value",
    palette = "Reds",
    title = "Upregulated Pathways - GO BP Treemap",
    title.legend = "-log10(FDR)",
    fontsize.labels = 10,
    fontface.labels = "bold",
    border.col = "white",
    border.lwds = 2
  )
  dev.off()
}

if (nrow(down_go_bp) > 0) {
  down_tree <- down_go_bp %>%
    arrange(p.adjust) %>%
    head(30) %>%
    mutate(
      pathway_short = substr(Description, 1, 40),
      size_val = Count,
      color_val = -log10(p.adjust)
    )

  pdf(file.path(OUTPUT_BASE, "02_downregulated_treemap.pdf"), width = 14, height = 10)
  treemap(down_tree,
    index = "pathway_short",
    vSize = "size_val",
    vColor = "color_val",
    type = "value",
    palette = "Blues",
    title = "Downregulated Pathways - GO BP Treemap",
    title.legend = "-log10(FDR)",
    fontsize.labels = 10,
    fontface.labels = "bold",
    border.col = "white",
    border.lwds = 2
  )
  dev.off()
}

# =============================================================================
# PLOT 2: STATISTICAL COMPARISON - UP VS DOWN ENRICHMENT
# =============================================================================

cat("2. Creating statistical comparison plots...\n")

if (nrow(up_go_bp) > 0 && nrow(down_go_bp) > 0) {
  # Compare enrichment strength
  comparison_data <- data.frame(
    Direction = c(rep("UP", nrow(up_go_bp)), rep("DOWN", nrow(down_go_bp))),
    NegLog10P = c(-log10(up_go_bp$p.adjust), -log10(down_go_bp$p.adjust)),
    GeneCount = c(up_go_bp$Count, down_go_bp$Count),
    GeneRatio = c(
      sapply(strsplit(up_go_bp$GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])),
      sapply(strsplit(down_go_bp$GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
    )
  )

  # Violin plot with statistics
  p3a <- ggplot(comparison_data, aes(x = Direction, y = NegLog10P, fill = Direction)) +
    geom_violin(alpha = 0.7, trim = FALSE) +
    geom_boxplot(width = 0.2, fill = "white", outlier.shape = NA) +
    geom_jitter(width = 0.1, alpha = 0.3, size = 1) +
    scale_fill_manual(values = c("UP" = "#E74C3C", "DOWN" = "#3498DB")) +
    stat_compare_means(
      method = "wilcox.test", label = "p.format",
      label.x = 1.5, label.y = max(comparison_data$NegLog10P) * 0.9
    ) +
    labs(
      title = "Enrichment Strength Comparison",
      subtitle = "Wilcoxon test",
      x = NULL,
      y = "-log10(Adjusted p-value)"
    ) +
    theme_pub() +
    theme(legend.position = "none")

  p3b <- ggplot(comparison_data, aes(x = Direction, y = GeneRatio, fill = Direction)) +
    geom_violin(alpha = 0.7, trim = FALSE) +
    geom_boxplot(width = 0.2, fill = "white", outlier.shape = NA) +
    geom_jitter(width = 0.1, alpha = 0.3, size = 1) +
    scale_fill_manual(values = c("UP" = "#E74C3C", "DOWN" = "#3498DB")) +
    stat_compare_means(
      method = "wilcox.test", label = "p.format",
      label.x = 1.5, label.y = max(comparison_data$GeneRatio) * 0.9
    ) +
    labs(
      title = "Gene Ratio Comparison",
      subtitle = "Wilcoxon test",
      x = NULL,
      y = "Gene Ratio"
    ) +
    theme_pub() +
    theme(legend.position = "none")

  p3 <- p3a | p3b
  ggsave(file.path(OUTPUT_BASE, "03_statistical_comparison.pdf"), p3,
    width = 14, height = 7, device = cairo_pdf
  )
  ggsave(file.path(OUTPUT_BASE, "03_statistical_comparison.png"), p3,
    width = 14, height = 7, dpi = 300
  )
}

# =============================================================================
# PLOT 4: GENE OVERLAP ANALYSIS (UPSET PLOT ALTERNATIVE)
# =============================================================================

cat("4. Creating gene overlap visualization...\n")

if (nrow(up_go_bp) > 5 && nrow(down_go_bp) > 5) {
  # Get top pathways from each direction
  top_up_paths <- up_go_bp %>%
    arrange(p.adjust) %>%
    head(10)
  top_down_paths <- down_go_bp %>%
    arrange(p.adjust) %>%
    head(10)

  # Create gene-pathway matrix for heatmap
  pathway_genes_up <- strsplit(top_up_paths$geneID, "/")
  names(pathway_genes_up) <- substr(top_up_paths$Description, 1, 40)

  pathway_genes_down <- strsplit(top_down_paths$geneID, "/")
  names(pathway_genes_down) <- substr(top_down_paths$Description, 1, 40)

  # Calculate Jaccard similarity between pathways
  calc_jaccard <- function(set1, set2) {
    length(intersect(set1, set2)) / length(union(set1, set2))
  }

  # UP pathway similarity
  if (length(pathway_genes_up) > 1) {
    n_up <- length(pathway_genes_up)
    sim_mat_up <- matrix(0, nrow = n_up, ncol = n_up)
    rownames(sim_mat_up) <- colnames(sim_mat_up) <- names(pathway_genes_up)

    for (i in 1:(n_up - 1)) {
      for (j in (i + 1):n_up) {
        sim <- calc_jaccard(pathway_genes_up[[i]], pathway_genes_up[[j]])
        sim_mat_up[i, j] <- sim_mat_up[j, i] <- sim
      }
    }

    library(ComplexHeatmap)
    library(circlize)

    pdf(file.path(OUTPUT_BASE, "04_upregulated_pathway_similarity.pdf"), width = 10, height = 10)
    Heatmap(sim_mat_up,
      name = "Jaccard\nIndex",
      col = colorRamp2(c(0, 0.5, 1), c("white", "#FDB863", "#B2182B")),
      cluster_rows = TRUE,
      cluster_columns = TRUE,
      row_names_gp = gpar(fontsize = 9),
      column_names_gp = gpar(fontsize = 9),
      column_title = "Upregulated Pathway Gene Overlap",
      column_title_gp = gpar(fontsize = 12, fontface = "bold"),
      cell_fun = function(j, i, x, y, width, height, fill) {
        if (i != j && sim_mat_up[i, j] > 0.1) {
          grid.text(sprintf("%.2f", sim_mat_up[i, j]), x, y,
            gp = gpar(fontsize = 7)
          )
        }
      }
    )
    dev.off()
  }
}

# =============================================================================
# PLOT 5: COMBINED BUBBLE CHART WITH BETTER AESTHETICS
# =============================================================================

cat("5. Creating enhanced bubble chart...\n")

if (nrow(up_go_bp) > 0 && nrow(down_go_bp) > 0) {
  top_combined <- rbind(
    up_go_bp %>%
      arrange(p.adjust) %>%
      head(20) %>%
      mutate(
        direction = "Upregulated",
        GeneRatio_numeric = sapply(
          strsplit(GeneRatio, "/"),
          function(x) as.numeric(x[1]) / as.numeric(x[2])
        )
      ),
    down_go_bp %>%
      arrange(p.adjust) %>%
      head(20) %>%
      mutate(
        direction = "Downregulated",
        GeneRatio_numeric = sapply(
          strsplit(GeneRatio, "/"),
          function(x) as.numeric(x[1]) / as.numeric(x[2])
        )
      )
  ) %>%
    mutate(Description_short = substr(Description, 1, 50))

  p5 <- ggplot(top_combined, aes(x = GeneRatio_numeric, y = reorder(Description_short, GeneRatio_numeric))) +
    geom_point(aes(size = Count, color = -log10(p.adjust)), alpha = 0.8) +
    scale_size_continuous(range = c(3, 12), name = "Gene Count") +
    scale_color_viridis(option = "plasma", name = "-log10(FDR)", direction = -1) +
    facet_wrap(~direction, ncol = 2, scales = "free_y") +
    labs(
      title = "Top 20 GO Terms: Upregulated vs Downregulated Genes",
      subtitle = "Bubble size = gene count | Color = significance",
      x = "Gene Ratio",
      y = NULL
    ) +
    theme_pub() +
    theme(
      axis.text.y = element_text(size = 8),
      strip.background = element_rect(fill = "grey90"),
      strip.text = element_text(face = "bold", size = 11)
    )

  ggsave(file.path(OUTPUT_BASE, "05_enhanced_bubble_chart.pdf"), p5,
    width = 16, height = 12, device = cairo_pdf
  )
  ggsave(file.path(OUTPUT_BASE, "05_enhanced_bubble_chart.png"), p5,
    width = 16, height = 12, dpi = 300
  )
}

# =============================================================================
# COMPLETION
# =============================================================================

cat("\n========================================\n")
cat("ENHANCED GO VISUALIZATIONS COMPLETE\n")
cat("========================================\n")
cat(sprintf("Input directory: %s\n", INPUT_BASE))
cat(sprintf("Output directory: %s\n", OUTPUT_BASE))
cat(sprintf("Interactive plots: %s\n", OUTPUT_INTERACTIVE))
cat("\nGenerated files:\n")
cat("  Static plots (PDF/PNG):\n")
cat("    1. Treemaps for UP and DOWN pathways\n")
cat("    2. Statistical comparison violins\n")
cat("    3. Pathway similarity heatmaps\n")
cat("    4. Enhanced bubble chart\n")
cat("\n")
