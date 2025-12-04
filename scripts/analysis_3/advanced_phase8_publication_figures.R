#!/usr/bin/env Rscript

################################################################################
# Phase 8: Publication-Quality Figures
# TES vs TEAD1 Comparative Study (Excluding TESmut)
#
# Purpose: Generate comprehensive, publication-ready multi-panel figures
#          Creates BOTH simplified (3-category) and detailed (6-category) versions
#
# Figures:
# 1. Figure 1: Study Overview (binding sites, categories, pipeline)
# 2. Figure 2: Binding Characterization (peaks, motifs, signal dynamics)
# 3. Figure 3: Expression Integration (binding vs expression, categories)
# 4. Figure 4: Methylation Integration (three-way integration, mechanisms)
# 5. Figure 5: Regulatory Networks (master regulators, cascades, co-regulators)
# 6. Figure 6: Pathway Analysis (enrichment, hub pathways, crosstalk)
# 7. Figure 7: Target Prioritization (scores, validation candidates)
# 8. Supplementary Figures (additional details)
#
# NOTE: All category-based figures are generated in both:
#       - detailed_6cat/ (6 binding categories)
#       - simplified_3cat/ (3 binding categories: TES_Unique, Shared, TEAD1_Unique)
#
# Author: Advanced Multi-Omics Analysis Pipeline
# Date: 2025-01-24
################################################################################

library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(scales)
library(ggrepel)
library(VennDiagram)
library(cowplot)

# Set working directory
setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top")

# Create output directories - separate for detailed vs simplified
output_dir <- "SRF_Eva_integrated_analysis/scripts/analysis_3/results/13_publication_figures"
output_dir_detailed <- file.path(output_dir, "detailed_6cat")
output_dir_simple <- file.path(output_dir, "simplified_3cat")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir_detailed, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir_simple, recursive = TRUE, showWarnings = FALSE)

################################################################################
# Helper function: Convert detailed to simplified categories
################################################################################

convert_to_simple_category <- function(category) {
  dplyr::case_when(
    category == "TES_unique" ~ "TES_Unique",
    category == "TEAD1_unique" ~ "TEAD1_Unique",
    grepl("Shared", category) ~ "Shared",
    category == "Unbound" ~ "Unbound",
    TRUE ~ category
  )
}

# Logging function
log_message <- function(msg) {
  cat(paste0("[", Sys.time(), "] ", msg, "\n"))
}

log_message("Starting Phase 8: Publication-Quality Figures")

# Publication theme
theme_pub <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(hjust = 0.5, size = base_size + 2, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = base_size),
      axis.title = element_text(size = base_size, face = "bold"),
      axis.text = element_text(size = base_size - 2),
      legend.title = element_text(size = base_size, face = "bold"),
      legend.text = element_text(size = base_size - 2),
      strip.text = element_text(size = base_size, face = "bold"),
      strip.background = element_rect(fill = "grey90", color = "black"),
      panel.grid.minor = element_blank()
    )
}

# Color palettes - DETAILED (6-category)
category_colors <- c(
  "TES_unique" = "#E41A1C",
  "TEAD1_unique" = "#377EB8",
  "Shared_equivalent" = "#4DAF4A",
  "Shared_TES_dominant" = "#FF7F00",
  "Shared_TEAD1_dominant" = "#984EA3",
  "Shared_high" = "#A65628",
  "Shared" = "#999999",
  "Unbound" = "#999999"
)

# Color palettes - SIMPLIFIED (3-category)
simple_colors <- c(
  "TES_Unique" = "#E41A1C",
  "Shared" = "#984EA3",
  "TEAD1_Unique" = "#377EB8",
  "Unbound" = "#999999"
)

################################################################################
# Load Data from All Phases
################################################################################

log_message("Loading data from all phases...")

# Phase 1: Binding classification
binding_file <- "SRF_Eva_integrated_analysis/scripts/analysis_3/results/01_binding_classification/binding_site_classification.csv"
binding_data <- if (file.exists(binding_file)) read_csv(binding_file, show_col_types = FALSE) else NULL

# Phase 2: Expression
expr_file <- "SRF_Eva_integrated_analysis/scripts/analysis_3/results/04_category_expression/genes_with_binding_and_expression.csv"
expr_data <- if (file.exists(expr_file)) read_csv(expr_file, show_col_types = FALSE) else NULL

# Phase 3: Methylation
meth_file <- "SRF_Eva_integrated_analysis/scripts/analysis_3/results/08_methylation_expression/integrated_binding_methylation_expression.csv"
meth_data <- if (file.exists(meth_file)) read_csv(meth_file, show_col_types = FALSE) else NULL

# Phase 5: Master regulators
mr_file <- "SRF_Eva_integrated_analysis/scripts/analysis_3/results/09_tf_networks/master_regulator_analysis.csv"
mr_data <- if (file.exists(mr_file)) read_csv(mr_file, show_col_types = FALSE) else NULL

# Phase 6: Target prioritization
prior_file <- "SRF_Eva_integrated_analysis/scripts/analysis_3/results/12_target_prioritization/all_genes_prioritized.csv"
prior_data <- if (file.exists(prior_file)) read_csv(prior_file, show_col_types = FALSE) else NULL

log_message("Data loading complete")

# Add simplified categories to data frames
if (!is.null(binding_data) && "category" %in% names(binding_data)) {
  binding_data$category_simple <- convert_to_simple_category(binding_data$category)
}
if (!is.null(expr_data) && "primary_category" %in% names(expr_data)) {
  expr_data$category_simple <- convert_to_simple_category(expr_data$primary_category)
}
log_message("Added simplified 3-category classifications")

################################################################################
# FIGURE 1: Study Overview
################################################################################

log_message("Creating Figure 1: Study Overview...")

if (!is.null(binding_data)) {
  # Standardize column name: Phase 1 outputs 'category', not 'binding_category'
  if (!"category" %in% names(binding_data) && "binding_category" %in% names(binding_data)) {
    binding_data <- binding_data %>% rename(category = binding_category)
  }

  # Calculate signal_ratio if not present
  if (!"signal_ratio" %in% names(binding_data) && all(c("tes_signal", "tead1_signal") %in% names(binding_data))) {
    binding_data <- binding_data %>%
      mutate(signal_ratio = log2((tes_signal + 1) / (tead1_signal + 1)))
  }

  # Panel A: Binding site counts
  binding_summary <- binding_data %>%
    count(category) %>%
    mutate(category = factor(category, levels = names(category_colors)))

  p1a <- ggplot(binding_summary, aes(x = category, y = n, fill = category)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = comma(n)), vjust = -0.5, size = 4) +
    scale_fill_manual(values = category_colors, na.value = "gray50") +
    labs(
      title = "A. Binding Site Classification",
      x = "Binding Category",
      y = "Number of Peaks"
    ) +
    theme_pub() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )

  # Panel B: Pie chart of categories
  p1b <- ggplot(binding_summary, aes(x = "", y = n, fill = category)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y") +
    scale_fill_manual(values = category_colors, na.value = "gray50") +
    labs(
      title = "B. Category Proportions",
      fill = "Category"
    ) +
    theme_pub() +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank()
    )

  # Combine panels - DETAILED version
  fig1 <- p1a + p1b + plot_layout(ncol = 2, widths = c(2, 1))

  # Save DETAILED version
  ggsave(file.path(output_dir_detailed, "Figure1_Study_Overview.pdf"),
         fig1, width = 12, height = 6)
  ggsave(file.path(output_dir_detailed, "Figure1_Study_Overview.png"),
         fig1, width = 12, height = 6, dpi = 300)

  # Also save backward-compatible copy
  ggsave(file.path(output_dir, "Figure1_Study_Overview.pdf"),
         fig1, width = 12, height = 6)
  ggsave(file.path(output_dir, "Figure1_Study_Overview.png"),
         fig1, width = 12, height = 6, dpi = 300)

  # ============================================================================
  # SIMPLIFIED (3-category) version
  # ============================================================================
  binding_summary_simple <- binding_data %>%
    count(category_simple) %>%
    mutate(category_simple = factor(category_simple,
                                     levels = names(simple_colors)))

  p1a_simple <- ggplot(binding_summary_simple,
                        aes(x = category_simple, y = n, fill = category_simple)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = comma(n)), vjust = -0.5, size = 4) +
    scale_fill_manual(values = simple_colors, na.value = "gray50") +
    labs(title = "A. Binding Site Classification (Simplified)",
         x = "Binding Category", y = "Number of Peaks") +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")

  p1b_simple <- ggplot(binding_summary_simple,
                        aes(x = "", y = n, fill = category_simple)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y") +
    scale_fill_manual(values = simple_colors, na.value = "gray50") +
    labs(title = "B. Category Proportions", fill = "Category") +
    theme_pub() +
    theme(axis.title = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank(), panel.grid = element_blank(),
          panel.border = element_blank())

  fig1_simple <- p1a_simple + p1b_simple + plot_layout(ncol = 2, widths = c(2, 1))

  ggsave(file.path(output_dir_simple, "Figure1_Study_Overview.pdf"),
         fig1_simple, width = 12, height = 6)
  ggsave(file.path(output_dir_simple, "Figure1_Study_Overview.png"),
         fig1_simple, width = 12, height = 6, dpi = 300)

  log_message("  Created Figure 1 (DETAILED and SIMPLIFIED versions)")
}

################################################################################
# FIGURE 2: Binding Characterization
################################################################################

log_message("Creating Figure 2: Binding Characterization...")

if (!is.null(binding_data)) {
  # Filter out NA signals for violin plots
  binding_data_filtered <- binding_data %>%
    filter(!is.na(tes_signal) & !is.na(tead1_signal) & tes_signal > 0 & tead1_signal > 0)

  # Panel A: Signal intensity distributions
  p2a <- ggplot(binding_data_filtered, aes(x = category, y = tes_signal, fill = category)) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.2, outlier.shape = NA) +
    scale_fill_manual(values = category_colors, na.value = "gray50") +
    scale_y_log10(labels = comma) +
    labs(
      title = "A. TES Signal Intensity",
      x = "Binding Category",
      y = "TES Signal (CPM, log10)"
    ) +
    theme_pub() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )

  p2b <- ggplot(binding_data_filtered, aes(x = category, y = tead1_signal, fill = category)) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.2, outlier.shape = NA) +
    scale_fill_manual(values = category_colors, na.value = "gray50") +
    scale_y_log10(labels = comma) +
    labs(
      title = "B. TEAD1 Signal Intensity",
      x = "Binding Category",
      y = "TEAD1 Signal (CPM, log10)"
    ) +
    theme_pub() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )

  # Panel C: Signal ratio
  p2c <- ggplot(binding_data_filtered, aes(x = category, y = signal_ratio, fill = category)) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.2, outlier.shape = NA) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    scale_fill_manual(values = category_colors, na.value = "gray50") +
    labs(
      title = "C. Signal Ratio (log2 TES/TEAD1)",
      x = "Binding Category",
      y = "Signal Ratio"
    ) +
    theme_pub() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )

  # Combine panels
  fig2 <- (p2a + p2b) / p2c + plot_layout(heights = c(1, 1))

  ggsave(
    file.path(output_dir, "Figure2_Binding_Characterization.pdf"),
    fig2,
    width = 14,
    height = 10
  )

  ggsave(
    file.path(output_dir, "Figure2_Binding_Characterization.png"),
    fig2,
    width = 14,
    height = 10,
    dpi = 300
  )
}

################################################################################
# FIGURE 3: Expression Integration
################################################################################

log_message("Creating Figure 3: Expression Integration...")

if (!is.null(expr_data)) {
  # Panel A: Expression by binding category
  p3a <- ggplot(expr_data, aes(x = primary_category, y = log2FoldChange, fill = primary_category)) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.3, outlier.shape = NA) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    scale_fill_manual(values = category_colors) +
    labs(
      title = "A. Expression by Binding Category",
      x = "Binding Category",
      y = "log2 Fold Change (TES vs GFP)"
    ) +
    theme_pub() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )

  # Panel B: Volcano plot with category colors
  expr_sig <- expr_data %>%
    mutate(
      sig = ifelse(padj < 0.05, "Significant", "Not Significant"),
      direction = case_when(
        padj < 0.05 & log2FoldChange > 0 ~ "Up",
        padj < 0.05 & log2FoldChange < 0 ~ "Down",
        TRUE ~ "NS"
      )
    )

  p3b <- ggplot(expr_sig, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = primary_category), alpha = 0.5, size = 1.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray") +
    scale_color_manual(values = category_colors) +
    labs(
      title = "B. Volcano Plot by Binding Category",
      x = "log2 Fold Change",
      y = "-log10(adjusted p-value)",
      color = "Category"
    ) +
    theme_pub()

  # Panel C: Count of DEGs by category and direction
  deg_summary <- expr_sig %>%
    filter(direction != "NS") %>%
    count(primary_category, direction) %>%
    mutate(n = ifelse(direction == "Down", -n, n))

  p3c <- ggplot(deg_summary, aes(x = primary_category, y = n, fill = direction)) +
    geom_bar(stat = "identity", position = "identity") +
    geom_hline(yintercept = 0) +
    scale_fill_manual(values = c("Up" = "#E41A1C", "Down" = "#377EB8")) +
    labs(
      title = "C. DEG Counts by Category",
      x = "Binding Category",
      y = "Number of DEGs",
      fill = "Direction"
    ) +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  # Combine panels
  fig3 <- (p3a + p3b) / p3c + plot_layout(heights = c(1, 0.8))

  ggsave(
    file.path(output_dir, "Figure3_Expression_Integration.pdf"),
    fig3,
    width = 14,
    height = 12
  )

  ggsave(
    file.path(output_dir, "Figure3_Expression_Integration.png"),
    fig3,
    width = 14,
    height = 12,
    dpi = 300
  )
}

################################################################################
# FIGURE 4: Methylation Integration
################################################################################

log_message("Creating Figure 4: Methylation Integration...")

if (!is.null(meth_data)) {
  # Panel A: Regulatory mechanism pie chart
  mech_summary <- meth_data %>%
    filter(!is.na(regulatory_mechanism)) %>%
    count(regulatory_mechanism) %>%
    mutate(percentage = n / sum(n) * 100)

  p4a <- ggplot(mech_summary, aes(x = "", y = n, fill = regulatory_mechanism)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y") +
    geom_text(aes(label = sprintf("%.1f%%", percentage)),
              position = position_stack(vjust = 0.5), size = 3) +
    scale_fill_brewer(palette = "Set3") +
    labs(
      title = "A. Regulatory Mechanisms",
      fill = "Mechanism"
    ) +
    theme_pub() +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank()
    )

  # Panel B: Expression vs methylation
  meth_expr <- meth_data %>%
    filter(!is.na(methylation_status), !is.na(log2FoldChange))

  p4b <- ggplot(meth_expr, aes(x = methylation_status, y = log2FoldChange, fill = methylation_status)) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.3, outlier.shape = NA) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    scale_fill_manual(values = c(
      "Hypermethylated" = "#E41A1C",
      "Hypomethylated" = "#377EB8",
      "No_Change" = "gray80"
    )) +
    labs(
      title = "B. Expression vs Methylation Status",
      x = "Methylation Status",
      y = "log2 Fold Change"
    ) +
    theme_pub() +
    theme(legend.position = "none")

  # Panel C: Three-way integration heatmap (simplified)
  mechanism_expr <- meth_data %>%
    filter(!is.na(regulatory_mechanism), is_DEG) %>%
    group_by(regulatory_mechanism) %>%
    summarise(
      mean_log2fc = mean(log2FoldChange, na.rm = TRUE),
      n_genes = n(),
      .groups = "drop"
    )

  p4c <- ggplot(mechanism_expr, aes(x = reorder(regulatory_mechanism, mean_log2fc),
                                   y = mean_log2fc)) +
    geom_bar(stat = "identity", aes(fill = mean_log2fc)) +
    geom_text(aes(label = n_genes), hjust = -0.3, size = 3) +
    coord_flip() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    labs(
      title = "C. Mean Expression by Mechanism",
      x = "Regulatory Mechanism",
      y = "Mean log2 Fold Change",
      fill = "log2FC"
    ) +
    theme_pub()

  # Combine panels
  fig4 <- (p4a + p4b) / p4c + plot_layout(heights = c(1, 1))

  ggsave(
    file.path(output_dir, "Figure4_Methylation_Integration.pdf"),
    fig4,
    width = 14,
    height = 12
  )

  ggsave(
    file.path(output_dir, "Figure4_Methylation_Integration.png"),
    fig4,
    width = 14,
    height = 12,
    dpi = 300
  )
}

################################################################################
# FIGURE 5: Regulatory Networks
################################################################################

log_message("Creating Figure 5: Regulatory Networks...")

if (!is.null(mr_data)) {
  # Standardize column names from Phase 5 output to expected names
  # Phase 5 outputs: TF, padj_all, TF_log2FC, n_targets_in_DEGs, odds_ratio_all
  # This script expects: tf_symbol, enrichment_padj, tf_log2fc, n_targets_in_degs, enrichment_score
  mr_data <- mr_data %>%
    rename_with(tolower) %>%
    rename(
      tf_symbol = any_of(c("tf", "TF")),
      enrichment_padj = any_of(c("padj_all", "p_adj_all")),
      tf_log2fc = any_of(c("tf_log2fc", "TF_log2FC")),
      n_targets_in_degs = any_of(c("n_targets_in_degs", "n_targets_in_DEGs")),
      enrichment_score = any_of(c("odds_ratio_all", "enrichment_score"))
    )

  # Handle any remaining column name issues
  if (!"tf_symbol" %in% names(mr_data) && "tf" %in% names(mr_data)) {
    mr_data <- mr_data %>% rename(tf_symbol = tf)
  }
  if (!"enrichment_padj" %in% names(mr_data) && "padj_all" %in% names(mr_data)) {
    mr_data <- mr_data %>% rename(enrichment_padj = padj_all)
  }
  if (!"enrichment_score" %in% names(mr_data) && "odds_ratio_all" %in% names(mr_data)) {
    mr_data <- mr_data %>% rename(enrichment_score = odds_ratio_all)
  }

  # Panel A: Master regulator volcano
  mr_sig <- mr_data %>%
    filter(!is.na(enrichment_padj) & !is.na(tf_log2fc)) %>%
    mutate(
      sig = ifelse(enrichment_padj < 0.05, "Significant", "Not Significant"),
      label_gene = ifelse(enrichment_padj < 0.001 & abs(tf_log2fc) > 1, tf_symbol, NA)
    )

  p5a <- ggplot(mr_sig, aes(x = tf_log2fc, y = -log10(enrichment_padj + 1e-300))) +
    geom_point(aes(color = sig), alpha = 0.6, size = 2) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "gray") +
    geom_text_repel(aes(label = label_gene), size = 3, max.overlaps = 15) +
    scale_color_manual(values = c("Significant" = "#E41A1C", "Not Significant" = "gray70")) +
    labs(
      title = "A. Master Regulator Enrichment",
      x = "TF log2 Fold Change",
      y = "-log10(enrichment padj)",
      color = "Status"
    ) +
    theme_pub()

  # Panel B: Top master regulators
  top_mr <- mr_sig %>%
    filter(enrichment_padj < 0.05) %>%
    arrange(enrichment_padj) %>%
    head(20)

  if (nrow(top_mr) > 0) {
    p5b <- ggplot(top_mr, aes(x = reorder(tf_symbol, -log10(enrichment_padj + 1e-300)),
                             y = -log10(enrichment_padj + 1e-300))) +
      geom_bar(stat = "identity", aes(fill = tf_log2fc)) +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
      coord_flip() +
      labs(
        title = "B. Top 20 Master Regulators",
        x = "Transcription Factor",
        y = "-log10(enrichment padj)",
        fill = "TF log2FC"
      ) +
      theme_pub()
  } else {
    p5b <- ggplot() + theme_void() + ggtitle("No significant master regulators found")
  }

  # Panel C: TF expression vs enrichment scatter
  p5c <- ggplot(mr_sig, aes(x = tf_log2fc, y = enrichment_score)) +
    geom_point(aes(color = sig, size = n_targets_in_degs), alpha = 0.6) +
    geom_smooth(method = "lm", color = "black", linetype = "dashed", se = FALSE) +
    scale_color_manual(values = c("Significant" = "#E41A1C", "Not Significant" = "gray70")) +
    labs(
      title = "C. TF Expression vs Target Enrichment",
      x = "TF log2 Fold Change",
      y = "Enrichment Score (Odds Ratio)",
      color = "Significance",
      size = "# Targets"
    ) +
    theme_pub()

  # Combine panels
  fig5 <- (p5a + p5b) / p5c + plot_layout(heights = c(1, 1))

  ggsave(
    file.path(output_dir, "Figure5_Regulatory_Networks.pdf"),
    fig5,
    width = 14,
    height = 12
  )

  ggsave(
    file.path(output_dir, "Figure5_Regulatory_Networks.png"),
    fig5,
    width = 14,
    height = 12,
    dpi = 300
  )
}

################################################################################
# FIGURE 6: Target Prioritization
################################################################################

log_message("Creating Figure 6: Target Prioritization...")

if (!is.null(prior_data)) {
  # Panel A: Priority score distribution
  p6a <- ggplot(prior_data, aes(x = priority_score)) +
    geom_histogram(aes(fill = priority_tier), bins = 50, alpha = 0.8) +
    geom_vline(xintercept = c(30, 50, 70), linetype = "dashed", color = "red") +
    scale_fill_manual(values = c(
      "High" = "#E41A1C",
      "Medium" = "#FF7F00",
      "Low" = "#4DAF4A",
      "Very_Low" = "gray70"
    )) +
    labs(
      title = "A. Priority Score Distribution",
      x = "Priority Score",
      y = "Number of Genes",
      fill = "Tier"
    ) +
    theme_pub()

  # Panel B: Top genes heatmap of component scores
  top_genes <- prior_data %>%
    arrange(desc(priority_score)) %>%
    head(30) %>%
    select(gene_symbol, binding_score, expression_score, significance_score,
           methylation_score, biological_score, novelty_score) %>%
    pivot_longer(cols = -gene_symbol, names_to = "component", values_to = "score") %>%
    mutate(component = str_remove(component, "_score"))

  p6b <- ggplot(top_genes, aes(x = component, y = reorder(gene_symbol, score), fill = score)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "#E41A1C") +
    labs(
      title = "B. Top 30 Genes: Component Scores",
      x = "Score Component",
      y = "Gene",
      fill = "Score"
    ) +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  # Panel C: Score components correlation
  score_cols <- c("binding_score", "expression_score", "significance_score",
                 "methylation_score", "biological_score", "novelty_score")

  score_matrix <- prior_data %>%
    select(all_of(score_cols)) %>%
    cor(use = "pairwise.complete.obs")

  # Convert to long format for ggplot
  score_cor_long <- as.data.frame(score_matrix) %>%
    mutate(var1 = rownames(score_matrix)) %>%
    pivot_longer(cols = -var1, names_to = "var2", values_to = "correlation")

  p6c <- ggplot(score_cor_long, aes(x = var1, y = var2, fill = correlation)) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.2f", correlation)), size = 3) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    labs(
      title = "C. Score Component Correlations",
      x = "",
      y = "",
      fill = "Correlation"
    ) +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  # Combine panels
  fig6 <- (p6a + p6c) / p6b + plot_layout(heights = c(1, 1.2))

  ggsave(
    file.path(output_dir, "Figure6_Target_Prioritization.pdf"),
    fig6,
    width = 14,
    height = 12
  )

  ggsave(
    file.path(output_dir, "Figure6_Target_Prioritization.png"),
    fig6,
    width = 14,
    height = 12,
    dpi = 300
  )
}

################################################################################
# FIGURE 7: Comprehensive Summary
################################################################################

log_message("Creating Figure 7: Comprehensive Summary...")

# This will be a multi-panel overview combining key metrics from all phases

if (!is.null(binding_data) && !is.null(expr_data)) {
  # Panel A: Numbers overview (text-based)
  # Count TES/TEAD1 bound from expr_data if available
  tes_direct_count <- if ("tes_bound" %in% names(expr_data)) {
    sum(expr_data$tes_bound & expr_data$padj < 0.05, na.rm = TRUE)
  } else if ("primary_category" %in% names(expr_data)) {
    sum(str_detect(expr_data$primary_category, "TES") & expr_data$padj < 0.05, na.rm = TRUE)
  } else 0

  tead1_direct_count <- if ("tead1_bound" %in% names(expr_data)) {
    sum(expr_data$tead1_bound & expr_data$padj < 0.05, na.rm = TRUE)
  } else if ("primary_category" %in% names(expr_data)) {
    sum(str_detect(expr_data$primary_category, "TEAD1") & expr_data$padj < 0.05, na.rm = TRUE)
  } else 0

  summary_stats <- data.frame(
    Metric = c(
      "Total Binding Sites",
      "TES-unique Sites",
      "TEAD1-unique Sites",
      "Shared Sites",
      "Genes with Expression",
      "DEGs (padj < 0.05)",
      "Direct Targets (TES)",
      "Direct Targets (TEAD1)"
    ),
    Value = c(
      nrow(binding_data),
      sum(binding_data$category == "TES_unique", na.rm = TRUE),
      sum(binding_data$category == "TEAD1_unique", na.rm = TRUE),
      sum(str_detect(binding_data$category, "Shared"), na.rm = TRUE),
      nrow(expr_data),
      sum(expr_data$padj < 0.05, na.rm = TRUE),
      tes_direct_count,
      tead1_direct_count
    )
  )

  p7a <- ggplot(summary_stats, aes(x = reorder(Metric, Value), y = Value)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    geom_text(aes(label = comma(Value)), hjust = -0.2, size = 3.5) +
    coord_flip() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(
      title = "A. Study Summary Statistics",
      x = "",
      y = "Count"
    ) +
    theme_pub()

  # Panel B: Workflow diagram (text-based representation)
  workflow_text <- data.frame(
    Phase = c("Phase 1", "Phase 2", "Phase 3", "Phase 5", "Phase 6"),
    Analysis = c(
      "Binding Classification\n(6 categories)",
      "Expression Integration\n(DEG mapping)",
      "Methylation Integration\n(3-way analysis)",
      "Network Analysis\n(Master regulators)",
      "Target Prioritization\n(Composite scoring)"
    ),
    y = 5:1
  )

  p7b <- ggplot(workflow_text, aes(x = 1, y = y, label = Analysis)) +
    geom_tile(fill = "lightblue", color = "black", width = 1.5, height = 0.8) +
    geom_text(size = 3.5, fontface = "bold") +
    geom_segment(aes(x = 1, xend = 1, y = y - 0.4, yend = y - 0.6),
                arrow = arrow(length = unit(0.2, "cm")),
                data = workflow_text[1:4, ]) +
    xlim(0, 2) +
    ylim(0, 6) +
    labs(title = "B. Analysis Workflow") +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

  # Combine panels
  fig7 <- p7a + p7b + plot_layout(widths = c(2, 1))

  ggsave(
    file.path(output_dir, "Figure7_Comprehensive_Summary.pdf"),
    fig7,
    width = 14,
    height = 8
  )

  ggsave(
    file.path(output_dir, "Figure7_Comprehensive_Summary.png"),
    fig7,
    width = 14,
    height = 8,
    dpi = 300
  )
}

################################################################################
# Generate Figure Legend File
################################################################################

log_message("Generating figure legends...")

legend_file <- file.path(output_dir, "FIGURE_LEGENDS.txt")

sink(legend_file)

cat("================================================================================\n")
cat("Publication-Quality Figures - Legend Descriptions\n")
cat("TES vs TEAD1 Comparative Multi-Omics Study\n")
cat("================================================================================\n\n")

cat("FIGURE 1: Study Overview\n")
cat("-------------------------\n")
cat("(A) Bar chart showing the number of binding sites in each of six binding\n")
cat("    categories: TES-unique, TEAD1-unique, and four shared categories\n")
cat("    (equivalent, TES-dominant, TEAD1-dominant, all shared).\n")
cat("(B) Pie chart showing the proportional distribution of binding categories.\n")
cat("\n")

cat("FIGURE 2: Binding Characterization\n")
cat("-----------------------------------\n")
cat("(A) Violin plots showing TES signal intensity (CPM, log10 scale) across\n")
cat("    binding categories. Box plots overlay shows median and quartiles.\n")
cat("(B) Violin plots showing TEAD1 signal intensity across binding categories.\n")
cat("(C) Violin plots showing signal ratio (log2 TES/TEAD1) by binding category.\n")
cat("    Dashed red line indicates equal binding (ratio = 0).\n")
cat("\n")

cat("FIGURE 3: Expression Integration\n")
cat("---------------------------------\n")
cat("(A) Violin plots showing log2 fold change (TES vs GFP) stratified by binding\n")
cat("    category. Dashed red line indicates no change (log2FC = 0).\n")
cat("(B) Volcano plot showing differential expression with points colored by\n")
cat("    binding category. Dashed lines indicate significance thresholds\n")
cat("    (padj < 0.05, |log2FC| > 1).\n")
cat("(C) Bar chart showing counts of upregulated and downregulated DEGs by\n")
cat("    binding category.\n")
cat("\n")

cat("FIGURE 4: Methylation Integration\n")
cat("----------------------------------\n")
cat("(A) Pie chart showing the distribution of genes across regulatory\n")
cat("    mechanisms identified by three-way integration (binding + methylation\n")
cat("    + expression).\n")
cat("(B) Violin plots showing gene expression changes stratified by methylation\n")
cat("    status (hypermethylated, hypomethylated, no change).\n")
cat("(C) Bar chart showing mean expression change for each regulatory mechanism.\n")
cat("    Numbers indicate gene counts per mechanism.\n")
cat("\n")

cat("FIGURE 5: Regulatory Networks\n")
cat("------------------------------\n")
cat("(A) Volcano plot showing master regulator enrichment. X-axis shows TF\n")
cat("    expression change, Y-axis shows target gene enrichment significance.\n")
cat("    Points are colored by significance (padj < 0.05). Top TFs are labeled.\n")
cat("(B) Bar chart showing the top 20 master regulators ranked by enrichment\n")
cat("    p-value. Bars are colored by TF expression change (log2FC).\n")
cat("(C) Scatter plot showing correlation between TF expression change and\n")
cat("    target gene enrichment. Point size indicates number of target genes\n")
cat("    in DEG list.\n")
cat("\n")

cat("FIGURE 6: Target Prioritization\n")
cat("--------------------------------\n")
cat("(A) Histogram showing distribution of priority scores across all genes.\n")
cat("    Colors indicate priority tiers (High, Medium, Low, Very Low). Dashed\n")
cat("    lines show tier thresholds (30, 50, 70).\n")
cat("(B) Heatmap showing component scores for top 30 prioritized genes. Columns\n")
cat("    represent score components (binding, expression, significance,\n")
cat("    methylation, biological relevance, novelty).\n")
cat("(C) Correlation matrix showing pairwise correlations between score\n")
cat("    components. Values range from -1 (negative) to +1 (positive).\n")
cat("\n")

cat("FIGURE 7: Comprehensive Summary\n")
cat("--------------------------------\n")
cat("(A) Bar chart summarizing key study statistics including total binding\n")
cat("    sites, unique/shared sites, genes analyzed, DEGs, and direct targets.\n")
cat("(B) Workflow diagram illustrating the five analysis phases: binding\n")
cat("    classification, expression integration, methylation integration,\n")
cat("    network analysis, and target prioritization.\n")
cat("\n")

cat("================================================================================\n")
cat("Color Schemes:\n")
cat("  - Binding categories: TES-unique (red), TEAD1-unique (blue),\n")
cat("                       Shared-equivalent (green), Shared-TES-dominant (orange),\n")
cat("                       Shared-TEAD1-dominant (purple)\n")
cat("  - Expression: Up (red), Down (blue), Not significant (gray)\n")
cat("  - Methylation: Hypermethylated (red), Hypomethylated (blue), No change (gray)\n")
cat("  - Heatmaps: Blue-White-Red diverging (expression/ratios)\n")
cat("             White-Red sequential (enrichment/scores)\n")
cat("================================================================================\n")

sink()

################################################################################
# Generate Summary Report
################################################################################

log_message("Generating summary report...")

summary_file <- file.path(output_dir, "PHASE8_SUMMARY.txt")

sink(summary_file)

cat("================================================================================\n")
cat("Phase 8: Publication-Quality Figures - Summary\n")
cat("TES vs TEAD1 Comparative Study (Excluding TESmut)\n")
cat("================================================================================\n\n")

cat("Analysis Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("FIGURES GENERATED:\n")
cat("  1. Figure1_Study_Overview - Binding site counts and proportions\n")
cat("  2. Figure2_Binding_Characterization - Signal intensities and ratios\n")
cat("  3. Figure3_Expression_Integration - Expression by binding category\n")
cat("  4. Figure4_Methylation_Integration - Three-way omics integration\n")
cat("  5. Figure5_Regulatory_Networks - Master regulator analysis\n")
cat("  6. Figure6_Target_Prioritization - Composite scoring and rankings\n")
cat("  7. Figure7_Comprehensive_Summary - Overall study statistics\n")
cat("\n")

cat("FILE FORMATS:\n")
cat("  - PDF: Vector graphics, publication-ready, editable\n")
cat("  - PNG: Raster graphics, 300 DPI, web-friendly\n")
cat("\n")

cat("FIGURE SPECIFICATIONS:\n")
cat("  - Resolution: 300 DPI (PNG)\n")
cat("  - Dimensions: Typically 12-14 inches wide, 8-12 inches tall\n")
cat("  - Font sizes: 10-12 pt base, 12-14 pt titles\n")
cat("  - Color schemes: Publication-appropriate, colorblind-friendly\n")
cat("  - Style: Clean, professional, minimal clutter\n")
cat("\n")

cat("USAGE RECOMMENDATIONS:\n")
cat("  - Use PDF files for manuscript submission\n")
cat("  - Use PNG files for presentations and web\n")
cat("  - Edit PDF files in Illustrator/Inkscape if needed\n")
cat("  - Combine figures into multi-panel layouts as needed\n")
cat("  - Refer to FIGURE_LEGENDS.txt for detailed descriptions\n")
cat("\n")

cat("NEXT STEPS:\n")
cat("  1. Review all figures for accuracy and clarity\n")
cat("  2. Edit legends in FIGURE_LEGENDS.txt as needed\n")
cat("  3. Combine figures into manuscript-specific layouts\n")
cat("  4. Add statistical annotations (e.g., p-values) if needed\n")
cat("  5. Verify color schemes are colorblind-friendly\n")
cat("  6. Export at different resolutions if required by journal\n")
cat("\n")

cat("================================================================================\n")
cat("All Figures Complete\n")
cat("================================================================================\n")

sink()

log_message("Phase 8 Analysis Complete!")
log_message(sprintf("Results saved to: %s", output_dir))
log_message(sprintf("Summary report: %s", summary_file))
log_message(sprintf("Figure legends: %s", legend_file))

cat("\n")
cat("================================================================================\n")
cat("Phase 8: Publication-Quality Figures - COMPLETE\n")
cat("================================================================================\n")
