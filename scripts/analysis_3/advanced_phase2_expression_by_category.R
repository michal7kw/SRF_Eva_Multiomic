#!/usr/bin/env Rscript

################################################################################
# Phase 2.1: Category-Specific Expression Analysis
#
# Purpose: Analyze gene expression changes for genes in each binding category
#          Creates BOTH simplified (3-category) and detailed (6-category) analyses
#
# Author: Advanced Multi-Omics Analysis Plan
# Date: 2025-01-24
################################################################################

message("=== Phase 2.1: Category-Specific Expression Analysis ===")
message("Start time: ", Sys.time())
message("NOTE: Creating both SIMPLIFIED (3-category) and DETAILED (6-category) analyses")

# Load required libraries
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(pheatmap)
  library(ComplexHeatmap)
  library(RColorBrewer)
  library(org.Hs.eg.db)
  library(ChIPseeker)
})

# Set working directory
setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top")

# Input files
BINDING_DATA <- "SRF_Eva_integrated_analysis/scripts/analysis_3/results/01_binding_classification/binding_classification_data.RData"
RNA_SEQ_FILE <- "SRF_Eva_RNA/results/05_deseq2/deseq2_results_TES_vs_GFP.txt"

# Output directories - separate for detailed vs simplified
OUTPUT_DIR <- "SRF_Eva_integrated_analysis/scripts/analysis_3/results/04_category_expression"
OUTPUT_DIR_DETAILED <- file.path(OUTPUT_DIR, "detailed_6cat")
OUTPUT_DIR_SIMPLE <- file.path(OUTPUT_DIR, "simplified_3cat")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTPUT_DIR_DETAILED, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTPUT_DIR_SIMPLE, recursive = TRUE, showWarnings = FALSE)

################################################################################
# Helper function: Convert detailed to simplified categories
################################################################################

convert_to_simple_category <- function(category) {
  case_when(
    category == "TES_unique" ~ "TES_Unique",
    category == "TEAD1_unique" ~ "TEAD1_Unique",
    grepl("Shared", category) ~ "Shared",
    category == "Unbound" ~ "Unbound",
    TRUE ~ category
  )
}

# Color palettes
DETAILED_COLORS <- c(
  "TES_unique" = "#E41A1C",
  "TEAD1_unique" = "#377EB8",
  "Shared_high" = "#984EA3",
  "Shared_TES_dominant" = "#FF7F00",
  "Shared_TEAD1_dominant" = "#4DAF4A",
  "Shared_equivalent" = "#A65628",
  "Unbound" = "#999999"
)

SIMPLE_COLORS <- c(
  "TES_Unique" = "#E41A1C",
  "Shared" = "#984EA3",
  "TEAD1_Unique" = "#377EB8",
  "Unbound" = "#999999"
)

################################################################################
# Step 1: Load binding and expression data
################################################################################

message("\n[Step 1] Loading data...")

# Load binding classification
if (!file.exists(BINDING_DATA)) {
  stop("Binding classification data not found! Run Phase 1.1 first.")
}
load(BINDING_DATA)
message("  Loaded binding classification: ", nrow(peaks_df), " peaks")

# Load RNA-seq results
rna_results <- read.delim(RNA_SEQ_FILE, stringsAsFactors = FALSE)
message("  Loaded RNA-seq results: ", nrow(rna_results), " genes")

# Filter for genes with adjusted p-value
rna_results <- rna_results[!is.na(rna_results$padj), ]
message("  Genes with padj values: ", nrow(rna_results))

################################################################################
# Step 2: Map peaks to genes
################################################################################

message("\n[Step 2] Mapping peaks to genes...")

# Function to assign peaks to genes based on proximity
assign_peaks_to_genes <- function(peaks, distance_threshold = 50000) {
  # For promoter peaks, use gene name from ChIPseeker annotation
  promoter_peaks <- peaks[grepl("Promoter", peaks$annotation), ]

  # For distal peaks, find genes within distance threshold
  distal_peaks <- peaks[!grepl("Promoter", peaks$annotation), ]

  # Create gene-peak mappings
  gene_mappings <- data.frame(
    gene_name = character(),
    peak_id = character(),
    category = character(),
    location = character(),
    distance = numeric(),
    tes_signal = numeric(),
    tead1_signal = numeric(),
    stringsAsFactors = FALSE
  )

  # Add promoter mappings
  if (nrow(promoter_peaks) > 0) {
    promoter_mappings <- data.frame(
      gene_name = promoter_peaks$gene_name,
      peak_id = promoter_peaks$peak_id,
      category = promoter_peaks$category,
      location = "Promoter",
      distance = promoter_peaks$distance_to_tss,
      tes_signal = promoter_peaks$tes_signal,
      tead1_signal = promoter_peaks$tead1_signal,
      stringsAsFactors = FALSE
    )
    gene_mappings <- rbind(gene_mappings, promoter_mappings)
  }

  # Add distal mappings within threshold
  if (nrow(distal_peaks) > 0) {
    distal_in_range <- distal_peaks[abs(distal_peaks$distance_to_tss) <= distance_threshold, ]

    if (nrow(distal_in_range) > 0) {
      distal_mappings <- data.frame(
        gene_name = distal_in_range$gene_name,
        peak_id = distal_in_range$peak_id,
        category = distal_in_range$category,
        location = ifelse(
          abs(distal_in_range$distance_to_tss) <= 10000,
          "Proximal_Enhancer",
          "Distal_Enhancer"
        ),
        distance = distal_in_range$distance_to_tss,
        tes_signal = distal_in_range$tes_signal,
        tead1_signal = distal_in_range$tead1_signal,
        stringsAsFactors = FALSE
      )
      gene_mappings <- rbind(gene_mappings, distal_mappings)
    }
  }

  # Remove NA gene names
  gene_mappings <- gene_mappings[!is.na(gene_mappings$gene_name) & gene_mappings$gene_name != "", ]

  return(gene_mappings)
}

# Assign peaks to genes
message("  Assigning peaks to genes (50kb threshold)...")
gene_peak_map <- assign_peaks_to_genes(peaks_df, distance_threshold = 50000)
message("  Created ", nrow(gene_peak_map), " peak-gene associations")

# For genes with multiple peaks, keep the closest peak per category
gene_peak_map <- gene_peak_map %>%
  group_by(gene_name, category) %>%
  arrange(abs(distance)) %>%
  slice(1) %>%
  ungroup()

message("  After filtering closest peaks: ", nrow(gene_peak_map), " associations")

################################################################################
# Step 3: Integrate with RNA-seq data
################################################################################

message("\n[Step 3] Integrating with RNA-seq data...")

# Match gene symbols
gene_peak_map <- gene_peak_map %>%
  left_join(
    rna_results %>% dplyr::select(gene_symbol, baseMean, log2FoldChange, padj),
    by = c("gene_name" = "gene_symbol")
  )

# Remove genes without expression data
gene_peak_map <- gene_peak_map[!is.na(gene_peak_map$baseMean), ]
message("  Genes with expression data: ", length(unique(gene_peak_map$gene_name)))

# Create gene-centric view (one row per gene with all binding info)
gene_categories <- gene_peak_map %>%
  group_by(gene_name) %>%
  summarise(
    n_peaks = n(),
    categories = paste(unique(category), collapse = ";"),
    primary_category = category[which.min(abs(distance))],
    closest_peak_distance = distance[which.min(abs(distance))],
    max_tes_signal = max(tes_signal, na.rm = TRUE),
    max_tead1_signal = max(tead1_signal, na.rm = TRUE),
    has_promoter_peak = any(location == "Promoter"),
    has_enhancer_peak = any(location != "Promoter"),
    baseMean = baseMean[1],
    log2FoldChange = log2FoldChange[1],
    padj = padj[1],
    .groups = "drop"
  )

# Add genes with no peaks (unbound category)
unbound_genes <- rna_results %>%
  dplyr::filter(!gene_symbol %in% gene_categories$gene_name) %>%
  mutate(
    gene_name = gene_symbol,
    n_peaks = 0,
    categories = "Unbound",
    primary_category = "Unbound",
    closest_peak_distance = NA,
    max_tes_signal = 0,
    max_tead1_signal = 0,
    has_promoter_peak = FALSE,
    has_enhancer_peak = FALSE
  ) %>%
  dplyr::select(gene_name, n_peaks, categories, primary_category, closest_peak_distance,
                max_tes_signal, max_tead1_signal, has_promoter_peak, has_enhancer_peak,
                baseMean, log2FoldChange, padj)

# Combine bound and unbound genes
all_genes <- bind_rows(gene_categories, unbound_genes)
message("  Total genes analyzed: ", nrow(all_genes))

# Add DE status
all_genes$is_DEG <- !is.na(all_genes$padj) & all_genes$padj < 0.05
all_genes$regulation <- ifelse(
  !all_genes$is_DEG, "Not_DE",
  ifelse(all_genes$log2FoldChange > 0, "Up", "Down")
)

# Add simplified category (3-category: TES_Unique, Shared, TEAD1_Unique, Unbound)
all_genes$primary_category_simple <- convert_to_simple_category(all_genes$primary_category)
message("  Added simplified 3-category classification")

################################################################################
# Step 4: Statistical comparisons
################################################################################

message("\n[Step 4] Performing statistical comparisons...")

# Get available categories from the data
available_categories <- unique(all_genes$primary_category)
message("  Available categories: ", paste(available_categories, collapse = ", "))

# Subset to available categories with sufficient data (>= 3 genes)
category_counts <- table(all_genes$primary_category)
valid_categories <- names(category_counts[category_counts >= 3])
message("  Categories with >= 3 genes: ", paste(valid_categories, collapse = ", "))

expr_data <- all_genes %>%
  dplyr::filter(primary_category %in% valid_categories)

# Wilcoxon tests comparing each category to Unbound (if Unbound exists)
ref_group <- if ("Unbound" %in% valid_categories) "Unbound" else valid_categories[1]
message("  Reference group for comparison: ", ref_group)

stat_results <- tryCatch({
  compare_means(
    log2FoldChange ~ primary_category,
    data = expr_data,
    method = "wilcox.test",
    ref.group = ref_group
  )
}, error = function(e) {
  message("  Warning: Statistical comparison failed - ", e$message)
  NULL
})

if (!is.null(stat_results)) {
  write.csv(stat_results,
            file.path(OUTPUT_DIR, "expression_statistical_comparison.csv"),
            row.names = FALSE)
  message("  Statistical results saved.")
}

# Chi-square test for proportions of up/down regulated genes
regulation_table <- table(expr_data$primary_category, expr_data$regulation)
chi_test <- chisq.test(regulation_table)

sink(file.path(OUTPUT_DIR, "regulation_proportions_test.txt"))
cat("Chi-square Test: Regulation Proportions by Category\n")
cat("===================================================\n\n")
print(regulation_table)
cat("\n")
print(chi_test)
sink()

################################################################################
# Step 5: Generate visualizations (BOTH detailed and simplified)
################################################################################

message("\n[Step 5] Generating visualizations...")

#' Helper function to create expression boxplot
create_expression_boxplot <- function(data, category_col, colors, output_file,
                                       title_suffix = "") {
  # Get category counts and valid categories
  category_counts <- table(data[[category_col]])
  valid_cats <- names(category_counts[category_counts >= 3])

  plot_data <- data %>%
    dplyr::filter(.data[[category_col]] %in% valid_cats)

  # Reference group
  ref_grp <- if ("Unbound" %in% valid_cats) "Unbound" else valid_cats[1]

  # Build comparisons
  comparisons_list <- list()
  for (cat in setdiff(valid_cats, ref_grp)) {
    comparisons_list <- c(comparisons_list, list(c(cat, ref_grp)))
  }

  pdf(output_file, width = 12, height = 8)
  p <- ggplot(plot_data,
              aes(x = .data[[category_col]], y = log2FoldChange,
                  fill = .data[[category_col]])) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.3, alpha = 0.5, outlier.shape = NA) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    scale_fill_manual(values = colors, na.value = "gray50") +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    ) +
    labs(
      title = paste("Gene Expression Changes by Binding Category", title_suffix),
      subtitle = "TES vs GFP (RNA-seq)",
      x = "Binding Category",
      y = "log2 Fold Change"
    )

  if (length(comparisons_list) > 0 && length(comparisons_list) <= 10) {
    p <- p + stat_compare_means(comparisons = comparisons_list,
                                 method = "wilcox.test")
  }
  print(p)
  dev.off()
}

#' Helper function to create cumulative distribution plot
create_ecdf_plot <- function(data, category_col, colors, output_file,
                              title_suffix = "") {
  category_counts <- table(data[[category_col]])
  valid_cats <- names(category_counts[category_counts >= 3])

  plot_data <- data %>%
    dplyr::filter(.data[[category_col]] %in% valid_cats)

  pdf(output_file, width = 10, height = 8)
  p <- ggplot(plot_data, aes(x = log2FoldChange, color = .data[[category_col]])) +
    stat_ecdf(geom = "step", size = 1) +
    scale_color_manual(values = colors, na.value = "gray50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    theme_classic() +
    labs(
      title = paste("Cumulative Distribution of Expression Changes", title_suffix),
      x = "log2 Fold Change (TES vs GFP)",
      y = "Cumulative Proportion",
      color = "Category"
    )
  print(p)
  dev.off()
}

# ============================================================================
# 5.1: DETAILED (6-category) plots
# ============================================================================
message("  Creating DETAILED (6-category) plots...")

create_expression_boxplot(
  all_genes, "primary_category", DETAILED_COLORS,
  file.path(OUTPUT_DIR_DETAILED, "expression_by_category_boxplots.pdf"),
  "(Detailed)"
)

create_ecdf_plot(
  all_genes, "primary_category", DETAILED_COLORS,
  file.path(OUTPUT_DIR_DETAILED, "expression_cumulative_distribution.pdf"),
  "(Detailed)"
)

# ============================================================================
# 5.2: SIMPLIFIED (3-category) plots
# ============================================================================
message("  Creating SIMPLIFIED (3-category) plots...")

create_expression_boxplot(
  all_genes, "primary_category_simple", SIMPLE_COLORS,
  file.path(OUTPUT_DIR_SIMPLE, "expression_by_category_boxplots.pdf"),
  "(Simplified)"
)

create_ecdf_plot(
  all_genes, "primary_category_simple", SIMPLE_COLORS,
  file.path(OUTPUT_DIR_SIMPLE, "expression_cumulative_distribution.pdf"),
  "(Simplified)"
)

# ============================================================================
# Also save copies to main OUTPUT_DIR for backward compatibility (use detailed)
# ============================================================================
message("  Creating backward-compatible copies in main output directory...")

create_expression_boxplot(
  all_genes, "primary_category", DETAILED_COLORS,
  file.path(OUTPUT_DIR, "expression_by_category_boxplots.pdf"),
  ""
)

create_ecdf_plot(
  all_genes, "primary_category", DETAILED_COLORS,
  file.path(OUTPUT_DIR, "expression_cumulative_distribution.pdf"),
  ""
)

# 5.3: Stacked bar chart of regulation proportions (BOTH versions)
#' Helper function to create regulation barplot
create_regulation_barplot <- function(data, category_col, output_file,
                                       title_suffix = "") {
  reg_counts <- data %>%
    dplyr::filter(is_DEG) %>%
    group_by(.data[[category_col]], regulation) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(.data[[category_col]]) %>%
    mutate(percentage = count / sum(count) * 100)

  pdf(output_file, width = 10, height = 6)
  p <- ggplot(reg_counts,
              aes(x = .data[[category_col]], y = count, fill = regulation)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_text(aes(label = count),
              position = position_dodge(width = 0.9), vjust = -0.5) +
    scale_fill_manual(values = c("Up" = "red", "Down" = "blue")) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = paste("Number of Up/Down-regulated DEGs by Category", title_suffix),
      x = "Binding Category",
      y = "Number of DEGs",
      fill = "Regulation"
    )
  print(p)
  dev.off()
}

# Create DETAILED version
create_regulation_barplot(
  all_genes, "primary_category",
  file.path(OUTPUT_DIR_DETAILED, "regulation_proportions_barplot.pdf"),
  "(Detailed)"
)

# Create SIMPLIFIED version
create_regulation_barplot(
  all_genes, "primary_category_simple",
  file.path(OUTPUT_DIR_SIMPLE, "regulation_proportions_barplot.pdf"),
  "(Simplified)"
)

# Backward compatible copy
create_regulation_barplot(
  all_genes, "primary_category",
  file.path(OUTPUT_DIR, "regulation_proportions_barplot.pdf"),
  ""
)

# 5.4: Scatter plot: binding signal vs expression change
pdf(file.path(OUTPUT_DIR, "binding_signal_vs_expression.pdf"), width = 12, height = 10)
p1 <- ggplot(expr_data %>% dplyr::filter(primary_category != "Unbound"),
             aes(x = log2(max_tes_signal + 1), y = log2FoldChange, color = primary_category)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_color_brewer(palette = "Set2") +
  theme_classic() +
  labs(
    title = "TES Signal vs Expression Change",
    x = "log2(TES Signal + 1)",
    y = "log2 Fold Change",
    color = "Category"
  )

p2 <- ggplot(expr_data %>% dplyr::filter(primary_category != "Unbound"),
             aes(x = log2(max_tead1_signal + 1), y = log2FoldChange, color = primary_category)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_color_brewer(palette = "Set2") +
  theme_classic() +
  labs(
    title = "TEAD1 Signal vs Expression Change",
    x = "log2(TEAD1 Signal + 1)",
    y = "log2 Fold Change",
    color = "Category"
  )

gridExtra::grid.arrange(p1, p2, ncol = 1)
dev.off()

# 5.5: Heatmap of top genes per category
top_genes_per_category <- expr_data %>%
  dplyr::filter(is_DEG) %>%
  group_by(primary_category) %>%
  arrange(padj) %>%
  slice_head(n = 20) %>%
  ungroup()

# Export for visualization
write.csv(top_genes_per_category,
          file.path(OUTPUT_DIR, "top20_DEGs_per_category.csv"),
          row.names = FALSE)

################################################################################
# Step 6: Export results
################################################################################

message("\n[Step 6] Exporting results...")

# Main output: gene classification with binding and expression
write.csv(all_genes,
          file.path(OUTPUT_DIR, "genes_with_binding_and_expression.csv"),
          row.names = FALSE)

# Peak-gene associations
write.csv(gene_peak_map,
          file.path(OUTPUT_DIR, "peak_gene_associations.csv"),
          row.names = FALSE)

# Summary statistics - DETAILED (6-category)
summary_stats <- all_genes %>%
  group_by(primary_category) %>%
  summarise(
    n_genes = n(),
    n_DEGs = sum(is_DEG),
    pct_DEG = round(sum(is_DEG) / n() * 100, 2),
    n_up = sum(regulation == "Up"),
    n_down = sum(regulation == "Down"),
    mean_log2FC = mean(log2FoldChange, na.rm = TRUE),
    median_log2FC = median(log2FoldChange, na.rm = TRUE),
    mean_baseMean = mean(baseMean, na.rm = TRUE),
    .groups = "drop"
  )

# Summary statistics - SIMPLIFIED (3-category)
summary_stats_simple <- all_genes %>%
  group_by(primary_category_simple) %>%
  summarise(
    n_genes = n(),
    n_DEGs = sum(is_DEG),
    pct_DEG = round(sum(is_DEG) / n() * 100, 2),
    n_up = sum(regulation == "Up"),
    n_down = sum(regulation == "Down"),
    mean_log2FC = mean(log2FoldChange, na.rm = TRUE),
    median_log2FC = median(log2FoldChange, na.rm = TRUE),
    mean_baseMean = mean(baseMean, na.rm = TRUE),
    .groups = "drop"
  )

# Save both detailed and simplified summaries
write.csv(summary_stats,
          file.path(OUTPUT_DIR, "category_summary_statistics.csv"),
          row.names = FALSE)
write.csv(summary_stats,
          file.path(OUTPUT_DIR_DETAILED, "category_summary_statistics.csv"),
          row.names = FALSE)
write.csv(summary_stats_simple,
          file.path(OUTPUT_DIR_SIMPLE, "category_summary_statistics.csv"),
          row.names = FALSE)

# Write summary report
sink(file.path(OUTPUT_DIR, "PHASE2_SUMMARY.txt"))
cat("=== Phase 2.1: Category-Specific Expression Analysis ===\n")
cat("Date:", as.character(Sys.time()), "\n\n")

cat("NOTE: Analysis creates BOTH detailed (6-category) and simplified (3-category)\n")
cat("      versions of all plots and statistics.\n\n")

cat("Output directories:\n")
cat("  - Detailed (6-category): ", OUTPUT_DIR_DETAILED, "\n")
cat("  - Simplified (3-category): ", OUTPUT_DIR_SIMPLE, "\n\n")

cat("=== DETAILED (6-category) Summary ===\n")
cat("Categories: TES_unique, TEAD1_unique, Shared_high, Shared_TES_dominant,\n")
cat("            Shared_TEAD1_dominant, Shared_equivalent, Unbound\n\n")
print(summary_stats)

cat("\n\n=== SIMPLIFIED (3-category) Summary ===\n")
cat("Categories: TES_Unique, Shared, TEAD1_Unique, Unbound\n\n")
print(summary_stats_simple)

cat("\n\nStatistical Comparison (vs Unbound):\n")
print(stat_results)
sink()

message("\n=== Analysis Complete ===")
message("Output directories:")
message("  Main: ", OUTPUT_DIR)
message("  Detailed (6-cat): ", OUTPUT_DIR_DETAILED)
message("  Simplified (3-cat): ", OUTPUT_DIR_SIMPLE)
message("End time: ", Sys.time())
