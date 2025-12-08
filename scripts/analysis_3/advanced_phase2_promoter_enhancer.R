#!/usr/bin/env Rscript

################################################################################
# Phase 2.2: Promoter vs Enhancer Binding Effects
#
# Purpose: Distinguish transcriptional effects of promoter vs distal binding
#          Creates BOTH simplified (3-category) and detailed (6-category) analyses
#
# Author: Advanced Multi-Omics Analysis Plan
# Date: 2025-01-24
################################################################################

message("=== Phase 2.2: Promoter vs Enhancer Binding Effects ===")
message("Start time: ", Sys.time())
message("NOTE: Creating both SIMPLIFIED (3-cat) and DETAILED (6-cat) analyses")

# Load required libraries
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(pheatmap)
  library(RColorBrewer)
})

setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top")

# Input files
PHASE2_1_DIR <- "SRF_Eva_integrated_analysis/scripts/analysis_3/results/04_category_expression"
BINDING_DATA <- "SRF_Eva_integrated_analysis/scripts/analysis_3/results/01_binding_classification/binding_classification_data.RData"

# Output directories - separate for detailed vs simplified
OUTPUT_DIR <- "SRF_Eva_integrated_analysis/scripts/analysis_3/results/05_promoter_enhancer"
OUTPUT_DIR_DETAILED <- file.path(OUTPUT_DIR, "detailed_6cat")
OUTPUT_DIR_SIMPLE <- file.path(OUTPUT_DIR, "simplified_3cat")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTPUT_DIR_DETAILED, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTPUT_DIR_SIMPLE, recursive = TRUE, showWarnings = FALSE)

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
# Step 1: Load data
################################################################################

message("\n[Step 1] Loading data...")

# Load binding classification
load(BINDING_DATA)

# Load gene-expression-binding data from Phase 2.1
gene_data <- read.csv(file.path(PHASE2_1_DIR, "genes_with_binding_and_expression.csv"))
peak_gene_data <- read.csv(file.path(PHASE2_1_DIR, "peak_gene_associations.csv"))

message("  Loaded ", nrow(gene_data), " genes")
message("  Loaded ", nrow(peak_gene_data), " peak-gene associations")

# Add simplified category to peak_gene_data
peak_gene_data$category_simple <- convert_to_simple_category(peak_gene_data$category)
message("  Added simplified 3-category classification")

################################################################################
# Step 2: Stratify by genomic location
################################################################################

message("\n[Step 2] Stratifying peaks by genomic location...")

# Define location types based on distance
peak_gene_data <- peak_gene_data %>%
  mutate(
    location_type = case_when(
      location == "Promoter" ~ "Promoter",
      abs(distance) <= 10000 ~ "Proximal_Enhancer",
      abs(distance) <= 50000 ~ "Distal_Enhancer",
      TRUE ~ "Very_Distal"
    ),
    location_detail = case_when(
      location == "Promoter" ~ "Promoter (±2kb)",
      abs(distance) <= 10000 ~ "Proximal (2-10kb)",
      abs(distance) <= 50000 ~ "Distal (10-50kb)",
      TRUE ~ "Very Distal (>50kb)"
    )
  )

# Count peaks by location
location_summary <- peak_gene_data %>%
  group_by(category, location_type) %>%
  summarise(n_peaks = n(), .groups = "drop")

write.csv(location_summary,
          file.path(OUTPUT_DIR, "peaks_by_location.csv"),
          row.names = FALSE)

################################################################################
# Step 3: Analyze expression by location type
################################################################################

message("\n[Step 3] Analyzing expression changes by location type...")

# For each gene, determine primary binding location
gene_location <- peak_gene_data %>%
  filter(!is.na(log2FoldChange)) %>%
  group_by(gene_name, category) %>%
  # Prioritize promoter peaks, then closest peak
  arrange(location_type != "Promoter", abs(distance)) %>%
  slice(1) %>%
  ungroup()

# Calculate statistics by location and category
location_expression <- gene_location %>%
  group_by(category, location_type) %>%
  summarise(
    n_genes = n(),
    n_DEGs = sum(padj < 0.05, na.rm = TRUE),
    pct_DEG = round(sum(padj < 0.05, na.rm = TRUE) / n() * 100, 2),
    mean_log2FC = mean(log2FoldChange, na.rm = TRUE),
    median_log2FC = median(log2FoldChange, na.rm = TRUE),
    sd_log2FC = sd(log2FoldChange, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(location_expression,
          file.path(OUTPUT_DIR, "expression_by_location.csv"),
          row.names = FALSE)

################################################################################
# Step 4: Distance decay analysis
################################################################################

message("\n[Step 4] Analyzing distance decay effect...")

# Filter for non-promoter peaks with valid expression data
distal_peaks <- peak_gene_data %>%
  filter(location_type != "Promoter",
         !is.na(log2FoldChange),
         abs(distance) <= 100000)  # Within 100kb

# Calculate correlation between distance and expression change
distance_cor_tes <- cor.test(
  abs(distal_peaks$distance[distal_peaks$category %in% c("TES_unique", "Shared_TES_dominant")]),
  abs(distal_peaks$log2FoldChange[distal_peaks$category %in% c("TES_unique", "Shared_TES_dominant")]),
  method = "spearman"
)

distance_cor_tead1 <- cor.test(
  abs(distal_peaks$distance[distal_peaks$category %in% c("TEAD1_unique", "Shared_TEAD1_dominant")]),
  abs(distal_peaks$log2FoldChange[distal_peaks$category %in% c("TEAD1_unique", "Shared_TEAD1_dominant")]),
  method = "spearman"
)

# Save correlation results
sink(file.path(OUTPUT_DIR, "distance_correlation_results.txt"))
cat("Distance-Expression Correlation Analysis\n")
cat("==========================================\n\n")
cat("TES-related peaks:\n")
print(distance_cor_tes)
cat("\n\nTEAD1-related peaks:\n")
print(distance_cor_tead1)
sink()

################################################################################
# Step 5: Multiple peak analysis
################################################################################

message("\n[Step 5] Analyzing genes with multiple peaks...")

# Count peaks per gene
peaks_per_gene <- peak_gene_data %>%
  filter(!is.na(log2FoldChange)) %>%
  group_by(gene_name) %>%
  summarise(
    n_total_peaks = n(),
    n_promoter_peaks = sum(location_type == "Promoter"),
    n_enhancer_peaks = sum(location_type != "Promoter"),
    categories = paste(unique(category), collapse = ";"),
    log2FoldChange = log2FoldChange[1],
    padj = padj[1],
    is_DEG = padj[1] < 0.05,
    .groups = "drop"
  )

# Compare single vs multiple peaks
multiple_peak_stats <- peaks_per_gene %>%
  mutate(peak_group = case_when(
    n_total_peaks == 1 ~ "Single Peak",
    n_total_peaks == 2 ~ "Two Peaks",
    n_total_peaks >= 3 ~ "Multiple Peaks (3+)"
  )) %>%
  group_by(peak_group) %>%
  summarise(
    n_genes = n(),
    n_DEGs = sum(is_DEG),
    pct_DEG = round(sum(is_DEG) / n() * 100, 2),
    mean_abs_log2FC = mean(abs(log2FoldChange), na.rm = TRUE),
    .groups = "drop"
  )

write.csv(peaks_per_gene,
          file.path(OUTPUT_DIR, "genes_with_peak_counts.csv"),
          row.names = FALSE)

write.csv(multiple_peak_stats,
          file.path(OUTPUT_DIR, "multiple_peak_statistics.csv"),
          row.names = FALSE)

# Test if multiple peaks increase effect size
wilcox_test_multiple <- wilcox.test(
  abs(peaks_per_gene$log2FoldChange[peaks_per_gene$n_total_peaks >= 3]),
  abs(peaks_per_gene$log2FoldChange[peaks_per_gene$n_total_peaks == 1])
)

################################################################################
# Step 6: Promoter vs Enhancer comparison
################################################################################

message("\n[Step 6] Comparing promoter vs enhancer effects...")

# For genes with both promoter and enhancer peaks
both_types <- peak_gene_data %>%
  filter(!is.na(log2FoldChange)) %>%
  group_by(gene_name) %>%
  summarise(
    has_promoter = any(location_type == "Promoter"),
    has_enhancer = any(location_type != "Promoter"),
    log2FoldChange = log2FoldChange[1],
    padj = padj[1],
    .groups = "drop"
  ) %>%
  mutate(
    binding_pattern = case_when(
      has_promoter & has_enhancer ~ "Both",
      has_promoter ~ "Promoter Only",
      has_enhancer ~ "Enhancer Only"
    )
  )

# Statistical comparison
binding_pattern_stats <- both_types %>%
  filter(!is.na(binding_pattern)) %>%
  group_by(binding_pattern) %>%
  summarise(
    n_genes = n(),
    n_DEGs = sum(padj < 0.05, na.rm = TRUE),
    pct_DEG = round(sum(padj < 0.05, na.rm = TRUE) / n() * 100, 2),
    mean_log2FC = mean(log2FoldChange, na.rm = TRUE),
    median_log2FC = median(log2FoldChange, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(binding_pattern_stats,
          file.path(OUTPUT_DIR, "promoter_enhancer_comparison.csv"),
          row.names = FALSE)

################################################################################
# Step 7: Visualizations (BOTH detailed and simplified)
################################################################################

message("\n[Step 7] Generating visualizations...")

# Add simplified category to gene_location
gene_location$category_simple <- convert_to_simple_category(gene_location$category)

#' Helper function to create location boxplot
create_location_boxplot <- function(data, category_col, colors, output_file,
                                     title_suffix = "") {
  pdf(output_file, width = 14, height = 8)
  p <- ggplot(data, aes(x = location_type, y = log2FoldChange,
                         fill = .data[[category_col]])) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.3, position = position_dodge(0.9),
                 alpha = 0.5, outlier.shape = NA) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    facet_wrap(as.formula(paste("~", category_col)), ncol = 3) +
    scale_fill_manual(values = colors, na.value = "gray50") +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      strip.text = element_text(face = "bold")
    ) +
    labs(
      title = paste("Gene Expression Changes by Peak Location", title_suffix),
      subtitle = "Stratified by binding category",
      x = "Peak Location Type",
      y = "log2 Fold Change (TES vs GFP)"
    )
  print(p)
  dev.off()
}

# Create DETAILED (6-category) version
create_location_boxplot(
  gene_location, "category", DETAILED_COLORS,
  file.path(OUTPUT_DIR_DETAILED, "expression_by_location_boxplots.pdf"),
  "(Detailed)"
)

# Create SIMPLIFIED (3-category) version
create_location_boxplot(
  gene_location, "category_simple", SIMPLE_COLORS,
  file.path(OUTPUT_DIR_SIMPLE, "expression_by_location_boxplots.pdf"),
  "(Simplified)"
)

# Backward compatible copy
create_location_boxplot(
  gene_location, "category", DETAILED_COLORS,
  file.path(OUTPUT_DIR, "expression_by_location_boxplots.pdf"),
  ""
)

# 7.2: Distance decay plot
pdf(file.path(OUTPUT_DIR, "distance_decay_plot.pdf"), width = 12, height = 10)
p1 <- ggplot(distal_peaks %>% filter(category %in% c("TES_unique", "Shared_TES_dominant")),
             aes(x = abs(distance)/1000, y = abs(log2FoldChange), color = category)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = TRUE) +
  scale_color_brewer(palette = "Reds") +
  scale_x_continuous(breaks = seq(0, 100, 20)) +
  theme_classic() +
  labs(
    title = "Distance Decay: TES-related Peaks",
    subtitle = paste0("Spearman rho = ", round(distance_cor_tes$estimate, 3),
                     ", p = ", format(distance_cor_tes$p.value, digits = 3)),
    x = "Distance to TSS (kb)",
    y = "|log2 Fold Change|",
    color = "Category"
  )

p2 <- ggplot(distal_peaks %>% filter(category %in% c("TEAD1_unique", "Shared_TEAD1_dominant")),
             aes(x = abs(distance)/1000, y = abs(log2FoldChange), color = category)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = TRUE) +
  scale_color_brewer(palette = "Blues") +
  scale_x_continuous(breaks = seq(0, 100, 20)) +
  theme_classic() +
  labs(
    title = "Distance Decay: TEAD1-related Peaks",
    subtitle = paste0("Spearman rho = ", round(distance_cor_tead1$estimate, 3),
                     ", p = ", format(distance_cor_tead1$p.value, digits = 3)),
    x = "Distance to TSS (kb)",
    y = "|log2 Fold Change|",
    color = "Category"
  )

gridExtra::grid.arrange(p1, p2, ncol = 1)
dev.off()

# 7.3: Multiple peak effects
pdf(file.path(OUTPUT_DIR, "multiple_peak_effects.pdf"), width = 10, height = 6)
ggplot(peaks_per_gene %>% filter(n_total_peaks <= 5),
       aes(x = factor(n_total_peaks), y = abs(log2FoldChange))) +
  geom_violin(fill = "skyblue", alpha = 0.7) +
  geom_boxplot(width = 0.3, alpha = 0.5) +
  stat_compare_means(comparisons = list(c("1", "2"), c("1", "3"), c("2", "3")),
                     method = "wilcox.test") +
  theme_classic() +
  labs(
    title = "Effect of Multiple Peaks on Gene Expression",
    subtitle = "Do genes with more peaks show stronger expression changes?",
    x = "Number of Peaks per Gene",
    y = "|log2 Fold Change|"
  )
dev.off()

# 7.4: Promoter vs Enhancer comparison
pdf(file.path(OUTPUT_DIR, "promoter_enhancer_effects.pdf"), width = 10, height = 6)
ggplot(both_types %>% filter(!is.na(binding_pattern)),
       aes(x = binding_pattern, y = log2FoldChange, fill = binding_pattern)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.3, alpha = 0.5, outlier.shape = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_fill_brewer(palette = "Set1") +
  stat_compare_means(comparisons = list(
    c("Promoter Only", "Enhancer Only"),
    c("Promoter Only", "Both"),
    c("Enhancer Only", "Both")
  ), method = "wilcox.test") +
  theme_classic() +
  theme(legend.position = "none") +
  labs(
    title = "Promoter vs Enhancer Binding Effects",
    x = "Binding Pattern",
    y = "log2 Fold Change (TES vs GFP)"
  )
dev.off()

# 7.5: DEG percentage by location
pdf(file.path(OUTPUT_DIR, "DEG_percentage_by_location.pdf"), width = 10, height = 6)
ggplot(location_expression, aes(x = location_type, y = pct_DEG, fill = category)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = paste0(round(pct_DEG, 1), "%")),
            position = position_dodge(width = 0.9), vjust = -0.5, size = 3) +
  scale_fill_brewer(palette = "Set2") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "DEG Percentage by Peak Location and Category",
    x = "Peak Location",
    y = "% of Genes that are DEGs",
    fill = "Category"
  )
dev.off()

################################################################################
# Step 8: Summary report
################################################################################

message("\n[Step 8] Creating summary report...")

sink(file.path(OUTPUT_DIR, "PHASE2_2_SUMMARY.txt"))
cat("=== Phase 2.2: Promoter vs Enhancer Binding Effects ===\n")
cat("Date:", as.character(Sys.time()), "\n\n")

cat("Peak Location Distribution:\n")
cat("============================\n")
print(location_summary)

cat("\n\nExpression Statistics by Location:\n")
cat("====================================\n")
print(location_expression)

cat("\n\nDistance-Expression Correlation:\n")
cat("==================================\n")
cat("TES-related peaks: rho =", round(distance_cor_tes$estimate, 3),
    ", p =", format(distance_cor_tes$p.value, scientific = TRUE), "\n")
cat("TEAD1-related peaks: rho =", round(distance_cor_tead1$estimate, 3),
    ", p =", format(distance_cor_tead1$p.value, scientific = TRUE), "\n")

cat("\n\nMultiple Peak Effects:\n")
cat("=======================\n")
print(multiple_peak_stats)
cat("\nWilcoxon test (3+ peaks vs 1 peak):\n")
cat("p-value =", format(wilcox_test_multiple$p.value, scientific = TRUE), "\n")

cat("\n\nPromoter vs Enhancer Binding:\n")
cat("==============================\n")
print(binding_pattern_stats)

sink()

message("\n=== Analysis Complete ===")
message("Output directory: ", OUTPUT_DIR)
message("End time: ", Sys.time())
