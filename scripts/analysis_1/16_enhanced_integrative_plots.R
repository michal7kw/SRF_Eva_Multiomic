#!/usr/bin/env Rscript
#
# ENHANCED INTEGRATIVE ANALYSIS VISUALIZATIONS
# Multi-omics integration plots (Cut&Tag + RNA-seq)
#
# Creates publication-quality figures showing:
# - Binding vs Expression relationships
# - Direct vs Indirect target classification
# - TES vs TEAD1 specificity
# - Pathway enrichment comparisons

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(ComplexHeatmap)
  library(circlize)
  library(ggsci)
  library(ggpubr)
  library(viridis)
  library(scales)
})

setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_1")

cat("=== ENHANCED INTEGRATIVE VISUALIZATIONS ===\n")
cat("Multi-omics analysis: Cut&Tag + RNA-seq\n\n")

# =============================================================================
# PATH CONFIGURATION
# =============================================================================
# Input paths from integrative analysis (script 10 for DEGs-only)
INPUT_BASE <- "output/10_final_integrative_analysis"
INPUT_DIRECT_TARGETS <- file.path(INPUT_BASE, "direct_targets")
INPUT_REGULATORY <- file.path(INPUT_BASE, "regulatory_networks")

# Output directory
OUTPUT_BASE <- "output/16_enhanced_integrative_plots"
dir.create(OUTPUT_BASE, showWarnings = FALSE, recursive = TRUE)

# Publication theme
theme_pub <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = rel(1.3), hjust = 0.5),
      plot.subtitle = element_text(size = rel(1.0), hjust = 0.5, margin = margin(b = 10)),
      axis.title = element_text(face = "bold", size = rel(1.1)),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      legend.position = "right"
    )
}

# =============================================================================
# LOAD DATA
# =============================================================================

cat("Loading integrative analysis results...\n")
cat(sprintf("Input directory: %s\n", INPUT_BASE))

# Load direct targets (from script 10 - DEGs only version)
tes_direct <- read.csv(file.path(INPUT_DIRECT_TARGETS, "TES_direct_targets.csv"))
tead1_direct <- read.csv(file.path(INPUT_DIRECT_TARGETS, "TEAD1_direct_targets.csv"))

# Load gene classification
gene_class <- read.csv(file.path(INPUT_REGULATORY, "gene_classification.csv"))

cat(sprintf("TES direct targets: %d\n", nrow(tes_direct)))
cat(sprintf("TEAD1 direct targets: %d\n", nrow(tead1_direct)))
cat(sprintf("Total classified genes: %d\n\n", nrow(gene_class)))

# =============================================================================
# PLOT 1: QUADRANT PLOT - BINDING VS EXPRESSION
# =============================================================================

cat("1. Creating binding vs expression quadrant plot...\n")

# Create combined dataset
# Note: All genes from script 10 are already significant DEGs, so is_significant = TRUE for all
combined_data <- gene_class %>%
  filter(!is.na(log2FoldChange)) %>%
  mutate(
    has_TES_peak = tes_bound,
    has_TEAD1_peak = tead1_bound,
    binding_category = case_when(
      has_TES_peak & has_TEAD1_peak ~ "Both",
      has_TES_peak ~ "TES only",
      has_TEAD1_peak ~ "TEAD1 only",
      TRUE ~ "No binding"
    ),
    is_significant = TRUE, # All genes from script 10 are significant DEGs
    quadrant = case_when(
      has_TES_peak & log2FoldChange > 1 ~ "TES bound + UP",
      has_TES_peak & log2FoldChange < -1 ~ "TES bound + DOWN",
      !has_TES_peak & log2FoldChange > 1 ~ "Unbound + UP",
      !has_TES_peak & log2FoldChange < -1 ~ "Unbound + DOWN",
      TRUE ~ "Other"
    )
  )

# Select genes to label (significant + extreme)
label_genes <- combined_data %>%
  filter(is_significant & (abs(log2FoldChange) > 2 |
    (has_TES_peak & abs(log2FoldChange) > 1.5))) %>%
  arrange(desc(abs(log2FoldChange))) %>%
  head(20)

combined_data$label <- ifelse(combined_data$gene_symbol %in% label_genes$gene_symbol,
  combined_data$gene_symbol, ""
)

p1 <- ggplot(combined_data, aes(x = log2FoldChange, y = -log10(padj))) +
  # Background points
  geom_point(
    data = subset(combined_data, binding_category == "No binding"),
    aes(color = binding_category), alpha = 0.3, size = 1.5
  ) +
  # Binding points
  geom_point(
    data = subset(combined_data, binding_category != "No binding"),
    aes(color = binding_category, size = binding_category),
    alpha = 0.7
  ) +
  # Labels
  geom_text_repel(aes(label = label),
    size = 3, max.overlaps = 20,
    box.padding = 0.5, segment.color = "grey50"
  ) +
  # Thresholds
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  # Scales
  scale_color_manual(
    values = c(
      "Both" = "#E41A1C",
      "TES only" = "#377EB8",
      "TEAD1 only" = "#4DAF4A",
      "No binding" = "grey70"
    ),
    name = "Binding"
  ) +
  scale_size_manual(
    values = c(
      "Both" = 3, "TES only" = 2.5,
      "TEAD1 only" = 2.5, "No binding" = 1
    ),
    guide = "none"
  ) +
  labs(
    title = "Integrative Analysis: TF Binding vs Gene Expression",
    subtitle = "Significant direct targets labeled (FDR < 0.05, |log2FC| > 1)",
    x = "log2 Fold Change (TES vs GFP)",
    y = "-log10(Adjusted p-value)"
  ) +
  theme_pub()

ggsave(file.path(OUTPUT_BASE, "01_quadrant_binding_expression.pdf"), p1,
  width = 14, height = 10, device = cairo_pdf
)
ggsave(file.path(OUTPUT_BASE, "01_quadrant_binding_expression.png"), p1,
  width = 14, height = 10, dpi = 300
)

# =============================================================================
# PLOT 2: OVERLAP BAR CHART - TES vs TEAD1
# =============================================================================

cat("2. Creating overlap bar chart for TES vs TEAD1...\n")

# Create gene sets
tes_genes <- unique(tes_direct$gene_symbol[!is.na(tes_direct$gene_symbol)])
tead1_genes <- unique(tead1_direct$gene_symbol[!is.na(tead1_direct$gene_symbol)])

# Calculate overlaps
overlap_genes <- intersect(tes_genes, tead1_genes)
tes_only <- setdiff(tes_genes, tead1_genes)
tead1_only <- setdiff(tead1_genes, tes_genes)

# Create summary data
overlap_data <- data.frame(
  Category = c("TES only", "Shared", "TEAD1 only"),
  Count = c(length(tes_only), length(overlap_genes), length(tead1_only)),
  Percent = c(
    length(tes_only) / length(tes_genes) * 100,
    length(overlap_genes) / length(tes_genes) * 100,
    length(tead1_only) / length(tead1_genes) * 100
  )
)
overlap_data$Category <- factor(overlap_data$Category, levels = c("TES only", "Shared", "TEAD1 only"))

cat(sprintf("  TES only: %d genes\n", length(tes_only)))
cat(sprintf("  Shared: %d genes\n", length(overlap_genes)))
cat(sprintf("  TEAD1 only: %d genes\n", length(tead1_only)))

# Create bar plot
p2 <- ggplot(overlap_data, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = paste0(Count, "\n(", round(Percent, 1), "%)")),
    vjust = -0.5, size = 5, fontface = "bold"
  ) +
  scale_fill_manual(values = c(
    "TES only" = "#E41A1C",
    "Shared" = "#984EA3",
    "TEAD1 only" = "#377EB8"
  )) +
  labs(
    title = "TES vs TEAD1 Direct Target Overlap",
    subtitle = sprintf(
      "Total: %d TES targets, %d TEAD1 targets",
      length(tes_genes), length(tead1_genes)
    ),
    x = NULL,
    y = "Number of Direct Target Genes"
  ) +
  theme_pub() +
  theme(legend.position = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

ggsave(file.path(OUTPUT_BASE, "02_overlap_TES_TEAD1.pdf"), p2,
  width = 10, height = 8, device = cairo_pdf
)
ggsave(file.path(OUTPUT_BASE, "02_overlap_TES_TEAD1.png"), p2,
  width = 10, height = 8, dpi = 300
)

# =============================================================================
# PLOT 3: HEATMAP - TOP DIRECT TARGETS
# =============================================================================

cat("3. Creating heatmap of top direct targets...\n")

# Get top TES and TEAD1 targets by significance and fold change
top_tes <- tes_direct %>%
  filter(abs(log2FoldChange) > 1) %>%
  arrange(padj, desc(abs(log2FoldChange))) %>%
  head(50)

top_tead1 <- tead1_direct %>%
  filter(abs(log2FoldChange) > 1) %>%
  arrange(padj, desc(abs(log2FoldChange))) %>%
  head(50)

# Combine and create matrix
top_genes <- unique(c(top_tes$gene_symbol, top_tead1$gene_symbol))
top_data <- combined_data %>%
  filter(gene_symbol %in% top_genes) %>%
  select(gene_symbol, log2FoldChange, has_TES_peak, has_TEAD1_peak, padj)

# Create annotation (simplified - no barplot)
ha_row <- rowAnnotation(
  "TES" = ifelse(top_data$has_TES_peak, "Bound", "Unbound"),
  "TEAD1" = ifelse(top_data$has_TEAD1_peak, "Bound", "Unbound"),
  col = list(
    "TES" = c("Bound" = "#E41A1C", "Unbound" = "grey90"),
    "TEAD1" = c("Bound" = "#377EB8", "Unbound" = "grey90")
  ),
  annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
  simple_anno_size = unit(0.4, "cm")
)

# Create main heatmap matrix (just log2FC for now)
mat <- matrix(top_data$log2FoldChange, ncol = 1)
rownames(mat) <- top_data$gene_symbol
colnames(mat) <- "log2FC"

pdf(file.path(OUTPUT_BASE, "03_top_direct_targets_heatmap.pdf"), width = 10, height = 14)
ht <- Heatmap(mat,
  name = "log2FC",
  col = colorRamp2(c(-4, 0, 4), c("#2166AC", "white", "#B2182B")),
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_dend = TRUE,
  right_annotation = ha_row,
  row_names_gp = gpar(fontsize = 7),
  column_names_gp = gpar(fontsize = 10, fontface = "bold"),
  column_title = "Top Direct Targets (TES + TEAD1)",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 9)
  )
)
draw(ht)
dev.off()

png(file.path(OUTPUT_BASE, "03_top_direct_targets_heatmap.png"),
  width = 10, height = 14, units = "in", res = 300
)
draw(ht)
dev.off()

# =============================================================================
# PLOT 4: COMPARATIVE BAR CHART - TARGET CATEGORIES
# =============================================================================

cat("4. Creating comparative target categories chart...\n")

# Count genes in each category
category_counts <- combined_data %>%
  mutate(
    category = case_when(
      has_TES_peak & has_TEAD1_peak & is_significant ~ "Shared targets",
      has_TES_peak & !has_TEAD1_peak & is_significant ~ "TES-specific",
      !has_TES_peak & has_TEAD1_peak & is_significant ~ "TEAD1-specific",
      (has_TES_peak | has_TEAD1_peak) & !is_significant ~ "Bound but not DE",
      !has_TES_peak & !has_TEAD1_peak & is_significant ~ "Indirect targets",
      TRUE ~ "Other"
    ),
    direction = ifelse(log2FoldChange > 0, "UP", "DOWN")
  ) %>%
  filter(category != "Other") %>%
  group_by(category, direction) %>%
  summarise(count = n(), .groups = "drop")

# ORIGINAL PLOT (raw counts - stacked bar)
p4 <- ggplot(category_counts, aes(x = reorder(category, count), y = count, fill = direction)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.5, alpha = 0.8) +
  geom_text(aes(label = count),
    position = position_stack(vjust = 0.5),
    fontface = "bold", size = 4
  ) +
  scale_fill_manual(
    values = c("UP" = "#E74C3C", "DOWN" = "#3498DB"),
    name = "Direction"
  ) +
  coord_flip() +
  labs(
    title = "Gene Regulatory Categories",
    subtitle = "Classification based on binding and expression",
    x = NULL,
    y = "Number of Genes"
  ) +
  theme_pub()

ggsave(file.path(OUTPUT_BASE, "04_target_categories.pdf"), p4,
  width = 12, height = 7, device = cairo_pdf
)
ggsave(file.path(OUTPUT_BASE, "04_target_categories.png"), p4,
  width = 12, height = 7, dpi = 300
)

# =============================================================================
# PLOT 4b: DIRECTIONAL BIAS CHART (PERCENTAGE-BASED)
# =============================================================================

cat("4b. Creating directional bias percentage chart...\n")

# Calculate percentages within each category
category_percentages <- category_counts %>%
  group_by(category) %>%
  mutate(
    total = sum(count),
    percentage = count / total * 100
  ) %>%
  ungroup()

# Print summary to log for verification
cat("\n=== DIRECTIONAL BIAS SUMMARY ===\n")
cat("Category-wise UP/DOWN percentages:\n")
category_percentages %>%
  arrange(category, direction) %>%
  mutate(summary = sprintf("  %s (%s): %d genes (%.1f%%)", category, direction, count, percentage)) %>%
  pull(summary) %>%
  cat(sep = "\n")

# Calculate UP percentage for each category for the bar chart
up_percentages <- category_percentages %>%
  filter(direction == "UP") %>%
  select(category, up_pct = percentage, up_count = count, total)

down_percentages <- category_percentages %>%
  filter(direction == "DOWN") %>%
  select(category, down_pct = percentage, down_count = count)

direction_summary <- up_percentages %>%
  left_join(down_percentages, by = "category") %>%
  mutate(
    bias = case_when(
      up_pct > 55 ~ "UP bias",
      up_pct < 45 ~ "DOWN bias",
      TRUE ~ "Balanced"
    ),
    label = sprintf("%.1f%% UP\n(n=%d)", up_pct, total)
  )

# Horizontal bar chart showing UP percentage with 50% reference line
p4b <- ggplot(direction_summary, aes(x = reorder(category, up_pct), y = up_pct, fill = bias)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.5, width = 0.7) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "black", linewidth = 1) +
  geom_text(aes(label = sprintf("%.1f%%", up_pct)),
            hjust = -0.2, fontface = "bold", size = 4) +
  geom_text(aes(y = 2, label = sprintf("n=%d", total)),
            hjust = 0, fontface = "plain", size = 3, color = "white") +
  scale_fill_manual(
    values = c("UP bias" = "#E74C3C", "DOWN bias" = "#3498DB", "Balanced" = "#95A5A6"),
    name = "Directional Bias"
  ) +
  coord_flip() +
  ylim(0, 70) +
  labs(
    title = "Directional Bias by Binding Category",
    subtitle = "Percentage of UP-regulated genes (dashed line = 50%)",
    x = NULL,
    y = "% UP-regulated"
  ) +
  theme_pub() +
  annotate("text", x = 0.5, y = 52, label = "50%", hjust = 0, fontface = "italic", size = 3)

ggsave(file.path(OUTPUT_BASE, "04b_directional_bias.pdf"), p4b,
  width = 12, height = 7, device = cairo_pdf
)
ggsave(file.path(OUTPUT_BASE, "04b_directional_bias.png"), p4b,
  width = 12, height = 7, dpi = 300
)

# =============================================================================
# PLOT 4c: DIVERGING BAR CHART (centered at 50%)
# =============================================================================

cat("4c. Creating diverging bar chart centered at 50%...\n")

# Create data for diverging bars (deviation from 50%)
diverging_data <- direction_summary %>%
  mutate(
    deviation = up_pct - 50,  # Positive = more UP, Negative = more DOWN
    category_label = sprintf("%s\n(n=%d)", category, total)
  )

p4c <- ggplot(diverging_data, aes(x = reorder(category, deviation), y = deviation, fill = deviation > 0)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.5, width = 0.7) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%+.1f%%", deviation)),
            hjust = ifelse(diverging_data$deviation > 0, -0.2, 1.2),
            fontface = "bold", size = 4) +
  scale_fill_manual(
    values = c("TRUE" = "#E74C3C", "FALSE" = "#3498DB"),
    labels = c("TRUE" = "More UP", "FALSE" = "More DOWN"),
    name = "Direction"
  ) +
  coord_flip() +
  ylim(-15, 15) +
  labs(
    title = "Deviation from Expected 50/50 Distribution",
    subtitle = "Positive = more upregulated, Negative = more downregulated",
    x = NULL,
    y = "Deviation from 50% (percentage points)"
  ) +
  theme_pub()

ggsave(file.path(OUTPUT_BASE, "04c_diverging_direction.pdf"), p4c,
  width = 12, height = 7, device = cairo_pdf
)
ggsave(file.path(OUTPUT_BASE, "04c_diverging_direction.png"), p4c,
  width = 12, height = 7, dpi = 300
)

# =============================================================================
# PLOT 4d: STACKED PROPORTION BAR (100% bars)
# =============================================================================

cat("4d. Creating 100% stacked proportion bar chart...\n")

# Add total to category_counts for ordering
category_counts_with_total <- category_counts %>%
  group_by(category) %>%
  mutate(total = sum(count)) %>%
  ungroup()

p4d <- ggplot(category_counts_with_total,
              aes(x = reorder(category, -total), y = count, fill = direction)) +
  geom_bar(stat = "identity", position = "fill", color = "black", linewidth = 0.5) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "white", linewidth = 1) +
  geom_text(aes(label = sprintf("%d\n(%.0f%%)", count, count/total*100)),
            position = position_fill(vjust = 0.5),
            fontface = "bold", size = 3.5, color = "white") +
  scale_fill_manual(
    values = c("UP" = "#E74C3C", "DOWN" = "#3498DB"),
    name = "Direction"
  ) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  coord_flip() +
  labs(
    title = "Gene Direction by Binding Category (Proportions)",
    subtitle = "Dashed line = 50% reference",
    x = NULL,
    y = "Proportion of Genes"
  ) +
  theme_pub()

ggsave(file.path(OUTPUT_BASE, "04d_stacked_proportions.pdf"), p4d,
  width = 12, height = 7, device = cairo_pdf
)
ggsave(file.path(OUTPUT_BASE, "04d_stacked_proportions.png"), p4d,
  width = 12, height = 7, dpi = 300
)

# =============================================================================
# PLOT 5: MULTI-PANEL SUMMARY
# =============================================================================

cat("5. Creating multi-panel integrative summary...\n")

# Panel A: Summary statistics
summary_stats <- data.frame(
  Category = c("Total DEGs", "TES direct", "TEAD1 direct", "Shared", "Indirect"),
  Count = c(
    sum(combined_data$is_significant),
    nrow(tes_direct),
    nrow(tead1_direct),
    sum(combined_data$has_TES_peak & combined_data$has_TEAD1_peak & combined_data$is_significant),
    sum(!combined_data$has_TES_peak & !combined_data$has_TEAD1_peak & combined_data$is_significant)
  )
)

pA <- ggplot(summary_stats, aes(x = reorder(Category, Count), y = Count, fill = Category)) +
  geom_bar(stat = "identity", color = "black", alpha = 0.8) +
  geom_text(aes(label = Count), hjust = -0.2, fontface = "bold", size = 4) +
  scale_fill_brewer(palette = "Set2") +
  coord_flip() +
  labs(title = "A) Overview", x = NULL, y = "Gene Count") +
  theme_pub() +
  theme(legend.position = "none") +
  ylim(0, max(summary_stats$Count) * 1.15)

# Panel B: Binding enrichment
binding_enrich <- combined_data %>%
  filter(is_significant) %>%
  mutate(binding = ifelse(has_TES_peak | has_TEAD1_peak, "Bound", "Unbound")) %>%
  group_by(binding, log2FoldChange > 0) %>%
  summarise(count = n(), .groups = "drop") %>%
  rename(upregulated = `log2FoldChange > 0`)

pB <- ggplot(binding_enrich, aes(x = binding, y = count, fill = upregulated)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", alpha = 0.8) +
  scale_fill_manual(
    values = c("TRUE" = "#E74C3C", "FALSE" = "#3498DB"),
    labels = c("Down", "Up"),
    name = "Direction"
  ) +
  labs(title = "B) Binding Enrichment", x = NULL, y = "Count") +
  theme_pub()

# Panel C: TES vs TEAD1 overlap (reuse from plot 2)
pC <- ggplot(overlap_data, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 0.7, color = "black", alpha = 0.8) +
  geom_text(aes(label = Count), vjust = -0.5, fontface = "bold", size = 4) +
  scale_fill_manual(values = c(
    "TES only" = "#E41A1C",
    "Shared" = "#984EA3",
    "TEAD1 only" = "#377EB8"
  )) +
  labs(
    title = "C) TES vs TEAD1 Overlap",
    x = NULL,
    y = "Direct Targets"
  ) +
  theme_pub() +
  theme(legend.position = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

# Panel D: FC distribution by category
fc_dist <- combined_data %>%
  filter(is_significant) %>%
  mutate(category = case_when(
    has_TES_peak & has_TEAD1_peak ~ "Both",
    has_TES_peak ~ "TES only",
    has_TEAD1_peak ~ "TEAD1 only",
    TRUE ~ "Indirect"
  ))

pD <- ggplot(fc_dist, aes(x = log2FoldChange, fill = category)) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(
    values = c(
      "Both" = "#E41A1C", "TES only" = "#377EB8",
      "TEAD1 only" = "#4DAF4A", "Indirect" = "#984EA3"
    ),
    name = "Category"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(
    title = "D) Expression Change Distribution",
    x = "log2 Fold Change",
    y = "Density"
  ) +
  theme_pub()

# Combine
combined_summary <- (pA | pB) / (pC | pD) +
  plot_annotation(
    title = "Integrative Analysis Summary: Cut&Tag + RNA-seq",
    subtitle = "TES/TEAD1 binding and transcriptional regulation",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5)
    )
  )

ggsave(file.path(OUTPUT_BASE, "05_multipanel_integrative_summary.pdf"), combined_summary,
  width = 16, height = 12, device = cairo_pdf
)
ggsave(file.path(OUTPUT_BASE, "05_multipanel_integrative_summary.png"), combined_summary,
  width = 16, height = 12, dpi = 300
)

# =============================================================================
# COMPLETION
# =============================================================================

cat("\n========================================\n")
cat("ENHANCED INTEGRATIVE PLOTS COMPLETE\n")
cat("========================================\n")
cat(sprintf("Input directory: %s\n", INPUT_BASE))
cat(sprintf("Output directory: %s\n", OUTPUT_BASE))
cat("\nGenerated files:\n")
cat("  1. 01_quadrant_binding_expression.pdf/png\n")
cat("  2. 02_overlap_TES_TEAD1.pdf/png\n")
cat("  3. 03_top_direct_targets_heatmap.pdf/png\n")
cat("  4. 04_target_categories.pdf/png (raw counts)\n")
cat("  4b. 04b_directional_bias.pdf/png (% UP with 50% reference)\n")
cat("  4c. 04c_diverging_direction.pdf/png (deviation from 50%)\n")
cat("  4d. 04d_stacked_proportions.pdf/png (100% stacked bars)\n")
cat("  5. 05_multipanel_integrative_summary.pdf/png\n")
cat("\n")
