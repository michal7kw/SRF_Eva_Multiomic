#!/usr/bin/env Rscript
#
# GSEA ANALYSIS: Cancer-Relevant Pathways
# Focused analysis of cell death, apoptosis, migration, and proliferation
# Comparing TES-only, TEAD1-only, and shared transcriptional targets

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(pheatmap)
  library(ComplexHeatmap)
  library(circlize)
  library(gridExtra)
  library(RColorBrewer)
})

# Set working directory and create output directories
setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_1")

cat("=== GSEA ANALYSIS: Cancer-Relevant Pathways ===\n")
cat("Analysis started:", as.character(Sys.time()), "\n\n")

# Create output directory
output_dir <- "output/17_gsea_cancer_pathways"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# PHASE 1: DATA LOADING
# =============================================================================

cat("=== PHASE 1: Loading Integrative Analysis Results ===\n")

# Load direct target classifications
cat("Loading direct target gene lists...\n")
tes_direct <- read.csv("output/11_final_integrative_analysis_all_genes/direct_targets/TES_direct_targets_all_genes.csv",
  stringsAsFactors = FALSE
)
tead1_direct <- read.csv("output/11_final_integrative_analysis_all_genes/direct_targets/TEAD1_direct_targets_all_genes.csv",
  stringsAsFactors = FALSE
)
tes_specific <- read.csv("output/11_final_integrative_analysis_all_genes/direct_targets/TES_specific_targets_all_genes.csv",
  stringsAsFactors = FALSE
)
tead1_specific <- read.csv("output/11_final_integrative_analysis_all_genes/direct_targets/TEAD1_specific_targets_all_genes.csv",
  stringsAsFactors = FALSE
)

# Calculate shared targets
shared_genes <- intersect(tes_direct$ensembl_id, tead1_direct$ensembl_id)
shared_targets <- tes_direct[tes_direct$ensembl_id %in% shared_genes, ]

cat(sprintf("✓ TES-only targets: %d genes\n", nrow(tes_specific)))
cat(sprintf("✓ TEAD1-only targets: %d genes\n", nrow(tead1_specific)))
cat(sprintf("✓ Shared TES/TEAD1 targets: %d genes\n\n", nrow(shared_targets)))

# Load GO enrichment results
cat("Loading GO enrichment pathway data...\n")
tes_go <- read.csv("output/11_final_integrative_analysis_all_genes/pathway_analysis/TES_direct_GO_enrichment_all_genes.csv",
  stringsAsFactors = FALSE
)
tead1_go <- read.csv("output/11_final_integrative_analysis_all_genes/pathway_analysis/TEAD1_direct_GO_enrichment_all_genes.csv",
  stringsAsFactors = FALSE
)
shared_go <- read.csv("output/11_final_integrative_analysis_all_genes/pathway_analysis/TES_TEAD1_shared_GO_enrichment_all_genes.csv",
  stringsAsFactors = FALSE
)

cat(sprintf("✓ TES direct pathways: %d\n", nrow(tes_go)))
cat(sprintf("✓ TEAD1 direct pathways: %d\n", nrow(tead1_go)))
cat(sprintf("✓ Shared pathways: %d\n\n", nrow(shared_go)))

# =============================================================================
# PHASE 2: FILTER FOR CANCER-RELEVANT PATHWAYS
# =============================================================================

cat("=== PHASE 2: Filtering for Cancer-Relevant Pathways ===\n")

# Define pathway keyword groups
cancer_keywords <- list(
  apoptosis = c(
    "apoptosis", "apoptotic", "cell death", "programmed cell death",
    "necrosis", "necrotic", "ferroptosis", "pyroptosis", "anoikis"
  ),
  migration = c(
    "migration", "migratory", "motility", "invasion", "invasive",
    "chemotaxis", "chemotactic", "cell movement", "locomotion"
  ),
  proliferation = c(
    "proliferation", "proliferative", "cell cycle", "mitosis",
    "mitotic", "cell division", "growth", "G1/S", "G2/M"
  )
)

# Function to filter pathways by keywords
filter_cancer_pathways <- function(go_results, keywords_list) {
  all_keywords <- unlist(keywords_list, use.names = FALSE)
  pattern <- paste(all_keywords, collapse = "|")

  filtered <- go_results[grep(pattern, go_results$Description, ignore.case = TRUE), ]

  # Add pathway category
  filtered$pathway_category <- NA
  for (category in names(keywords_list)) {
    category_pattern <- paste(keywords_list[[category]], collapse = "|")
    filtered$pathway_category[grep(category_pattern, filtered$Description, ignore.case = TRUE)] <- category
  }

  return(filtered)
}

# Filter pathways for each target group
cat("Filtering pathways by cancer-relevant keywords...\n")
tes_cancer <- filter_cancer_pathways(tes_go, cancer_keywords)
tead1_cancer <- filter_cancer_pathways(tead1_go, cancer_keywords)
shared_cancer <- filter_cancer_pathways(shared_go, cancer_keywords)

cat(sprintf(
  "✓ TES cancer pathways: %d (%d apoptosis, %d migration, %d proliferation)\n",
  nrow(tes_cancer),
  sum(tes_cancer$pathway_category == "apoptosis", na.rm = TRUE),
  sum(tes_cancer$pathway_category == "migration", na.rm = TRUE),
  sum(tes_cancer$pathway_category == "proliferation", na.rm = TRUE)
))
cat(sprintf(
  "✓ TEAD1 cancer pathways: %d (%d apoptosis, %d migration, %d proliferation)\n",
  nrow(tead1_cancer),
  sum(tead1_cancer$pathway_category == "apoptosis", na.rm = TRUE),
  sum(tead1_cancer$pathway_category == "migration", na.rm = TRUE),
  sum(tead1_cancer$pathway_category == "proliferation", na.rm = TRUE)
))
cat(sprintf(
  "✓ Shared cancer pathways: %d (%d apoptosis, %d migration, %d proliferation)\n\n",
  nrow(shared_cancer),
  sum(shared_cancer$pathway_category == "apoptosis", na.rm = TRUE),
  sum(shared_cancer$pathway_category == "migration", na.rm = TRUE),
  sum(shared_cancer$pathway_category == "proliferation", na.rm = TRUE)
))

# =============================================================================
# PHASE 3: PREPARE DATA FOR VISUALIZATION
# =============================================================================

cat("=== PHASE 3: Preparing Comparative Datasets ===\n")

# Add source annotation
tes_cancer$target_group <- "TES-only"
tead1_cancer$target_group <- "TEAD1-only"
shared_cancer$target_group <- "Shared"

# Combine all cancer pathways
all_cancer_pathways <- rbind(
  tes_cancer[, c(
    "ID", "Description", "GeneRatio", "pvalue", "p.adjust",
    "qvalue", "Count", "pathway_category", "target_group"
  )],
  tead1_cancer[, c(
    "ID", "Description", "GeneRatio", "pvalue", "p.adjust",
    "qvalue", "Count", "pathway_category", "target_group"
  )],
  shared_cancer[, c(
    "ID", "Description", "GeneRatio", "pvalue", "p.adjust",
    "qvalue", "Count", "pathway_category", "target_group"
  )]
)

# Parse GeneRatio for plotting
all_cancer_pathways$GeneRatio_numeric <- sapply(
  strsplit(all_cancer_pathways$GeneRatio, "/"),
  function(x) as.numeric(x[1]) / as.numeric(x[2])
)

cat(sprintf("✓ Total cancer-relevant pathway entries: %d\n", nrow(all_cancer_pathways)))
cat(sprintf(
  "✓ Unique pathways across all groups: %d\n\n",
  length(unique(all_cancer_pathways$ID))
))

# =============================================================================
# PHASE 4: EXPORT FILTERED RESULTS
# =============================================================================

cat("=== PHASE 4: Exporting Filtered Pathway Tables ===\n")

write.csv(tes_cancer,
  file.path(output_dir, "TES_only_cancer_pathways.csv"),
  row.names = FALSE
)
write.csv(tead1_cancer,
  file.path(output_dir, "TEAD1_only_cancer_pathways.csv"),
  row.names = FALSE
)
write.csv(shared_cancer,
  file.path(output_dir, "Shared_cancer_pathways.csv"),
  row.names = FALSE
)
write.csv(all_cancer_pathways,
  file.path(output_dir, "All_cancer_pathways_combined.csv"),
  row.names = FALSE
)

cat("✓ Pathway tables exported\n\n")

# =============================================================================
# PHASE 5: VISUALIZATION - PATHWAY COMPARISON
# =============================================================================

cat("=== PHASE 5: Creating Visualizations ===\n")

# Plot 1: Pathway category counts by target group
cat("Creating pathway category count barplot...\n")
category_counts <- all_cancer_pathways %>%
  filter(!is.na(pathway_category)) %>%
  group_by(target_group, pathway_category) %>%
  summarise(count = n(), .groups = "drop")

pdf(file.path(output_dir, "01_pathway_category_counts.pdf"), width = 12, height = 8)
p1 <- ggplot(category_counts, aes(x = pathway_category, y = count, fill = target_group)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", linewidth = 0.3) +
  labs(
    title = "Cancer-Relevant Pathway Enrichment by Target Group",
    subtitle = "Cell Death/Apoptosis, Migration, and Proliferation Pathways",
    x = "Pathway Category",
    y = "Number of Enriched Pathways",
    fill = "Target Group"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),
    axis.text.y = element_text(size = 12),
    legend.position = "right"
  ) +
  scale_fill_brewer(palette = "Set2")
print(p1)
dev.off()

# Plot 2: Top pathways dotplot (top 10 per category)
cat("Creating top pathways dotplot...\n")
top_pathways <- all_cancer_pathways %>%
  filter(!is.na(pathway_category)) %>%
  group_by(target_group, pathway_category) %>%
  arrange(p.adjust) %>%
  slice_head(n = 10) %>%
  ungroup()

pdf(file.path(output_dir, "02_top_pathways_dotplot.pdf"), width = 16, height = 12)
p2 <- ggplot(top_pathways, aes(
  x = target_group, y = Description,
  color = p.adjust, size = GeneRatio_numeric
)) +
  geom_point() +
  facet_wrap(~pathway_category, scales = "free_y", ncol = 1) +
  scale_color_gradient(low = "red", high = "blue", name = "Adjusted p-value") +
  scale_size_continuous(range = c(3, 10), name = "Gene Ratio") +
  labs(
    title = "Top 10 Cancer Pathways per Category and Target Group",
    x = "Target Group",
    y = "Pathway Description"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 9),
    strip.text = element_text(face = "bold", size = 12),
    strip.background = element_rect(fill = "lightgray")
  )
print(p2)
dev.off()

# Plot 3: Heatmap of pathway presence across target groups
cat("Creating pathway presence heatmap...\n")

# Get top pathways (by significance) from each category
top_pathways_for_heatmap <- all_cancer_pathways %>%
  filter(!is.na(pathway_category)) %>%
  group_by(pathway_category) %>%
  arrange(p.adjust) %>%
  slice_head(n = 20) %>%
  ungroup()

# Create presence/absence matrix
pathway_matrix <- top_pathways_for_heatmap %>%
  select(Description, target_group, p.adjust) %>%
  tidyr::pivot_wider(
    names_from = target_group, values_from = p.adjust,
    values_fill = 1
  ) %>%
  as.data.frame()

rownames(pathway_matrix) <- pathway_matrix$Description
pathway_matrix <- pathway_matrix[, -1]

# Convert p-values to -log10 for better visualization
pathway_matrix_log <- as.matrix(-log10(pathway_matrix + 1e-100))

# Create annotation for pathway categories
pathway_annotation <- top_pathways_for_heatmap %>%
  select(Description, pathway_category) %>%
  distinct() %>%
  as.data.frame()
rownames(pathway_annotation) <- pathway_annotation$Description
pathway_annotation <- pathway_annotation[rownames(pathway_matrix_log), , drop = FALSE]

pdf(file.path(output_dir, "03_pathway_heatmap.pdf"), width = 10, height = 14)
pheatmap(pathway_matrix_log,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("white", "yellow", "orange", "red", "darkred"))(100),
  border_color = "gray80",
  fontsize_row = 8,
  fontsize_col = 12,
  annotation_row = pathway_annotation[, "pathway_category", drop = FALSE],
  annotation_colors = list(pathway_category = c(
    apoptosis = "#E41A1C",
    migration = "#377EB8",
    proliferation = "#4DAF4A"
  )),
  main = "Cancer Pathway Enrichment Heatmap\n(-log10 adjusted p-value)",
  show_rownames = TRUE,
  show_colnames = TRUE
)
dev.off()

# Plot 4: Volcano-style plots showing enrichment strength
cat("Creating enrichment strength comparison...\n")

pdf(file.path(output_dir, "04_enrichment_strength_by_category.pdf"), width = 14, height = 10)
p4 <- ggplot(
  all_cancer_pathways %>% filter(!is.na(pathway_category)),
  aes(x = GeneRatio_numeric, y = -log10(p.adjust), color = target_group)
) +
  geom_point(size = 3, alpha = 0.7) +
  facet_wrap(~pathway_category, scales = "free") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40") +
  labs(
    title = "Pathway Enrichment Strength by Category",
    subtitle = "Gene Ratio vs Statistical Significance",
    x = "Gene Ratio (proportion of genes in pathway)",
    y = "-log10(Adjusted p-value)",
    color = "Target Group"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    strip.text = element_text(face = "bold", size = 12),
    strip.background = element_rect(fill = "lightgray")
  ) +
  scale_color_brewer(palette = "Set1")
print(p4)
dev.off()

# Plot 5: Gene count distribution
cat("Creating gene count distribution plot...\n")

pdf(file.path(output_dir, "05_gene_count_distribution.pdf"), width = 12, height = 8)
p5 <- ggplot(
  all_cancer_pathways %>% filter(!is.na(pathway_category)),
  aes(x = Count, fill = target_group)
) +
  geom_histogram(bins = 30, alpha = 0.7, position = "identity", color = "black", linewidth = 0.3) +
  facet_wrap(~pathway_category, scales = "free_y") +
  labs(
    title = "Distribution of Gene Counts in Cancer Pathways",
    x = "Number of Genes in Pathway",
    y = "Frequency",
    fill = "Target Group"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    strip.text = element_text(face = "bold", size = 12),
    strip.background = element_rect(fill = "lightgray")
  ) +
  scale_fill_brewer(palette = "Set2")
print(p5)
dev.off()

cat("✓ All visualizations created\n\n")

# =============================================================================
# PHASE 6: SUMMARY STATISTICS
# =============================================================================

cat("=== PHASE 6: Generating Summary Report ===\n")

# Create summary statistics
summary_stats <- data.frame(
  Metric = c(
    "TES-only targets (genes)",
    "TEAD1-only targets (genes)",
    "Shared TES/TEAD1 targets (genes)",
    "",
    "TES-only cancer pathways",
    "  - Apoptosis pathways",
    "  - Migration pathways",
    "  - Proliferation pathways",
    "",
    "TEAD1-only cancer pathways",
    "  - Apoptosis pathways",
    "  - Migration pathways",
    "  - Proliferation pathways",
    "",
    "Shared cancer pathways",
    "  - Apoptosis pathways",
    "  - Migration pathways",
    "  - Proliferation pathways"
  ),
  Value = c(
    nrow(tes_specific),
    nrow(tead1_specific),
    nrow(shared_targets),
    "",
    nrow(tes_cancer),
    sum(tes_cancer$pathway_category == "apoptosis", na.rm = TRUE),
    sum(tes_cancer$pathway_category == "migration", na.rm = TRUE),
    sum(tes_cancer$pathway_category == "proliferation", na.rm = TRUE),
    "",
    nrow(tead1_cancer),
    sum(tead1_cancer$pathway_category == "apoptosis", na.rm = TRUE),
    sum(tead1_cancer$pathway_category == "migration", na.rm = TRUE),
    sum(tead1_cancer$pathway_category == "proliferation", na.rm = TRUE),
    "",
    nrow(shared_cancer),
    sum(shared_cancer$pathway_category == "apoptosis", na.rm = TRUE),
    sum(shared_cancer$pathway_category == "migration", na.rm = TRUE),
    sum(shared_cancer$pathway_category == "proliferation", na.rm = TRUE)
  )
)

# Write summary to file
summary_file <- file.path(output_dir, "GSEA_CANCER_PATHWAYS_SUMMARY.txt")
cat("GSEA CANCER PATHWAYS ANALYSIS SUMMARY\n", file = summary_file)
cat("=====================================\n\n", file = summary_file, append = TRUE)
cat(paste("Generated:", Sys.time(), "\n\n"), file = summary_file, append = TRUE)
cat("TARGET GROUPS AND PATHWAY ENRICHMENT\n", file = summary_file, append = TRUE)
cat("------------------------------------\n\n", file = summary_file, append = TRUE)

write.table(summary_stats,
  file = summary_file, append = TRUE,
  row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t"
)

cat("\n\nKEY FINDINGS:\n", file = summary_file, append = TRUE)
cat("-------------\n\n", file = summary_file, append = TRUE)

# Calculate percentages
tes_apop_pct <- round(100 * sum(tes_cancer$pathway_category == "apoptosis", na.rm = TRUE) / nrow(tes_cancer), 1)
tes_migr_pct <- round(100 * sum(tes_cancer$pathway_category == "migration", na.rm = TRUE) / nrow(tes_cancer), 1)
tes_prol_pct <- round(100 * sum(tes_cancer$pathway_category == "proliferation", na.rm = TRUE) / nrow(tes_cancer), 1)

tead1_apop_pct <- round(100 * sum(tead1_cancer$pathway_category == "apoptosis", na.rm = TRUE) / nrow(tead1_cancer), 1)
tead1_migr_pct <- round(100 * sum(tead1_cancer$pathway_category == "migration", na.rm = TRUE) / nrow(tead1_cancer), 1)
tead1_prol_pct <- round(100 * sum(tead1_cancer$pathway_category == "proliferation", na.rm = TRUE) / nrow(tead1_cancer), 1)

cat(sprintf("TES-only pathway distribution:\n"), file = summary_file, append = TRUE)
cat(
  sprintf(
    "  Apoptosis: %.1f%%, Migration: %.1f%%, Proliferation: %.1f%%\n\n",
    tes_apop_pct, tes_migr_pct, tes_prol_pct
  ),
  file = summary_file, append = TRUE
)

cat(sprintf("TEAD1-only pathway distribution:\n"), file = summary_file, append = TRUE)
cat(
  sprintf(
    "  Apoptosis: %.1f%%, Migration: %.1f%%, Proliferation: %.1f%%\n\n",
    tead1_apop_pct, tead1_migr_pct, tead1_prol_pct
  ),
  file = summary_file, append = TRUE
)

cat("\n\nOUTPUT FILES:\n", file = summary_file, append = TRUE)
cat("-------------\n", file = summary_file, append = TRUE)
cat("CSV Tables:\n", file = summary_file, append = TRUE)
cat("  - TES_only_cancer_pathways.csv\n", file = summary_file, append = TRUE)
cat("  - TEAD1_only_cancer_pathways.csv\n", file = summary_file, append = TRUE)
cat("  - Shared_cancer_pathways.csv\n", file = summary_file, append = TRUE)
cat("  - All_cancer_pathways_combined.csv\n\n", file = summary_file, append = TRUE)
cat("Visualizations:\n", file = summary_file, append = TRUE)
cat("  - 01_pathway_category_counts.pdf\n", file = summary_file, append = TRUE)
cat("  - 02_top_pathways_dotplot.pdf\n", file = summary_file, append = TRUE)
cat("  - 03_pathway_heatmap.pdf\n", file = summary_file, append = TRUE)
cat("  - 04_enrichment_strength_by_category.pdf\n", file = summary_file, append = TRUE)
cat("  - 05_gene_count_distribution.pdf\n", file = summary_file, append = TRUE)

cat("✓ Summary report saved\n\n")

# =============================================================================
# COMPLETION
# =============================================================================

cat("========================================\n")
cat("GSEA CANCER PATHWAYS ANALYSIS COMPLETE\n")
cat("========================================\n")
cat("Completed:", as.character(Sys.time()), "\n")
cat(sprintf("Output directory: %s\n", output_dir))
cat("\nAll results exported successfully!\n")
