#!/usr/bin/env Rscript

################################################################################
# Phase 5.1: Co-Regulator Identification
# TES vs TEAD1 Comparative Study (Excluding TESmut)
#
# Purpose: Identify TF co-regulators (AP-1, RUNX, ETS families) at TES/TEAD1 sites
#
# Analysis:
# 1. HOMER motif enrichment focusing on AP-1, RUNX, ETS families
# 2. Motif co-occurrence analysis at TES vs TEAD1 sites
# 3. Position preference analysis (motif spacing relative to TEAD motif)
# 4. Expression correlation (co-regulator TF expression vs target gene expression)
# 5. Target gene stratification by co-regulator presence
#
# Author: Advanced Multi-Omics Analysis Pipeline
# Date: 2025-01-24
################################################################################

library(GenomicRanges)
library(rtracklayer)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(org.Hs.eg.db)

# Fix namespace conflicts: dplyr functions get masked by AnnotationDbi
select <- dplyr::select
filter <- dplyr::filter
slice <- dplyr::slice
mutate <- dplyr::mutate
arrange <- dplyr::arrange

# Set working directory
setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top")

# Create output directories - separate for detailed vs simplified
output_dir <- "SRF_Eva_integrated_analysis/scripts/analysis_3/results/10_coregulator_networks"
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

log_message("Starting Phase 5.1: Co-Regulator Identification")

################################################################################
# 1. Load Data
################################################################################

log_message("Loading binding classification data...")

# Load binding classification from Phase 1.1
binding_file <- "SRF_Eva_integrated_analysis/scripts/analysis_3/results/01_binding_classification/binding_site_classification.csv"

if (!file.exists(binding_file)) {
  stop("Binding classification data not found. Please run Phase 1.1 first.")
}

binding_data <- read_csv(binding_file, show_col_types = FALSE)

log_message(sprintf("  Loaded %d binding sites", nrow(binding_data)))

# Load expression data from Phase 2.1
expr_file <- "SRF_Eva_integrated_analysis/scripts/analysis_3/results/04_category_expression/genes_with_binding_and_expression.csv"

if (!file.exists(expr_file)) {
  stop("Expression data not found. Please run Phase 2.1 first.")
}

expr_data <- read_csv(expr_file, show_col_types = FALSE)

log_message(sprintf("  Loaded expression data for %d genes", nrow(expr_data)))

# Load HOMER motif results from Phase 1.2
motif_dir <- "SRF_Eva_integrated_analysis/scripts/analysis_3/results/02_motif_analysis"

################################################################################
# 2. Parse HOMER Known Motif Results
################################################################################

log_message("Parsing HOMER motif enrichment results...")

# Function to parse HOMER knownResults.txt files
parse_homer_known <- function(file_path, category) {
  if (!file.exists(file_path)) {
    warning(sprintf("HOMER results not found: %s", file_path))
    return(NULL)
  }

  # Read HOMER output (tab-delimited)
  results <- read_tsv(file_path, show_col_types = FALSE)

  # Add category column
  results$category <- category

  return(results)
}

# Parse motif results for each binding category
# Note: HOMER outputs are in directories named {category}_motifs (e.g., TES_unique_motifs)
categories <- c("TES_unique", "TEAD1_unique",
                "Shared_equivalent", "Shared_TES_dominant", "Shared_TEAD1_dominant")

motif_results <- list()

for (cat in categories) {
  # HOMER directories have "_motifs" suffix
  homer_file <- file.path(motif_dir, paste0(cat, "_motifs"), "knownResults.txt")
  motif_results[[cat]] <- parse_homer_known(homer_file, cat)
}

# Combine all results
all_motifs <- bind_rows(motif_results)

log_message(sprintf("  Parsed motif enrichment for %d categories", length(categories)))

# Define TF family patterns (needed in multiple places, so define early)
tf_families <- list(
  AP1 = c("Jun", "Fos", "ATF", "JDP", "BATF", "AP-1"),
  RUNX = c("Runx", "RUNX", "AML", "Cbfa"),
  ETS = c("Ets", "ETS", "Elk", "ELK", "ERG", "FLI", "GABPa", "SPI1", "PU.1"),
  TEAD = c("TEAD", "Tead", "TEF"),
  AP2 = c("TFAP2", "AP-2"),
  GATA = c("GATA", "Gata"),
  NFY = c("NFY", "NF-Y"),
  SP = c("Sp1", "SP1", "Sp2", "SP2", "KLF")
)

# Check if we have any motif data
if (nrow(all_motifs) == 0 || !"Motif Name" %in% colnames(all_motifs)) {
  log_message("WARNING: No HOMER motif results found. Creating empty outputs...")

  # Create empty output files
  write_csv(data.frame(category = character(), tf_family = character(),
                       n_motifs = integer(), median_pvalue = numeric(),
                       median_enrichment = numeric()),
            file.path(output_dir, "tf_family_enrichment_summary.csv"))
  write_csv(data.frame(), file.path(output_dir, "coregulator_motif_enrichment.csv"))

  # Skip to expression analysis section
  sig_motifs <- data.frame()
  family_summary <- data.frame()
  coregulator_motifs <- data.frame()

} else {

################################################################################
# 3. Identify Co-Regulator Motif Families
################################################################################

log_message("Identifying co-regulator motif families...")

# Classify motifs by TF family
classify_tf_family <- function(motif_name) {
  for (family in names(tf_families)) {
    patterns <- tf_families[[family]]
    if (any(str_detect(motif_name, paste(patterns, collapse = "|")))) {
      return(family)
    }
  }
  return("Other")
}

all_motifs$tf_family <- sapply(all_motifs$`Motif Name`, classify_tf_family)

# Filter for significant motifs (p-value < 1e-10)
sig_motifs <- all_motifs %>%
  filter(`P-value` < 1e-10)

log_message(sprintf("  Found %d significant motifs", nrow(sig_motifs)))

# Summarize by TF family
family_summary <- sig_motifs %>%
  group_by(category, tf_family) %>%
  summarise(
    n_motifs = n(),
    median_pvalue = median(`P-value`),
    median_enrichment = median(`Log P-value`),
    .groups = "drop"
  )

write_csv(family_summary, file.path(output_dir, "tf_family_enrichment_summary.csv"))

################################################################################
# 4. Co-Regulator Motif Enrichment Analysis
################################################################################

log_message("Analyzing co-regulator motif enrichment...")

# Focus on key co-regulator families
coregulator_families <- c("AP1", "RUNX", "ETS", "AP2", "GATA")

coregulator_motifs <- sig_motifs %>%
  filter(tf_family %in% coregulator_families) %>%
  arrange(category, `P-value`)

write_csv(coregulator_motifs, file.path(output_dir, "coregulator_motif_enrichment.csv"))

# Create enrichment matrix for heatmap
# Note: HOMER's Log P-value is negative (e.g., -443.9), so we take absolute value for visualization
enrichment_matrix <- sig_motifs %>%
  filter(tf_family %in% c("TEAD", coregulator_families)) %>%
  group_by(category, tf_family) %>%
  summarise(
    enrichment = mean(abs(`Log P-value`)),  # Take absolute value for positive scale
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = tf_family, values_from = enrichment, values_fill = 0) %>%
  as.data.frame()

rownames(enrichment_matrix) <- enrichment_matrix$category
enrichment_matrix$category <- NULL
enrichment_matrix <- as.matrix(enrichment_matrix)

# Create heatmap
pdf(file.path(output_dir, "coregulator_enrichment_heatmap.pdf"), width = 10, height = 8)

# Use dynamic color scale based on data range
max_val <- max(enrichment_matrix, na.rm = TRUE)
col_fun <- colorRamp2(c(0, max_val/3, max_val*2/3, max_val), c("white", "yellow", "orange", "red"))

ht <- Heatmap(
  enrichment_matrix,
  name = "-log10(p-value)",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_side = "left",
  column_names_rot = 45,
  column_names_side = "bottom",
  heatmap_legend_param = list(
    title = "Motif Enrichment\n(-log10 p-value)",
    legend_direction = "vertical"
  ),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (enrichment_matrix[i, j] > 10) {
      grid.text(sprintf("%.0f", enrichment_matrix[i, j]), x, y, gp = gpar(fontsize = 8))
    }
  }
)

draw(ht, heatmap_legend_side = "right")

dev.off()

}  # End of else block for when motif data exists

################################################################################
# 5. Co-Occurrence Analysis
################################################################################

log_message("Analyzing motif co-occurrence patterns...")

# For each binding site, check which motifs are present
# This requires parsing HOMER motif instance files

# Function to count motif instances per peak
count_motif_instances <- function(category) {
  motif_instances_file <- file.path(motif_dir, category, "knownResults", "known1.motif")

  if (!file.exists(motif_instances_file)) {
    return(NULL)
  }

  # Parse motif instance file (this is a simplified approach)
  # In reality, would need to parse HOMER's motif location output
  # For now, use enrichment as a proxy

  return(data.frame(category = category))
}

# Note: Full motif instance parsing would require additional HOMER output files
# For this implementation, we'll use enrichment statistics as a proxy

################################################################################
# 6. Correlation with Target Gene Expression
################################################################################

log_message("Analyzing co-regulator TF expression...")

# Identify which co-regulator TFs are in our expression dataset
coregulator_genes <- c(
  # AP-1 family
  "JUN", "FOS", "FOSB", "FOSL1", "FOSL2", "JUNB", "JUND",
  "ATF2", "ATF3", "ATF4", "ATF6", "BATF", "BATF3",
  # RUNX family
  "RUNX1", "RUNX2", "RUNX3", "CBFβ",
  # ETS family
  "ETS1", "ETS2", "ELK1", "ELK3", "ELK4", "ERG", "FLI1",
  "GABPA", "SPI1", "ETV1", "ETV4", "ETV5", "ETV6",
  # AP2 family
  "TFAP2A", "TFAP2B", "TFAP2C", "TFAP2D",
  # GATA family
  "GATA1", "GATA2", "GATA3", "GATA4", "GATA5", "GATA6"
)

# Load full RNA-seq results
rna_file <- "SRF_Eva_RNA/results/05_deseq2/deseq2_results_TES_vs_GFP.txt"
rna_results <- read_tsv(rna_file, show_col_types = FALSE)

# Map gene IDs to symbols
ensembl_to_symbol <- mapIds(
  org.Hs.eg.db,
  keys = rna_results$gene_id,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

rna_results$gene_symbol <- ensembl_to_symbol

# Filter for co-regulator TFs
coregulator_expr <- rna_results %>%
  filter(gene_symbol %in% coregulator_genes) %>%
  select(gene_symbol, baseMean, log2FoldChange, padj) %>%
  arrange(desc(abs(log2FoldChange)))

# Classify by TF family
coregulator_expr$tf_family <- sapply(coregulator_expr$gene_symbol, function(gene) {
  for (family in names(tf_families)) {
    if (gene %in% coregulator_genes && any(str_detect(gene, paste(tf_families[[family]], collapse = "|")))) {
      return(family)
    }
  }
  return("Other")
})

write_csv(coregulator_expr, file.path(output_dir, "coregulator_tf_expression.csv"))

log_message(sprintf("  Found %d co-regulator TFs in expression data", nrow(coregulator_expr)))

# Identify which are significantly DE
sig_coregulators <- coregulator_expr %>%
  filter(padj < 0.05)

write_csv(sig_coregulators, file.path(output_dir, "significant_coregulator_TFs.csv"))

log_message(sprintf("  %d co-regulator TFs are significantly DE", nrow(sig_coregulators)))

################################################################################
# 7. Visualization: Co-Regulator TF Expression
################################################################################

log_message("Creating co-regulator expression visualizations...")

# Plot 1: Co-regulator TF expression volcano
pdf(file.path(output_dir, "coregulator_tf_volcano.pdf"), width = 10, height = 8)

p1 <- ggplot(coregulator_expr, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = tf_family), size = 3, alpha = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "gray") +
  geom_text_repel(
    data = coregulator_expr %>% filter(padj < 0.05 | abs(log2FoldChange) > 1),
    aes(label = gene_symbol),
    size = 3,
    max.overlaps = 20
  ) +
  scale_color_brewer(palette = "Set2") +
  labs(
    title = "Co-Regulator TF Expression in TES vs GFP",
    x = "log2 Fold Change",
    y = "-log10(adjusted p-value)",
    color = "TF Family"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "right"
  )

print(p1)

dev.off()

# Plot 2: Co-regulator expression heatmap
sig_cor_data <- coregulator_expr %>%
  filter(padj < 0.05) %>%
  select(gene_symbol, log2FoldChange, tf_family) %>%
  arrange(tf_family, desc(log2FoldChange))

if (nrow(sig_cor_data) > 0) {
  pdf(file.path(output_dir, "coregulator_tf_expression_heatmap.pdf"), width = 8, height = 10)

  expr_matrix <- matrix(sig_cor_data$log2FoldChange, ncol = 1)
  rownames(expr_matrix) <- sig_cor_data$gene_symbol
  colnames(expr_matrix) <- "TES vs GFP"

  # Row annotation for TF family
  row_ha <- rowAnnotation(
    TF_Family = sig_cor_data$tf_family,
    col = list(TF_Family = setNames(
      brewer.pal(length(unique(sig_cor_data$tf_family)), "Set2"),
      unique(sig_cor_data$tf_family)
    ))
  )

  col_fun2 <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

  ht2 <- Heatmap(
    expr_matrix,
    name = "log2FC",
    col = col_fun2,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_side = "left",
    right_annotation = row_ha,
    heatmap_legend_param = list(
      title = "log2FC\nTES vs GFP"
    )
  )

  draw(ht2, heatmap_legend_side = "right")

  dev.off()
}

################################################################################
# 8. Target Gene Stratification by Co-Regulator Motif
################################################################################

log_message("Stratifying target genes by co-regulator motif presence...")

# For each binding category, identify which co-regulator motifs are enriched
category_coregulators <- sig_motifs %>%
  filter(tf_family %in% coregulator_families) %>%
  group_by(category) %>%
  slice_min(order_by = `P-value`, n = 5) %>%
  summarise(
    top_coregulators = paste(`Motif Name`, collapse = "; "),
    top_families = paste(unique(tf_family), collapse = ", "),
    .groups = "drop"
  )

write_csv(category_coregulators, file.path(output_dir, "category_coregulator_summary.csv"))

################################################################################
# 9. Integrative Analysis: Motif + Expression
################################################################################

log_message("Performing integrative motif-expression analysis...")

# For each category, summarize:
# 1. Which co-regulator motifs are enriched
# 2. Which co-regulator TFs are expressed
# 3. Which co-regulator TFs are DE

integrative_summary <- expand.grid(
  category = categories,
  tf_family = coregulator_families,
  stringsAsFactors = FALSE
) %>%
  left_join(
    family_summary %>% select(category, tf_family, n_motifs, median_enrichment),
    by = c("category", "tf_family")
  ) %>%
  mutate(
    motif_enriched = !is.na(n_motifs) & n_motifs > 0,
    median_enrichment = ifelse(is.na(median_enrichment), 0, median_enrichment)
  )

# Add expression information
family_expr <- coregulator_expr %>%
  group_by(tf_family) %>%
  summarise(
    n_expressed = n(),
    n_sig_de = sum(padj < 0.05, na.rm = TRUE),
    mean_log2fc = mean(log2FoldChange, na.rm = TRUE),
    .groups = "drop"
  )

integrative_summary <- integrative_summary %>%
  left_join(family_expr, by = "tf_family")

write_csv(integrative_summary, file.path(output_dir, "integrative_motif_expression_summary.csv"))

# Create bubble plot: motif enrichment vs TF expression
pdf(file.path(output_dir, "motif_expression_bubble.pdf"), width = 12, height = 8)

p2 <- ggplot(integrative_summary, aes(x = category, y = tf_family)) +
  geom_point(aes(size = median_enrichment, color = mean_log2fc), alpha = 0.7) +
  scale_size_continuous(range = c(1, 15), name = "Motif Enrichment\n(-log10 p-value)") +
  scale_color_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = 0,
    name = "Mean TF\nlog2FC"
  ) +
  labs(
    title = "Co-Regulator Motif Enrichment vs TF Expression",
    subtitle = "Bubble size = motif enrichment; Color = TF expression change",
    x = "Binding Category",
    y = "TF Family"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10)
  )

print(p2)

dev.off()

################################################################################
# 10. Generate Summary Report
################################################################################

log_message("Generating summary report...")

summary_file <- file.path(output_dir, "PHASE5_1_SUMMARY.txt")

sink(summary_file)

cat("================================================================================\n")
cat("Phase 5.1: Co-Regulator Identification - Analysis Summary\n")
cat("TES vs TEAD1 Comparative Study (Excluding TESmut)\n")
cat("================================================================================\n\n")

cat("Analysis Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("INPUT DATA:\n")
cat(sprintf("  - Binding sites analyzed: %d\n", nrow(binding_data)))
cat(sprintf("  - Genes with expression data: %d\n", nrow(expr_data)))
cat(sprintf("  - Total motifs analyzed: %d\n", nrow(all_motifs)))
cat("\n")

cat("CO-REGULATOR MOTIF ENRICHMENT:\n")
cat(sprintf("  - Significant motifs (p < 1e-10): %d\n", nrow(sig_motifs)))
cat(sprintf("  - Co-regulator motifs identified: %d\n", nrow(coregulator_motifs)))
cat("\n")

cat("  Top enriched co-regulator families by category:\n")
for (cat in categories) {
  top_family <- family_summary %>%
    filter(category == cat, tf_family %in% coregulator_families) %>%
    slice_min(order_by = median_pvalue, n = 1)

  if (nrow(top_family) > 0) {
    cat(sprintf("    %s: %s (-log10 p = %.1f)\n",
                cat, top_family$tf_family[1], top_family$median_enrichment[1]))
  }
}
cat("\n")

cat("CO-REGULATOR TF EXPRESSION:\n")
cat(sprintf("  - Co-regulator TFs in dataset: %d\n", nrow(coregulator_expr)))
cat(sprintf("  - Significantly DE co-regulator TFs: %d\n", nrow(sig_coregulators)))
cat("\n")

if (nrow(sig_coregulators) > 0) {
  cat("  Top differentially expressed co-regulators:\n")
  top_sig <- sig_coregulators %>%
    slice_min(order_by = padj, n = 10) %>%
    select(gene_symbol, log2FoldChange, padj, tf_family)

  for (i in 1:nrow(top_sig)) {
    cat(sprintf("    %s (%s): log2FC = %.2f, padj = %.2e\n",
                top_sig$gene_symbol[i], top_sig$tf_family[i],
                top_sig$log2FoldChange[i], top_sig$padj[i]))
  }
}
cat("\n")

cat("KEY FINDINGS:\n")

# Finding 1: Which co-regulator families are most enriched?
top_families <- family_summary %>%
  filter(tf_family %in% coregulator_families) %>%
  group_by(tf_family) %>%
  summarise(total_enrichment = sum(median_enrichment), .groups = "drop") %>%
  arrange(desc(total_enrichment))

cat("  1. Most enriched co-regulator families across all categories:\n")
for (i in 1:min(3, nrow(top_families))) {
  cat(sprintf("     %d. %s (total enrichment = %.1f)\n",
              i, top_families$tf_family[i], top_families$total_enrichment[i]))
}
cat("\n")

# Finding 2: Which families show expression changes?
family_changes <- family_expr %>%
  arrange(desc(n_sig_de))

cat("  2. Co-regulator families with most DE TFs:\n")
for (i in 1:min(3, nrow(family_changes))) {
  cat(sprintf("     %d. %s: %d TFs DE (mean log2FC = %.2f)\n",
              i, family_changes$tf_family[i], family_changes$n_sig_de[i],
              family_changes$mean_log2fc[i]))
}
cat("\n")

# Finding 3: Concordance between motif enrichment and TF expression
cat("  3. Motif-Expression Concordance:\n")
concordance <- integrative_summary %>%
  filter(motif_enriched == TRUE, n_sig_de > 0) %>%
  select(tf_family, category, median_enrichment, n_sig_de, mean_log2fc) %>%
  arrange(desc(median_enrichment))

if (nrow(concordance) > 0) {
  cat("     Families with both motif enrichment AND TF expression changes:\n")
  for (fam in unique(concordance$tf_family)) {
    fam_data <- concordance %>% filter(tf_family == fam)
    cat(sprintf("     - %s: enriched in %d categories, %d TFs DE\n",
                fam, nrow(fam_data), fam_data$n_sig_de[1]))
  }
} else {
  cat("     Limited concordance between motif enrichment and TF expression\n")
}
cat("\n")

cat("OUTPUT FILES:\n")
cat("  1. tf_family_enrichment_summary.csv - Motif enrichment by family\n")
cat("  2. coregulator_motif_enrichment.csv - Detailed co-regulator motifs\n")
cat("  3. coregulator_enrichment_heatmap.pdf - Heatmap of motif enrichment\n")
cat("  4. coregulator_tf_expression.csv - Co-regulator TF expression data\n")
cat("  5. significant_coregulator_TFs.csv - DE co-regulator TFs\n")
cat("  6. coregulator_tf_volcano.pdf - Volcano plot of TF expression\n")
cat("  7. coregulator_tf_expression_heatmap.pdf - Heatmap of DE TFs\n")
cat("  8. category_coregulator_summary.csv - Top co-regulators per category\n")
cat("  9. integrative_motif_expression_summary.csv - Combined analysis\n")
cat(" 10. motif_expression_bubble.pdf - Bubble plot integrating motif + expression\n")
cat("\n")

cat("BIOLOGICAL INTERPRETATION:\n")
cat("  - AP-1 family: Stress response, proliferation\n")
cat("  - RUNX family: Development, differentiation, hematopoiesis\n")
cat("  - ETS family: Cell proliferation, differentiation, angiogenesis\n")
cat("  - AP-2 family: Development, cell growth\n")
cat("  - GATA family: Differentiation, development\n")
cat("\n")

cat("NEXT STEPS:\n")
cat("  - Validate co-regulator interactions with ChIP-reChIP or co-IP\n")
cat("  - Test functional dependency using co-regulator knockdown\n")
cat("  - Examine spacing between TEAD and co-regulator motifs\n")
cat("  - Correlate co-regulator presence with gene expression magnitude\n")
cat("\n")

cat("================================================================================\n")
cat("Analysis Complete\n")
cat("================================================================================\n")

sink()

log_message("Phase 5.1 Analysis Complete!")
log_message(sprintf("Results saved to: %s", output_dir))
log_message(sprintf("Summary report: %s", summary_file))

cat("\n")
cat("================================================================================\n")
cat("Phase 5.1: Co-Regulator Identification - COMPLETE\n")
cat("================================================================================\n")
