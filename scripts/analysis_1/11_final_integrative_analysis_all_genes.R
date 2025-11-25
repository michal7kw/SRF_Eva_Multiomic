#!/usr/bin/env Rscript
#
# FINAL INTEGRATIVE ANALYSIS: Cut&Tag + RNA-seq (ALL GENES VERSION)
# Comprehensive analysis of TES/TEAD1 transcriptional regulatory networks
#

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(GenomicRanges)
  library(ChIPseeker)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(VennDiagram)
  library(pheatmap)
  library(ComplexHeatmap)
  library(enrichplot)
  library(readr)
  library(stringr)
  library(gridExtra)
})

# Set working directory and create output directories
setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_1")

# =============================================================================
# OUTPUT PATH CONFIGURATION
# =============================================================================
# Define all output paths as variables for easy modification
DATA_DIR <- "../../data"
OUTPUT_BASE <- "output/11_final_integrative_analysis_all_genes"
OUTPUT_DIRECT_TARGETS <- file.path(OUTPUT_BASE, "direct_targets")
OUTPUT_PATHWAY <- file.path(OUTPUT_BASE, "pathway_analysis")
OUTPUT_REGULATORY <- file.path(OUTPUT_BASE, "regulatory_networks")
OUTPUT_PLOTS <- file.path(OUTPUT_BASE, "plots")
LOG_DIR <- "logs"

# Input paths
RNA_SEQ_RESULTS <- "../../../SRF_Eva_RNA/results/05_deseq2/deseq2_results_TES_vs_GFP.txt"
RNA_SEQ_SIGNIFICANT <- "../../../SRF_Eva_RNA/results/05_deseq2/significant_genes_TES_vs_GFP.txt"
TES_PEAKS_FILE <- "../../../SRF_Eva_CUTandTAG/results/07_analysis_narrow/TES_peaks_annotated.csv"
TEAD1_PEAKS_FILE <- "../../../SRF_Eva_CUTandTAG/results/07_analysis_narrow/TEAD1_peaks_annotated.csv"

# Ensure all output directories exist
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(OUTPUT_DIRECT_TARGETS, showWarnings = FALSE, recursive = TRUE)
dir.create(OUTPUT_PATHWAY, showWarnings = FALSE, recursive = TRUE)
dir.create(OUTPUT_REGULATORY, showWarnings = FALSE, recursive = TRUE)
dir.create(OUTPUT_PLOTS, showWarnings = FALSE, recursive = TRUE)
dir.create(LOG_DIR, showWarnings = FALSE, recursive = TRUE)

cat("=== FINAL INTEGRATIVE ANALYSIS: Cut&Tag + RNA-seq (ALL GENES) ===\n")
cat("Analysis started:", as.character(Sys.time()), "\n\n")

# =============================================================================
# PHASE 1: DATA LOADING AND PROCESSING
# =============================================================================

cat("=== PHASE 1: Data Loading and Processing ===\n")

# Load RNA-seq differential expression data (ALL GENES)
cat("Loading RNA-seq differential expression data...\n")
rna_results <- read.delim(RNA_SEQ_RESULTS,
  stringsAsFactors = FALSE
)
rna_significant <- read.delim(RNA_SEQ_SIGNIFICANT,
  stringsAsFactors = FALSE, header = FALSE
)
colnames(rna_significant) <- "gene_symbol"

cat(sprintf(
  "✓ RNA-seq data loaded: %d total genes, %d significant DE genes\n",
  nrow(rna_results), nrow(rna_significant)
))
cat("✓ Using ALL genes for comprehensive analysis\n")

# Load Cut&Tag peak annotation data
cat("Loading Cut&Tag peak annotation data...\n")
tes_peaks <- read.csv(TES_PEAKS_FILE,
  stringsAsFactors = FALSE
)
tead1_peaks <- read.csv(TEAD1_PEAKS_FILE,
  stringsAsFactors = FALSE
)

cat(sprintf(
  "✓ Cut&Tag data loaded: %d TES peaks, %d TEAD1 peaks\n",
  nrow(tes_peaks), nrow(tead1_peaks)
))

# Process gene IDs for integration
cat("Processing gene IDs for integration...\n")

# Clean Ensembl IDs (remove version numbers)
rna_results$ensembl_id <- gsub("\\..*", "", rna_results$gene_id)
rna_results$is_significant <- !is.na(rna_results$padj) & rna_results$padj < 0.05

# Gene ID conversion functions
convert_ensembl_to_symbol <- function(ensembl_ids) {
  symbols <- mapIds(org.Hs.eg.db,
    keys = ensembl_ids,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  return(symbols)
}

convert_entrez_to_ensembl <- function(entrez_ids) {
  ensembl <- mapIds(org.Hs.eg.db,
    keys = as.character(entrez_ids),
    column = "ENSEMBL",
    keytype = "ENTREZID",
    multiVals = "first"
  )
  return(ensembl)
}

# Perform gene ID conversions
cat("Converting gene IDs...\n")
rna_results$gene_symbol <- convert_ensembl_to_symbol(rna_results$ensembl_id)
tes_peaks$ensembl_id <- convert_entrez_to_ensembl(tes_peaks$geneId)
tead1_peaks$ensembl_id <- convert_entrez_to_ensembl(tead1_peaks$geneId)

# Report conversion success
cat(sprintf(
  "✓ RNA-seq: %d genes with symbols mapped (%.1f%%)\n",
  sum(!is.na(rna_results$gene_symbol)),
  100 * sum(!is.na(rna_results$gene_symbol)) / nrow(rna_results)
))
cat(sprintf(
  "✓ TES peaks: %d unique Entrez genes, %d converted to Ensembl (%.1f%%)\n",
  length(unique(tes_peaks$geneId)),
  sum(!is.na(tes_peaks$ensembl_id)),
  100 * sum(!is.na(tes_peaks$ensembl_id)) / nrow(tes_peaks)
))
cat(sprintf(
  "✓ TEAD1 peaks: %d unique Entrez genes, %d converted to Ensembl (%.1f%%)\n",
  length(unique(tead1_peaks$geneId)),
  sum(!is.na(tead1_peaks$ensembl_id)),
  100 * sum(!is.na(tead1_peaks$ensembl_id)) / nrow(tead1_peaks)
))

# Save processed data
save(rna_results, rna_significant, tes_peaks, tead1_peaks,
  file = file.path(DATA_DIR, "integrated_data_all_genes.RData")
)
cat(sprintf("✓ Processed data saved to %s/integrated_data_all_genes.RData\n\n", DATA_DIR))

# =============================================================================
# PHASE 2: PEAK-GENE ASSOCIATION ANALYSIS
# =============================================================================

cat("=== PHASE 2: Peak-Gene Association Analysis ===\n")

# Advanced peak-gene mapping function
get_peak_genes <- function(peaks, distance_threshold = 50000, tf_name = "TF") {
  valid_peaks <- peaks[!is.na(peaks$ensembl_id), ]

  # Classify peaks by genomic location
  promoter_peaks <- valid_peaks[grepl("Promoter", valid_peaks$annotation), ]
  enhancer_peaks <- valid_peaks[abs(valid_peaks$distanceToTSS) <= distance_threshold &
    !grepl("Promoter", valid_peaks$annotation), ]

  # Get unique genes
  promoter_genes <- unique(promoter_peaks$ensembl_id)
  enhancer_genes <- unique(enhancer_peaks$ensembl_id)
  output <- unique(c(promoter_genes, enhancer_genes))

  # Report statistics
  cat(sprintf("%s peak analysis:\n", tf_name))
  cat(sprintf("  Valid peaks with Ensembl ID: %d\n", nrow(valid_peaks)))
  cat(sprintf("  Promoter peaks: %d (genes: %d)\n", nrow(promoter_peaks), length(promoter_genes)))
  cat(sprintf(
    "  Enhancer peaks (<%dkb): %d (genes: %d)\n",
    distance_threshold / 1000, nrow(enhancer_peaks), length(enhancer_genes)
  ))
  cat(sprintf("  Total bound genes: %d\n\n", length(output)))

  return(list(
    promoter_genes = promoter_genes,
    enhancer_genes = enhancer_genes,
    output = output,
    promoter_peaks = promoter_peaks,
    enhancer_peaks = enhancer_peaks,
    valid_peaks = valid_peaks
  ))
}

# Analyze TES and TEAD1 peak-gene associations
tes_gene_mapping <- get_peak_genes(tes_peaks, tf_name = "TES")
tead1_gene_mapping <- get_peak_genes(tead1_peaks, tf_name = "TEAD1")

# =============================================================================
# PHASE 3: DIRECT TARGET IDENTIFICATION
# =============================================================================

cat("=== PHASE 3: Direct Target Identification ===\n")

# Get bound gene lists
tes_bound_genes <- tes_gene_mapping$output
tead1_bound_genes <- tead1_gene_mapping$output

# Identify bound genes in RNA-seq data
rna_results$tes_bound <- rna_results$ensembl_id %in% tes_bound_genes
rna_results$tead1_bound <- rna_results$ensembl_id %in% tead1_bound_genes

# Identify direct targets (bound + differentially expressed)
tes_direct_targets <- rna_results[rna_results$tes_bound & rna_results$is_significant, ]
tead1_direct_targets <- rna_results[rna_results$tead1_bound & rna_results$is_significant, ]

# Calculate overlaps and specific targets
shared_targets <- intersect(tes_direct_targets$ensembl_id, tead1_direct_targets$ensembl_id)
tes_specific <- tes_direct_targets[!tes_direct_targets$ensembl_id %in% tead1_direct_targets$ensembl_id, ]
tead1_specific <- tead1_direct_targets[!tead1_direct_targets$ensembl_id %in% tes_direct_targets$ensembl_id, ]

# Report direct target statistics
cat("Direct targets identified:\n")
cat(sprintf("  TES: %d genes (bound + DE)\n", nrow(tes_direct_targets)))
cat(sprintf("  TEAD1: %d genes (bound + DE)\n", nrow(tead1_direct_targets)))
cat(sprintf("  Shared TES/TEAD1: %d genes\n", length(shared_targets)))
cat(sprintf("  TES-specific: %d genes\n", nrow(tes_specific)))
cat(sprintf("  TEAD1-specific: %d genes\n\n", nrow(tead1_specific)))

# Activation vs repression analysis
analyze_regulation_mode <- function(targets, tf_name) {
  upregulated <- sum(targets$log2FoldChange > 0)
  downregulated <- sum(targets$log2FoldChange < 0)
  total <- nrow(targets)

  cat(sprintf("%s regulation mode:\n", tf_name))
  cat(sprintf("  Upregulated: %d (%.1f%%)\n", upregulated, 100 * upregulated / total))
  cat(sprintf("  Downregulated: %d (%.1f%%)\n", downregulated, 100 * downregulated / total))
  cat(sprintf("  Net effect: %.1f%% activation\n\n", 100 * upregulated / total))

  return(c(up = upregulated, down = downregulated, total = total))
}

tes_regulation <- analyze_regulation_mode(tes_direct_targets, "TES")
tead1_regulation <- analyze_regulation_mode(tead1_direct_targets, "TEAD1")

# Save direct target lists
write.csv(tes_direct_targets, file.path(OUTPUT_DIRECT_TARGETS, "TES_direct_targets_all_genes.csv"), row.names = FALSE)
write.csv(tead1_direct_targets, file.path(OUTPUT_DIRECT_TARGETS, "TEAD1_direct_targets_all_genes.csv"), row.names = FALSE)
write.csv(tes_specific, file.path(OUTPUT_DIRECT_TARGETS, "TES_specific_targets_all_genes.csv"), row.names = FALSE)
write.csv(tead1_specific, file.path(OUTPUT_DIRECT_TARGETS, "TEAD1_specific_targets_all_genes.csv"), row.names = FALSE)
cat(sprintf("✓ Direct target lists saved to %s/\n\n", OUTPUT_DIRECT_TARGETS))

# =============================================================================
# PHASE 4: COMPREHENSIVE GENE CLASSIFICATION
# =============================================================================

cat("=== PHASE 4: Comprehensive Gene Classification ===\n")

# Create comprehensive classification table
classification <- data.frame(
  gene_id = rna_results$ensembl_id,
  gene_symbol = rna_results$gene_symbol,
  log2FoldChange = rna_results$log2FoldChange,
  padj = rna_results$padj,
  tes_bound = rna_results$tes_bound,
  tead1_bound = rna_results$tead1_bound,
  is_significant = rna_results$is_significant,
  stringsAsFactors = FALSE
)

# Assign regulatory classes
classification$regulatory_class <- "Unbound"
classification$regulatory_class[classification$tes_bound & classification$is_significant] <- "TES_direct"
classification$regulatory_class[classification$tead1_bound & classification$is_significant] <- "TEAD1_direct"
classification$regulatory_class[classification$tes_bound & classification$tead1_bound & classification$is_significant] <- "TES_TEAD1_shared"
classification$regulatory_class[!classification$tes_bound & !classification$tead1_bound & classification$is_significant] <- "Indirect"

# Report classification statistics
class_summary <- table(classification$regulatory_class)
cat("Gene classification summary:\n")
print(class_summary)
cat("\n")

# Save classification table
write.csv(classification, file.path(OUTPUT_REGULATORY, "gene_classification_all_genes.csv"), row.names = FALSE)
cat(sprintf("✓ Gene classification table saved to %s/\n\n", OUTPUT_REGULATORY))

# =============================================================================
# PHASE 5: PATHWAY ENRICHMENT ANALYSIS
# =============================================================================

cat("=== PHASE 5: Pathway Enrichment Analysis ===\n")

# Pathway enrichment function
perform_enrichment <- function(gene_list, background, analysis_name, output_dir = OUTPUT_PATHWAY) {
  if (length(gene_list) < 10) {
    cat(sprintf("Skipping %s: too few genes (%d)\n", analysis_name, length(gene_list)))
    return(NULL)
  }

  # Remove genes not in background
  gene_list <- gene_list[gene_list %in% background]
  if (length(gene_list) < 10) {
    cat(sprintf("Skipping %s: too few genes after background filtering (%d)\n", analysis_name, length(gene_list)))
    return(NULL)
  }

  cat(sprintf("Running GO enrichment for %s (%d genes)...\n", analysis_name, length(gene_list)))

  ego <- enrichGO(
    gene = gene_list,
    universe = background,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05,
    keyType = "ENSEMBL"
  )

  if (nrow(ego) > 0) {
    cat(sprintf("✓ %s: %d enriched pathways found\n", analysis_name, nrow(ego)))
    write.csv(ego@result, file.path(output_dir, sprintf("%s_GO_enrichment_all_genes.csv", analysis_name)), row.names = FALSE)
    return(ego)
  } else {
    cat(sprintf("✗ %s: no enriched pathways found\n", analysis_name))
    return(NULL)
  }
}

# Define background gene set (all genes with symbols)
background_genes <- rna_results$ensembl_id[!is.na(rna_results$gene_symbol)]

cat(sprintf("Background for pathway analysis: %d genes with gene symbols\n", length(background_genes)))

# Run pathway enrichment for different gene sets
tes_go <- perform_enrichment(tes_direct_targets$ensembl_id, background_genes, "TES_direct")
tead1_go <- perform_enrichment(tead1_direct_targets$ensembl_id, background_genes, "TEAD1_direct")
shared_go <- perform_enrichment(shared_targets, background_genes, "TES_TEAD1_shared")
indirect_genes <- classification$gene_id[classification$regulatory_class == "Indirect"]
indirect_go <- perform_enrichment(indirect_genes, background_genes, "Indirect")

cat("✓ Pathway enrichment analysis completed\n\n")

# =============================================================================
# PHASE 6: COMPREHENSIVE VISUALIZATION (SEPARATE PLOTS)
# =============================================================================

cat("=== PHASE 6: Comprehensive Visualization ===\n")

# Plot 1: Gene classification pie chart
cat("Creating gene classification pie chart...\n")
class_counts <- table(classification$regulatory_class)
pie_data <- data.frame(
  class = names(class_counts),
  count = as.numeric(class_counts),
  percentage = round(100 * as.numeric(class_counts) / sum(class_counts), 1)
)

pdf(file.path(OUTPUT_PLOTS, "01_gene_classification_pie_all_genes.pdf"), width = 10, height = 8)
p1 <- ggplot(pie_data, aes(x = "", y = count, fill = class)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  labs(
    title = "Gene Classification: Binding vs Expression (All Genes)",
    subtitle = sprintf("Total genes analyzed: %d", nrow(classification)),
    fill = "Regulatory Class"
  ) +
  theme_void() +
  theme(
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(size = 14, hjust = 0.5),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold")
  ) +
  scale_fill_brewer(type = "qual", palette = "Set2") +
  geom_text(aes(label = paste0(count, "\n(", percentage, "%)")),
    position = position_stack(vjust = 0.5), size = 4, fontface = "bold"
  )
print(p1)
dev.off()

# Plot 2: Expression changes histograms by regulatory class
cat("Creating expression histograms by regulatory class...\n")
sig_genes <- classification[classification$is_significant, ]
pdf(file.path(OUTPUT_PLOTS, "02_expression_histograms_all_genes.pdf"), width = 14, height = 10)
p2 <- ggplot(sig_genes, aes(x = log2FoldChange, fill = regulatory_class)) +
  geom_histogram(bins = 40, alpha = 0.8, color = "white") +
  facet_wrap(~regulatory_class, scales = "free_y", ncol = 2) +
  labs(
    title = "Expression Changes by Regulatory Class (All Genes)",
    subtitle = "Distribution of Log2 Fold Changes (TES vs GFP)",
    x = "Log2 Fold Change (TES vs GFP)",
    y = "Number of Genes"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14),
    strip.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 12),
    legend.position = "none"
  ) +
  scale_fill_brewer(type = "qual", palette = "Set2") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", alpha = 0.7)
print(p2)
dev.off()

# Plot 3: Boxplot comparison
cat("Creating boxplot comparison...\n")
boxplot_data <- sig_genes[sig_genes$regulatory_class %in% c("TES_direct", "TEAD1_direct", "TES_TEAD1_shared", "Indirect"), ]
pdf(file.path(OUTPUT_PLOTS, "03_expression_boxplots_all_genes.pdf"), width = 12, height = 8)
p3 <- ggplot(boxplot_data, aes(x = regulatory_class, y = log2FoldChange, fill = regulatory_class)) +
  geom_boxplot(alpha = 0.8, outlier.alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 1) +
  labs(
    title = "Expression Changes by Regulatory Mechanism (All Genes)",
    subtitle = "Comparison of direct vs indirect regulation",
    x = "Regulatory Class",
    y = "Log2 Fold Change (TES vs GFP)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.title = element_text(size = 12),
    legend.position = "none"
  ) +
  scale_fill_brewer(type = "qual", palette = "Set2") +
  stat_summary(fun = median, geom = "point", shape = 23, size = 3, fill = "white")
print(p3)
dev.off()

# Plot 4: Venn diagram
cat("Creating Venn diagram...\n")
venn_list <- list(
  "TES Direct" = tes_direct_targets$ensembl_id,
  "TEAD1 Direct" = tead1_direct_targets$ensembl_id,
  "All DE Genes" = rna_results$ensembl_id[rna_results$is_significant]
)

pdf(file.path(OUTPUT_PLOTS, "04_target_overlap_venn_all_genes.pdf"), width = 10, height = 8)
grid.newpage()
venn_plot <- venn.diagram(
  x = venn_list,
  category.names = names(venn_list),
  filename = NULL,
  output = TRUE,
  col = "transparent",
  fill = c("#E31A1C", "#1F78B4", "#33A02C"),
  alpha = 0.6,
  cex = 1.8,
  fontfamily = "serif",
  cat.cex = 1.4,
  cat.fontfamily = "serif",
  cat.col = c("#E31A1C", "#1F78B4", "#33A02C"),
  main = "Overlap: Direct Targets vs All DE Genes",
  main.cex = 1.8,
  main.fontfamily = "serif"
)
grid.draw(venn_plot)
dev.off()

# Plot 5: Summary bar chart
cat("Creating summary bar chart...\n")
summary_data <- data.frame(
  Category = c("TES Direct", "TEAD1 Direct", "Shared", "TES Specific", "TEAD1 Specific", "Indirect", "Unbound"),
  Count = c(
    nrow(tes_direct_targets), nrow(tead1_direct_targets), length(shared_targets),
    nrow(tes_specific), nrow(tead1_specific), sum(classification$regulatory_class == "Indirect"),
    sum(classification$regulatory_class == "Unbound")
  ),
  Type = c("Direct", "Direct", "Direct", "Specific", "Specific", "Indirect", "Unbound")
)

pdf(file.path(OUTPUT_PLOTS, "05_summary_barchart_all_genes.pdf"), width = 12, height = 8)
p5 <- ggplot(summary_data, aes(x = reorder(Category, Count), y = Count, fill = Type)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  labs(
    title = "Summary of Target Gene Categories (All Genes)",
    subtitle = "Number of genes in each regulatory class",
    x = "Regulatory Category",
    y = "Number of Genes"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold")
  ) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  geom_text(aes(label = Count), vjust = -0.3, size = 4, fontface = "bold") +
  coord_flip()
print(p5)
dev.off()

# Create pathway enrichment plots
if (!is.null(tes_go) && nrow(tes_go) > 0) {
  pdf(file.path(OUTPUT_PLOTS, "TES_pathway_enrichment_all_genes.pdf"), width = 14, height = 10)
  print(dotplot(tes_go, showCategory = 20) +
    ggtitle("TES Direct Targets: GO Biological Process Enrichment (All Genes)") +
    theme(plot.title = element_text(size = 16)))
  dev.off()
}

if (!is.null(tead1_go) && nrow(tead1_go) > 0) {
  pdf(file.path(OUTPUT_PLOTS, "TEAD1_pathway_enrichment_all_genes.pdf"), width = 14, height = 10)
  print(dotplot(tead1_go, showCategory = 20) +
    ggtitle("TEAD1 Direct Targets: GO Biological Process Enrichment (All Genes)") +
    theme(plot.title = element_text(size = 16)))
  dev.off()
}

if (!is.null(shared_go) && nrow(shared_go) > 0) {
  pdf(file.path(OUTPUT_PLOTS, "Shared_pathway_enrichment_all_genes.pdf"), width = 14, height = 10)
  print(dotplot(shared_go, showCategory = 20) +
    ggtitle("TES/TEAD1 Shared Targets: GO Biological Process Enrichment (All Genes)") +
    theme(plot.title = element_text(size = 16)))
  dev.off()
}

cat(sprintf("✓ Individual visualization plots created in %s/:\n", OUTPUT_PLOTS))
cat("  - 01_gene_classification_pie_all_genes.pdf\n")
cat("  - 02_expression_histograms_all_genes.pdf\n")
cat("  - 03_expression_boxplots_all_genes.pdf\n")
cat("  - 04_target_overlap_venn_all_genes.pdf\n")
cat("  - 05_summary_barchart_all_genes.pdf\n")
cat("  - TES_pathway_enrichment_all_genes.pdf (if enrichment found)\n")
cat("  - TEAD1_pathway_enrichment_all_genes.pdf (if enrichment found)\n")
cat("  - Shared_pathway_enrichment_all_genes.pdf (if enrichment found)\n\n")

# =============================================================================
# PHASE 7: COMPREHENSIVE SUMMARY REPORT
# =============================================================================

cat("=== PHASE 7: Final Summary Report ===\n")

# Compile comprehensive summary statistics
summary_stats <- list(
  analysis_timestamp = as.character(Sys.time()),

  # Data overview
  total_genes_analyzed = nrow(rna_results),
  genes_with_symbols = sum(!is.na(rna_results$gene_symbol)),
  significant_de_genes = sum(rna_results$is_significant),

  # Peak statistics
  tes_total_peaks = nrow(tes_peaks),
  tes_peaks_with_ensembl = sum(!is.na(tes_peaks$ensembl_id)),
  tes_bound_genes = length(tes_bound_genes),
  tead1_total_peaks = nrow(tead1_peaks),
  tead1_peaks_with_ensembl = sum(!is.na(tead1_peaks$ensembl_id)),
  tead1_bound_genes = length(tead1_bound_genes),

  # Direct targets
  tes_direct_targets = nrow(tes_direct_targets),
  tead1_direct_targets = nrow(tead1_direct_targets),
  shared_targets = length(shared_targets),
  tes_specific_targets = nrow(tes_specific),
  tead1_specific_targets = nrow(tead1_specific),

  # Regulatory modes
  tes_upregulated = tes_regulation["up"],
  tes_downregulated = tes_regulation["down"],
  tead1_upregulated = tead1_regulation["up"],
  tead1_downregulated = tead1_regulation["down"],

  # Classification
  indirect_targets = sum(classification$regulatory_class == "Indirect"),
  unbound_genes = sum(classification$regulatory_class == "Unbound"),

  # Pathway enrichment
  tes_enriched_pathways = ifelse(!is.null(tes_go), nrow(tes_go), 0),
  tead1_enriched_pathways = ifelse(!is.null(tead1_go), nrow(tead1_go), 0),
  shared_enriched_pathways = ifelse(!is.null(shared_go), nrow(shared_go), 0)
)

# Print summary
cat("FINAL ANALYSIS SUMMARY (ALL GENES):\n")
cat("===================================\n")
for (name in names(summary_stats)) {
  cat(sprintf("%-25s: %s\n", name, summary_stats[[name]]))
}

# Save detailed summary
writeLines(c(
  "INTEGRATIVE ANALYSIS SUMMARY REPORT (ALL GENES)",
  "===============================================",
  paste("Generated:", summary_stats$analysis_timestamp),
  "",
  "DATA OVERVIEW:",
  sprintf("  Total genes analyzed: %d", summary_stats$total_genes_analyzed),
  sprintf(
    "  Genes with symbols: %d (%.1f%%)", summary_stats$genes_with_symbols,
    100 * summary_stats$genes_with_symbols / summary_stats$total_genes_analyzed
  ),
  sprintf(
    "  Significant DE genes: %d (%.1f%%)", summary_stats$significant_de_genes,
    100 * summary_stats$significant_de_genes / summary_stats$total_genes_analyzed
  ),
  "",
  "BINDING ANALYSIS:",
  sprintf(
    "  TES peaks: %d total, %d with Ensembl ID, %d bound genes",
    summary_stats$tes_total_peaks, summary_stats$tes_peaks_with_ensembl, summary_stats$tes_bound_genes
  ),
  sprintf(
    "  TEAD1 peaks: %d total, %d with Ensembl ID, %d bound genes",
    summary_stats$tead1_total_peaks, summary_stats$tead1_peaks_with_ensembl, summary_stats$tead1_bound_genes
  ),
  "",
  "DIRECT TARGETS:",
  sprintf("  TES direct targets: %d", summary_stats$tes_direct_targets),
  sprintf("  TEAD1 direct targets: %d", summary_stats$tead1_direct_targets),
  sprintf("  Shared targets: %d", summary_stats$shared_targets),
  sprintf("  TES-specific: %d", summary_stats$tes_specific_targets),
  sprintf("  TEAD1-specific: %d", summary_stats$tead1_specific_targets),
  "",
  "REGULATORY MODES:",
  sprintf(
    "  TES: %d upregulated, %d downregulated (%.1f%% activation)",
    summary_stats$tes_upregulated, summary_stats$tes_downregulated,
    100 * summary_stats$tes_upregulated / (summary_stats$tes_upregulated + summary_stats$tes_downregulated)
  ),
  sprintf(
    "  TEAD1: %d upregulated, %d downregulated (%.1f%% activation)",
    summary_stats$tead1_upregulated, summary_stats$tead1_downregulated,
    100 * summary_stats$tead1_upregulated / (summary_stats$tead1_upregulated + summary_stats$tead1_downregulated)
  ),
  "",
  "GENE CLASSIFICATION:",
  sprintf("  Direct targets: %d", summary_stats$tes_direct_targets + summary_stats$tead1_direct_targets - summary_stats$shared_targets),
  sprintf("  Indirect targets: %d", summary_stats$indirect_targets),
  sprintf("  Unbound genes: %d", summary_stats$unbound_genes),
  "",
  "PATHWAY ENRICHMENT:",
  sprintf("  TES pathways: %d", summary_stats$tes_enriched_pathways),
  sprintf("  TEAD1 pathways: %d", summary_stats$tead1_enriched_pathways),
  sprintf("  Shared pathways: %d", summary_stats$shared_enriched_pathways)
), file.path(OUTPUT_BASE, "FINAL_ANALYSIS_SUMMARY_ALL_GENES.txt"))

# Save R object with all results
save(rna_results, rna_significant, tes_peaks, tead1_peaks,
  tes_gene_mapping, tead1_gene_mapping,
  tes_direct_targets, tead1_direct_targets,
  classification, summary_stats,
  tes_go, tead1_go, shared_go,
  file = file.path(OUTPUT_BASE, "FINAL_ANALYSIS_RESULTS_ALL_GENES.RData")
)

cat(sprintf("\n✓ Complete analysis results saved to %s/FINAL_ANALYSIS_RESULTS_ALL_GENES.RData\n", OUTPUT_BASE))
cat(sprintf("✓ Summary report saved to %s/FINAL_ANALYSIS_SUMMARY_ALL_GENES.txt\n", OUTPUT_BASE))

cat("\n", rep("=", 80), "\n")
cat("INTEGRATIVE ANALYSIS (ALL GENES) COMPLETED SUCCESSFULLY!\n")
cat("Analysis finished:", as.character(Sys.time()), "\n")
cat(sprintf("Results available in: %s/\n", OUTPUT_BASE))
cat(sprintf("Plots available in: %s/\n", OUTPUT_PLOTS))
cat(rep("=", 80), "\n")
