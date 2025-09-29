#!/usr/bin/env Rscript
#
# Create Individual Plots from Existing Analysis Results
# This script loads the completed analysis data and generates clean, separate plots
#

library(ggplot2)
library(VennDiagram)
library(dplyr)

setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/integrative_analysis")

# Load the completed analysis results
cat("Loading analysis results...\n")
load("data/integrated_data.RData")

# Recreate the classification (quick version)
rna_results$ensembl_id <- gsub("\\..*", "", rna_results$gene_id)
rna_results$is_significant <- !is.na(rna_results$padj) & rna_results$padj < 0.05

# Load existing analysis if available
if(file.exists("results/FINAL_ANALYSIS_RESULTS.RData")) {
  load("results/FINAL_ANALYSIS_RESULTS.RData")
  cat("Using existing complete analysis results\n")
} else {
  cat("Recreating classification from basic data...\n")
  # Quick recreation of key variables
  library(org.Hs.eg.db)

  convert_entrez_to_ensembl <- function(entrez_ids) {
    ensembl <- mapIds(org.Hs.eg.db,
                     keys = as.character(entrez_ids),
                     column = "ENSEMBL",
                     keytype = "ENTREZID",
                     multiVals = "first")
    return(ensembl)
  }

  tes_peaks$ensembl_id <- convert_entrez_to_ensembl(tes_peaks$geneId)
  tead1_peaks$ensembl_id <- convert_entrez_to_ensembl(tead1_peaks$geneId)

  tes_bound_genes <- unique(tes_peaks$ensembl_id[!is.na(tes_peaks$ensembl_id)])
  tead1_bound_genes <- unique(tead1_peaks$ensembl_id[!is.na(tead1_peaks$ensembl_id)])

  rna_results$tes_bound <- rna_results$ensembl_id %in% tes_bound_genes
  rna_results$tead1_bound <- rna_results$ensembl_id %in% tead1_bound_genes

  tes_direct_targets <- rna_results[rna_results$tes_bound & rna_results$is_significant, ]
  tead1_direct_targets <- rna_results[rna_results$tead1_bound & rna_results$is_significant, ]

  shared_targets <- intersect(tes_direct_targets$ensembl_id, tead1_direct_targets$ensembl_id)
  tes_specific <- tes_direct_targets[!tes_direct_targets$ensembl_id %in% tead1_direct_targets$ensembl_id, ]
  tead1_specific <- tead1_direct_targets[!tead1_direct_targets$ensembl_id %in% tes_direct_targets$ensembl_id, ]

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

  classification$regulatory_class <- "Unbound"
  classification$regulatory_class[classification$tes_bound & classification$is_significant] <- "TES_direct"
  classification$regulatory_class[classification$tead1_bound & classification$is_significant] <- "TEAD1_direct"
  classification$regulatory_class[classification$tes_bound & classification$tead1_bound & classification$is_significant] <- "TES_TEAD1_shared"
  classification$regulatory_class[!classification$tes_bound & !classification$tead1_bound & classification$is_significant] <- "Indirect"
}

cat("Creating individual plots...\n")

# Plot 1: Gene classification pie chart
cat("1. Creating gene classification pie chart...\n")
class_counts <- table(classification$regulatory_class)
pie_data <- data.frame(
  class = names(class_counts),
  count = as.numeric(class_counts),
  percentage = round(100 * as.numeric(class_counts) / sum(class_counts), 1)
)

pdf("plots/01_gene_classification_pie.pdf", width = 10, height = 8)
p1 <- ggplot(pie_data, aes(x = "", y = count, fill = class)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  labs(title = "Gene Classification: Binding vs Expression",
       subtitle = sprintf("Total genes analyzed: %d", nrow(classification)),
       fill = "Regulatory Class") +
  theme_void() +
  theme(plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(size = 14, hjust = 0.5),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold")) +
  scale_fill_brewer(type = "qual", palette = "Set2") +
  geom_text(aes(label = paste0(count, "\n(", percentage, "%)")),
            position = position_stack(vjust = 0.5), size = 4, fontface = "bold")
print(p1)
dev.off()

# Plot 2: Expression changes histograms by regulatory class
cat("2. Creating expression histograms by regulatory class...\n")
sig_genes <- classification[classification$is_significant, ]
pdf("plots/02_expression_histograms.pdf", width = 14, height = 10)
p2 <- ggplot(sig_genes, aes(x = log2FoldChange, fill = regulatory_class)) +
  geom_histogram(bins = 40, alpha = 0.8, color = "white") +
  facet_wrap(~regulatory_class, scales = "free_y", ncol = 2) +
  labs(title = "Expression Changes by Regulatory Class",
       subtitle = "Distribution of Log2 Fold Changes (TES vs GFP)",
       x = "Log2 Fold Change (TES vs GFP)",
       y = "Number of Genes") +
  theme_minimal() +
  theme(plot.title = element_text(size = 18, face = "bold"),
        plot.subtitle = element_text(size = 14),
        strip.text = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 12),
        legend.position = "none") +
  scale_fill_brewer(type = "qual", palette = "Set2") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", alpha = 0.7)
print(p2)
dev.off()

# Plot 3: Boxplot comparison
cat("3. Creating boxplot comparison...\n")
boxplot_data <- sig_genes[sig_genes$regulatory_class %in% c("TES_direct", "TEAD1_direct", "TES_TEAD1_shared", "Indirect"), ]
pdf("plots/03_expression_boxplots.pdf", width = 12, height = 8)
p3 <- ggplot(boxplot_data, aes(x = regulatory_class, y = log2FoldChange, fill = regulatory_class)) +
  geom_boxplot(alpha = 0.8, outlier.alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1) +
  labs(title = "Expression Changes by Regulatory Mechanism",
       subtitle = "Comparison of direct vs indirect regulation",
       x = "Regulatory Class",
       y = "Log2 Fold Change (TES vs GFP)") +
  theme_minimal() +
  theme(plot.title = element_text(size = 18, face = "bold"),
        plot.subtitle = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.title = element_text(size = 12),
        legend.position = "none") +
  scale_fill_brewer(type = "qual", palette = "Set2") +
  stat_summary(fun = median, geom = "point", shape = 23, size = 3, fill = "white")
print(p3)
dev.off()

# Plot 4: Venn diagram
cat("4. Creating Venn diagram...\n")
venn_list <- list(
  "TES Direct" = tes_direct_targets$ensembl_id,
  "TEAD1 Direct" = tead1_direct_targets$ensembl_id,
  "All DE Genes" = rna_results$ensembl_id[rna_results$is_significant]
)

pdf("plots/04_target_overlap_venn.pdf", width = 10, height = 8)
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
cat("5. Creating summary bar chart...\n")
summary_data <- data.frame(
  Category = c("TES Direct", "TEAD1 Direct", "Shared", "TES Specific", "TEAD1 Specific", "Indirect"),
  Count = c(nrow(tes_direct_targets), nrow(tead1_direct_targets), length(shared_targets),
            nrow(tes_specific), nrow(tead1_specific), sum(classification$regulatory_class == "Indirect")),
  Type = c("Direct", "Direct", "Direct", "Specific", "Specific", "Indirect")
)

pdf("plots/05_summary_barchart.pdf", width = 12, height = 8)
p5 <- ggplot(summary_data, aes(x = reorder(Category, Count), y = Count, fill = Type)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  labs(title = "Summary of Target Gene Categories",
       subtitle = "Number of genes in each regulatory class",
       x = "Regulatory Category",
       y = "Number of Genes") +
  theme_minimal() +
  theme(plot.title = element_text(size = 18, face = "bold"),
        plot.subtitle = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.title = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold")) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  geom_text(aes(label = Count), vjust = -0.3, size = 4, fontface = "bold") +
  coord_flip()
print(p5)
dev.off()

cat("\n✓ All individual plots created successfully!\n")
cat("Generated plots:\n")
cat("  - 01_gene_classification_pie.pdf\n")
cat("  - 02_expression_histograms.pdf\n")
cat("  - 03_expression_boxplots.pdf\n")
cat("  - 04_target_overlap_venn.pdf\n")
cat("  - 05_summary_barchart.pdf\n")
cat("\nAll plots saved in: plots/ directory\n")