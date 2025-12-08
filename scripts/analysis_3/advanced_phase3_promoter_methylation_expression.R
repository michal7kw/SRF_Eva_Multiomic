#!/usr/bin/env Rscript

################################################################################
# Phase 3.2: Promoter Methylation and Gene Silencing
#
# Purpose: Three-way integration of binding + methylation + expression
#          Creates BOTH simplified (3-category) and detailed (6-category) analyses
#
# Author: Advanced Multi-Omics Analysis Plan
# Date: 2025-01-24
################################################################################

message("=== Phase 3.2: Promoter Methylation and Gene Silencing ===")
message("Start time: ", Sys.time())
message("NOTE: Creating both SIMPLIFIED (3-cat) and DETAILED (6-cat) analyses")

# Load required libraries
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(pheatmap)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(org.Hs.eg.db)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(gridExtra)  # For grid.arrange
})

# Fix namespace conflicts: dplyr functions get masked by AnnotationDbi
# Note: For TxDb queries, use AnnotationDbi::select() explicitly
select <- dplyr::select
filter <- dplyr::filter
slice <- dplyr::slice
mutate <- dplyr::mutate
arrange <- dplyr::arrange

setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top")

# Input files
BINDING_DATA <- "SRF_Eva_integrated_analysis/scripts/analysis_3/results/01_binding_classification/binding_classification_data.RData"
GENE_EXPR_BINDING <- "SRF_Eva_integrated_analysis/scripts/analysis_3/results/04_category_expression/genes_with_binding_and_expression.csv"
MEDIP_DMRS <- "meDIP/results/07_differential_MEDIPS/TES_vs_GFP_DMRs_FDR05_FC2.csv"

# Output directories - separate for detailed vs simplified
OUTPUT_DIR <- "SRF_Eva_integrated_analysis/scripts/analysis_3/results/08_methylation_expression"
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

################################################################################
# Step 1: Load all data
################################################################################

message("\n[Step 1] Loading data...")

# Load binding classification
load(BINDING_DATA)
message("  Loaded binding peaks: ", nrow(peaks_df))

# Load gene-level binding and expression
gene_data <- read.csv(GENE_EXPR_BINDING)
message("  Loaded genes: ", nrow(gene_data))

# Load meDIP DMRs
if (!file.exists(MEDIP_DMRS)) {
  stop("meDIP DMRs file not found: ", MEDIP_DMRS)
}
dmrs <- read.csv(MEDIP_DMRS)
message("  Loaded DMRs: ", nrow(dmrs))

# Convert DMRs to GRanges
# Note: Column names vary between DMR files - handle both formats
if ("edgeR.logFC" %in% colnames(dmrs)) {
  logfc_col <- dmrs$edgeR.logFC
  fdr_col <- dmrs$edgeR.adj.p.val
} else {
  logfc_col <- dmrs$logFC
  fdr_col <- dmrs$FDR
}

# Handle mean signal columns (may be named differently)
if ("mean_TES" %in% colnames(dmrs)) {
  mean_tes_col <- dmrs$mean_TES
  mean_gfp_col <- dmrs$mean_GFP
} else if ("mean_group1" %in% colnames(dmrs)) {
  mean_tes_col <- dmrs$mean_group1
  mean_gfp_col <- dmrs$mean_group2
} else {
  mean_tes_col <- rep(NA, nrow(dmrs))
  mean_gfp_col <- rep(NA, nrow(dmrs))
}

dmrs_gr <- GRanges(
  seqnames = dmrs$chr,
  ranges = IRanges(start = dmrs$start, end = dmrs$stop),
  logFC = logfc_col,
  FDR = fdr_col,
  mean_TES = mean_tes_col,
  mean_GFP = mean_gfp_col
)

# Classify DMRs
dmrs_gr$methylation_change <- ifelse(dmrs_gr$logFC > 0, "Hypermethylated", "Hypomethylated")

################################################################################
# Step 2: Get promoter regions for all genes
################################################################################

message("\n[Step 2] Defining promoter regions...")

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Get promoters (±3kb from TSS)
promoters <- promoters(txdb, upstream = 3000, downstream = 3000)

# Get gene symbols using the correct approach:
# 1. First get transcript-to-gene mapping from TxDb
# 2. Then map Entrez IDs to symbols
message("  Mapping transcripts to genes...")

# Get transcript to gene mapping (use AnnotationDbi::select for TxDb objects)
tx_to_gene <- AnnotationDbi::select(txdb, keys = promoters$tx_name,
                                    columns = c("TXNAME", "GENEID"),
                                    keytype = "TXNAME")

# Map Entrez IDs (GENEID) to gene symbols
entrez_ids <- tx_to_gene$GENEID[match(promoters$tx_name, tx_to_gene$TXNAME)]

gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = entrez_ids[!is.na(entrez_ids)],
                       column = "SYMBOL",
                       keytype = "ENTREZID",
                       multiVals = "first")

# Create a named vector for easy lookup
symbol_lookup <- setNames(gene_symbols, names(gene_symbols))
promoters$gene_symbol <- symbol_lookup[entrez_ids]

# Remove NA gene symbols
promoters <- promoters[!is.na(promoters$gene_symbol)]

message("  Defined promoters for ", length(unique(promoters$gene_symbol)), " genes")

################################################################################
# Step 3: Find DMRs at gene promoters
################################################################################

message("\n[Step 3] Finding DMRs at gene promoters...")

# Overlap DMRs with promoters
dmr_promoter_overlaps <- findOverlaps(promoters, dmrs_gr)

# Create promoter-DMR mapping
promoter_dmr_map <- data.frame(
  gene_symbol = promoters$gene_symbol[queryHits(dmr_promoter_overlaps)],
  dmr_chr = as.character(seqnames(dmrs_gr[subjectHits(dmr_promoter_overlaps)])),
  dmr_start = start(dmrs_gr[subjectHits(dmr_promoter_overlaps)]),
  dmr_end = end(dmrs_gr[subjectHits(dmr_promoter_overlaps)]),
  dmr_logFC = dmrs_gr$logFC[subjectHits(dmr_promoter_overlaps)],
  dmr_FDR = dmrs_gr$FDR[subjectHits(dmr_promoter_overlaps)],
  methylation_change = dmrs_gr$methylation_change[subjectHits(dmr_promoter_overlaps)]
)

message("  Found ", nrow(promoter_dmr_map), " promoter-DMR associations")

# Summarize at gene level (if multiple DMRs, keep one with strongest effect)
gene_promoter_meth <- promoter_dmr_map %>%
  group_by(gene_symbol) %>%
  arrange(desc(abs(dmr_logFC))) %>%
  slice(1) %>%
  ungroup() %>%
  select(gene_symbol, dmr_logFC, dmr_FDR, methylation_change)

message("  Genes with promoter DMRs: ", nrow(gene_promoter_meth))

################################################################################
# Step 4: Three-way integration
################################################################################

message("\n[Step 4] Three-way integration: Binding + Methylation + Expression...")

# Merge all data
integrated_data <- gene_data %>%
  left_join(gene_promoter_meth, by = c("gene_name" = "gene_symbol"))

# Add methylation status
integrated_data$has_promoter_dmr <- !is.na(integrated_data$dmr_logFC)
integrated_data$methylation_status <- ifelse(
  is.na(integrated_data$methylation_change), "No_Change",
  as.character(integrated_data$methylation_change)
)

# Classification into mechanistic categories
integrated_data <- integrated_data %>%
  mutate(
    # Determine if bound by TES or TEAD1
    tes_bound = primary_category %in% c("TES_unique", "Shared_high", "Shared_TES_dominant", "Shared_equivalent"),
    tead1_bound = primary_category %in% c("TEAD1_unique", "Shared_high", "Shared_TEAD1_dominant", "Shared_equivalent"),

    # Mechanistic classification
    regulatory_mechanism = case_when(
      # Category A: Direct epigenetic silencing
      tes_bound & methylation_status == "Hypermethylated" & is_DEG & log2FoldChange < 0 ~
        "Direct_Epigenetic_Silencing",

      # Category B: TES binding + methylation but no DE (buffered)
      tes_bound & methylation_status == "Hypermethylated" & !is_DEG ~
        "Methylated_No_Effect",

      # Category C: TES binding + DE but no methylation (KRAB-mediated?)
      tes_bound & methylation_status == "No_Change" & is_DEG & log2FoldChange < 0 ~
        "Non_Methylation_Repression",

      # Category D: Methylation + DE but not bound (indirect)
      !tes_bound & methylation_status == "Hypermethylated" & is_DEG & log2FoldChange < 0 ~
        "Indirect_Methylation",

      # TES-bound and activated
      tes_bound & is_DEG & log2FoldChange > 0 ~
        "TES_Activation",

      # TEAD1-specific effects
      tead1_bound & !tes_bound & is_DEG ~
        "TEAD1_Specific",

      # Other bound genes
      (tes_bound | tead1_bound) & is_DEG ~
        "Other_Direct",

      # Indirect effects
      is_DEG & primary_category == "Unbound" ~
        "Indirect",

      # Not DE
      TRUE ~ "Not_DE"
    )
  )

message("  Integrated data for ", nrow(integrated_data), " genes")

# Count genes in each category
mechanism_counts <- table(integrated_data$regulatory_mechanism)
print(mechanism_counts)

################################################################################
# Step 5: Correlation analysis
################################################################################

message("\n[Step 5] Analyzing promoter methylation-expression correlation...")

# Filter for genes with promoter DMRs
genes_with_dmr <- integrated_data %>%
  filter(has_promoter_dmr & !is.na(log2FoldChange))

message("  Genes with both promoter DMR and expression data: ", nrow(genes_with_dmr))

# Correlation test
if (nrow(genes_with_dmr) >= 10) {
  cor_test <- cor.test(
    genes_with_dmr$dmr_logFC,
    genes_with_dmr$log2FoldChange,
    method = "spearman"
  )

  message("  Spearman correlation: rho = ", round(cor_test$estimate, 3),
          ", p = ", format(cor_test$p.value, scientific = TRUE))
} else {
  cor_test <- NULL
  message("  WARNING: Not enough genes for correlation analysis")
}

# Separate analysis by binding status
genes_tes_bound_dmr <- genes_with_dmr %>% filter(tes_bound)
genes_unbound_dmr <- genes_with_dmr %>% filter(!tes_bound & !tead1_bound)

if (nrow(genes_tes_bound_dmr) >= 10) {
  cor_test_tes <- cor.test(
    genes_tes_bound_dmr$dmr_logFC,
    genes_tes_bound_dmr$log2FoldChange,
    method = "spearman"
  )
} else {
  cor_test_tes <- NULL
}

if (nrow(genes_unbound_dmr) >= 10) {
  cor_test_unbound <- cor.test(
    genes_unbound_dmr$dmr_logFC,
    genes_unbound_dmr$log2FoldChange,
    method = "spearman"
  )
} else {
  cor_test_unbound <- NULL
}

################################################################################
# Step 6: Export results
################################################################################

message("\n[Step 6] Exporting results...")

write.csv(integrated_data,
          file.path(OUTPUT_DIR, "integrated_binding_methylation_expression.csv"),
          row.names = FALSE)

write.csv(promoter_dmr_map,
          file.path(OUTPUT_DIR, "promoter_dmr_associations.csv"),
          row.names = FALSE)

# Summary statistics
mechanism_stats <- integrated_data %>%
  group_by(regulatory_mechanism) %>%
  summarise(
    n_genes = n(),
    n_with_dmr = sum(has_promoter_dmr),
    mean_log2FC = mean(log2FoldChange, na.rm = TRUE),
    mean_meth_FC = mean(dmr_logFC, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(mechanism_stats,
          file.path(OUTPUT_DIR, "mechanism_summary_statistics.csv"),
          row.names = FALSE)

################################################################################
# Step 7: Visualizations
################################################################################

message("\n[Step 7] Generating visualizations...")

# 7.1: Methylation-Expression correlation
pdf(file.path(OUTPUT_DIR, "methylation_expression_correlation.pdf"), width = 12, height = 10)

p1 <- ggplot(genes_with_dmr, aes(x = dmr_logFC, y = log2FoldChange)) +
  geom_point(aes(color = tes_bound), alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "darkred") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue"),
                     labels = c("TRUE" = "TES-bound", "FALSE" = "Not bound")) +
  theme_classic() +
  labs(
    title = "Promoter Methylation vs Gene Expression",
    subtitle = ifelse(!is.null(cor_test),
                     paste0("Spearman rho = ", round(cor_test$estimate, 3),
                           ", p = ", format(cor_test$p.value, digits = 3)),
                     ""),
    x = "Promoter Methylation log2FC (TES vs GFP)",
    y = "Gene Expression log2FC (TES vs GFP)",
    color = "Binding Status"
  )

# Separate plots for bound vs unbound
if (!is.null(cor_test_tes) && nrow(genes_tes_bound_dmr) >= 10) {
  p2 <- ggplot(genes_tes_bound_dmr, aes(x = dmr_logFC, y = log2FoldChange)) +
    geom_point(color = "red", alpha = 0.6, size = 2) +
    geom_smooth(method = "lm", se = TRUE, color = "darkred") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    theme_classic() +
    labs(
      title = "TES-Bound Genes Only",
      subtitle = paste0("Spearman rho = ", round(cor_test_tes$estimate, 3),
                       ", p = ", format(cor_test_tes$p.value, digits = 3)),
      x = "Promoter Methylation log2FC",
      y = "Gene Expression log2FC"
    )
} else {
  p2 <- ggplot() + theme_void() + ggtitle("Insufficient data for TES-bound genes")
}

gridExtra::grid.arrange(p1, p2, ncol = 1)
dev.off()

# 7.2: Three-way integration heatmap
# Select top genes from each mechanistic category
top_genes_per_mechanism <- integrated_data %>%
  filter(regulatory_mechanism %in% c("Direct_Epigenetic_Silencing",
                                     "Non_Methylation_Repression",
                                     "Indirect_Methylation",
                                     "TES_Activation",
                                     "TEAD1_Specific")) %>%
  group_by(regulatory_mechanism) %>%
  arrange(padj) %>%
  slice_head(n = 10) %>%
  ungroup()

if (nrow(top_genes_per_mechanism) > 0) {
  # Create matrix for heatmap
  heatmap_data <- top_genes_per_mechanism %>%
    select(gene_name, log2FoldChange, dmr_logFC, max_tes_signal, max_tead1_signal) %>%
    as.data.frame()

  rownames(heatmap_data) <- heatmap_data$gene_name
  heatmap_matrix <- as.matrix(heatmap_data[, -1])

  # Scale columns
  heatmap_matrix_scaled <- scale(heatmap_matrix)

  # Annotation
  ha <- rowAnnotation(
    Mechanism = top_genes_per_mechanism$regulatory_mechanism,
    col = list(Mechanism = c(
      "Direct_Epigenetic_Silencing" = "#D73027",
      "Non_Methylation_Repression" = "#FC8D59",
      "Indirect_Methylation" = "#FEE090",
      "TES_Activation" = "#91BFDB",
      "TEAD1_Specific" = "#4575B4"
    ))
  )

  pdf(file.path(OUTPUT_DIR, "three_way_integration_heatmap.pdf"), width = 10, height = 12)
  Heatmap(
    heatmap_matrix_scaled,
    name = "Z-score",
    column_labels = c("Expression\nlog2FC", "Methylation\nlog2FC",
                     "TES\nSignal", "TEAD1\nSignal"),
    row_names_gp = gpar(fontsize = 8),
    right_annotation = ha,
    col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_row_dend = TRUE,
    column_title = "Multi-Omics Profile of Top Genes by Mechanism",
    heatmap_legend_param = list(title = "Z-score")
  )
  dev.off()
}

# 7.3: Mechanistic classification pie chart
pdf(file.path(OUTPUT_DIR, "mechanistic_classification_pie.pdf"), width = 10, height = 8)
mechanism_df <- as.data.frame(mechanism_counts)
colnames(mechanism_df) <- c("Mechanism", "Count")

# Filter for key mechanisms
mechanism_df <- mechanism_df %>%
  filter(Mechanism != "Not_DE" | Count > 100) %>%
  mutate(Percentage = round(Count / sum(Count) * 100, 1))

ggplot(mechanism_df, aes(x = "", y = Count, fill = Mechanism)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  geom_text(aes(label = paste0(Mechanism, "\n", Count, " (", Percentage, "%)")),
            position = position_stack(vjust = 0.5), size = 3) +
  scale_fill_brewer(palette = "Set3") +
  theme_void() +
  theme(legend.position = "bottom") +
  labs(title = "Mechanistic Classification of Genes")
dev.off()

# 7.4: Comparison of mechanisms
pdf(file.path(OUTPUT_DIR, "mechanism_comparison_boxplots.pdf"), width = 14, height = 10)
key_mechanisms <- integrated_data %>%
  filter(regulatory_mechanism %in% c("Direct_Epigenetic_Silencing",
                                     "Methylated_No_Effect",
                                     "Non_Methylation_Repression",
                                     "Indirect_Methylation",
                                     "TES_Activation",
                                     "TEAD1_Specific"))

p1 <- ggplot(key_mechanisms, aes(x = regulatory_mechanism, y = log2FoldChange, fill = regulatory_mechanism)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.2, alpha = 0.5, outlier.shape = NA) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_brewer(palette = "Set2") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  labs(title = "Expression Changes", x = "", y = "log2 Fold Change")

p2 <- ggplot(key_mechanisms %>% filter(has_promoter_dmr),
             aes(x = regulatory_mechanism, y = dmr_logFC, fill = regulatory_mechanism)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.2, alpha = 0.5, outlier.shape = NA) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_brewer(palette = "Set2") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  labs(title = "Promoter Methylation Changes", x = "Mechanism", y = "Methylation log2FC")

gridExtra::grid.arrange(p1, p2, ncol = 1)
dev.off()

################################################################################
# Step 8: Summary report
################################################################################

message("\n[Step 8] Creating summary report...")

sink(file.path(OUTPUT_DIR, "PHASE3_2_SUMMARY.txt"))
cat("=== Phase 3.2: Promoter Methylation and Gene Silencing ===\n")
cat("Date:", as.character(Sys.time()), "\n\n")

cat("Genes with Promoter DMRs:", nrow(gene_promoter_meth), "\n")
cat("Hypermethylated promoters:", sum(gene_promoter_meth$methylation_change == "Hypermethylated"), "\n")
cat("Hypomethylated promoters:", sum(gene_promoter_meth$methylation_change == "Hypomethylated"), "\n\n")

cat("Mechanistic Classification:\n")
cat("===========================\n")
print(mechanism_counts)
cat("\n")

cat("Mechanism Statistics:\n")
cat("=====================\n")
print(mechanism_stats)
cat("\n")

if (!is.null(cor_test)) {
  cat("Overall Methylation-Expression Correlation:\n")
  cat("============================================\n")
  print(cor_test)
  cat("\n")
}

if (!is.null(cor_test_tes)) {
  cat("TES-Bound Genes Methylation-Expression Correlation:\n")
  cat("====================================================\n")
  print(cor_test_tes)
  cat("\n")
}

sink()

message("\n=== Analysis Complete ===")
message("Output directory: ", OUTPUT_DIR)
message("End time: ", Sys.time())
