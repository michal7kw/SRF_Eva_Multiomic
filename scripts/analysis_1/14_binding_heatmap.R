#!/usr/bin/env Rscript

# Cut&Tag Binding Heatmap for TES Target Genes
# Creates a clustered heatmap showing mean promoter binding signal
# Similar to RNA-seq expression heatmap format

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(ComplexHeatmap)
  library(circlize)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(dplyr)
  library(grid)
})

# ===================== Configuration =====================
BASE_DIR <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top"
SCRIPT_DIR <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_1"
setwd(SCRIPT_DIR)
BIGWIG_DIR <- file.path(BASE_DIR, "SRF_Eva_CUTandTAG/results/06_bigwig")
OUTPUT_DIR <- "output/14_binding_heatmap"
DEG_FILE <- file.path(BASE_DIR, "SRF_Eva_integrated_analysis/data/TES_degs.txt")
DESEQ_FILE <- file.path(BASE_DIR, "SRF_Eva_RNA/results/05_deseq2/deseq2_results_TES_vs_GFP.txt")

# Create output directory
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Promoter region definition
PROMOTER_UPSTREAM <- 3000
PROMOTER_DOWNSTREAM <- 3000

cat("Starting Cut&Tag binding heatmap analysis...\n")
cat("Output directory:", OUTPUT_DIR, "\n\n")

# ===================== Load Data =====================

cat("Loading TES DEGs list...\n")
tes_degs <- read.table(DEG_FILE, header = FALSE, stringsAsFactors = FALSE)$V1
cat("Found", length(tes_degs), "DEGs\n")

cat("\nLoading DESeq2 results for expression data...\n")
deseq_results <- read.table(DESEQ_FILE, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Extract gene symbols
if ("gene_symbol" %in% colnames(deseq_results)) {
  deseq_results$symbol <- deseq_results$gene_symbol
} else if (grepl("ENSG", deseq_results[1, 1])) {
  ensembl_ids <- gsub("\\..*", "", deseq_results[, 1])
  deseq_results$symbol <- mapIds(org.Hs.eg.db,
    keys = ensembl_ids,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
}

# Filter for TES DEGs
deseq_degs <- deseq_results[deseq_results$symbol %in% tes_degs &
  !is.na(deseq_results$padj) &
  deseq_results$padj < 0.05, ]
cat("Matched", nrow(deseq_degs), "DEGs with DESeq2 results\n")

# Classify as UP or DOWN regulated
deseq_degs$regulation <- ifelse(deseq_degs$log2FoldChange > 0, "UP", "DOWN")
cat("UP-regulated:", sum(deseq_degs$regulation == "UP"), "\n")
cat("DOWN-regulated:", sum(deseq_degs$regulation == "DOWN"), "\n")

# ===================== Get Gene Coordinates =====================

cat("\nRetrieving gene coordinates from TxDb...\n")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
genes <- genes(txdb)

# Map gene symbols to Entrez IDs
symbol_to_entrez <- mapIds(org.Hs.eg.db,
  keys = deseq_degs$symbol,
  column = "ENTREZID",
  keytype = "SYMBOL",
  multiVals = "first"
)

deg_entrez <- na.omit(symbol_to_entrez)
target_genes <- genes[genes$gene_id %in% deg_entrez]

# Add gene symbols
target_genes$symbol <- mapIds(org.Hs.eg.db,
  keys = target_genes$gene_id,
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)

# Add regulation direction and log2FC
target_genes$regulation <- deseq_degs$regulation[match(target_genes$symbol, deseq_degs$symbol)]
target_genes$log2FoldChange <- deseq_degs$log2FoldChange[match(target_genes$symbol, deseq_degs$symbol)]

cat("Found coordinates for", length(target_genes), "genes\n")

# Define promoter regions
promoters_gr <- promoters(target_genes,
  upstream = PROMOTER_UPSTREAM,
  downstream = PROMOTER_DOWNSTREAM
)
mcols(promoters_gr) <- mcols(target_genes)

cat("Promoter regions created:", length(promoters_gr), "\n")

# ===================== Load BigWig Files and Calculate Mean Signal =====================

cat("\nLoading BigWig files and calculating mean promoter signal...\n")

# Define samples to use
samples <- list(
  TES_1 = file.path(BIGWIG_DIR, "TES-1_CPM.bw"),
  TES_2 = file.path(BIGWIG_DIR, "TES-2_CPM.bw"),
  TES_3 = file.path(BIGWIG_DIR, "TES-3_CPM.bw"),
  TEAD1_1 = file.path(BIGWIG_DIR, "TEAD1-1_CPM.bw"),
  TEAD1_2 = file.path(BIGWIG_DIR, "TEAD1-2_CPM.bw"),
  TEAD1_3 = file.path(BIGWIG_DIR, "TEAD1-3_CPM.bw")
)

# Check which files exist
samples <- samples[file.exists(unlist(samples))]
cat("Found", length(samples), "sample files\n")

# Function to calculate mean signal in promoter regions
calculate_promoter_signal <- function(bw_file, regions) {
  cat("  Processing:", basename(bw_file), "\n")

  # Import BigWig
  coverage <- import(bw_file, format = "BigWig")

  # Fix chromosome naming if needed
  if (length(intersect(seqlevels(coverage), seqlevels(regions))) == 0) {
    if (any(grepl("^chr", seqlevels(regions))) && !any(grepl("^chr", seqlevels(coverage)))) {
      seqlevels(coverage) <- paste0("chr", seqlevels(coverage))
    }
  }

  # Subset to common chromosomes
  common_chrs <- intersect(seqlevels(coverage), seqlevels(regions))
  coverage <- keepSeqlevels(coverage, common_chrs, pruning.mode = "coarse")
  regions <- keepSeqlevels(regions, common_chrs, pruning.mode = "coarse")

  # Find overlaps and calculate mean signal
  overlaps <- findOverlaps(regions, coverage)

  # Calculate weighted mean signal for each promoter
  mean_signals <- sapply(seq_along(regions), function(i) {
    hits_idx <- which(queryHits(overlaps) == i)
    if (length(hits_idx) == 0) {
      return(0)
    }

    # Get overlapping coverage regions
    subject_idx <- subjectHits(overlaps)[hits_idx]
    cov_regions <- coverage[subject_idx]

    # Calculate overlap widths for each pair
    region_rep <- rep(regions[i], length(cov_regions))
    overlap_ranges <- pintersect(region_rep, cov_regions)
    overlap_widths <- width(overlap_ranges)

    # Weighted mean by overlap width
    weighted.mean(cov_regions$score, overlap_widths)
  })

  return(mean_signals)
}

# Calculate signals for all samples
signal_matrix <- sapply(samples, function(bw) {
  calculate_promoter_signal(bw, promoters_gr)
})

# Add gene symbols as rownames
rownames(signal_matrix) <- promoters_gr$symbol

cat("\nSignal matrix created:\n")
cat("  Dimensions:", nrow(signal_matrix), "genes x", ncol(signal_matrix), "samples\n")
cat("  Value range:", min(signal_matrix), "to", max(signal_matrix), "\n")
cat("  Mean signal:", mean(signal_matrix), "\n")
cat("  Median signal:", median(signal_matrix), "\n")
cat("  Genes with zero signal in all samples:", sum(rowSums(signal_matrix) == 0), "\n")
cat("  Sample correlations (Pearson):\n")
cor_matrix <- cor(signal_matrix, method = "pearson")
print(round(cor_matrix, 3))

# ===================== Prepare Data for Heatmap =====================

cat("\nPreparing data for heatmap...\n")

# Check for genes with zero variance (will cause issues with z-score)
row_vars <- apply(signal_matrix, 1, var)
zero_var_genes <- sum(row_vars == 0 | is.na(row_vars))
if (zero_var_genes > 0) {
  cat("WARNING:", zero_var_genes, "genes have zero variance across samples\n")
  cat("These genes will be set to zero after scaling\n")
}

# Z-score normalization (center and scale by row)
# This shows relative binding strength: each gene's pattern across samples
signal_matrix_scaled <- t(scale(t(signal_matrix)))

cat("After scaling:\n")
cat(
  "  Range:", min(signal_matrix_scaled, na.rm = TRUE), "to",
  max(signal_matrix_scaled, na.rm = TRUE), "\n"
)
cat("  Mean:", mean(signal_matrix_scaled, na.rm = TRUE), "\n")
cat("  SD:", sd(signal_matrix_scaled, na.rm = TRUE), "\n")

# Handle any NA or Inf values (genes with no variance)
# These occur when a gene has identical signal across all samples
signal_matrix_scaled[is.na(signal_matrix_scaled)] <- 0
signal_matrix_scaled[is.infinite(signal_matrix_scaled)] <- 0

cat("  Genes with zero variance:", sum(row_vars == 0, na.rm = TRUE), "\n")

# Create annotation data frame
annotation_df <- data.frame(
  gene = promoters_gr$symbol,
  regulation = promoters_gr$regulation,
  log2FC = promoters_gr$log2FoldChange,
  stringsAsFactors = FALSE
)

# ===================== Create Heatmap =====================

cat("\nGenerating heatmap...\n")

# Color scheme: blue (low) -> white (mean) -> red (high)
col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

# Column annotation (condition grouping)
sample_condition <- factor(c("TES", "TES", "TES", "TEAD1", "TEAD1", "TEAD1"),
  levels = c("TES", "TEAD1")
)
sample_number <- factor(c("1", "2", "3", "1", "2", "3"))

ha_column <- HeatmapAnnotation(
  condition = sample_condition,
  samples = sample_number,
  col = list(
    condition = c("TES" = "#FFA500", "TEAD1" = "#FFD700"),
    samples = c("1" = "#4D4D4D", "2" = "#808080", "3" = "#B3B3B3")
  ),
  annotation_name_side = "left",
  show_annotation_name = TRUE,
  annotation_legend_param = list(
    condition = list(title = "condition"),
    samples = list(title = "samples")
  )
)

# Row annotation (DEG status)
ha_row <- rowAnnotation(
  DEGs = annotation_df$regulation,
  col = list(DEGs = c("DOWN" = "blue", "UP" = "red")),
  show_annotation_name = TRUE,
  annotation_legend_param = list(
    DEGs = list(title = "DEGs")
  ),
  width = unit(5, "mm")
)

# Create main heatmap
ht <- Heatmap(
  signal_matrix_scaled,
  name = "Z-score\nBinding Signal",
  col = col_fun,

  # Clustering
  cluster_rows = TRUE,
  cluster_columns = FALSE, # Keep sample order
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "complete",

  # Row (gene) settings
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 6),
  row_names_side = "right",
  row_dend_width = unit(15, "mm"),

  # Column (sample) settings
  show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 10),
  column_names_rot = 0,
  column_names_centered = TRUE,

  # Split by DEG status
  row_split = annotation_df$regulation,
  row_title_rot = 0,
  row_title_gp = gpar(fontsize = 10),
  row_gap = unit(2, "mm"),

  # Annotations
  top_annotation = ha_column,

  # Display options
  use_raster = FALSE,
  heatmap_legend_param = list(
    title = "Z-score\nBinding",
    at = c(-2, -1, 0, 1, 2),
    labels = c("-2", "-1", "0", "1", "2"),
    legend_height = unit(4, "cm")
  ),

  # Size
  width = unit(8, "cm"),
  height = unit(20, "cm")
)

# Save heatmap
pdf(file.path(OUTPUT_DIR, "promoter_binding_heatmap.pdf"),
  width = 10, height = 14
)
draw(ha_row + ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

cat("  Saved:", file.path(OUTPUT_DIR, "promoter_binding_heatmap.pdf"), "\n")

# Also create a version with log2FC annotation
ha_row_extended <- rowAnnotation(
  DEGs = annotation_df$regulation,
  log2FC = anno_barplot(
    annotation_df$log2FC,
    gp = gpar(fill = ifelse(annotation_df$regulation == "UP", "red", "blue")),
    width = unit(2, "cm")
  ),
  col = list(DEGs = c("DOWN" = "blue", "UP" = "red")),
  show_annotation_name = TRUE,
  annotation_legend_param = list(
    DEGs = list(title = "DEGs")
  )
)

pdf(file.path(OUTPUT_DIR, "promoter_binding_heatmap_with_fc.pdf"),
  width = 12, height = 14
)
draw(ha_row_extended + ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

cat("  Saved:", file.path(OUTPUT_DIR, "promoter_binding_heatmap_with_fc.pdf"), "\n")

# ===================== Save Data =====================

cat("\nSaving binding signal matrix...\n")

# Create output data frame with raw signals
output_df_raw <- data.frame(
  gene = rownames(signal_matrix),
  regulation = annotation_df$regulation,
  log2FoldChange = annotation_df$log2FC,
  signal_matrix,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Create data frame with scaled signals
output_df_scaled <- data.frame(
  gene = rownames(signal_matrix_scaled),
  regulation = annotation_df$regulation,
  log2FoldChange = annotation_df$log2FC,
  signal_matrix_scaled,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Rename scaled columns to avoid confusion
colnames(output_df_scaled)[4:ncol(output_df_scaled)] <-
  paste0(colnames(output_df_scaled)[4:ncol(output_df_scaled)], "_zscore")

# Combine both raw and scaled
output_df <- cbind(output_df_raw, output_df_scaled[, 4:ncol(output_df_scaled)])

write.table(output_df,
  file.path(OUTPUT_DIR, "promoter_binding_signals.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cat("  Saved:", file.path(OUTPUT_DIR, "promoter_binding_signals.txt"), "\n")

# ===================== Summary Statistics =====================

cat("\n=== Analysis Summary ===\n")
cat("Total genes analyzed:", nrow(signal_matrix), "\n")
cat("UP-regulated:", sum(annotation_df$regulation == "UP"), "\n")
cat("DOWN-regulated:", sum(annotation_df$regulation == "DOWN"), "\n")

cat("\n--- Mean promoter binding signal (CPM) ---\n")
for (i in 1:ncol(signal_matrix)) {
  cat(sprintf(
    "  %s: %.2f (SD: %.2f)\n",
    colnames(signal_matrix)[i],
    mean(signal_matrix[, i]),
    sd(signal_matrix[, i])
  ))
}

cat("\n--- Replicate concordance ---\n")
tes_cor <- cor(signal_matrix[, 1:3])
tead1_cor <- cor(signal_matrix[, 4:6])
cat(sprintf("  TES replicates (mean r): %.3f\n", mean(tes_cor[upper.tri(tes_cor)])))
cat(sprintf("  TEAD1 replicates (mean r): %.3f\n", mean(tead1_cor[upper.tri(tead1_cor)])))

cat("\n--- Top 10 genes with highest TES binding ---\n")
tes_mean <- rowMeans(signal_matrix[, 1:3])
top_tes <- head(sort(tes_mean, decreasing = TRUE), 10)
for (i in 1:length(top_tes)) {
  gene <- names(top_tes)[i]
  reg <- annotation_df$regulation[annotation_df$gene == gene]
  fc <- annotation_df$log2FC[annotation_df$gene == gene]
  cat(sprintf(
    "  %d. %s (%s, log2FC=%.2f): %.2f CPM\n",
    i, gene, reg, fc, top_tes[i]
  ))
}

cat("\n--- Top 10 genes with highest TEAD1 binding ---\n")
tead1_mean <- rowMeans(signal_matrix[, 4:6])
top_tead1 <- head(sort(tead1_mean, decreasing = TRUE), 10)
for (i in 1:length(top_tead1)) {
  gene <- names(top_tead1)[i]
  reg <- annotation_df$regulation[annotation_df$gene == gene]
  fc <- annotation_df$log2FC[annotation_df$gene == gene]
  cat(sprintf(
    "  %d. %s (%s, log2FC=%.2f): %.2f CPM\n",
    i, gene, reg, fc, top_tead1[i]
  ))
}

cat("\n--- Binding vs Expression correlation ---\n")
tes_mean_all <- rowMeans(signal_matrix[, 1:3])
tead1_mean_all <- rowMeans(signal_matrix[, 4:6])
tes_expr_cor <- cor(tes_mean_all, annotation_df$log2FC, method = "spearman")
tead1_expr_cor <- cor(tead1_mean_all, annotation_df$log2FC, method = "spearman")
cat(sprintf("  TES binding vs log2FC (Spearman rho): %.3f\n", tes_expr_cor))
cat(sprintf("  TEAD1 binding vs log2FC (Spearman rho): %.3f\n", tead1_expr_cor))

cat("\n--- Differential binding between conditions ---\n")
up_genes <- annotation_df$regulation == "UP"
down_genes <- annotation_df$regulation == "DOWN"
cat(sprintf("  Mean TES binding in UP genes: %.2f CPM\n", mean(tes_mean_all[up_genes])))
cat(sprintf("  Mean TES binding in DOWN genes: %.2f CPM\n", mean(tes_mean_all[down_genes])))
cat(sprintf("  Mean TEAD1 binding in UP genes: %.2f CPM\n", mean(tead1_mean_all[up_genes])))
cat(sprintf("  Mean TEAD1 binding in DOWN genes: %.2f CPM\n", mean(tead1_mean_all[down_genes])))

cat("\n==================================================\n")
cat("Analysis complete!\n")
cat("Results saved to:", OUTPUT_DIR, "\n")
cat("==================================================\n")
