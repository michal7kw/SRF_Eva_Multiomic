#!/usr/bin/env Rscript

# Cut&Tag Density Plot Analysis for TES Target Genes
# This script generates heatmaps and profile plots showing Cut&Tag signal
# enrichment around TES differentially expressed genes

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(EnrichedHeatmap)
  library(ComplexHeatmap)
  library(circlize)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(DESeq2)
})

# ===================== Configuration =====================
BASE_DIR <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top"
SCRIPT_DIR <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_1"
setwd(SCRIPT_DIR)
BIGWIG_DIR <- file.path(BASE_DIR, "SRF_Eva_CUTandTAG/results/06_bigwig")
OUTPUT_DIR <- "output/15_cutandtag_density_plot"
DEG_FILE <- file.path(BASE_DIR, "SRF_Eva_integrated_analysis/data/TES_degs.txt")
DESEQ_FILE <- file.path(BASE_DIR, "SRF_Eva_RNA/results/05_deseq2/deseq2_results_TES_vs_GFP.txt")

# Create output directory
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Analysis parameters
UPSTREAM <- 3000
DOWNSTREAM <- 3000
WINDOW_SIZE <- 50

cat("Starting Cut&Tag density plot analysis...\n")
cat("Output directory:", OUTPUT_DIR, "\n\n")

# ===================== Load Data =====================

cat("Loading TES DEGs list...\n")
tes_degs <- read.table(DEG_FILE, header = FALSE, stringsAsFactors = FALSE)$V1
cat("Found", length(tes_degs), "DEGs\n")

cat("\nLoading DESeq2 results for expression data...\n")
deseq_results <- read.table(DESEQ_FILE, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# Extract gene symbols from the rownames or ID column
if ("gene_symbol" %in% colnames(deseq_results)) {
  deseq_results$symbol <- deseq_results$gene_symbol
} else if (grepl("ENSG", deseq_results[1, 1])) {
  # Map Ensembl IDs to gene symbols
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

# Get all genes
genes <- genes(txdb)

cat("Total genes in TxDb:", length(genes), "\n")
cat("Chromosome names in TxDb:", paste(head(seqlevels(genes), 10), collapse = ", "), "...\n")

# Map gene symbols to Entrez IDs
symbol_to_entrez <- mapIds(org.Hs.eg.db,
  keys = deseq_degs$symbol,
  column = "ENTREZID",
  keytype = "SYMBOL",
  multiVals = "first"
)

# Filter genes for our DEGs
deg_entrez <- na.omit(symbol_to_entrez)
cat("Mapped", length(deg_entrez), "DEGs to Entrez IDs\n")

target_genes <- genes[genes$gene_id %in% deg_entrez]

# Add gene symbols
target_genes$symbol <- mapIds(org.Hs.eg.db,
  keys = target_genes$gene_id,
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)

# Add regulation direction
target_genes$regulation <- deseq_degs$regulation[match(target_genes$symbol, deseq_degs$symbol)]

cat("Found coordinates for", length(target_genes), "genes\n")
cat("Chromosome distribution:", table(seqnames(target_genes)), "\n")

# Define TSS regions
tss <- promoters(target_genes, upstream = 0, downstream = 1)
mcols(tss) <- mcols(target_genes)

cat("TSS regions created:", length(tss), "\n")
cat("Example TSS regions:\n")
print(head(tss, 3))

# ===================== Load BigWig Files =====================

cat("\nLoading BigWig files...\n")
bigwig_files <- list(
  TES_rep1 = file.path(BIGWIG_DIR, "TES-1_CPM.bw"),
  TES_rep2 = file.path(BIGWIG_DIR, "TES-2_CPM.bw"),
  TES_rep3 = file.path(BIGWIG_DIR, "TES-3_CPM.bw"),
  TES_combined = file.path(BIGWIG_DIR, "TES_comb.bw"),
  TEAD1_combined = file.path(BIGWIG_DIR, "TEAD1_comb.bw")
)

# Check which files exist
bigwig_files <- bigwig_files[file.exists(unlist(bigwig_files))]
cat("Found", length(bigwig_files), "BigWig files\n")

# ===================== Create Signal Matrices =====================

cat("\nCreating signal matrices around TSS...\n")

# Function to create normalized matrix
create_signal_matrix <- function(bw_file, regions, upstream, downstream, window_size) {
  cat("  Processing:", basename(bw_file), "\n")

  # Import BigWig
  coverage <- import(bw_file, format = "BigWig")

  # Check coverage statistics
  cat("    Coverage regions:", length(coverage), "\n")
  cat("    Coverage chromosomes:", paste(head(unique(seqnames(coverage)), 10), collapse = ", "), "\n")
  cat(
    "    Score range:", min(coverage$score, na.rm = TRUE), "to",
    max(coverage$score, na.rm = TRUE), "\n"
  )
  cat("    Mean score:", mean(coverage$score, na.rm = TRUE), "\n")

  # Check for chromosome naming mismatch
  cat("    Gene chromosomes:", paste(head(unique(seqnames(regions)), 10), collapse = ", "), "\n")

  # Find common chromosomes
  common_chrs <- intersect(seqlevels(coverage), seqlevels(regions))
  cat("    Common chromosomes:", length(common_chrs), "\n")

  if (length(common_chrs) == 0) {
    cat("    WARNING: No common chromosomes! Attempting to fix chromosome naming...\n")
    # Try to fix chromosome naming (chr1 vs 1)
    if (any(grepl("^chr", seqlevels(coverage))) && !any(grepl("^chr", seqlevels(regions)))) {
      seqlevels(regions) <- paste0("chr", seqlevels(regions))
      cat("    Added 'chr' prefix to gene regions\n")
    } else if (!any(grepl("^chr", seqlevels(coverage))) && any(grepl("^chr", seqlevels(regions)))) {
      seqlevels(coverage) <- paste0("chr", seqlevels(coverage))
      cat("    Added 'chr' prefix to coverage\n")
    }
  }

  # Subset to common chromosomes
  coverage <- keepSeqlevels(coverage, intersect(seqlevels(coverage), seqlevels(regions)),
    pruning.mode = "coarse"
  )
  regions <- keepSeqlevels(regions, intersect(seqlevels(coverage), seqlevels(regions)),
    pruning.mode = "coarse"
  )

  cat("    Regions after filtering:", length(regions), "\n")

  # Create matrix - DO NOT smooth, keep raw signal
  mat <- normalizeToMatrix(coverage, regions,
    value_column = "score",
    extend = c(upstream, downstream),
    mean_mode = "absolute",
    w = window_size,
    smooth = FALSE,
    empty_value = NA,
    keep = c(0, 0.99)
  )

  cat(
    "    Matrix range before cleanup:", min(mat, na.rm = TRUE), "to",
    max(mat, na.rm = TRUE), "\n"
  )
  cat("    NA count before cleanup:", sum(is.na(mat)), "\n")

  # Replace NA with 0 only after matrix creation
  mat[is.na(mat)] <- 0
  mat[is.nan(mat)] <- 0
  mat[is.infinite(mat)] <- 0

  cat("    Matrix range after cleanup:", min(mat), "to", max(mat), "\n")

  return(mat)
}

# Create matrices for all samples
signal_matrices <- lapply(bigwig_files, function(bw) {
  create_signal_matrix(bw, tss, UPSTREAM, DOWNSTREAM, WINDOW_SIZE)
})

# ===================== Generate Heatmaps =====================

cat("\nGenerating heatmaps...\n")

# Determine appropriate color scale based on actual data range
all_values <- unlist(lapply(signal_matrices, function(m) as.vector(m)))
all_values <- all_values[!is.na(all_values) & !is.infinite(all_values)]

data_min <- quantile(all_values, 0.01, na.rm = TRUE)
data_median <- quantile(all_values, 0.5, na.rm = TRUE)
data_max <- quantile(all_values, 0.99, na.rm = TRUE)

cat("\nData statistics:\n")
cat("  Min (1%):", data_min, "\n")
cat("  Median:", data_median, "\n")
cat("  Max (99%):", data_max, "\n")
cat("  Non-zero values:", sum(all_values > 0), "/", length(all_values), "\n")

# Color schemes
col_up <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
col_down <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

# Use dynamic range for signal, ensuring valid range
if (data_max > data_min && data_max > 0) {
  # Use actual data range
  col_signal <- colorRamp2(
    c(0, data_median, data_max),
    c("white", "yellow", "red")
  )
  cat("  Color scale: 0 ->", data_median, "->", data_max, "\n")
} else if (data_max == 0) {
  # All zeros - something is wrong, use minimal scale
  cat("  WARNING: All values are zero! Using minimal scale.\n")
  col_signal <- colorRamp2(c(0, 0.1, 1), c("white", "yellow", "red"))
} else {
  # Fallback
  col_signal <- colorRamp2(c(0, 0.5, 1), c("white", "yellow", "red"))
}

# Split genes by regulation
up_genes <- tss[tss$regulation == "UP"]
down_genes <- tss[tss$regulation == "DOWN"]

# Function to create heatmap for one sample
create_heatmap <- function(mat, sample_name, gene_regions) {
  # Order by mean signal
  gene_order <- order(rowMeans(mat), decreasing = TRUE)

  ht <- EnrichedHeatmap(mat[gene_order, ],
    col = col_signal,
    name = "Signal",
    column_title = sample_name,
    top_annotation = HeatmapAnnotation(
      lines = anno_enriched(gp = gpar(col = 2:4))
    ),
    axis_name = c(paste0("-", UPSTREAM / 1000, "kb"), "TSS", paste0("+", DOWNSTREAM / 1000, "kb")),
    use_raster = TRUE,
    raster_quality = 2
  )

  return(ht)
}

# Generate heatmaps for key samples
for (sample_name in c("TES_combined", "TEAD1_combined")) {
  if (sample_name %in% names(signal_matrices)) {
    cat("  Creating heatmap for", sample_name, "\n")

    mat_full <- signal_matrices[[sample_name]]

    # Validate matrix dimensions and values
    cat("    Matrix dimensions:", nrow(mat_full), "x", ncol(mat_full), "\n")
    cat("    Value range:", min(mat_full, na.rm = TRUE), "to", max(mat_full, na.rm = TRUE), "\n")
    cat("    NA count:", sum(is.na(mat_full)), "\n")

    # Skip if matrix is invalid
    if (nrow(mat_full) == 0 || all(is.na(mat_full))) {
      cat("    WARNING: Invalid matrix, skipping...\n")
      next
    }

    # Split by regulation and order by signal strength within each group
    up_idx <- which(tss$regulation == "UP")
    down_idx <- which(tss$regulation == "DOWN")

    if (length(up_idx) == 0 || length(down_idx) == 0) {
      cat("    WARNING: Missing UP or DOWN genes, skipping split...\n")
      next
    }

    # Calculate mean signal for ordering
    up_means <- rowMeans(mat_full[up_idx, ], na.rm = TRUE)
    down_means <- rowMeans(mat_full[down_idx, ], na.rm = TRUE)

    # Create ordering within each group
    up_order <- up_idx[order(up_means, decreasing = TRUE)]
    down_order <- down_idx[order(down_means, decreasing = TRUE)]

    # Combine orders
    row_order <- c(up_order, down_order)

    # Create combined heatmap with split
    pdf(file.path(OUTPUT_DIR, paste0(sample_name, "_heatmap_split.pdf")),
      width = 8, height = 10
    )

    # Set heatmap options to suppress messages and handle edge cases
    ht_opt$message <- FALSE

    ht <- EnrichedHeatmap(mat_full,
      col = col_signal,
      name = "CPM",
      column_title = paste(sample_name, "at TES DEGs"),
      split = tss$regulation,
      row_order = row_order,
      top_annotation = HeatmapAnnotation(
        enriched = anno_enriched(gp = gpar(col = c("red", "blue")))
      ),
      axis_name = c(
        paste0("-", UPSTREAM / 1000, "kb"), "TSS",
        paste0("+", DOWNSTREAM / 1000, "kb")
      ),
      use_raster = TRUE,
      raster_quality = 2,
      row_title_rot = 0,
      cluster_rows = FALSE,
      show_row_names = FALSE,
      heatmap_height = unit(0.8, "npc")
    ) # Set explicit height

    # Add gene annotation
    gene_anno <- rowAnnotation(
      Regulation = tss$regulation,
      col = list(Regulation = c("UP" = "red", "DOWN" = "blue")),
      show_annotation_name = TRUE,
      width = unit(5, "mm")
    )

    tryCatch(
      {
        draw(ht + gene_anno, heatmap_legend_side = "right")
      },
      error = function(e) {
        cat("    ERROR drawing heatmap:", e$message, "\n")
        cat("    Trying simplified version without annotation...\n")
        draw(ht, heatmap_legend_side = "right")
      }
    )

    dev.off()

    cat("    Saved:", file.path(OUTPUT_DIR, paste0(sample_name, "_heatmap_split.pdf")), "\n")
  }
}

# ===================== Generate Profile Plots =====================

cat("\nGenerating profile plots...\n")

# Function to calculate mean profile
calculate_profile <- function(mat) {
  profile <- colMeans(mat, na.rm = TRUE)
  return(profile)
}

# Calculate profiles
profiles <- lapply(signal_matrices, calculate_profile)

# Create comparison plot
pdf(file.path(OUTPUT_DIR, "profile_comparison.pdf"), width = 10, height = 6)

par(mfrow = c(1, 2))

# Plot 1: All samples combined
x_coords <- seq(-UPSTREAM, DOWNSTREAM, length.out = ncol(signal_matrices[[1]]))
plot(NULL,
  xlim = c(-UPSTREAM, DOWNSTREAM), ylim = c(0, max(unlist(profiles), na.rm = TRUE)),
  xlab = "Distance from TSS (bp)", ylab = "Mean CPM Signal",
  main = "Cut&Tag Signal at TES DEGs (All Genes)"
)
abline(v = 0, lty = 2, col = "gray")

colors <- rainbow(length(profiles))
for (i in seq_along(profiles)) {
  lines(x_coords, profiles[[i]], col = colors[i], lwd = 2)
}
legend("topright", legend = names(profiles), col = colors, lwd = 2, cex = 0.7)

# Plot 2: Separated by UP/DOWN regulation
mat_tes <- signal_matrices[["TES_combined"]]
profile_up <- colMeans(mat_tes[tss$regulation == "UP", ], na.rm = TRUE)
profile_down <- colMeans(mat_tes[tss$regulation == "DOWN", ], na.rm = TRUE)

plot(x_coords, profile_up,
  type = "l", col = "red", lwd = 2,
  xlim = c(-UPSTREAM, DOWNSTREAM),
  ylim = c(0, max(c(profile_up, profile_down), na.rm = TRUE)),
  xlab = "Distance from TSS (bp)", ylab = "Mean CPM Signal",
  main = "TES Signal: UP vs DOWN Regulated Genes"
)
lines(x_coords, profile_down, col = "blue", lwd = 2)
abline(v = 0, lty = 2, col = "gray")
legend("topright",
  legend = c("UP-regulated", "DOWN-regulated"),
  col = c("red", "blue"), lwd = 2
)

dev.off()

cat("  Saved:", file.path(OUTPUT_DIR, "profile_comparison.pdf"), "\n")

# ===================== Generate Combined Heatmap with All Samples =====================

cat("\nGenerating comprehensive multi-sample heatmap...\n")

pdf(file.path(OUTPUT_DIR, "all_samples_comprehensive_heatmap.pdf"),
  width = 14, height = 10
)

# Set heatmap options
ht_opt$message <- FALSE

# Create list of heatmaps
ht_list <- NULL

for (sample_name in c("TES_combined", "TEAD1_combined")) {
  if (sample_name %in% names(signal_matrices)) {
    cat("  Adding", sample_name, "to combined heatmap...\n")

    ht <- EnrichedHeatmap(signal_matrices[[sample_name]],
      col = col_signal,
      name = paste0(sample_name, "\nCPM"),
      column_title = sample_name,
      split = tss$regulation,
      axis_name = c("-3kb", "TSS", "+3kb"),
      use_raster = TRUE,
      raster_quality = 2,
      cluster_rows = FALSE,
      show_row_names = FALSE,
      row_title_rot = 0,
      width = unit(3, "cm"),
      heatmap_height = unit(0.75, "npc")
    )

    if (is.null(ht_list)) {
      ht_list <- ht
    } else {
      ht_list <- ht_list + ht
    }
  }
}

# Add annotation
gene_anno <- rowAnnotation(
  DEG = tss$regulation,
  col = list(DEG = c("UP" = "red", "DOWN" = "blue")),
  show_annotation_name = TRUE,
  width = unit(5, "mm")
)

tryCatch(
  {
    draw(gene_anno + ht_list,
      heatmap_legend_side = "right",
      padding = unit(c(2, 10, 2, 2), "mm")
    )
  },
  error = function(e) {
    cat("  ERROR drawing combined heatmap:", e$message, "\n")
    cat("  Trying without annotation...\n")
    draw(ht_list,
      heatmap_legend_side = "right",
      padding = unit(c(2, 10, 2, 2), "mm")
    )
  }
)

dev.off()

cat("  Saved:", file.path(OUTPUT_DIR, "all_samples_comprehensive_heatmap.pdf"), "\n")

# ===================== Summary Statistics =====================

cat("\n=== Analysis Summary ===\n")
cat("Total DEGs analyzed:", length(tss), "\n")
cat("UP-regulated genes:", sum(tss$regulation == "UP"), "\n")
cat("DOWN-regulated genes:", sum(tss$regulation == "DOWN"), "\n")
cat("\nMean signal at TSS:\n")
for (sample_name in names(signal_matrices)) {
  tss_col <- ncol(signal_matrices[[sample_name]]) / 2
  mean_signal <- mean(signal_matrices[[sample_name]][, tss_col], na.rm = TRUE)
  cat(sprintf("  %s: %.2f CPM\n", sample_name, mean_signal))
}

cat("\nAnalysis complete!\n")
cat("Results saved to:", OUTPUT_DIR, "\n")
