#!/usr/bin/env Rscript

################################################################################
# Cut&Tag Density Heatmap Analysis for TES Target Genes
#
# This script generates density plots and heatmaps showing Cut&Tag signal
# enrichment around TES differentially expressed genes
################################################################################

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
  library(readr)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(EnrichedHeatmap)
})

# Set up paths
base_dir <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top"
deg_file <- file.path(base_dir, "SRF_Eva_integrated_analysis/data/TES_degs.txt")
bigwig_dir <- file.path(base_dir, "SRF_Eva_CUTandTAG/results/06_bigwig")
output_dir <- file.path(base_dir, "SRF_Eva_integrated_analysis/plots")
deseq_file <- file.path(base_dir, "SRF_Eva_RNA/results/05_deseq2/deseq2_results_TES_vs_GFP.txt")

# Create output directory if needed
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

cat("Loading data...\n")

# Load TES DEGs
deg_genes <- read_lines(deg_file)
cat(sprintf("Loaded %d TES target genes\n", length(deg_genes)))

# Load full DESeq2 results to get expression values
deseq_results <- read_tsv(deseq_file, show_col_types = FALSE)

# Get gene annotations from TxDb
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Get genes and convert to GRanges
genes_gr <- genes(txdb)

# Convert gene IDs to symbols
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = names(genes_gr),
                       column = "SYMBOL",
                       keytype = "ENTREZID",
                       multiVals = "first")

genes_gr$symbol <- gene_symbols

# Filter to only TES DEGs that we can find
genes_gr_filtered <- genes_gr[genes_gr$symbol %in% deg_genes]
genes_gr_filtered <- genes_gr_filtered[!is.na(genes_gr_filtered$symbol)]

cat(sprintf("Mapped %d/%d genes to genomic locations\n",
            length(genes_gr_filtered), length(deg_genes)))

# Get TSS positions
tss_gr <- promoters(genes_gr_filtered, upstream = 0, downstream = 1)

# Load BigWig files
cat("Loading BigWig files...\n")

bigwig_files <- list(
  TES_1 = file.path(bigwig_dir, "TES-1_CPM.bw"),
  TES_2 = file.path(bigwig_dir, "TES-2_CPM.bw"),
  TES_3 = file.path(bigwig_dir, "TES-3_CPM.bw"),
  GFP_1 = file.path(bigwig_dir, "IggRb_CPM.bw"),
  GFP_2 = file.path(bigwig_dir, "IggMs_CPM.bw")
)

# Define window around TSS
window_size <- 3000  # 3kb upstream and downstream

# Function to create normalized matrix from BigWig
create_signal_matrix <- function(bw_file, regions, upstream = 3000, downstream = 3000) {
  cat(sprintf("Processing %s...\n", basename(bw_file)))

  # Import BigWig
  bw <- import(bw_file, format = "BigWig")

  # Create windows around regions
  mat <- normalizeToMatrix(
    bw,
    regions,
    value_column = "score",
    extend = c(upstream, downstream),
    mean_mode = "w0",
    w = 50,  # 50bp bins
    background = 0,
    smooth = FALSE
  )

  return(mat)
}

# Create matrices for all samples
cat("\nGenerating signal matrices...\n")
signal_matrices <- list()

for (sample_name in names(bigwig_files)) {
  if (file.exists(bigwig_files[[sample_name]])) {
    signal_matrices[[sample_name]] <- create_signal_matrix(
      bigwig_files[[sample_name]],
      tss_gr,
      upstream = window_size,
      downstream = window_size
    )
  } else {
    cat(sprintf("Warning: BigWig file not found: %s\n", bigwig_files[[sample_name]]))
  }
}

# Merge expression data
expression_values <- deseq_results %>%
  filter(gene_name %in% deg_genes) %>%
  dplyr::select(gene_name, log2FoldChange, padj) %>%
  mutate(regulation = ifelse(log2FoldChange > 0, "UP", "DOWN"))

# Match expression to our TSS regions
tss_df <- as.data.frame(tss_gr)
tss_df$gene_name <- tss_df$symbol

expression_matched <- left_join(
  tss_df %>% dplyr::select(gene_name),
  expression_values,
  by = "gene_name"
)

cat("\nGenerating heatmaps...\n")

# Create color functions
col_fun_signal <- colorRamp2(c(0, 1, 2, 3, 5), c("white", "lightblue", "blue", "darkblue", "black"))
col_fun_expr <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

# Sort genes by expression level
gene_order <- order(expression_matched$log2FoldChange, decreasing = TRUE, na.last = TRUE)

# --- Heatmap 1: TES replicates vs GFP controls ---
pdf(file.path(output_dir, "TES_DEG_cutntag_heatmap.pdf"), width = 12, height = 16)

# Create heatmap list
ht_list <- NULL

# Add expression annotation
ha_expr <- rowAnnotation(
  log2FC = expression_matched$log2FoldChange[gene_order],
  regulation = expression_matched$regulation[gene_order],
  col = list(
    log2FC = col_fun_expr,
    regulation = c("UP" = "red", "DOWN" = "blue")
  ),
  annotation_name_side = "top",
  width = unit(1, "cm")
)

# Add TES replicates
for (i in 1:3) {
  sample_name <- paste0("TES_", i)
  if (!is.null(signal_matrices[[sample_name]])) {
    ht <- EnrichedHeatmap(
      signal_matrices[[sample_name]][gene_order, ],
      col = col_fun_signal,
      name = sample_name,
      column_title = sample_name,
      axis_name = c("-3kb", "TSS", "+3kb"),
      axis_name_gp = gpar(fontsize = 8),
      use_raster = TRUE,
      raster_quality = 2,
      top_annotation = HeatmapAnnotation(
        lines = anno_enriched(
          gp = gpar(col = "red", lwd = 2)
        )
      ),
      width = unit(3, "cm")
    )

    if (is.null(ht_list)) {
      ht_list <- ht
    } else {
      ht_list <- ht_list + ht
    }
  }
}

# Add GFP controls
for (i in 1:2) {
  sample_name <- paste0("GFP_", i)
  if (!is.null(signal_matrices[[sample_name]])) {
    ht <- EnrichedHeatmap(
      signal_matrices[[sample_name]][gene_order, ],
      col = col_fun_signal,
      name = sample_name,
      column_title = sample_name,
      axis_name = c("-3kb", "TSS", "+3kb"),
      axis_name_gp = gpar(fontsize = 8),
      use_raster = TRUE,
      raster_quality = 2,
      top_annotation = HeatmapAnnotation(
        lines = anno_enriched(
          gp = gpar(col = "blue", lwd = 2)
        )
      ),
      width = unit(3, "cm")
    )

    ht_list <- ht_list + ht
  }
}

# Add expression annotation
ht_list <- ht_list + ha_expr

# Draw heatmap
draw(ht_list,
     column_title = "Cut&Tag Signal at TES Target Genes (Sorted by Expression)",
     heatmap_legend_side = "bottom")

dev.off()

# --- Heatmap 2: Combined samples (TES, TESmut, TEAD1) ---
cat("Creating combined sample heatmap...\n")

combined_bigwig <- list(
  TES = file.path(bigwig_dir, "TES_comb.bw"),
  TESmut = file.path(bigwig_dir, "TESmut_comb.bw"),
  TEAD1 = file.path(bigwig_dir, "TEAD1_comb.bw")
)

combined_matrices <- list()
for (sample_name in names(combined_bigwig)) {
  if (file.exists(combined_bigwig[[sample_name]])) {
    combined_matrices[[sample_name]] <- create_signal_matrix(
      combined_bigwig[[sample_name]],
      tss_gr,
      upstream = window_size,
      downstream = window_size
    )
  }
}

pdf(file.path(output_dir, "TES_DEG_cutntag_combined_heatmap.pdf"), width = 10, height = 16)

ht_list2 <- NULL

# Add combined heatmaps
for (sample_name in names(combined_matrices)) {
  if (!is.null(combined_matrices[[sample_name]])) {
    ht <- EnrichedHeatmap(
      combined_matrices[[sample_name]][gene_order, ],
      col = col_fun_signal,
      name = sample_name,
      column_title = sample_name,
      axis_name = c("-3kb", "TSS", "+3kb"),
      axis_name_gp = gpar(fontsize = 8),
      use_raster = TRUE,
      raster_quality = 2,
      top_annotation = HeatmapAnnotation(
        lines = anno_enriched(
          gp = gpar(col = c("TES" = "red", "TESmut" = "orange", "TEAD1" = "green")[sample_name], lwd = 2)
        )
      ),
      width = unit(4, "cm")
    )

    if (is.null(ht_list2)) {
      ht_list2 <- ht
    } else {
      ht_list2 <- ht_list2 + ht
    }
  }
}

# Add expression annotation
ht_list2 <- ht_list2 + ha_expr

# Draw heatmap
draw(ht_list2,
     column_title = "Cut&Tag Signal at TES Target Genes - Combined Samples",
     heatmap_legend_side = "bottom")

dev.off()

# --- Profile plots ---
cat("Creating profile plots...\n")

pdf(file.path(output_dir, "TES_DEG_cutntag_profiles.pdf"), width = 10, height = 8)

par(mfrow = c(2, 2))

# Function to calculate mean profile
calc_profile <- function(mat) {
  colMeans(mat, na.rm = TRUE)
}

# Plot TES replicates
plot(NULL, xlim = c(-window_size, window_size),
     ylim = c(0, max(sapply(signal_matrices[1:3], function(m) max(calc_profile(m), na.rm = TRUE)))),
     xlab = "Distance from TSS (bp)",
     ylab = "Mean Cut&Tag Signal",
     main = "TES Replicates")
abline(v = 0, col = "gray", lty = 2)

colors <- c("red", "darkred", "firebrick")
for (i in 1:3) {
  sample_name <- paste0("TES_", i)
  if (!is.null(signal_matrices[[sample_name]])) {
    profile <- calc_profile(signal_matrices[[sample_name]])
    x_coords <- seq(-window_size, window_size, length.out = length(profile))
    lines(x_coords, profile, col = colors[i], lwd = 2)
  }
}
legend("topright", legend = c("TES-1", "TES-2", "TES-3"),
       col = colors, lwd = 2, bty = "n")

# Plot combined samples
if (all(sapply(combined_matrices, function(x) !is.null(x)))) {
  plot(NULL, xlim = c(-window_size, window_size),
       ylim = c(0, max(sapply(combined_matrices, function(m) max(calc_profile(m), na.rm = TRUE)))),
       xlab = "Distance from TSS (bp)",
       ylab = "Mean Cut&Tag Signal",
       main = "Combined Samples")
  abline(v = 0, col = "gray", lty = 2)

  sample_colors <- c("TES" = "red", "TESmut" = "orange", "TEAD1" = "green")
  for (sample_name in names(combined_matrices)) {
    profile <- calc_profile(combined_matrices[[sample_name]])
    x_coords <- seq(-window_size, window_size, length.out = length(profile))
    lines(x_coords, profile, col = sample_colors[sample_name], lwd = 2)
  }
  legend("topright", legend = names(combined_matrices),
         col = sample_colors[names(combined_matrices)], lwd = 2, bty = "n")
}

# Plot separated by up/down regulation
up_genes <- which(expression_matched$regulation == "UP" & !is.na(expression_matched$regulation))
down_genes <- which(expression_matched$regulation == "DOWN" & !is.na(expression_matched$regulation))

if (length(up_genes) > 0 && !is.null(combined_matrices$TES)) {
  profile_up <- colMeans(combined_matrices$TES[up_genes, ], na.rm = TRUE)
  profile_down <- colMeans(combined_matrices$TES[down_genes, ], na.rm = TRUE)

  plot(NULL, xlim = c(-window_size, window_size),
       ylim = c(0, max(c(profile_up, profile_down), na.rm = TRUE)),
       xlab = "Distance from TSS (bp)",
       ylab = "Mean TES Cut&Tag Signal",
       main = "TES Signal by Gene Regulation")
  abline(v = 0, col = "gray", lty = 2)

  x_coords <- seq(-window_size, window_size, length.out = length(profile_up))
  lines(x_coords, profile_up, col = "red", lwd = 2)
  lines(x_coords, profile_down, col = "blue", lwd = 2)
  legend("topright", legend = c(sprintf("Upregulated (n=%d)", length(up_genes)),
                                 sprintf("Downregulated (n=%d)", length(down_genes))),
         col = c("red", "blue"), lwd = 2, bty = "n")
}

dev.off()

# Create summary statistics
cat("\n=== Analysis Summary ===\n")
cat(sprintf("Total TES DEGs analyzed: %d\n", length(deg_genes)))
cat(sprintf("Genes mapped to genome: %d\n", length(genes_gr_filtered)))
cat(sprintf("Upregulated genes: %d\n", sum(expression_matched$regulation == "UP", na.rm = TRUE)))
cat(sprintf("Downregulated genes: %d\n", sum(expression_matched$regulation == "DOWN", na.rm = TRUE)))

# Save summary to file
summary_file <- file.path(output_dir, "cutntag_density_summary.txt")
sink(summary_file)
cat("Cut&Tag Density Analysis Summary\n")
cat("================================\n\n")
cat(sprintf("Analysis date: %s\n", Sys.Date()))
cat(sprintf("Total TES DEGs analyzed: %d\n", length(deg_genes)))
cat(sprintf("Genes mapped to genome: %d\n", length(genes_gr_filtered)))
cat(sprintf("Upregulated genes: %d\n", sum(expression_matched$regulation == "UP", na.rm = TRUE)))
cat(sprintf("Downregulated genes: %d\n", sum(expression_matched$regulation == "DOWN", na.rm = TRUE)))
cat(sprintf("Window size: +/- %d bp around TSS\n", window_size))
cat("\nOutput files:\n")
cat("  - TES_DEG_cutntag_heatmap.pdf: Individual replicate heatmaps\n")
cat("  - TES_DEG_cutntag_combined_heatmap.pdf: Combined sample heatmaps\n")
cat("  - TES_DEG_cutntag_profiles.pdf: Profile plots\n")
sink()

cat("\nAnalysis complete! Output files saved to:\n")
cat(sprintf("  %s\n", output_dir))
cat("\nGenerated files:\n")
cat("  - TES_DEG_cutntag_heatmap.pdf\n")
cat("  - TES_DEG_cutntag_combined_heatmap.pdf\n")
cat("  - TES_DEG_cutntag_profiles.pdf\n")
cat("  - cutntag_density_summary.txt\n")
