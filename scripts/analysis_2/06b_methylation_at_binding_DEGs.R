#!/usr/bin/env Rscript

# 06b_methylation_at_binding_DEGs.R
# Modified version of 06_methylation_at_regulated_genes.R
# Filters for DEGs with TES or TEAD1 binding (all peaks combined)
#
# Key difference: Only analyzes DEGs that have TES OR TEAD1 binding peaks
# at their promoters (within 2kb of TSS)

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(gridExtra)
})

# ===================== Configuration =====================

BASE_DIR <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top"
SCRIPT_DIR <- file.path(BASE_DIR, "SRF_Eva_integrated_analysis/scripts/analysis_2")

# Input files
RNA_SEQ_FILE <- file.path(BASE_DIR, "SRF_Eva_RNA/results/05_deseq2/deseq2_results_TES_vs_GFP.txt")
DMR_FILE <- file.path(BASE_DIR, "meDIP/results/07_differential_MEDIPS/TES_vs_GFP_DMRs_FDR05.csv")

# Peak files - TES and TEAD1 consensus peaks
TES_PEAKS_FILE <- file.path(BASE_DIR, "SRF_Eva_CUTandTAG/results/11_combined_replicates_narrow/peaks/TES_consensus_peaks.bed")
TEAD1_PEAKS_FILE <- file.path(BASE_DIR, "SRF_Eva_CUTandTAG/results/11_combined_replicates_narrow/peaks/TEAD1_consensus_peaks.bed")

# meDIP BigWig files for signal quantification
BIGWIG_DIR <- file.path(BASE_DIR, "meDIP/results/05_bigwig")
TES_BIGWIGS <- c(
  file.path(BIGWIG_DIR, "TES-1-IP_RPKM.bw"),
  file.path(BIGWIG_DIR, "TES-2-IP_RPKM.bw")
)
GFP_BIGWIGS <- c(
  file.path(BIGWIG_DIR, "GFP-1-IP_RPKM.bw"),
  file.path(BIGWIG_DIR, "GFP-2-IP_RPKM.bw")
)

# Output
OUTPUT_DIR <- file.path(SCRIPT_DIR, "results/06b_methylation_at_binding_DEGs")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Thresholds
PADJ_THRESHOLD <- 0.05
PROMOTER_UPSTREAM <- 2000
PROMOTER_DOWNSTREAM <- 500

# Gene body definition: TSS to TES (transcription end site)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# ===================== Functions =====================

# Get gene body regions (TSS to TES)
get_gene_bodies <- function(txdb) {
  cat("Extracting gene body regions...\n")

  # Get genes
  genes <- genes(txdb)

  # Add gene symbols
  symbols <- mapIds(org.Hs.eg.db,
                    keys = genes$gene_id,
                    column = "SYMBOL",
                    keytype = "ENTREZID",
                    multiVals = "first")
  genes$symbol <- symbols

  # Add Ensembl IDs for matching with RNA-seq
  ensembl <- mapIds(org.Hs.eg.db,
                    keys = genes$gene_id,
                    column = "ENSEMBL",
                    keytype = "ENTREZID",
                    multiVals = "first")
  genes$ensembl_id <- ensembl

  return(genes)
}

# Load and combine TES + TEAD1 peaks
load_combined_peaks <- function(tes_file, tead1_file) {
  cat("Loading TES and TEAD1 peaks...\n")

  # Load TES peaks - custom format: chr, start, end, count, samples
  tes_df <- read.delim(tes_file, header = FALSE, stringsAsFactors = FALSE)
  colnames(tes_df) <- c("chr", "start", "end", "count", "samples")
  tes_peaks <- GRanges(
    seqnames = paste0("chr", tes_df$chr),
    ranges = IRanges(start = tes_df$start, end = tes_df$end)
  )
  cat(sprintf("  TES peaks: %d\n", length(tes_peaks)))

  # Load TEAD1 peaks - custom format: chr, start, end, count, samples
  tead1_df <- read.delim(tead1_file, header = FALSE, stringsAsFactors = FALSE)
  colnames(tead1_df) <- c("chr", "start", "end", "count", "samples")
  tead1_peaks <- GRanges(
    seqnames = paste0("chr", tead1_df$chr),
    ranges = IRanges(start = tead1_df$start, end = tead1_df$end)
  )
  cat(sprintf("  TEAD1 peaks: %d\n", length(tead1_peaks)))

  # Combine (union)
  all_peaks <- c(tes_peaks, tead1_peaks)

  # Reduce overlapping peaks
  combined_peaks <- reduce(all_peaks)
  cat(sprintf("  Combined (merged) peaks: %d\n", length(combined_peaks)))

  return(combined_peaks)
}

# Get promoter regions for genes
get_promoters <- function(txdb, upstream = 2000, downstream = 500) {
  # Get transcripts
  txs <- transcripts(txdb)

  # Get promoters
  proms <- promoters(txs, upstream = upstream, downstream = downstream)

  return(proms)
}

# Extract mean methylation signal from BigWig files at given regions
extract_bigwig_signal <- function(regions, bigwig_files, name = "signal") {
  cat(sprintf("  Extracting %s from %d BigWig files...\n", name, length(bigwig_files)))

  signals <- matrix(NA, nrow = length(regions), ncol = length(bigwig_files))

  for (i in seq_along(bigwig_files)) {
    bw_file <- bigwig_files[i]
    cat(sprintf("    Processing: %s\n", basename(bw_file)))

    if (!file.exists(bw_file)) {
      warning(paste("BigWig not found:", bw_file))
      next
    }

    tryCatch({
      # Import BigWig
      bw <- import(bw_file, format = "BigWig", as = "RleList")
      bw_chroms <- names(bw)

      # Determine chromosome naming convention
      bw_has_chr <- any(grepl("^chr", bw_chroms))
      regions_has_chr <- any(grepl("^chr", seqlevels(regions)))

      # Extract signal for each region
      for (j in seq_along(regions)) {
        chr <- as.character(seqnames(regions[j]))

        # Map chromosome name if needed
        if (bw_has_chr && !grepl("^chr", chr)) {
          target_chr <- paste0("chr", chr)
        } else if (!bw_has_chr && grepl("^chr", chr)) {
          target_chr <- gsub("^chr", "", chr)
        } else {
          target_chr <- chr
        }

        if (!target_chr %in% bw_chroms) next

        tryCatch({
          st <- start(regions[j])
          en <- end(regions[j])
          chr_len <- length(bw[[target_chr]])

          st <- max(1, st)
          en <- min(chr_len, en)

          if (st <= en && en > 0) {
            region_signal <- as.numeric(bw[[target_chr]][st:en])
            signals[j, i] <- mean(region_signal, na.rm = TRUE)
          }
        }, error = function(e) NULL)
      }

    }, error = function(e) {
      warning(paste("Error reading BigWig:", e$message))
    })
  }

  # Return mean across replicates
  result <- rowMeans(signals, na.rm = TRUE)
  result[is.nan(result)] <- NA

  n_valid <- sum(!is.na(result))
  cat(sprintf("    Valid signals: %d/%d (%.1f%%)\n",
              n_valid, length(result), 100 * n_valid / length(result)))

  return(result)
}

# ===================== Main Analysis =====================

cat("=== Module 6b: Methylation at Binding DEGs ===\n\n")
cat("MODIFIED VERSION: Only DEGs with TES or TEAD1 binding\n\n")
cat(sprintf("Date: %s\n", Sys.time()))
cat(sprintf("Output: %s\n\n", OUTPUT_DIR))

# 1. Load RNA-seq differential expression results
cat("1. Loading RNA-seq differential expression results...\n")

if (!file.exists(RNA_SEQ_FILE)) {
  stop(paste("RNA-seq file not found:", RNA_SEQ_FILE))
}

rna_df <- read.delim(RNA_SEQ_FILE, stringsAsFactors = FALSE)
cat(sprintf("   Total genes in RNA-seq: %d\n", nrow(rna_df)))

# Clean gene IDs (remove version numbers)
rna_df$ensembl_clean <- gsub("\\..*", "", rna_df$gene_id)

# Identify DEGs
rna_df$is_deg <- !is.na(rna_df$padj) & rna_df$padj < PADJ_THRESHOLD
rna_df$direction <- ifelse(rna_df$log2FoldChange > 0, "UP", "DOWN")

n_degs <- sum(rna_df$is_deg, na.rm = TRUE)
n_up <- sum(rna_df$is_deg & rna_df$direction == "UP", na.rm = TRUE)
n_down <- sum(rna_df$is_deg & rna_df$direction == "DOWN", na.rm = TRUE)

cat(sprintf("   Significant DEGs (padj < %.2f): %d\n", PADJ_THRESHOLD, n_degs))
cat(sprintf("   - Upregulated: %d\n", n_up))
cat(sprintf("   - Downregulated: %d\n", n_down))

# 2. Load TES and TEAD1 peaks (combined)
cat("\n2. Loading and combining TES + TEAD1 binding peaks...\n")

if (!file.exists(TES_PEAKS_FILE)) {
  stop(paste("TES peaks file not found:", TES_PEAKS_FILE))
}
if (!file.exists(TEAD1_PEAKS_FILE)) {
  stop(paste("TEAD1 peaks file not found:", TEAD1_PEAKS_FILE))
}

combined_peaks <- load_combined_peaks(TES_PEAKS_FILE, TEAD1_PEAKS_FILE)

# 3. Get gene body regions and identify bound genes
cat("\n3. Getting gene body regions and identifying bound genes...\n")

gene_bodies <- get_gene_bodies(txdb)
cat(sprintf("   Total genes with coordinates: %d\n", length(gene_bodies)))

# Filter to standard chromosomes
standard_chr <- paste0("chr", c(1:22, "X", "Y"))
gene_bodies <- keepSeqlevels(gene_bodies, standard_chr[standard_chr %in% seqlevels(gene_bodies)],
                              pruning.mode = "coarse")

# Get promoter regions for each gene
gene_promoters <- promoters(gene_bodies, upstream = PROMOTER_UPSTREAM, downstream = PROMOTER_DOWNSTREAM)

# Find genes with peaks at promoters
cat("\n   Finding genes with binding peaks at promoters...\n")
promoter_hits <- findOverlaps(gene_promoters, combined_peaks)
bound_gene_idx <- unique(queryHits(promoter_hits))
cat(sprintf("   Genes with TES or TEAD1 binding at promoter: %d\n", length(bound_gene_idx)))

# Mark which genes have binding
gene_bodies$has_binding <- seq_along(gene_bodies) %in% bound_gene_idx

# Match genes to RNA-seq data
gene_bodies$in_rnaseq <- gene_bodies$ensembl_id %in% rna_df$ensembl_clean
gene_bodies$is_deg <- gene_bodies$ensembl_id %in% rna_df$ensembl_clean[rna_df$is_deg]

# Add expression direction
rna_deg <- rna_df[rna_df$is_deg, ]
deg_direction <- setNames(rna_deg$direction, rna_deg$ensembl_clean)
gene_bodies$expr_direction <- deg_direction[gene_bodies$ensembl_id]
gene_bodies$expr_direction[is.na(gene_bodies$expr_direction)] <- "Not_DEG"

cat(sprintf("   Genes matched to RNA-seq: %d\n", sum(gene_bodies$in_rnaseq)))
cat(sprintf("   DEGs with coordinates: %d\n", sum(gene_bodies$is_deg)))

# KEY FILTER: DEGs with binding
gene_bodies$is_binding_deg <- gene_bodies$is_deg & gene_bodies$has_binding
n_binding_degs <- sum(gene_bodies$is_binding_deg)
cat(sprintf("\n   *** DEGs WITH TES/TEAD1 BINDING: %d ***\n", n_binding_degs))

# Count by direction
n_bound_up <- sum(gene_bodies$is_binding_deg & gene_bodies$expr_direction == "UP")
n_bound_down <- sum(gene_bodies$is_binding_deg & gene_bodies$expr_direction == "DOWN")
cat(sprintf("   - Upregulated with binding: %d\n", n_bound_up))
cat(sprintf("   - Downregulated with binding: %d\n", n_bound_down))

# 4. Load DMRs (for reference)
cat("\n4. Loading DMRs...\n")

if (!file.exists(DMR_FILE)) {
  stop(paste("DMR file not found:", DMR_FILE))
}

dmr_df <- read.csv(DMR_FILE, stringsAsFactors = FALSE)
cat(sprintf("   Total DMRs: %d\n", nrow(dmr_df)))

# Create GRanges for DMRs
dmr_chr <- dmr_df$chr
if (!grepl("^chr", dmr_chr[1])) {
  dmr_chr <- paste0("chr", dmr_chr)
}

dmrs <- GRanges(
  seqnames = dmr_chr,
  ranges = IRanges(start = dmr_df$start, end = dmr_df$stop),
  logFC = dmr_df$logFC,
  FDR = dmr_df$FDR
)

# 5. BigWig Methylation Signal Analysis at Binding DEGs
cat("\n5. Analyzing methylation signal from BigWig files (BINDING DEGs ONLY)...\n")

# Check if BigWig files exist
tes_bw_exist <- all(file.exists(TES_BIGWIGS))
gfp_bw_exist <- all(file.exists(GFP_BIGWIGS))

if (tes_bw_exist && gfp_bw_exist) {

  # Get BINDING DEG gene bodies for signal extraction
  binding_deg_gene_bodies <- gene_bodies[gene_bodies$is_binding_deg]
  cat(sprintf("  Binding DEG gene bodies to analyze: %d\n", length(binding_deg_gene_bodies)))

  # Extract methylation signal from TES and GFP BigWigs
  cat("\n  Extracting TES methylation signal...\n")
  binding_deg_gene_bodies$tes_meth_signal <- extract_bigwig_signal(
    binding_deg_gene_bodies, TES_BIGWIGS, "TES methylation"
  )

  cat("\n  Extracting GFP methylation signal...\n")
  binding_deg_gene_bodies$gfp_meth_signal <- extract_bigwig_signal(
    binding_deg_gene_bodies, GFP_BIGWIGS, "GFP methylation"
  )

  # Calculate methylation difference (TES - GFP)
  binding_deg_gene_bodies$meth_diff <- binding_deg_gene_bodies$tes_meth_signal - binding_deg_gene_bodies$gfp_meth_signal

  # Create data frame for analysis
  meth_signal_df <- data.frame(
    ensembl_id = binding_deg_gene_bodies$ensembl_id,
    symbol = binding_deg_gene_bodies$symbol,
    expr_direction = binding_deg_gene_bodies$expr_direction,
    tes_meth = binding_deg_gene_bodies$tes_meth_signal,
    gfp_meth = binding_deg_gene_bodies$gfp_meth_signal,
    meth_diff = binding_deg_gene_bodies$meth_diff,
    stringsAsFactors = FALSE
  )

  # Filter for valid signals
  meth_signal_df <- meth_signal_df[!is.na(meth_signal_df$meth_diff), ]
  cat(sprintf("\n  Binding DEGs with valid methylation signal: %d\n", nrow(meth_signal_df)))

  # 5.1 Summary statistics by expression direction
  cat("\n  5.1 Methylation signal summary by expression direction:\n")

  meth_summary <- meth_signal_df %>%
    group_by(expr_direction) %>%
    summarise(
      n_genes = n(),
      mean_tes_meth = mean(tes_meth, na.rm = TRUE),
      median_tes_meth = median(tes_meth, na.rm = TRUE),
      mean_gfp_meth = mean(gfp_meth, na.rm = TRUE),
      median_gfp_meth = median(gfp_meth, na.rm = TRUE),
      mean_meth_diff = mean(meth_diff, na.rm = TRUE),
      median_meth_diff = median(meth_diff, na.rm = TRUE),
      .groups = "drop"
    )

  print(meth_summary)

  # 5.2 Statistical tests
  cat("\n  5.2 Statistical tests:\n")

  up_meth <- meth_signal_df$meth_diff[meth_signal_df$expr_direction == "UP"]
  down_meth <- meth_signal_df$meth_diff[meth_signal_df$expr_direction == "DOWN"]

  wilcox_test <- NULL
  ttest_up <- NULL
  ttest_down <- NULL

  # Wilcoxon test comparing UP vs DOWN
  if (length(up_meth) > 10 && length(down_meth) > 10) {
    wilcox_test <- wilcox.test(up_meth, down_meth)
    cat(sprintf("    Wilcoxon test (UP vs DOWN methylation difference):\n"))
    cat(sprintf("      p-value: %.4e\n", wilcox_test$p.value))
    cat(sprintf("      UP genes mean meth diff: %.4f\n", mean(up_meth, na.rm = TRUE)))
    cat(sprintf("      DOWN genes mean meth diff: %.4f\n", mean(down_meth, na.rm = TRUE)))

    # One-sample t-test: Is methylation difference significantly different from 0?
    ttest_up <- t.test(up_meth)
    ttest_down <- t.test(down_meth)
    cat(sprintf("\n    One-sample t-test (meth diff vs 0):\n"))
    cat(sprintf("      UP genes: p = %.4e (mean = %.4f)\n",
                ttest_up$p.value, ttest_up$estimate))
    cat(sprintf("      DOWN genes: p = %.4e (mean = %.4f)\n",
                ttest_down$p.value, ttest_down$estimate))
  } else {
    cat("    Not enough samples for statistical tests (need >10 in each group)\n")
    cat(sprintf("    UP genes: %d, DOWN genes: %d\n", length(up_meth), length(down_meth)))
  }

  # 5.3 Visualizations
  cat("\n  5.3 Generating methylation signal plots...\n")

  # Plot 1: Methylation signal boxplot by expression direction
  pdf(file.path(OUTPUT_DIR, "01_Methylation_signal_by_direction.pdf"), width = 12, height = 6)

  # Panel A: TES methylation
  p1a <- ggplot(meth_signal_df, aes(x = expr_direction, y = tes_meth, fill = expr_direction)) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.2, outlier.shape = NA) +
    scale_fill_manual(values = c("UP" = "#4DAF4A", "DOWN" = "#984EA3")) +
    theme_minimal(base_size = 14) +
    labs(
      title = "TES Condition Methylation",
      x = "Expression Direction",
      y = "Mean Methylation Signal (RPKM)"
    ) +
    theme(legend.position = "none")

  # Panel B: GFP methylation
  p1b <- ggplot(meth_signal_df, aes(x = expr_direction, y = gfp_meth, fill = expr_direction)) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.2, outlier.shape = NA) +
    scale_fill_manual(values = c("UP" = "#4DAF4A", "DOWN" = "#984EA3")) +
    theme_minimal(base_size = 14) +
    labs(
      title = "GFP Condition Methylation",
      x = "Expression Direction",
      y = "Mean Methylation Signal (RPKM)"
    ) +
    theme(legend.position = "none")

  # Combine panels
  gridExtra::grid.arrange(p1a, p1b, ncol = 2,
                          top = "Methylation Signal at Gene Bodies of BINDING DEGs")
  dev.off()

  # Plot 2: Methylation difference (TES - GFP) by expression direction
  pdf(file.path(OUTPUT_DIR, "02_Methylation_difference_by_direction.pdf"), width = 10, height = 6)

  wilcox_label <- if (!is.null(wilcox_test)) {
    sprintf("Wilcoxon p = %.2e", wilcox_test$p.value)
  } else {
    ""
  }

  p2 <- ggplot(meth_signal_df, aes(x = expr_direction, y = meth_diff, fill = expr_direction)) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.2, outlier.shape = NA) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    scale_fill_manual(values = c("UP" = "#4DAF4A", "DOWN" = "#984EA3")) +
    theme_minimal(base_size = 14) +
    labs(
      title = "Methylation Change (TES - GFP) at Gene Bodies",
      subtitle = paste("BINDING DEGs ONLY (TES or TEAD1 bound) |", wilcox_label),
      x = "Expression Direction",
      y = "Methylation Difference (TES - GFP)"
    ) +
    theme(legend.position = "none") +
    coord_cartesian(ylim = quantile(meth_signal_df$meth_diff, c(0.01, 0.99), na.rm = TRUE))
  print(p2)

  dev.off()

  # Plot 3: Scatter plot of expression vs methylation change
  pdf(file.path(OUTPUT_DIR, "03_Expression_vs_Methylation_scatter.pdf"), width = 10, height = 8)

  # Add log2FC to meth_signal_df
  rna_fc <- setNames(rna_df$log2FoldChange, rna_df$ensembl_clean)
  meth_signal_df$log2FC <- rna_fc[meth_signal_df$ensembl_id]

  # Calculate correlation
  cor_test <- cor.test(meth_signal_df$log2FC, meth_signal_df$meth_diff,
                       method = "spearman", use = "complete.obs")

  p3 <- ggplot(meth_signal_df, aes(x = log2FC, y = meth_diff, color = expr_direction)) +
    geom_point(alpha = 0.5, size = 2) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
    geom_hline(yintercept = 0, linetype = "dotted", color = "gray50") +
    geom_vline(xintercept = 0, linetype = "dotted", color = "gray50") +
    scale_color_manual(values = c("UP" = "#4DAF4A", "DOWN" = "#984EA3")) +
    theme_minimal(base_size = 14) +
    labs(
      title = "Expression Change vs Methylation Change (BINDING DEGs)",
      subtitle = sprintf("Spearman rho = %.3f, p = %.2e", cor_test$estimate, cor_test$p.value),
      x = "log2 Fold Change (Expression)",
      y = "Methylation Difference (TES - GFP)",
      color = "Direction"
    ) +
    theme(legend.position = "right")
  print(p3)

  dev.off()

  # Plot 4: Density plot of methylation difference
  pdf(file.path(OUTPUT_DIR, "04_Methylation_difference_density.pdf"), width = 10, height = 6)

  p4 <- ggplot(meth_signal_df, aes(x = meth_diff, fill = expr_direction, color = expr_direction)) +
    geom_density(alpha = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray30") +
    scale_fill_manual(values = c("UP" = "#4DAF4A", "DOWN" = "#984EA3")) +
    scale_color_manual(values = c("UP" = "#4DAF4A", "DOWN" = "#984EA3")) +
    theme_minimal(base_size = 14) +
    labs(
      title = "Distribution of Methylation Change at BINDING DEG Gene Bodies",
      x = "Methylation Difference (TES - GFP)",
      y = "Density",
      fill = "Expression", color = "Expression"
    ) +
    theme(legend.position = "right") +
    coord_cartesian(xlim = quantile(meth_signal_df$meth_diff, c(0.01, 0.99), na.rm = TRUE))
  print(p4)

  dev.off()

  # Save methylation signal data
  write.csv(meth_signal_df, file.path(OUTPUT_DIR, "Binding_DEG_methylation_signal.csv"), row.names = FALSE)
  write.csv(meth_summary, file.path(OUTPUT_DIR, "Methylation_signal_summary.csv"), row.names = FALSE)

  bigwig_analysis_done <- TRUE

} else {
  cat("  WARNING: BigWig files not found. Skipping signal analysis.\n")
  cat(sprintf("    TES BigWigs exist: %s\n", tes_bw_exist))
  cat(sprintf("    GFP BigWigs exist: %s\n", gfp_bw_exist))
  bigwig_analysis_done <- FALSE
}

# 6. Save final results and generate summary
cat("\n6. Generating summary report...\n")

# Build summary report
summary_text <- sprintf("
=== Methylation at BINDING DEGs Analysis ===
Date: %s

MODIFICATION: This analysis filters for DEGs with TES OR TEAD1 binding
at promoters (within %d bp upstream, %d bp downstream of TSS)

Input Data:
  RNA-seq: %s
  DMRs: %s
  TES peaks: %s
  TEAD1 peaks: %s

Gene Counts:
  Total DEGs (padj < %.2f): %d
    - Upregulated: %d
    - Downregulated: %d

  Genes with TES/TEAD1 binding at promoter: %d

  *** BINDING DEGs (DEG + Bound): %d ***
    - Upregulated with binding: %d
    - Downregulated with binding: %d

DMR Statistics:
  Total DMRs: %d
",
Sys.time(),
PROMOTER_UPSTREAM, PROMOTER_DOWNSTREAM,
RNA_SEQ_FILE,
DMR_FILE,
TES_PEAKS_FILE,
TEAD1_PEAKS_FILE,
PADJ_THRESHOLD,
n_degs, n_up, n_down,
length(bound_gene_idx),
n_binding_degs,
n_bound_up, n_bound_down,
length(dmrs)
)

# Add BigWig analysis section if it was performed
if (bigwig_analysis_done && exists("meth_signal_df") && exists("meth_summary")) {
  bigwig_summary <- sprintf("
BigWig Methylation Signal Analysis (BINDING DEGs ONLY):
  Binding DEGs with valid signal: %d

  Methylation by Expression Direction:
    UP-regulated genes:
      Mean TES methylation: %.4f
      Mean GFP methylation: %.4f
      Mean difference (TES-GFP): %.4f
      N genes: %d
    DOWN-regulated genes:
      Mean TES methylation: %.4f
      Mean GFP methylation: %.4f
      Mean difference (TES-GFP): %.4f
      N genes: %d
",
    nrow(meth_signal_df),
    meth_summary$mean_tes_meth[meth_summary$expr_direction == "UP"],
    meth_summary$mean_gfp_meth[meth_summary$expr_direction == "UP"],
    meth_summary$mean_meth_diff[meth_summary$expr_direction == "UP"],
    meth_summary$n_genes[meth_summary$expr_direction == "UP"],
    meth_summary$mean_tes_meth[meth_summary$expr_direction == "DOWN"],
    meth_summary$mean_gfp_meth[meth_summary$expr_direction == "DOWN"],
    meth_summary$mean_meth_diff[meth_summary$expr_direction == "DOWN"],
    meth_summary$n_genes[meth_summary$expr_direction == "DOWN"]
  )

  # Add statistical test results if they were computed
  if (!is.null(wilcox_test)) {
    bigwig_summary <- paste0(bigwig_summary, sprintf("
  Statistical Tests:
    Wilcoxon test (UP vs DOWN meth diff):
      p-value: %.4e
    Spearman correlation (expression vs methylation):
      rho: %.4f
      p-value: %.4e
",
      wilcox_test$p.value,
      if(exists("cor_test")) cor_test$estimate else NA,
      if(exists("cor_test")) cor_test$p.value else NA
    ))
  }

  summary_text <- paste0(summary_text, bigwig_summary)
}

# Add output files section
output_files_text <- "
Output Files:
  - Binding_DEG_methylation_signal.csv: Methylation signal per binding DEG
  - Methylation_signal_summary.csv: Summary by direction
  - 01_Methylation_signal_by_direction.pdf: Signal boxplots
  - 02_Methylation_difference_by_direction.pdf: TES-GFP difference
  - 03_Expression_vs_Methylation_scatter.pdf: Correlation plot
  - 04_Methylation_difference_density.pdf: Density distributions
  - ANALYSIS_SUMMARY.txt: This file
"

summary_text <- paste0(summary_text, output_files_text)

writeLines(summary_text, file.path(OUTPUT_DIR, "ANALYSIS_SUMMARY.txt"))
cat(summary_text)

cat("\n=== Module 6b Complete ===\n")
cat(sprintf("Results saved to: %s\n", OUTPUT_DIR))
