#!/usr/bin/env Rscript

# 06_methylation_at_regulated_genes.R
# Examine whether TES-induced methylation changes occur at gene bodies
# of transcriptionally regulated genes (indirect methylation hypothesis)
#
# Tests the hypothesis that DNMT3A/3L domains in TES cause methylation
# changes at regulated genes rather than at TES binding sites

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(gridExtra)  # For multi-panel plots
})

# ===================== Configuration =====================

BASE_DIR <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top"
SCRIPT_DIR <- file.path(BASE_DIR, "SRF_Eva_integrated_analysis/scripts/analysis_2")

# Input files
RNA_SEQ_FILE <- file.path(BASE_DIR, "SRF_Eva_RNA/results/05_deseq2/deseq2_results_TES_vs_GFP.txt")
DMR_FILE <- file.path(BASE_DIR, "meDIP/results/07_differential_MEDIPS/TES_vs_GFP_DMRs_FDR05.csv")

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
OUTPUT_DIR <- file.path(SCRIPT_DIR, "results/06_methylation_at_regulated_genes")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Thresholds
PADJ_THRESHOLD <- 0.05
LOG2FC_THRESHOLD <- 0  # Any significant change

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

# Extend gene bodies to include promoter (upstream region)
get_extended_gene_regions <- function(genes, upstream = 2000, downstream = 0) {
  # Extend upstream of TSS and downstream of TES
  promoters_and_bodies <- promoters(genes, upstream = upstream, downstream = width(genes) + downstream)
  promoters_and_bodies$gene_id <- genes$gene_id
  promoters_and_bodies$symbol <- genes$symbol
  promoters_and_bodies$ensembl_id <- genes$ensembl_id
  return(promoters_and_bodies)
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

cat("=== Module 6: Methylation at Regulated Genes ===\n\n")
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

# 2. Load DMRs
cat("\n2. Loading DMRs...\n")

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

# Classify DMRs by direction
dmrs$dmr_direction <- ifelse(dmrs$logFC > 0, "Hypermethylated", "Hypomethylated")
cat(sprintf("   - Hypermethylated (TES > GFP): %d\n", sum(dmrs$dmr_direction == "Hypermethylated")))
cat(sprintf("   - Hypomethylated (TES < GFP): %d\n", sum(dmrs$dmr_direction == "Hypomethylated")))

# 3. Get gene body regions
cat("\n3. Getting gene body regions...\n")

gene_bodies <- get_gene_bodies(txdb)
cat(sprintf("   Total genes with coordinates: %d\n", length(gene_bodies)))

# Filter to standard chromosomes
standard_chr <- paste0("chr", c(1:22, "X", "Y"))
gene_bodies <- keepSeqlevels(gene_bodies, standard_chr[standard_chr %in% seqlevels(gene_bodies)],
                              pruning.mode = "coarse")

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

# 4. Find DMRs overlapping gene bodies
cat("\n4. Finding DMRs at gene bodies...\n")

# Find overlaps
hits <- findOverlaps(gene_bodies, dmrs)

# Count DMRs per gene
dmr_counts <- data.frame(
  gene_idx = seq_along(gene_bodies),
  ensembl_id = gene_bodies$ensembl_id,
  symbol = gene_bodies$symbol,
  is_deg = gene_bodies$is_deg,
  expr_direction = gene_bodies$expr_direction,
  gene_width = width(gene_bodies)
)

# Count overlapping DMRs
dmr_per_gene <- as.data.frame(table(queryHits(hits)))
colnames(dmr_per_gene) <- c("gene_idx", "n_dmrs")
dmr_per_gene$gene_idx <- as.integer(as.character(dmr_per_gene$gene_idx))

dmr_counts <- merge(dmr_counts, dmr_per_gene, by = "gene_idx", all.x = TRUE)
dmr_counts$n_dmrs[is.na(dmr_counts$n_dmrs)] <- 0
dmr_counts$has_dmr <- dmr_counts$n_dmrs > 0

# Normalize by gene length (DMRs per Mb)
dmr_counts$dmr_density <- dmr_counts$n_dmrs / (dmr_counts$gene_width / 1e6)

# 5. Compare DEGs vs non-DEGs
cat("\n5. Comparing DMR enrichment at DEGs vs non-DEGs...\n")

# Summary statistics
deg_stats <- dmr_counts %>%
  filter(!is.na(is_deg)) %>%
  group_by(is_deg) %>%
  summarise(
    n_genes = n(),
    genes_with_dmr = sum(has_dmr),
    pct_with_dmr = mean(has_dmr) * 100,
    mean_dmr_count = mean(n_dmrs),
    mean_dmr_density = mean(dmr_density, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n   DMR enrichment summary:\n")
print(deg_stats)

# Statistical test: Are DEGs more likely to have DMRs?
contingency <- matrix(c(
  sum(dmr_counts$is_deg & dmr_counts$has_dmr, na.rm = TRUE),
  sum(dmr_counts$is_deg & !dmr_counts$has_dmr, na.rm = TRUE),
  sum(!dmr_counts$is_deg & dmr_counts$has_dmr, na.rm = TRUE),
  sum(!dmr_counts$is_deg & !dmr_counts$has_dmr, na.rm = TRUE)
), nrow = 2, byrow = TRUE,
dimnames = list(c("DEG", "Non-DEG"), c("Has_DMR", "No_DMR")))

cat("\n   Contingency table:\n")
print(contingency)

fisher_test <- fisher.test(contingency)
cat(sprintf("\n   Fisher's exact test p-value: %.2e\n", fisher_test$p.value))
cat(sprintf("   Odds ratio: %.2f\n", fisher_test$estimate))

# 6. Direction analysis: Do methylation changes correlate with expression changes?
cat("\n6. Analyzing correlation between methylation and expression direction...\n")

# For DEGs with DMRs, get the predominant DMR direction
deg_with_dmr_idx <- queryHits(hits)[gene_bodies$is_deg[queryHits(hits)]]
dmr_at_deg_idx <- subjectHits(hits)[gene_bodies$is_deg[queryHits(hits)]]

if (length(deg_with_dmr_idx) > 0) {
  direction_df <- data.frame(
    gene_idx = deg_with_dmr_idx,
    ensembl_id = gene_bodies$ensembl_id[deg_with_dmr_idx],
    symbol = gene_bodies$symbol[deg_with_dmr_idx],
    expr_direction = gene_bodies$expr_direction[deg_with_dmr_idx],
    dmr_logFC = dmrs$logFC[dmr_at_deg_idx],
    dmr_direction = dmrs$dmr_direction[dmr_at_deg_idx]
  )

  # Summarize by gene (take mean logFC if multiple DMRs)
  direction_summary <- direction_df %>%
    group_by(ensembl_id, symbol, expr_direction) %>%
    summarise(
      n_dmrs = n(),
      mean_dmr_logFC = mean(dmr_logFC),
      predominant_dmr_direction = ifelse(mean(dmr_logFC) > 0, "Hypermethylated", "Hypomethylated"),
      .groups = "drop"
    )

  # Cross-tabulation
  direction_table <- table(direction_summary$expr_direction, direction_summary$predominant_dmr_direction)
  cat("\n   Expression vs Methylation direction:\n")
  print(direction_table)

  # Expected pattern if methylation represses:
  # Hypermethylated -> Downregulated
  # Hypomethylated -> Upregulated
  if (nrow(direction_table) >= 2 && ncol(direction_table) >= 2) {
    concordant <- direction_table["DOWN", "Hypermethylated"] + direction_table["UP", "Hypomethylated"]
    discordant <- direction_table["UP", "Hypermethylated"] + direction_table["DOWN", "Hypomethylated"]
    total <- concordant + discordant

    cat(sprintf("\n   Concordant (Hyper->Down OR Hypo->Up): %d (%.1f%%)\n",
                concordant, 100*concordant/total))
    cat(sprintf("   Discordant: %d (%.1f%%)\n", discordant, 100*discordant/total))

    # Binomial test for concordance
    binom_test <- binom.test(concordant, total, p = 0.5)
    cat(sprintf("   Binomial test p-value: %.3f\n", binom_test$p.value))
  }

  # Save direction summary
  write.csv(direction_summary,
            file.path(OUTPUT_DIR, "DEGs_with_DMRs_direction.csv"),
            row.names = FALSE)
}

# 7. Separate analysis by expression direction
cat("\n7. Analyzing UP vs DOWN regulated genes separately...\n")

direction_stats <- dmr_counts %>%
  filter(is_deg) %>%
  group_by(expr_direction) %>%
  summarise(
    n_genes = n(),
    genes_with_dmr = sum(has_dmr),
    pct_with_dmr = mean(has_dmr) * 100,
    mean_dmr_count = mean(n_dmrs),
    .groups = "drop"
  )

cat("\n   DMR enrichment by expression direction:\n")
print(direction_stats)

# 8. Visualization
cat("\n8. Generating visualizations...\n")

# Plot 1: DMR enrichment at DEGs vs non-DEGs
pdf(file.path(OUTPUT_DIR, "01_DMR_enrichment_DEG_vs_nonDEG.pdf"), width = 10, height = 6)

p1_data <- dmr_counts %>%
  filter(!is.na(is_deg)) %>%
  mutate(gene_type = ifelse(is_deg, "DEG", "Non-DEG"))

p1 <- ggplot(p1_data, aes(x = gene_type, y = n_dmrs)) +
  geom_violin(aes(fill = gene_type), alpha = 0.7) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  scale_fill_manual(values = c("DEG" = "#E41A1C", "Non-DEG" = "#377EB8")) +
  theme_minimal(base_size = 14) +
  labs(
    title = "DMR Count at Gene Bodies: DEGs vs Non-DEGs",
    subtitle = sprintf("Fisher's test p = %.2e, OR = %.2f", fisher_test$p.value, fisher_test$estimate),
    x = "",
    y = "Number of DMRs per Gene"
  ) +
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(0, quantile(p1_data$n_dmrs, 0.99)))
print(p1)

dev.off()

# Plot 2: Percentage of genes with DMRs
pdf(file.path(OUTPUT_DIR, "02_Pct_genes_with_DMRs.pdf"), width = 8, height = 6)

p2_data <- deg_stats %>%
  mutate(gene_type = ifelse(is_deg, "DEG", "Non-DEG"))

p2 <- ggplot(p2_data, aes(x = gene_type, y = pct_with_dmr, fill = gene_type)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%\n(%d/%d)", pct_with_dmr, genes_with_dmr, n_genes)),
            vjust = -0.3, size = 4) +
  scale_fill_manual(values = c("DEG" = "#E41A1C", "Non-DEG" = "#377EB8")) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Percentage of Genes Containing DMRs",
    subtitle = "Gene bodies (TSS to TES)",
    x = "",
    y = "% of Genes with ≥1 DMR"
  ) +
  theme(legend.position = "none") +
  ylim(0, max(p2_data$pct_with_dmr) * 1.3)
print(p2)

dev.off()

# Plot 3: Expression direction breakdown
pdf(file.path(OUTPUT_DIR, "03_DMR_by_expression_direction.pdf"), width = 10, height = 6)

p3_data <- dmr_counts %>%
  filter(is_deg) %>%
  group_by(expr_direction) %>%
  summarise(
    pct_with_dmr = mean(has_dmr) * 100,
    n_genes = n(),
    genes_with_dmr = sum(has_dmr),
    .groups = "drop"
  )

p3 <- ggplot(p3_data, aes(x = expr_direction, y = pct_with_dmr, fill = expr_direction)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%\n(%d/%d)", pct_with_dmr, genes_with_dmr, n_genes)),
            vjust = -0.3, size = 4) +
  scale_fill_manual(values = c("UP" = "#4DAF4A", "DOWN" = "#984EA3")) +
  theme_minimal(base_size = 14) +
  labs(
    title = "DMR Presence at Upregulated vs Downregulated Genes",
    x = "Expression Direction (TES vs GFP)",
    y = "% of DEGs with ≥1 DMR in Gene Body"
  ) +
  theme(legend.position = "none") +
  ylim(0, max(p3_data$pct_with_dmr) * 1.3)
print(p3)

dev.off()

# Plot 4: Methylation-Expression correlation heatmap
if (exists("direction_summary") && nrow(direction_summary) > 0) {
  pdf(file.path(OUTPUT_DIR, "04_Methylation_Expression_correlation.pdf"), width = 8, height = 6)

  p4_data <- direction_summary %>%
    count(expr_direction, predominant_dmr_direction) %>%
    complete(expr_direction, predominant_dmr_direction, fill = list(n = 0))

  p4 <- ggplot(p4_data, aes(x = predominant_dmr_direction, y = expr_direction, fill = n)) +
    geom_tile(color = "white", size = 1) +
    geom_text(aes(label = n), size = 8, fontface = "bold") +
    scale_fill_gradient(low = "white", high = "#E41A1C", name = "Gene Count") +
    theme_minimal(base_size = 14) +
    labs(
      title = "Methylation vs Expression Direction at DEGs",
      subtitle = "Expected if methylation represses: Hyper→Down, Hypo→Up",
      x = "Methylation Change Direction",
      y = "Expression Change Direction"
    ) +
    theme(panel.grid = element_blank())
  print(p4)

  dev.off()
}

# Plot 5: DMR density (DMRs per Mb of gene)
pdf(file.path(OUTPUT_DIR, "05_DMR_density_comparison.pdf"), width = 10, height = 6)

p5_data <- dmr_counts %>%
  filter(!is.na(is_deg), dmr_density < Inf) %>%
  mutate(gene_type = ifelse(is_deg, "DEG", "Non-DEG"))

p5 <- ggplot(p5_data, aes(x = gene_type, y = dmr_density, fill = gene_type)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  scale_fill_manual(values = c("DEG" = "#E41A1C", "Non-DEG" = "#377EB8")) +
  scale_y_log10() +
  theme_minimal(base_size = 14) +
  labs(
    title = "DMR Density at Gene Bodies",
    subtitle = "Normalized by gene length",
    x = "",
    y = "DMRs per Megabase (log10)"
  ) +
  theme(legend.position = "none")
print(p5)

dev.off()

# 9. Save DMR results
cat("\n9. Saving DMR results...\n")

write.csv(dmr_counts, file.path(OUTPUT_DIR, "Gene_DMR_counts.csv"), row.names = FALSE)
write.csv(deg_stats, file.path(OUTPUT_DIR, "DEG_vs_nonDEG_DMR_stats.csv"), row.names = FALSE)

# =====================================================================
# 10. BigWig Methylation Signal Analysis at UP vs DOWN Regulated Genes
# =====================================================================

cat("\n10. Analyzing methylation signal from BigWig files...\n")

# Check if BigWig files exist
tes_bw_exist <- all(file.exists(TES_BIGWIGS))
gfp_bw_exist <- all(file.exists(GFP_BIGWIGS))

if (tes_bw_exist && gfp_bw_exist) {

  # Get DEG gene bodies for signal extraction
  deg_gene_bodies <- gene_bodies[gene_bodies$is_deg]
  cat(sprintf("  DEG gene bodies to analyze: %d\n", length(deg_gene_bodies)))

  # Extract methylation signal from TES and GFP BigWigs
  cat("\n  Extracting TES methylation signal...\n")
  deg_gene_bodies$tes_meth_signal <- extract_bigwig_signal(
    deg_gene_bodies, TES_BIGWIGS, "TES methylation"
  )

  cat("\n  Extracting GFP methylation signal...\n")
  deg_gene_bodies$gfp_meth_signal <- extract_bigwig_signal(
    deg_gene_bodies, GFP_BIGWIGS, "GFP methylation"
  )

  # Calculate methylation difference (TES - GFP)
  deg_gene_bodies$meth_diff <- deg_gene_bodies$tes_meth_signal - deg_gene_bodies$gfp_meth_signal

  # Create data frame for analysis
  meth_signal_df <- data.frame(
    ensembl_id = deg_gene_bodies$ensembl_id,
    symbol = deg_gene_bodies$symbol,
    expr_direction = deg_gene_bodies$expr_direction,
    tes_meth = deg_gene_bodies$tes_meth_signal,
    gfp_meth = deg_gene_bodies$gfp_meth_signal,
    meth_diff = deg_gene_bodies$meth_diff,
    stringsAsFactors = FALSE
  )

  # Filter for valid signals
  meth_signal_df <- meth_signal_df[!is.na(meth_signal_df$meth_diff), ]
  cat(sprintf("\n  DEGs with valid methylation signal: %d\n", nrow(meth_signal_df)))

  # 10.1 Summary statistics by expression direction
  cat("\n  10.1 Methylation signal summary by expression direction:\n")

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

  # 10.2 Statistical tests
  cat("\n  10.2 Statistical tests:\n")

  up_meth <- meth_signal_df$meth_diff[meth_signal_df$expr_direction == "UP"]
  down_meth <- meth_signal_df$meth_diff[meth_signal_df$expr_direction == "DOWN"]

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
  }

  # 10.3 Visualizations
  cat("\n  10.3 Generating methylation signal plots...\n")

  # Plot 6: Methylation signal boxplot by expression direction
  pdf(file.path(OUTPUT_DIR, "06_Methylation_signal_by_direction.pdf"), width = 12, height = 6)

  # Panel A: TES methylation
  p6a <- ggplot(meth_signal_df, aes(x = expr_direction, y = tes_meth, fill = expr_direction)) +
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
  p6b <- ggplot(meth_signal_df, aes(x = expr_direction, y = gfp_meth, fill = expr_direction)) +
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
  gridExtra::grid.arrange(p6a, p6b, ncol = 2,
                          top = "Methylation Signal at Gene Bodies of DEGs")
  dev.off()

  # Plot 7: Methylation difference (TES - GFP) by expression direction
  pdf(file.path(OUTPUT_DIR, "07_Methylation_difference_by_direction.pdf"), width = 10, height = 6)

  wilcox_label <- if (exists("wilcox_test")) {
    sprintf("Wilcoxon p = %.2e", wilcox_test$p.value)
  } else {
    ""
  }

  p7 <- ggplot(meth_signal_df, aes(x = expr_direction, y = meth_diff, fill = expr_direction)) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.2, outlier.shape = NA) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    scale_fill_manual(values = c("UP" = "#4DAF4A", "DOWN" = "#984EA3")) +
    theme_minimal(base_size = 14) +
    labs(
      title = "Methylation Change (TES - GFP) at Gene Bodies",
      subtitle = wilcox_label,
      x = "Expression Direction",
      y = "Methylation Difference (TES - GFP)"
    ) +
    theme(legend.position = "none") +
    coord_cartesian(ylim = quantile(meth_signal_df$meth_diff, c(0.01, 0.99), na.rm = TRUE))
  print(p7)

  dev.off()

  # Plot 8: Scatter plot of expression vs methylation change
  pdf(file.path(OUTPUT_DIR, "08_Expression_vs_Methylation_scatter.pdf"), width = 10, height = 8)

  # Add log2FC to meth_signal_df
  rna_fc <- setNames(rna_df$log2FoldChange, rna_df$ensembl_clean)
  meth_signal_df$log2FC <- rna_fc[meth_signal_df$ensembl_id]

  # Calculate correlation
  cor_test <- cor.test(meth_signal_df$log2FC, meth_signal_df$meth_diff,
                       method = "spearman", use = "complete.obs")

  p8 <- ggplot(meth_signal_df, aes(x = log2FC, y = meth_diff, color = expr_direction)) +
    geom_point(alpha = 0.5, size = 2) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
    geom_hline(yintercept = 0, linetype = "dotted", color = "gray50") +
    geom_vline(xintercept = 0, linetype = "dotted", color = "gray50") +
    scale_color_manual(values = c("UP" = "#4DAF4A", "DOWN" = "#984EA3")) +
    theme_minimal(base_size = 14) +
    labs(
      title = "Expression Change vs Methylation Change",
      subtitle = sprintf("Spearman rho = %.3f, p = %.2e", cor_test$estimate, cor_test$p.value),
      x = "log2 Fold Change (Expression)",
      y = "Methylation Difference (TES - GFP)",
      color = "Direction"
    ) +
    theme(legend.position = "right")
  print(p8)

  dev.off()

  # Plot 9: Density plot of methylation difference
  pdf(file.path(OUTPUT_DIR, "09_Methylation_difference_density.pdf"), width = 10, height = 6)

  p9 <- ggplot(meth_signal_df, aes(x = meth_diff, fill = expr_direction, color = expr_direction)) +
    geom_density(alpha = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray30") +
    scale_fill_manual(values = c("UP" = "#4DAF4A", "DOWN" = "#984EA3")) +
    scale_color_manual(values = c("UP" = "#4DAF4A", "DOWN" = "#984EA3")) +
    theme_minimal(base_size = 14) +
    labs(
      title = "Distribution of Methylation Change at DEG Gene Bodies",
      x = "Methylation Difference (TES - GFP)",
      y = "Density",
      fill = "Expression", color = "Expression"
    ) +
    theme(legend.position = "right") +
    coord_cartesian(xlim = quantile(meth_signal_df$meth_diff, c(0.01, 0.99), na.rm = TRUE))
  print(p9)

  dev.off()

  # Save methylation signal data
  write.csv(meth_signal_df, file.path(OUTPUT_DIR, "DEG_methylation_signal.csv"), row.names = FALSE)
  write.csv(meth_summary, file.path(OUTPUT_DIR, "Methylation_signal_summary.csv"), row.names = FALSE)

  # Store results for summary
  bigwig_analysis_done <- TRUE

} else {
  cat("  WARNING: BigWig files not found. Skipping signal analysis.\n")
  cat(sprintf("    TES BigWigs exist: %s\n", tes_bw_exist))
  cat(sprintf("    GFP BigWigs exist: %s\n", gfp_bw_exist))
  bigwig_analysis_done <- FALSE
}

# 11. Save final results and generate summary

# Build summary report
summary_text <- sprintf("
=== Methylation at Regulated Genes Analysis ===
Date: %s

Input Data:
  RNA-seq: %s
  DMRs: %s

Gene Counts:
  Total genes analyzed: %d
  DEGs (padj < %.2f): %d
    - Upregulated: %d
    - Downregulated: %d

DMR Statistics:
  Total DMRs: %d
    - Hypermethylated: %d
    - Hypomethylated: %d

DMR Enrichment at DEGs:
  DEGs with ≥1 DMR: %d / %d (%.1f%%)
  Non-DEGs with ≥1 DMR: %d / %d (%.1f%%)

  Fisher's exact test:
    p-value: %.2e
    Odds ratio: %.2f

  Interpretation: %s
",
Sys.time(),
RNA_SEQ_FILE,
DMR_FILE,
nrow(dmr_counts),
PADJ_THRESHOLD,
n_degs, n_up, n_down,
length(dmrs),
sum(dmrs$dmr_direction == "Hypermethylated"),
sum(dmrs$dmr_direction == "Hypomethylated"),
deg_stats$genes_with_dmr[deg_stats$is_deg == TRUE],
deg_stats$n_genes[deg_stats$is_deg == TRUE],
deg_stats$pct_with_dmr[deg_stats$is_deg == TRUE],
deg_stats$genes_with_dmr[deg_stats$is_deg == FALSE],
deg_stats$n_genes[deg_stats$is_deg == FALSE],
deg_stats$pct_with_dmr[deg_stats$is_deg == FALSE],
fisher_test$p.value,
fisher_test$estimate,
ifelse(fisher_test$p.value < 0.05 & fisher_test$estimate > 1,
       "DEGs are ENRICHED for DMRs - supports indirect methylation hypothesis",
       ifelse(fisher_test$p.value < 0.05 & fisher_test$estimate < 1,
              "DEGs are DEPLETED for DMRs",
              "No significant enrichment/depletion"))
)

# Add BigWig analysis section if it was performed
if (bigwig_analysis_done && exists("meth_signal_df") && exists("meth_summary")) {
  bigwig_summary <- sprintf("
BigWig Methylation Signal Analysis:
  DEGs with valid signal: %d

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
  if (exists("wilcox_test")) {
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
  DMR Analysis:
    - Gene_DMR_counts.csv: All genes with DMR counts
    - DEG_vs_nonDEG_DMR_stats.csv: Summary statistics
    - DEGs_with_DMRs_direction.csv: Direction correlation analysis
    - 01-05_*.pdf: DMR visualization plots
"

if (bigwig_analysis_done) {
  output_files_text <- paste0(output_files_text, "
  BigWig Signal Analysis:
    - DEG_methylation_signal.csv: Methylation signal per DEG
    - Methylation_signal_summary.csv: Summary by direction
    - 06_Methylation_signal_by_direction.pdf: Signal boxplots
    - 07_Methylation_difference_by_direction.pdf: TES-GFP difference
    - 08_Expression_vs_Methylation_scatter.pdf: Correlation plot
    - 09_Methylation_difference_density.pdf: Density distributions
")
}

summary_text <- paste0(summary_text, output_files_text)

writeLines(summary_text, file.path(OUTPUT_DIR, "ANALYSIS_SUMMARY.txt"))
cat(summary_text)

cat("\n=== Module 6 Complete ===\n")
cat(sprintf("Results saved to: %s\n", OUTPUT_DIR))
