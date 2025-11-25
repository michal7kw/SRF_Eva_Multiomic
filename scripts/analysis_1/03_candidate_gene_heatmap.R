#!/usr/bin/env Rscript
#
# CANDIDATE GENE HEATMAP: TES_degs.txt genes vs random control
# Addresses Todo #6: Heatmap for specific gene list vs random genes
# Compares TES/TEAD1 binding at 223 candidate genes

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(GenomicRanges)
  library(rtracklayer)
  library(GenomicFeatures)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
})

# Fix namespace conflicts: GenomicFeatures loads AnnotationDbi which masks dplyr functions
# Also clusterProfiler masks stats::filter
# Explicitly bind dplyr functions to avoid conflicts
select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate
arrange <- dplyr::arrange

setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_1")

cat("=== CANDIDATE GENE HEATMAP ANALYSIS ===\n")
cat("Comparing TES/TEAD1 binding at 223 candidate genes vs random controls\n")
cat("Analysis started:", as.character(Sys.time()), "\n\n")

# Create output directory
output_dir <- "output/03_candidate_gene_heatmap"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# PHASE 1: LOAD CANDIDATE GENES
# =============================================================================

cat("=== PHASE 1: Loading Candidate Gene List ===\n")

# Load candidate genes
candidate_genes <- read.table("../../data/TES_degs.txt",
  stringsAsFactors = FALSE,
  header = FALSE
)
colnames(candidate_genes) <- "gene_symbol"

cat(sprintf("✓ Loaded %d candidate genes from TES_degs.txt\n", nrow(candidate_genes)))
cat(sprintf("  First 10: %s\n\n", paste(head(candidate_genes$gene_symbol, 10), collapse = ", ")))

# Load RNA-seq results to get expression info
rna_results <- read.delim("../../../SRF_Eva_RNA/results/05_deseq2/deseq2_results_TES_vs_GFP.txt",
  stringsAsFactors = FALSE
)

# Convert Ensembl to gene symbols
library(org.Hs.eg.db)
rna_results$ensembl_id <- gsub("\\..*", "", rna_results$gene_id)
rna_results$gene_symbol <- mapIds(org.Hs.eg.db,
  keys = rna_results$ensembl_id,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

# Match candidate genes to RNA-seq data
candidate_with_data <- rna_results[rna_results$gene_symbol %in% candidate_genes$gene_symbol, ]

cat(sprintf("✓ Matched %d candidate genes to RNA-seq data\n", nrow(candidate_with_data)))

# Check significance
candidate_sig <- candidate_with_data[!is.na(candidate_with_data$padj) &
  candidate_with_data$padj < 0.05, ]
cat(sprintf("  - Significant DEGs: %d\n", nrow(candidate_sig)))
cat(sprintf("  - Non-significant: %d\n\n", nrow(candidate_with_data) - nrow(candidate_sig)))

# =============================================================================
# PHASE 2: SELECT MATCHED RANDOM CONTROLS
# =============================================================================

cat("=== PHASE 2: Selecting Matched Random Controls ===\n")

# Match random genes by:
# 1. Expression level (baseMean)
# 2. Chromosome distribution
# 3. Gene length (approximate by using similar genes)

# Load GTF for gene information
gtf_file <- "/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/annotation/gencode.v44.annotation.gtf"

# Try to use makeTxDbFromGFF (preferred method - creates proper TxDb object)
# Falls back to rtracklayer if txdbmaker is not installed
annotation_result <- tryCatch(
  {
    cat("Loading gene annotations using makeTxDbFromGFF (preferred method)...\n")
    txdb <- makeTxDbFromGFF(gtf_file, format = "gtf")
    list(genes_gr = genes(txdb), txdb = txdb)
  },
  error = function(e) {
    cat("Note: txdbmaker not available, using rtracklayer fallback method\n")
    cat("  To use preferred method, install: BiocManager::install('txdbmaker')\n")
    gtf <- rtracklayer::import(gtf_file)
    genes_subset <- gtf[gtf$type == "gene"]
    names(genes_subset) <- genes_subset$gene_id
    list(genes_gr = genes_subset, txdb = NULL)
  }
)

genes_gr <- annotation_result$genes_gr
txdb <- annotation_result$txdb

cat(sprintf("Loaded %d gene annotations\n", length(genes_gr)))

# Function to select matched random genes
select_matched_random <- function(target_genes_data, all_genes_data, n_controls = 1, seed = 42) {
  set.seed(seed)

  random_genes <- data.frame()

  for (i in 1:nrow(target_genes_data)) {
    target_gene <- target_genes_data[i, ]

    # Find genes with similar expression level
    expression_range <- c(target_gene$baseMean * 0.5, target_gene$baseMean * 2)

    similar_genes <- all_genes_data %>%
      filter(baseMean >= expression_range[1] & baseMean <= expression_range[2]) %>%
      filter(!gene_symbol %in% target_genes_data$gene_symbol) %>% # Exclude candidates
      filter(!gene_symbol %in% random_genes$gene_symbol) # Avoid duplicates

    if (nrow(similar_genes) > 0) {
      # Sample random gene
      selected <- similar_genes[sample(nrow(similar_genes), min(n_controls, nrow(similar_genes))), ]
      random_genes <- rbind(random_genes, selected)
    }
  }

  return(random_genes)
}

cat("Selecting matched random control genes...\n")
random_controls <- select_matched_random(candidate_with_data,
  rna_results[!is.na(rna_results$gene_symbol), ],
  n_controls = 1
)

cat(sprintf("✓ Selected %d random control genes\n", nrow(random_controls)))
cat(sprintf(
  "  Expression match: %.1f%% within 2-fold range\n\n",
  100 * sum(abs(log2(random_controls$baseMean / candidate_with_data$baseMean[1:nrow(random_controls)])) < 1) /
    nrow(random_controls)
))

# =============================================================================
# PHASE 3: EXTRACT GENE COORDINATES AND CREATE BED FILES
# =============================================================================

cat("=== PHASE 3: Creating BED Files for deepTools ===\n")

# Function to create BED file from gene list
create_gene_bed <- function(gene_data, output_file, txdb_obj = NULL, genes_granges = NULL) {
  # Clean Ensembl IDs
  gene_data$ensembl_id_clean <- gsub("\\..*", "", gene_data$ensembl_id)

  # Get gene coordinates - use either txdb or genes_granges
  if (!is.null(txdb_obj)) {
    genes_gr_all <- genes(txdb_obj)
  } else if (!is.null(genes_granges)) {
    genes_gr_all <- genes_granges
  } else {
    cat("Error: Neither txdb_obj nor genes_granges provided\n")
    return(NULL)
  }

  # Clean names to match Ensembl IDs (remove version suffixes like .16 from ENSG00000000003.16)
  names(genes_gr_all) <- gsub("\\..*", "", names(genes_gr_all))

  matched_genes <- genes_gr_all[names(genes_gr_all) %in% gene_data$ensembl_id_clean]

  if (length(matched_genes) == 0) {
    cat("Warning: No genes matched for", basename(output_file), "\n")
    cat("  Searched for", length(gene_data$ensembl_id_clean), "Ensembl IDs\n")
    cat("  Available gene names:", length(unique(names(genes_gr_all))), "\n")
    cat("  Example searched IDs:", paste(head(gene_data$ensembl_id_clean, 3), collapse = ", "), "\n")
    cat("  Example available IDs:", paste(head(unique(names(genes_gr_all)), 3), collapse = ", "), "\n")
    return(NULL)
  }

  # Extend to promoter region (TSS +/- 2kb)
  promoters_gr <- promoters(matched_genes, upstream = 2000, downstream = 2000)

  # Create BED format
  bed_df <- data.frame(
    chr = as.character(seqnames(promoters_gr)),
    start = start(promoters_gr),
    end = end(promoters_gr),
    name = names(promoters_gr),
    score = 0,
    strand = as.character(strand(promoters_gr))
  )

  # Sort
  bed_df <- bed_df %>% arrange(chr, start)

  # Write BED
  write.table(bed_df, output_file,
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
  )

  cat(sprintf("✓ Created %s with %d regions\n", basename(output_file), nrow(bed_df)))
  return(nrow(bed_df))
}

# Create BED files - use txdb if available, otherwise use genes_gr
if (!is.null(txdb)) {
  create_gene_bed(candidate_with_data,
    file.path(output_dir, "candidate_genes.bed"),
    txdb_obj = txdb
  )

  create_gene_bed(random_controls,
    file.path(output_dir, "random_control_genes.bed"),
    txdb_obj = txdb
  )
} else {
  create_gene_bed(candidate_with_data,
    file.path(output_dir, "candidate_genes.bed"),
    genes_granges = genes_gr
  )

  create_gene_bed(random_controls,
    file.path(output_dir, "random_control_genes.bed"),
    genes_granges = genes_gr
  )
}

cat("\n")

# =============================================================================
# PHASE 4: LOAD PEAK DATA FOR CANDIDATE GENES
# =============================================================================

cat("=== PHASE 4: Analyzing Peak Coverage at Candidate Genes ===\n")

# Load peak annotations
tes_peaks <- read.csv("../../../SRF_Eva_CUTandTAG/results/07_analysis_narrow/TES_peaks_annotated.csv",
  stringsAsFactors = FALSE
)
tead1_peaks <- read.csv("../../../SRF_Eva_CUTandTAG/results/07_analysis_narrow/TEAD1_peaks_annotated.csv",
  stringsAsFactors = FALSE
)

# Convert Entrez to Ensembl for matching
tes_peaks$ensembl_id <- mapIds(org.Hs.eg.db,
  keys = as.character(tes_peaks$geneId),
  column = "ENSEMBL",
  keytype = "ENTREZID",
  multiVals = "first"
)

tead1_peaks$ensembl_id <- mapIds(org.Hs.eg.db,
  keys = as.character(tead1_peaks$geneId),
  column = "ENSEMBL",
  keytype = "ENTREZID",
  multiVals = "first"
)

# Count peaks near candidate genes
candidate_with_data$tes_peaks <- sapply(candidate_with_data$ensembl_id, function(eid) {
  sum(tes_peaks$ensembl_id == eid, na.rm = TRUE)
})

candidate_with_data$tead1_peaks <- sapply(candidate_with_data$ensembl_id, function(eid) {
  sum(tead1_peaks$ensembl_id == eid, na.rm = TRUE)
})

random_controls$tes_peaks <- sapply(random_controls$ensembl_id, function(eid) {
  sum(tes_peaks$ensembl_id == eid, na.rm = TRUE)
})

random_controls$tead1_peaks <- sapply(random_controls$ensembl_id, function(eid) {
  sum(tead1_peaks$ensembl_id == eid, na.rm = TRUE)
})

cat(sprintf("Candidate genes:\n"))
cat(sprintf(
  "  - With TES peaks: %d (%.1f%%)\n",
  sum(candidate_with_data$tes_peaks > 0),
  100 * sum(candidate_with_data$tes_peaks > 0) / nrow(candidate_with_data)
))
cat(sprintf(
  "  - With TEAD1 peaks: %d (%.1f%%)\n",
  sum(candidate_with_data$tead1_peaks > 0),
  100 * sum(candidate_with_data$tead1_peaks > 0) / nrow(candidate_with_data)
))

cat(sprintf("\nRandom controls:\n"))
cat(sprintf(
  "  - With TES peaks: %d (%.1f%%)\n",
  sum(random_controls$tes_peaks > 0),
  100 * sum(random_controls$tes_peaks > 0) / nrow(random_controls)
))
cat(sprintf(
  "  - With TEAD1 peaks: %d (%.1f%%)\n\n",
  sum(random_controls$tead1_peaks > 0),
  100 * sum(random_controls$tead1_peaks > 0) / nrow(random_controls)
))

# =============================================================================
# PHASE 5: EXPORT RESULTS
# =============================================================================

cat("=== PHASE 5: Exporting Results ===\n")

write.csv(candidate_with_data, file.path(output_dir, "candidate_genes_with_peaks.csv"), row.names = FALSE)
write.csv(random_controls, file.path(output_dir, "random_control_genes_with_peaks.csv"), row.names = FALSE)

cat("✓ Gene lists with peak data exported\n\n")

# =============================================================================
# PHASE 6: STATISTICAL COMPARISON
# =============================================================================

cat("=== PHASE 6: Statistical Comparison ===\n")

# Chi-square test for binding enrichment
tes_contingency <- matrix(c(
  sum(candidate_with_data$tes_peaks > 0),
  sum(candidate_with_data$tes_peaks == 0),
  sum(random_controls$tes_peaks > 0),
  sum(random_controls$tes_peaks == 0)
), nrow = 2, byrow = TRUE)

tead1_contingency <- matrix(c(
  sum(candidate_with_data$tead1_peaks > 0),
  sum(candidate_with_data$tead1_peaks == 0),
  sum(random_controls$tead1_peaks > 0),
  sum(random_controls$tead1_peaks == 0)
), nrow = 2, byrow = TRUE)

tes_test <- chisq.test(tes_contingency)
tead1_test <- chisq.test(tead1_contingency)

cat("Chi-square tests for binding enrichment:\n")
cat(sprintf("  TES: X² = %.2f, p = %.2e\n", tes_test$statistic, tes_test$p.value))
cat(sprintf("  TEAD1: X² = %.2f, p = %.2e\n\n", tead1_test$statistic, tead1_test$p.value))

# =============================================================================
# PHASE 7: VISUALIZATIONS
# =============================================================================

cat("=== PHASE 7: Creating Visualizations ===\n")

# Plot 1: Binding comparison barplot
cat("Creating binding comparison plot...\n")

binding_data <- data.frame(
  group = rep(c("Candidate", "Random"), each = 2),
  tf = rep(c("TES", "TEAD1"), 2),
  percentage = c(
    100 * sum(candidate_with_data$tes_peaks > 0) / nrow(candidate_with_data),
    100 * sum(candidate_with_data$tead1_peaks > 0) / nrow(candidate_with_data),
    100 * sum(random_controls$tes_peaks > 0) / nrow(random_controls),
    100 * sum(random_controls$tead1_peaks > 0) / nrow(random_controls)
  )
)

pdf(file.path(output_dir, "01_binding_comparison.pdf"), width = 10, height = 7)
p1 <- ggplot(binding_data, aes(x = tf, y = percentage, fill = group)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", linewidth = 0.3) +
  labs(
    title = "TF Binding at Candidate vs Random Genes",
    subtitle = sprintf(
      "%d candidate genes vs %d matched controls",
      nrow(candidate_with_data), nrow(random_controls)
    ),
    x = "Transcription Factor",
    y = "Percentage of Genes with Binding",
    fill = "Gene Set"
  ) +
  scale_fill_manual(values = c("Candidate" = "#E31A1C", "Random" = "#1F78B4")) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5)
  ) +
  geom_text(aes(label = sprintf("%.1f%%", percentage)),
    position = position_dodge(width = 0.9), vjust = -0.5, size = 4
  )
print(p1)
dev.off()

# Plot 2: Peak count distribution
cat("Creating peak count distribution plot...\n")

peak_counts <- rbind(
  data.frame(gene_set = "Candidate", tf = "TES", count = candidate_with_data$tes_peaks),
  data.frame(gene_set = "Candidate", tf = "TEAD1", count = candidate_with_data$tead1_peaks),
  data.frame(gene_set = "Random", tf = "TES", count = random_controls$tes_peaks),
  data.frame(gene_set = "Random", tf = "TEAD1", count = random_controls$tead1_peaks)
)

pdf(file.path(output_dir, "02_peak_count_distribution.pdf"), width = 12, height = 7)
p2 <- ggplot(peak_counts, aes(x = count, fill = gene_set)) +
  geom_histogram(bins = 20, alpha = 0.7, position = "identity", color = "black", linewidth = 0.3) +
  facet_wrap(~tf, scales = "free_y") +
  labs(
    title = "Peak Count Distribution",
    x = "Number of Peaks per Gene",
    y = "Frequency",
    fill = "Gene Set"
  ) +
  scale_fill_manual(values = c("Candidate" = "#E31A1C", "Random" = "#1F78B4")) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5))
print(p2)
dev.off()

# Plot 3: Expression vs binding scatter
cat("Creating expression vs binding scatter plot...\n")

candidate_plot_data <- candidate_with_data %>%
  mutate(gene_set = "Candidate") %>%
  select(gene_symbol, log2FoldChange, tes_peaks, tead1_peaks, gene_set)

random_plot_data <- random_controls %>%
  mutate(gene_set = "Random") %>%
  select(gene_symbol, log2FoldChange, tes_peaks, tead1_peaks, gene_set)

combined_data <- rbind(candidate_plot_data, random_plot_data)

pdf(file.path(output_dir, "03_expression_vs_binding.pdf"), width = 14, height = 6)
par(mfrow = c(1, 2))

# TES binding vs expression
ggplot(combined_data, aes(x = log2FoldChange, y = tes_peaks, color = gene_set)) +
  geom_point(alpha = 0.6, size = 2) +
  labs(
    title = "TES Binding vs Expression Change",
    x = "Log2 Fold Change (TES vs GFP)",
    y = "Number of TES Peaks",
    color = "Gene Set"
  ) +
  scale_color_manual(values = c("Candidate" = "#E31A1C", "Random" = "#1F78B4")) +
  theme_bw(base_size = 12) -> p3a

# TEAD1 binding vs expression
ggplot(combined_data, aes(x = log2FoldChange, y = tead1_peaks, color = gene_set)) +
  geom_point(alpha = 0.6, size = 2) +
  labs(
    title = "TEAD1 Binding vs Expression Change",
    x = "Log2 Fold Change (TES vs GFP)",
    y = "Number of TEAD1 Peaks",
    color = "Gene Set"
  ) +
  scale_color_manual(values = c("Candidate" = "#E31A1C", "Random" = "#1F78B4")) +
  theme_bw(base_size = 12) -> p3b

gridExtra::grid.arrange(p3a, p3b, ncol = 2)
dev.off()

cat("✓ All visualizations created\n\n")

# =============================================================================
# PHASE 8: SUMMARY REPORT
# =============================================================================

cat("=== PHASE 8: Generating Summary Report ===\n")

summary_file <- file.path(output_dir, "CANDIDATE_GENE_ANALYSIS_SUMMARY.txt")
cat("CANDIDATE GENE HEATMAP ANALYSIS SUMMARY\n", file = summary_file)
cat("========================================\n\n", file = summary_file, append = TRUE)
cat(paste("Generated:", Sys.time(), "\n\n"), file = summary_file, append = TRUE)

cat("CANDIDATE GENES:\n", file = summary_file, append = TRUE)
cat(sprintf("  Total loaded: %d\n", nrow(candidate_genes)), file = summary_file, append = TRUE)
cat(sprintf("  Matched to RNA-seq: %d\n", nrow(candidate_with_data)), file = summary_file, append = TRUE)
cat(sprintf("  Significant DEGs: %d\n\n", nrow(candidate_sig)), file = summary_file, append = TRUE)

cat("RANDOM CONTROLS:\n", file = summary_file, append = TRUE)
cat(sprintf("  Total selected: %d\n", nrow(random_controls)), file = summary_file, append = TRUE)
cat("  Matching criteria: Expression level (2-fold range)\n\n", file = summary_file, append = TRUE)

cat("TF BINDING COMPARISON:\n", file = summary_file, append = TRUE)
cat(sprintf("Candidate genes:\n"), file = summary_file, append = TRUE)
cat(
  sprintf(
    "  - TES binding: %d genes (%.1f%%)\n",
    sum(candidate_with_data$tes_peaks > 0),
    100 * sum(candidate_with_data$tes_peaks > 0) / nrow(candidate_with_data)
  ),
  file = summary_file, append = TRUE
)
cat(
  sprintf(
    "  - TEAD1 binding: %d genes (%.1f%%)\n",
    sum(candidate_with_data$tead1_peaks > 0),
    100 * sum(candidate_with_data$tead1_peaks > 0) / nrow(candidate_with_data)
  ),
  file = summary_file, append = TRUE
)

cat(sprintf("\nRandom controls:\n"), file = summary_file, append = TRUE)
cat(
  sprintf(
    "  - TES binding: %d genes (%.1f%%)\n",
    sum(random_controls$tes_peaks > 0),
    100 * sum(random_controls$tes_peaks > 0) / nrow(random_controls)
  ),
  file = summary_file, append = TRUE
)
cat(
  sprintf(
    "  - TEAD1 binding: %d genes (%.1f%%)\n\n",
    sum(random_controls$tead1_peaks > 0),
    100 * sum(random_controls$tead1_peaks > 0) / nrow(random_controls)
  ),
  file = summary_file, append = TRUE
)

cat("STATISTICAL TESTS:\n", file = summary_file, append = TRUE)
cat(sprintf("  TES enrichment: X² = %.2f, p = %.2e\n", tes_test$statistic, tes_test$p.value),
  file = summary_file, append = TRUE
)
cat(sprintf("  TEAD1 enrichment: X² = %.2f, p = %.2e\n\n", tead1_test$statistic, tead1_test$p.value),
  file = summary_file, append = TRUE
)

if (tes_test$p.value < 0.05) {
  cat("→ TES binding is significantly enriched at candidate genes!\n", file = summary_file, append = TRUE)
}
if (tead1_test$p.value < 0.05) {
  cat("→ TEAD1 binding is significantly enriched at candidate genes!\n", file = summary_file, append = TRUE)
}

cat("\n\nNEXT STEPS FOR DEEPTOOLS HEATMAP:\n", file = summary_file, append = TRUE)
cat("----------------------------------\n", file = summary_file, append = TRUE)
cat("BED files created:\n", file = summary_file, append = TRUE)
cat("  - candidate_genes.bed\n", file = summary_file, append = TRUE)
cat("  - random_control_genes.bed\n\n", file = summary_file, append = TRUE)

cat("Run deepTools to generate heatmap:\n", file = summary_file, append = TRUE)
cat("  cd ../../../SRF_Eva_CUTandTAG\n", file = summary_file, append = TRUE)
cat("  See generate_candidate_heatmap.sh script\n", file = summary_file, append = TRUE)

cat("\n✓ Summary report saved\n\n")

cat("========================================\n")
cat("CANDIDATE GENE ANALYSIS COMPLETE\n")
cat("========================================\n")
cat("Completed:", as.character(Sys.time()), "\n")
cat(sprintf("Output directory: %s\n\n", output_dir))
cat("Key files:\n")
cat("  - candidate_genes_with_peaks.csv\n")
cat("  - random_control_genes_with_peaks.csv\n")
cat("  - candidate_genes.bed (for deepTools)\n")
cat("  - random_control_genes.bed (for deepTools)\n")
cat("  - 01_binding_comparison.pdf\n")
cat("  - 02_peak_count_distribution.pdf\n")
cat("  - 03_expression_vs_binding.pdf\n")
cat("  - CANDIDATE_GENE_ANALYSIS_SUMMARY.txt\n\n")
cat("To generate deepTools heatmap, use the BED files created above.\n")
