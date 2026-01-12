#!/usr/bin/env Rscript

# 04_gene_regulatory_logic_combined.R
# Module 4 VARIANT: Gene Regulatory Logic - ALL PEAKS COMBINED
#
# This version combines all peak categories (Shared, TES_Unique, TEAD1_Unique)
# into a single group to compare Promoter vs Enhancer/Distal expression changes.

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(dplyr)
})

# ===================== Configuration =====================

# Base directories - use absolute paths throughout
BASE_DIR <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top"
SCRIPT_DIR <- file.path(BASE_DIR, "SRF_Eva_integrated_analysis/scripts/analysis_2")

# Input: Peak Categories from Module 1
PEAK_DIR <- file.path(SCRIPT_DIR, "results/01_peak_classification")
SHARED_PEAKS <- file.path(PEAK_DIR, "Shared_TES_TEAD1.bed")
TES_UNIQUE <- file.path(PEAK_DIR, "TES_Unique.bed")
TEAD1_UNIQUE <- file.path(PEAK_DIR, "TEAD1_Unique.bed")

# Input: RNA-seq Data
RNA_SEQ_FILE <- file.path(BASE_DIR,
                          "SRF_Eva_RNA/results/05_deseq2/deseq2_results_TES_vs_GFP.txt")

# Output
OUTPUT_DIR <- file.path(SCRIPT_DIR, "results/04_gene_regulatory_logic")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Parameters
PADJ_THRESHOLD <- 0.05  # Significance threshold for DEGs

# Annotation Database
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# ===================== Functions =====================

load_peaks <- function(bed_file) {
  if (!file.exists(bed_file)) {
    stop(paste("File not found:", bed_file))
  }
  import(bed_file, format = "BED")
}

annotate_peaks <- function(peaks) {
  if (length(peaks) == 0) {
    warning("No peaks provided")
    return(data.frame())
  }

  anno <- annotatePeak(peaks,
                       TxDb = txdb,
                       tssRegion = c(-2000, 2000),
                       verbose = FALSE)
  df <- as.data.frame(anno)

  # Classify as Promoter or Enhancer/Distal
  df$Regulatory_Type <- ifelse(grepl("Promoter", df$annotation),
                               "Promoter", "Enhancer/Distal")

  df
}

# ===================== Main Analysis =====================

cat("=== Module 4 (Combined): Gene Regulatory Logic - All Peaks ===\n\n")

# 1. Load and Combine All Peaks
cat("Loading data...\n")

# Check if peak files exist
if (!file.exists(SHARED_PEAKS)) {
  stop(paste("Peak file not found:", SHARED_PEAKS,
             "\nPlease run Module 1 (01_peak_classification.R) first."))
}

peaks_shared <- load_peaks(SHARED_PEAKS)
peaks_tes <- load_peaks(TES_UNIQUE)
peaks_tead1 <- load_peaks(TEAD1_UNIQUE)

cat(sprintf("  Shared peaks: %d\n", length(peaks_shared)))
cat(sprintf("  TES_Unique peaks: %d\n", length(peaks_tes)))
cat(sprintf("  TEAD1_Unique peaks: %d\n", length(peaks_tead1)))

# COMBINE ALL PEAKS INTO ONE SET
all_peaks <- c(peaks_shared, peaks_tes, peaks_tead1)
cat(sprintf("  TOTAL combined peaks: %d\n", length(all_peaks)))

# Load RNA-seq data
if (!file.exists(RNA_SEQ_FILE)) {
  stop(paste("RNA-seq file not found:", RNA_SEQ_FILE))
}

rna_data <- read.table(RNA_SEQ_FILE, header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE)
cat(sprintf("  RNA-seq genes: %d\n", nrow(rna_data)))

# 2. Annotate All Peaks Together (Promoter vs Enhancer)
cat("\nAnnotating all peaks together...\n")

combined_df <- annotate_peaks(all_peaks)

cat(sprintf("  Total annotated peaks: %d\n", nrow(combined_df)))
cat(sprintf("  Promoter peaks: %d\n", sum(combined_df$Regulatory_Type == "Promoter")))
cat(sprintf("  Enhancer/Distal peaks: %d\n", sum(combined_df$Regulatory_Type == "Enhancer/Distal")))

# 3. Link to Expression
cat("\nLinking to gene expression...\n")

# Map Entrez IDs (from ChIPseeker) to Ensembl IDs (for RNA-seq merge)
ensembl_ids <- mapIds(org.Hs.eg.db,
                      keys = as.character(combined_df$geneId),
                      column = "ENSEMBL",
                      keytype = "ENTREZID",
                      multiVals = "first")
combined_df$ensembl_id <- ensembl_ids

cat(sprintf("  Mapped %d out of %d peaks to Ensembl IDs\n",
            sum(!is.na(combined_df$ensembl_id)),
            nrow(combined_df)))

# Prepare RNA-seq data for merge
if (!"gene_id" %in% colnames(rna_data)) {
  rna_data$gene_id <- rownames(rna_data)
}

# Clean version numbers if present (ENSG00000xxx.1 -> ENSG00000xxx)
rna_data$clean_id <- gsub("\\..*", "", rna_data$gene_id)

# Merge peaks with expression
merged_df <- merge(combined_df, rna_data,
                   by.x = "ensembl_id", by.y = "clean_id",
                   all.x = TRUE)

# Add significance status
merged_df$is_significant <- !is.na(merged_df$padj) &
                            merged_df$padj < PADJ_THRESHOLD

# Add direction
merged_df$Direction <- ifelse(is.na(merged_df$log2FoldChange), "Unknown",
                              ifelse(merged_df$log2FoldChange > 0, "Up", "Down"))

# Save merged data (all peaks combined)
write.csv(merged_df, file.path(OUTPUT_DIR, "Peaks_Linked_to_Expression_Combined.csv"),
          row.names = FALSE)

cat(sprintf("  Total peaks with expression data: %d\n",
            sum(!is.na(merged_df$log2FoldChange))))
cat(sprintf("  Peaks linked to significant DEGs (padj < %.2f): %d\n",
            PADJ_THRESHOLD, sum(merged_df$is_significant, na.rm = TRUE)))

# 4. Generate summary statistics
cat("\n=== Expression Summary by Regulatory Type (All Peaks Combined) ===\n")

summary_stats <- merged_df %>%
  filter(!is.na(log2FoldChange)) %>%
  group_by(Regulatory_Type) %>%
  summarise(
    n_peaks = n(),
    n_significant = sum(is_significant, na.rm = TRUE),
    pct_significant = 100 * n_significant / n_peaks,
    mean_log2FC = mean(log2FoldChange, na.rm = TRUE),
    median_log2FC = median(log2FoldChange, na.rm = TRUE),
    n_up = sum(Direction == "Up", na.rm = TRUE),
    n_down = sum(Direction == "Down", na.rm = TRUE),
    .groups = "drop"
  )

print(summary_stats)
write.csv(summary_stats, file.path(OUTPUT_DIR, "Expression_Summary_Stats_Combined.csv"),
          row.names = FALSE)

# 5. Statistical Test: Promoter vs Enhancer/Distal
cat("\n=== Statistical Comparison ===\n")
promoter_fc <- merged_df$log2FoldChange[merged_df$Regulatory_Type == "Promoter" &
                                         !is.na(merged_df$log2FoldChange)]
enhancer_fc <- merged_df$log2FoldChange[merged_df$Regulatory_Type == "Enhancer/Distal" &
                                         !is.na(merged_df$log2FoldChange)]

wilcox_test <- wilcox.test(promoter_fc, enhancer_fc)
cat(sprintf("Wilcoxon test (Promoter vs Enhancer/Distal):\n"))
cat(sprintf("  Promoter median log2FC: %.4f (n=%d)\n", median(promoter_fc), length(promoter_fc)))
cat(sprintf("  Enhancer/Distal median log2FC: %.4f (n=%d)\n", median(enhancer_fc), length(enhancer_fc)))
cat(sprintf("  p-value: %.2e\n", wilcox_test$p.value))

# 6. Visualization
cat("\nGenerating expression plots...\n")

# Filter for peaks with valid expression data
plot_data_all <- merged_df[!is.na(merged_df$log2FoldChange), ]
# Filter for significant DEGs only
plot_data_sig <- merged_df[merged_df$is_significant == TRUE, ]

cat(sprintf("  Plotting all genes: %d peaks\n", nrow(plot_data_all)))
cat(sprintf("  Plotting significant DEGs: %d peaks\n", nrow(plot_data_sig)))

# --- Plot 1: Promoter vs Enhancer (All genes) ---
pdf(file.path(OUTPUT_DIR, "Expression_Promoter_vs_Enhancer_AllPeaks.pdf"),
    width = 8, height = 6)
p1 <- ggplot(plot_data_all, aes(x = Regulatory_Type, y = log2FoldChange,
                                fill = Regulatory_Type)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.1, size = 0.5) +
  coord_cartesian(ylim = c(-2, 2)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  theme_minimal(base_size = 14) +
  labs(title = "Expression Changes: Promoter vs Enhancer Binding",
       subtitle = sprintf("All peaks combined (n=%d), Wilcoxon p=%.2e",
                         nrow(plot_data_all), wilcox_test$p.value),
       y = "Log2 Fold Change (TES vs GFP)",
       x = "Regulatory Type") +
  scale_fill_manual(values = c("Promoter" = "#E41A1C",
                               "Enhancer/Distal" = "#377EB8")) +
  theme(legend.position = "none")
print(p1)
dev.off()

# --- Plot 2: Promoter vs Enhancer (Significant DEGs only) ---
if (nrow(plot_data_sig) > 10) {
  # Stats for significant only
  promoter_fc_sig <- plot_data_sig$log2FoldChange[plot_data_sig$Regulatory_Type == "Promoter"]
  enhancer_fc_sig <- plot_data_sig$log2FoldChange[plot_data_sig$Regulatory_Type == "Enhancer/Distal"]

  if (length(promoter_fc_sig) > 0 && length(enhancer_fc_sig) > 0) {
    wilcox_sig <- wilcox.test(promoter_fc_sig, enhancer_fc_sig)
    pval_sig <- wilcox_sig$p.value
  } else {
    pval_sig <- NA
  }

  pdf(file.path(OUTPUT_DIR, "Expression_Promoter_vs_Enhancer_AllPeaks_Significant.pdf"),
      width = 8, height = 6)
  p2 <- ggplot(plot_data_sig, aes(x = Regulatory_Type, y = log2FoldChange,
                                  fill = Regulatory_Type)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
    coord_cartesian(ylim = quantile(plot_data_sig$log2FoldChange,
                                    c(0.01, 0.99), na.rm = TRUE)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    theme_minimal(base_size = 14) +
    labs(title = "Expression: Promoter vs Enhancer (Significant DEGs)",
         subtitle = sprintf("padj < %.2f (n=%d), Wilcoxon p=%.2e",
                           PADJ_THRESHOLD, nrow(plot_data_sig),
                           ifelse(is.na(pval_sig), 1, pval_sig)),
         y = "Log2 Fold Change (TES vs GFP)",
         x = "Regulatory Type") +
    scale_fill_manual(values = c("Promoter" = "#E41A1C",
                                 "Enhancer/Distal" = "#377EB8")) +
    theme(legend.position = "none")
  print(p2)
  dev.off()
} else {
  cat("  Too few significant DEGs to plot separately\n")
}

# --- Plot 3: Direction barplot ---
direction_counts <- plot_data_sig %>%
  filter(Direction %in% c("Up", "Down")) %>%
  group_by(Regulatory_Type, Direction) %>%
  summarise(Count = n(), .groups = "drop")

if (nrow(direction_counts) > 0) {
  pdf(file.path(OUTPUT_DIR, "Expression_Direction_Barplot_AllPeaks.pdf"),
      width = 8, height = 6)
  p3 <- ggplot(direction_counts, aes(x = Regulatory_Type, y = Count, fill = Direction)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal(base_size = 14) +
    labs(title = "Direction of Expression Changes by Regulatory Type",
         subtitle = "Significant DEGs only, all peaks combined",
         y = "Number of Bound Genes",
         x = "Regulatory Type") +
    scale_fill_manual(values = c("Up" = "#D73027", "Down" = "#4575B4"))
  print(p3)
  dev.off()
}

# --- Plot 4: Violin plot for better visualization ---
pdf(file.path(OUTPUT_DIR, "Expression_Promoter_vs_Enhancer_Violin_AllPeaks.pdf"),
    width = 8, height = 6)
p4 <- ggplot(plot_data_all, aes(x = Regulatory_Type, y = log2FoldChange,
                                fill = Regulatory_Type)) +
  geom_violin(trim = TRUE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  coord_cartesian(ylim = c(-3, 3)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  theme_minimal(base_size = 14) +
  labs(title = "Expression Changes: Promoter vs Enhancer Binding",
       subtitle = sprintf("All peaks combined (n=%d)", nrow(plot_data_all)),
       y = "Log2 Fold Change (TES vs GFP)",
       x = "Regulatory Type") +
  scale_fill_manual(values = c("Promoter" = "#E41A1C",
                               "Enhancer/Distal" = "#377EB8")) +
  theme(legend.position = "none")
print(p4)
dev.off()

cat("\nAnalysis complete. Results saved to", OUTPUT_DIR, "\n")
cat("\nNew output files:\n")
cat("  - Expression_Promoter_vs_Enhancer_AllPeaks.pdf\n")
cat("  - Expression_Promoter_vs_Enhancer_AllPeaks_Significant.pdf\n")
cat("  - Expression_Direction_Barplot_AllPeaks.pdf\n")
cat("  - Expression_Promoter_vs_Enhancer_Violin_AllPeaks.pdf\n")
cat("  - Peaks_Linked_to_Expression_Combined.csv\n")
cat("  - Expression_Summary_Stats_Combined.csv\n")
