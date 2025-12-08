#!/usr/bin/env Rscript

# 04_gene_regulatory_logic.R
# Module 4: Gene Regulatory Logic
# Part of the SRF_Eva_integrated_analysis pipeline
#
# FIXED: Absolute paths, significance filtering for DEGs

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

annotate_and_assign <- function(peaks, category) {
  if (length(peaks) == 0) {
    warning(paste("No peaks for category:", category))
    return(data.frame())
  }

  anno <- annotatePeak(peaks,
                       TxDb = txdb,
                       tssRegion = c(-2000, 2000),
                       verbose = FALSE)
  df <- as.data.frame(anno)
  df$Category <- category

  # Classify as Promoter or Enhancer/Distal
  df$Regulatory_Type <- ifelse(grepl("Promoter", df$annotation),
                               "Promoter", "Enhancer/Distal")

  df
}

# ===================== Main Analysis =====================

cat("=== Module 4: Gene Regulatory Logic ===\n\n")

# 1. Load Data
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

# Load RNA-seq data
if (!file.exists(RNA_SEQ_FILE)) {
  stop(paste("RNA-seq file not found:", RNA_SEQ_FILE))
}

rna_data <- read.table(RNA_SEQ_FILE, header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE)
cat(sprintf("  RNA-seq genes: %d\n", nrow(rna_data)))

# 2. Annotate Peaks (Promoter vs Enhancer)
cat("\nAnnotating peaks...\n")

df_shared <- annotate_and_assign(peaks_shared, "Shared")
df_tes <- annotate_and_assign(peaks_tes, "TES_Unique")
df_tead1 <- annotate_and_assign(peaks_tead1, "TEAD1_Unique")

combined_df <- rbind(df_shared, df_tes, df_tead1)

cat(sprintf("  Total annotated peaks: %d\n", nrow(combined_df)))

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

# Save merged data (all peaks)
write.csv(merged_df, file.path(OUTPUT_DIR, "Peaks_Linked_to_Expression.csv"),
          row.names = FALSE)

cat(sprintf("  Total peaks with expression data: %d\n",
            sum(!is.na(merged_df$log2FoldChange))))
cat(sprintf("  Peaks linked to significant DEGs (padj < %.2f): %d\n",
            PADJ_THRESHOLD, sum(merged_df$is_significant, na.rm = TRUE)))

# 4. Generate summary statistics
cat("\n=== Expression Summary by Binding Category ===\n")

summary_stats <- merged_df %>%
  filter(!is.na(log2FoldChange)) %>%
  group_by(Category) %>%
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
write.csv(summary_stats, file.path(OUTPUT_DIR, "Expression_Summary_Stats.csv"),
          row.names = FALSE)

# 5. Visualization
cat("\nGenerating expression plots...\n")

# Filter for peaks with valid expression data
plot_data_all <- merged_df[!is.na(merged_df$log2FoldChange), ]
# Filter for significant DEGs only
plot_data_sig <- merged_df[merged_df$is_significant == TRUE, ]

cat(sprintf("  Plotting all genes: %d peaks\n", nrow(plot_data_all)))
cat(sprintf("  Plotting significant DEGs: %d peaks\n", nrow(plot_data_sig)))

# --- Plot 1: All genes by Binding Category ---
pdf(file.path(OUTPUT_DIR, "Expression_by_Binding_Category.pdf"),
    width = 8, height = 6)
p1 <- ggplot(plot_data_all, aes(x = Category, y = log2FoldChange,
                                fill = Category)) +
  geom_boxplot(outlier.shape = NA) +
  coord_cartesian(ylim = c(-2, 2)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  theme_minimal(base_size = 14) +
  labs(title = "Gene Expression Changes by Binding Category",
       subtitle = "All genes (not filtered for significance)",
       y = "Log2 Fold Change (TES vs GFP)",
       x = "") +
  scale_fill_manual(values = c("Shared" = "#984EA3",
                               "TES_Unique" = "#E41A1C",
                               "TEAD1_Unique" = "#377EB8")) +
  theme(legend.position = "none")
print(p1)
dev.off()

# --- Plot 2: Significant DEGs only ---
if (nrow(plot_data_sig) > 10) {
  pdf(file.path(OUTPUT_DIR, "Expression_Significant_DEGs_Only.pdf"),
      width = 8, height = 6)
  p2 <- ggplot(plot_data_sig, aes(x = Category, y = log2FoldChange,
                                  fill = Category)) +
    geom_boxplot(outlier.shape = NA) +
    coord_cartesian(ylim = quantile(plot_data_sig$log2FoldChange,
                                    c(0.01, 0.99), na.rm = TRUE)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    theme_minimal(base_size = 14) +
    labs(title = "Expression of Significant DEGs by Binding Category",
         subtitle = sprintf("Filtered for padj < %.2f", PADJ_THRESHOLD),
         y = "Log2 Fold Change (TES vs GFP)",
         x = "") +
    scale_fill_manual(values = c("Shared" = "#984EA3",
                                 "TES_Unique" = "#E41A1C",
                                 "TEAD1_Unique" = "#377EB8")) +
    theme(legend.position = "none")
  print(p2)
  dev.off()
} else {
  cat("  Too few significant DEGs to plot separately\n")
}

# --- Plot 3: Promoter vs Enhancer (All genes) ---
pdf(file.path(OUTPUT_DIR, "Expression_Promoter_vs_Enhancer.pdf"),
    width = 10, height = 6)
p3 <- ggplot(plot_data_all, aes(x = Category, y = log2FoldChange,
                                fill = Category)) +
  geom_boxplot(outlier.shape = NA) +
  coord_cartesian(ylim = c(-2, 2)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~Regulatory_Type) +
  theme_minimal(base_size = 14) +
  labs(title = "Expression Changes: Promoter vs Enhancer Binding",
       subtitle = "All genes",
       y = "Log2 Fold Change (TES vs GFP)",
       x = "") +
  scale_fill_manual(values = c("Shared" = "#984EA3",
                               "TES_Unique" = "#E41A1C",
                               "TEAD1_Unique" = "#377EB8")) +
  theme(legend.position = "none")
print(p3)
dev.off()

# --- Plot 4: Promoter vs Enhancer (Significant DEGs) ---
if (nrow(plot_data_sig) > 10) {
  pdf(file.path(OUTPUT_DIR, "Expression_Promoter_vs_Enhancer_Significant.pdf"),
      width = 10, height = 6)
  p4 <- ggplot(plot_data_sig, aes(x = Category, y = log2FoldChange,
                                  fill = Category)) +
    geom_boxplot(outlier.shape = NA) +
    coord_cartesian(ylim = quantile(plot_data_sig$log2FoldChange,
                                    c(0.01, 0.99), na.rm = TRUE)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    facet_wrap(~Regulatory_Type) +
    theme_minimal(base_size = 14) +
    labs(title = "Expression: Promoter vs Enhancer (Significant DEGs)",
         subtitle = sprintf("Filtered for padj < %.2f", PADJ_THRESHOLD),
         y = "Log2 Fold Change (TES vs GFP)",
         x = "") +
    scale_fill_manual(values = c("Shared" = "#984EA3",
                                 "TES_Unique" = "#E41A1C",
                                 "TEAD1_Unique" = "#377EB8")) +
    theme(legend.position = "none")
  print(p4)
  dev.off()
}

# --- Plot 5: Direction barplot ---
direction_counts <- plot_data_sig %>%
  filter(Direction %in% c("Up", "Down")) %>%
  group_by(Category, Direction) %>%
  summarise(Count = n(), .groups = "drop")

if (nrow(direction_counts) > 0) {
  pdf(file.path(OUTPUT_DIR, "Expression_Direction_Barplot.pdf"),
      width = 8, height = 6)
  p5 <- ggplot(direction_counts, aes(x = Category, y = Count, fill = Direction)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal(base_size = 14) +
    labs(title = "Direction of Expression Changes (Significant DEGs)",
         y = "Number of Bound Genes",
         x = "") +
    scale_fill_manual(values = c("Up" = "#D73027", "Down" = "#4575B4"))
  print(p5)
  dev.off()
}

cat("\nAnalysis complete. Results saved to", OUTPUT_DIR, "\n")
