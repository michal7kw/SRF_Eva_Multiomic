#!/usr/bin/env Rscript

################################################################################
# Phase 5.2: Master Regulator Analysis
#
# Purpose: Identify key transcription factors controlling TES-responsive genes
#          Creates BOTH simplified (3-category) and detailed (6-category) analyses
#
# Author: Advanced Multi-Omics Analysis Plan
# Date: 2025-01-24
################################################################################

message("=== Phase 5.2: Master Regulator Analysis ===")
message("Start time: ", Sys.time())
message("NOTE: Creating both SIMPLIFIED (3-cat) and DETAILED (6-cat) analyses")

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(msigdbr)
})

# Fix namespace conflicts: dplyr functions get masked by AnnotationDbi
select <- dplyr::select
filter <- dplyr::filter
slice <- dplyr::slice
mutate <- dplyr::mutate
arrange <- dplyr::arrange

setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top")

# Input files
RNA_SEQ_FILE <- "SRF_Eva_RNA/results/05_deseq2/deseq2_results_TES_vs_GFP.txt"
GENE_EXPR_BINDING <- "SRF_Eva_integrated_analysis/scripts/analysis_3/results/04_category_expression/genes_with_binding_and_expression.csv"
BINDING_DATA <- "SRF_Eva_integrated_analysis/scripts/analysis_3/results/01_binding_classification/binding_classification_data.RData"

# Output directories - separate for detailed vs simplified
OUTPUT_DIR <- "SRF_Eva_integrated_analysis/scripts/analysis_3/results/09_tf_networks"
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
# Step 1: Load data
################################################################################

message("\n[Step 1] Loading data...")

# Load RNA-seq results
rna_results <- read.delim(RNA_SEQ_FILE, stringsAsFactors = FALSE)
message("  Loaded RNA-seq: ", nrow(rna_results), " genes")

# Load binding classification
load(BINDING_DATA)

# Load gene-binding data
gene_data <- read.csv(GENE_EXPR_BINDING)
message("  Loaded gene-binding data: ", nrow(gene_data), " genes")

################################################################################
# Step 2: Identify transcription factors among DEGs
################################################################################

message("\n[Step 2] Identifying transcription factors...")

# Get TF annotation from GO Molecular Function
tf_go_terms <- c(
  "GO:0003700", # DNA-binding transcription factor activity
  "GO:0000981", # DNA-binding transcription factor activity, RNA polymerase II-specific
  "GO:0001228"  # DNA-binding transcription factor activity, RNA polymerase II transcription regulatory region sequence-specific binding
)

# Get genes annotated as TFs
tf_annotations <- lapply(tf_go_terms, function(go_term) {
  genes <- get(go_term, org.Hs.egGO2ALLEGS)
  symbols <- mapIds(org.Hs.eg.db, keys = genes, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
  symbols[!is.na(symbols)]
})

all_tfs <- unique(unlist(tf_annotations))
message("  Found ", length(all_tfs), " transcription factors")

# Add TF annotation to RNA-seq data
rna_results$is_TF <- rna_results$gene_symbol %in% all_tfs

# Get DEG TFs
deg_tfs <- rna_results %>%
  filter(is_TF, !is.na(padj), padj < 0.05) %>%
  arrange(padj)

message("  Differentially expressed TFs: ", nrow(deg_tfs))

# Check which TFs are TES/TEAD1 targets
deg_tfs <- deg_tfs %>%
  left_join(gene_data %>% dplyr::select(gene_name, primary_category, has_promoter_peak),
            by = c("gene_symbol" = "gene_name"))

deg_tfs$is_direct_target <- deg_tfs$primary_category %in% c(
  "TES_unique", "TEAD1_unique", "Shared_high",
  "Shared_TES_dominant", "Shared_TEAD1_dominant", "Shared_equivalent"
)

message("  TFs that are direct TES/TEAD1 targets: ", sum(deg_tfs$is_direct_target, na.rm = TRUE))

################################################################################
# Step 3: Get TF target gene sets from MSigDB
################################################################################

message("\n[Step 3] Loading TF target gene sets from MSigDB...")

# Get TF target gene sets (C3 collection)
msigdb_tfs <- msigdbr(species = "Homo sapiens", category = "C3", subcategory = "TFT:GTRD")

message("  Loaded ", length(unique(msigdb_tfs$gs_name)), " TF target gene sets")

# Filter for our DEG TFs
deg_tf_symbols <- deg_tfs$gene_symbol

# Find gene sets related to our DEG TFs
relevant_tf_sets <- msigdb_tfs %>%
  mutate(tf_name = toupper(gsub("_.*", "", gs_name))) %>%
  filter(tf_name %in% toupper(deg_tf_symbols))

message("  Found gene sets for ", length(unique(relevant_tf_sets$gs_name)), " of our DEG TFs")

################################################################################
# Step 4: Test for master regulators
################################################################################

message("\n[Step 4] Testing for master regulator enrichment...")

# Get all DEGs (up and down separately)
all_degs <- rna_results %>%
  filter(!is.na(padj), padj < 0.05)

up_degs <- all_degs %>% filter(log2FoldChange > 0) %>% pull(gene_symbol)
down_degs <- all_degs %>% filter(log2FoldChange < 0) %>% pull(gene_symbol)
all_deg_symbols <- all_degs$gene_symbol

message("  Upregulated DEGs: ", length(up_degs))
message("  Downregulated DEGs: ", length(down_degs))

# For each TF, test if its targets are enriched in DEGs
test_tf_enrichment <- function(tf_targets, deg_list, all_genes) {
  # Fisher's exact test
  in_targets_and_deg <- sum(deg_list %in% tf_targets)
  in_targets_not_deg <- sum(all_genes %in% tf_targets) - in_targets_and_deg
  not_in_targets_deg <- length(deg_list) - in_targets_and_deg
  not_in_targets_not_deg <- length(all_genes) - in_targets_and_deg - in_targets_not_deg - not_in_targets_deg

  contingency_table <- matrix(c(
    in_targets_and_deg, in_targets_not_deg,
    not_in_targets_deg, not_in_targets_not_deg
  ), nrow = 2)

  test_result <- fisher.test(contingency_table, alternative = "greater")

  return(list(
    n_targets = length(tf_targets),
    n_overlap = in_targets_and_deg,
    pct_overlap = round(in_targets_and_deg / length(deg_list) * 100, 2),
    p_value = test_result$p.value,
    odds_ratio = test_result$estimate
  ))
}

# Test each TF
all_genes <- rna_results$gene_symbol[!is.na(rna_results$baseMean)]

master_regulator_results <- lapply(unique(relevant_tf_sets$gs_name), function(gs) {
  tf_targets <- relevant_tf_sets %>%
    filter(gs_name == gs) %>%
    pull(gene_symbol)

  # Extract TF name
  tf_name <- gsub("_.*", "", gs)

  # Get TF's own expression
  tf_expr <- rna_results %>%
    filter(toupper(gene_symbol) == toupper(tf_name)) %>%
    slice(1)

  # Test enrichment
  enrich_all <- test_tf_enrichment(tf_targets, all_deg_symbols, all_genes)
  enrich_up <- test_tf_enrichment(tf_targets, up_degs, all_genes)
  enrich_down <- test_tf_enrichment(tf_targets, down_degs, all_genes)

  data.frame(
    TF = tf_name,
    n_targets_total = enrich_all$n_targets,
    n_targets_in_DEGs = enrich_all$n_overlap,
    p_value_all = enrich_all$p_value,
    odds_ratio_all = enrich_all$odds_ratio,
    n_targets_up = enrich_up$n_overlap,
    p_value_up = enrich_up$p_value,
    n_targets_down = enrich_down$n_overlap,
    p_value_down = enrich_down$p_value,
    TF_log2FC = ifelse(nrow(tf_expr) > 0, tf_expr$log2FoldChange, NA),
    TF_padj = ifelse(nrow(tf_expr) > 0, tf_expr$padj, NA),
    stringsAsFactors = FALSE
  )
})

master_regulator_df <- do.call(rbind, master_regulator_results)

# Adjust p-values
master_regulator_df$padj_all <- p.adjust(master_regulator_df$p_value_all, method = "BH")
master_regulator_df$padj_up <- p.adjust(master_regulator_df$p_value_up, method = "BH")
master_regulator_df$padj_down <- p.adjust(master_regulator_df$p_value_down, method = "BH")

# Rank by significance
master_regulator_df <- master_regulator_df %>%
  arrange(padj_all) %>%
  mutate(
    is_significant = padj_all < 0.05,
    direction_bias = case_when(
      padj_up < 0.05 & padj_down >= 0.05 ~ "Up",
      padj_down < 0.05 & padj_up >= 0.05 ~ "Down",
      padj_up < 0.05 & padj_down < 0.05 ~ "Both",
      TRUE ~ "None"
    )
  )

message("  Significant master regulators (padj < 0.05): ",
        sum(master_regulator_df$is_significant))

################################################################################
# Step 5: Identify which master regulators are TES/TEAD1 targets
################################################################################

message("\n[Step 5] Checking if master regulators are TES/TEAD1 targets...")

master_regulator_df <- master_regulator_df %>%
  left_join(deg_tfs %>% dplyr::select(gene_symbol, primary_category, is_direct_target),
            by = c("TF" = "gene_symbol"))

master_regulator_df$is_TES_TEAD1_target <- !is.na(master_regulator_df$is_direct_target) &
                                            master_regulator_df$is_direct_target

message("  Master regulators that are direct targets: ",
        sum(master_regulator_df$is_TES_TEAD1_target, na.rm = TRUE))

################################################################################
# Step 6: Export results
################################################################################

message("\n[Step 6] Exporting results...")

write.csv(master_regulator_df,
          file.path(OUTPUT_DIR, "master_regulator_analysis.csv"),
          row.names = FALSE)

write.csv(deg_tfs,
          file.path(OUTPUT_DIR, "differentially_expressed_TFs.csv"),
          row.names = FALSE)

# Top master regulators
top_mrs <- master_regulator_df %>%
  filter(is_significant) %>%
  head(20)

write.csv(top_mrs,
          file.path(OUTPUT_DIR, "top_20_master_regulators.csv"),
          row.names = FALSE)

################################################################################
# Step 7: Visualizations
################################################################################

message("\n[Step 7] Generating visualizations...")

# 7.1: Master regulator volcano plot
pdf(file.path(OUTPUT_DIR, "master_regulator_volcano.pdf"), width = 10, height = 8)
ggplot(master_regulator_df,
       aes(x = log2(odds_ratio_all), y = -log10(padj_all))) +
  geom_point(aes(color = is_TES_TEAD1_target, size = n_targets_in_DEGs), alpha = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_text_repel(data = master_regulator_df %>% filter(padj_all < 0.001),
                  aes(label = TF), size = 3, max.overlaps = 20) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray"),
                     labels = c("TRUE" = "TES/TEAD1 Target", "FALSE" = "Not Target"),
                     name = "Direct Target") +
  scale_size_continuous(name = "# Targets\nin DEGs") +
  theme_classic() +
  labs(
    title = "Master Regulator Enrichment Analysis",
    subtitle = "TF targets enriched in TES-responsive genes",
    x = "log2(Odds Ratio)",
    y = "-log10(Adjusted P-value)"
  )
dev.off()

# 7.2: Top master regulators bar plot
if (nrow(top_mrs) > 0) {
  pdf(file.path(OUTPUT_DIR, "top_master_regulators_barplot.pdf"), width = 12, height = 8)
  ggplot(top_mrs, aes(x = reorder(TF, -log10(padj_all)), y = -log10(padj_all),
                      fill = direction_bias)) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    scale_fill_manual(values = c("Up" = "#E41A1C", "Down" = "#377EB8",
                                 "Both" = "#984EA3", "None" = "gray"),
                      name = "Target\nDirection") +
    coord_flip() +
    theme_classic() +
    labs(
      title = "Top 20 Master Regulators",
      subtitle = "Transcription factors with targets enriched in TES-responsive genes",
      x = "Transcription Factor",
      y = "-log10(Adjusted P-value)"
    )
  dev.off()
}

# 7.3: TF expression vs target enrichment
pdf(file.path(OUTPUT_DIR, "TF_expression_vs_enrichment.pdf"), width = 10, height = 8)
mr_with_expr <- master_regulator_df %>%
  filter(!is.na(TF_log2FC), is_significant)

if (nrow(mr_with_expr) > 0) {
  ggplot(mr_with_expr,
         aes(x = TF_log2FC, y = -log10(padj_all))) +
    geom_point(aes(color = is_TES_TEAD1_target, size = n_targets_in_DEGs), alpha = 0.6) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_text_repel(aes(label = TF), size = 3, max.overlaps = 15) +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue"),
                       labels = c("TRUE" = "TES/TEAD1 Target", "FALSE" = "Not Target")) +
    theme_classic() +
    labs(
      title = "Master Regulator Expression vs Target Enrichment",
      x = "TF Expression log2FC (TES vs GFP)",
      y = "Target Enrichment -log10(padj)",
      color = "Direct Target",
      size = "# Targets\nin DEGs"
    )
}
dev.off()

# 7.4: Heatmap of top MRs and their expression
if (nrow(top_mrs) > 0) {
  mr_matrix <- top_mrs %>%
    dplyr::select(TF, TF_log2FC, n_targets_in_DEGs, odds_ratio_all) %>%
    filter(!is.na(TF_log2FC)) %>%
    as.data.frame()

  if (nrow(mr_matrix) > 2) {
    rownames(mr_matrix) <- mr_matrix$TF
    mr_matrix_scaled <- scale(mr_matrix[, -1])

    pdf(file.path(OUTPUT_DIR, "master_regulator_heatmap.pdf"), width = 8, height = 10)
    pheatmap(
      mr_matrix_scaled,
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      color = colorRampPalette(c("blue", "white", "red"))(100),
      main = "Top Master Regulators",
      labels_col = c("TF Expression\nlog2FC", "# Targets\nin DEGs", "Odds Ratio"),
      fontsize_row = 8
    )
    dev.off()
  }
}

################################################################################
# Step 8: Summary report
################################################################################

message("\n[Step 8] Creating summary report...")

sink(file.path(OUTPUT_DIR, "PHASE5_2_SUMMARY.txt"))
cat("=== Phase 5.2: Master Regulator Analysis ===\n")
cat("Date:", as.character(Sys.time()), "\n\n")

cat("Total transcription factors identified:", length(all_tfs), "\n")
cat("Differentially expressed TFs:", nrow(deg_tfs), "\n")
cat("  - Direct TES/TEAD1 targets:", sum(deg_tfs$is_direct_target, na.rm = TRUE), "\n\n")

cat("Master Regulator Enrichment:\n")
cat("============================\n")
cat("TFs tested:", nrow(master_regulator_df), "\n")
cat("Significant master regulators (padj < 0.05):", sum(master_regulator_df$is_significant), "\n")
cat("  - With upregulated targets:", sum(master_regulator_df$direction_bias == "Up"), "\n")
cat("  - With downregulated targets:", sum(master_regulator_df$direction_bias == "Down"), "\n")
cat("  - With both:", sum(master_regulator_df$direction_bias == "Both"), "\n\n")

cat("Master Regulators that are TES/TEAD1 Targets:\n")
cat("==============================================\n")
cat(sum(master_regulator_df$is_TES_TEAD1_target, na.rm = TRUE), "master regulators are direct targets\n\n")

cat("Top 10 Master Regulators:\n")
cat("=========================\n")
print(head(master_regulator_df %>%
           dplyr::select(TF, n_targets_in_DEGs, padj_all, odds_ratio_all, TF_log2FC, direction_bias),
           10))

sink()

message("\n=== Analysis Complete ===")
message("Output directory: ", OUTPUT_DIR)
message("End time: ", Sys.time())
