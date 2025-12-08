#!/usr/bin/env Rscript
# Load all data required for enrichment analyses

cat("=================================================================\n")
cat("Loading Data for Enrichment Analysis\n")
cat("=================================================================\n\n")

# Source helper functions
source("SRF_Eva_integrated_analysis/focused_enrichment_analysis/scripts/00_helper_functions.R")

# Load libraries
load_libraries()

# Set working directory
setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top")

################################################################################
# LOAD DATA
################################################################################
cat("\n=== Loading Data ===\n")

# RNA-seq differential expression results
cat("Loading RNA-seq DEG data...\n")
deseq_results <- read.table(
  "SRF_Eva_RNA/results/05_deseq2/deseq2_results_TES_vs_GFP.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

# Integrative analysis results
cat("Loading integrative analysis results...\n")
tes_direct <- read_csv("SRF_Eva_integrated_analysis/output/results/10_direct_targets/TES_direct_targets_all_genes.csv")
tead1_direct <- read_csv("SRF_Eva_integrated_analysis/output/results/10_direct_targets/TEAD1_direct_targets_all_genes.csv")
tes_specific <- read_csv("SRF_Eva_integrated_analysis/output/results/10_direct_targets/TES_specific_targets_all_genes.csv")

# Cut&Tag annotated peaks
cat("Loading Cut&Tag peak annotations...\n")
tes_peaks_annotated <- read_csv("SRF_Eva_CUTandTAG/results/07_analysis_narrow/TES_peaks_annotated.csv")
tead1_peaks_annotated <- read_csv("SRF_Eva_CUTandTAG/results/07_analysis_narrow/TEAD1_peaks_annotated.csv")

# DiffBind results
cat("Loading DiffBind results...\n")
diffbind_tes_vs_tesmut <- read_csv("SRF_Eva_CUTandTAG/results/07_analysis_narrow/DiffBind_TES_vs_TESmut.csv")

# MSigDB gene sets
cat("Loading MSigDB gene sets...\n")
msigdb_hallmark <- msigdbr(species = "Homo sapiens", collection = "H")
msigdb_c2 <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME")

cat("Data loading complete!\n")
cat("  DEGs: ", nrow(deseq_results), "\n")
cat("  TES direct targets: ", nrow(tes_direct), "\n")
cat("  TEAD1 direct targets: ", nrow(tead1_direct), "\n")
cat("  TES peaks annotated: ", nrow(tes_peaks_annotated), "\n")

################################################################################
# PREPARE COMMON VARIABLES
################################################################################

# Background: all genes with reasonable expression
all_expressed_genes <- deseq_results$gene_symbol[!is.na(deseq_results$baseMean) &
  deseq_results$baseMean > 10 &
  !is.na(deseq_results$gene_symbol)]

# DEGs
degs <- deseq_results %>%
  filter(!is.na(padj) & padj < 0.05 & !is.na(gene_symbol))

cat("\n=== Background and DEG Summary ===\n")
cat("Background (expressed genes):", length(all_expressed_genes), "\n")
cat("DEGs (padj < 0.05):", nrow(degs), "\n")

################################################################################
# SAVE WORKSPACE
################################################################################

cat("\n=== Saving workspace ===\n")
save.image("SRF_Eva_integrated_analysis/focused_enrichment_analysis/data/loaded_data.RData")
cat("Workspace saved to: SRF_Eva_integrated_analysis/focused_enrichment_analysis/data/loaded_data.RData\n")

cat("\n=================================================================\n")
cat("Data loading complete!\n")
cat("=================================================================\n")
