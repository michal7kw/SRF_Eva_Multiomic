#!/usr/bin/env Rscript

cat("=================================================================\n")
cat("APPROACH 4: High-Confidence Peaks (Top 50% by signal)\n")
cat("=================================================================\n")

# Load workspace and helper functions
source("SRF_Eva_integrated_analysis/focused_enrichment_analysis/scripts/00_helper_functions.R")
load_libraries()
setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top")
load("SRF_Eva_integrated_analysis/focused_enrichment_analysis/data/loaded_data.RData")

# Output directory
base_dir <- "SRF_Eva_integrated_analysis/focused_enrichment_analysis"
approach4_dir <- file.path(base_dir, "approach4_high_confidence")

# Create directories
dir.create(file.path(approach4_dir, "results"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(approach4_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(approach4_dir, "gene_lists"), recursive = TRUE, showWarnings = FALSE)

################################################################################
# ANALYSIS
################################################################################

# Filter for high-confidence peaks (top 50% by qValue or fold enrichment)
# Note: In narrowPeak format from ChIPseeker:
#  V7 = signalValue (fold enrichment)
#  V8 = pValue (-log10)
#  V9 = qValue (-log10)

# Check which columns are available
if ("V9" %in% colnames(tes_peaks_annotated)) {
  # Use qValue (-log10, so higher is better)
  cat("Using V9 (qValue -log10) for filtering...\n")
  median_qval <- median(tes_peaks_annotated$V9, na.rm = TRUE)
  tes_highconf_peaks <- tes_peaks_annotated %>%
    filter(V9 >= median_qval)
} else if ("V7" %in% colnames(tes_peaks_annotated)) {
  # Use signalValue (fold enrichment)
  cat("Using V7 (signalValue) for filtering...\n")
  median_signal <- median(tes_peaks_annotated$V7, na.rm = TRUE)
  tes_highconf_peaks <- tes_peaks_annotated %>%
    filter(V7 >= median_signal)
} else if ("V8" %in% colnames(tes_peaks_annotated)) {
  # Use pValue (-log10, so higher is better)
  cat("Using V8 (pValue -log10) for filtering...\n")
  median_pval <- median(tes_peaks_annotated$V8, na.rm = TRUE)
  tes_highconf_peaks <- tes_peaks_annotated %>%
    filter(V8 >= median_pval)
} else {
  cat("WARNING: No suitable quality columns found, using all peaks\n")
  tes_highconf_peaks <- tes_peaks_annotated
}

cat("Total TES peaks:", nrow(tes_peaks_annotated), "\n")
cat("High-confidence peaks:", nrow(tes_highconf_peaks), "\n")

# Extract genes
tes_highconf_genes <- extract_genes_from_peaks(tes_highconf_peaks)
cat("Unique genes with high-conf peaks:", length(tes_highconf_genes), "\n")

# Intersect with DEGs
tes_highconf_degs <- intersect(tes_highconf_genes, degs$gene_symbol)
cat("High-conf genes that are DEGs:", length(tes_highconf_degs), "\n")

# Save gene lists
write.table(tes_highconf_genes,
            file.path(approach4_dir, "gene_lists", "TES_highconf_peak_genes.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(tes_highconf_degs,
            file.path(approach4_dir, "gene_lists", "TES_highconf_DEGs.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# GO enrichment
cat("Running GO enrichment on high-confidence DEGs...\n")
ego_approach4 <- perform_GO_enrichment(
  gene_list = tes_highconf_degs,
  background_genes = all_expressed_genes,
  ont = "BP"
)

results_approach4 <- save_enrichment_results(ego_approach4, "approach4_TES_highconf_DEGs", approach4_dir)

# Extract migration terms
if (!is.null(results_approach4)) {
  migration_approach4 <- extract_migration_terms(results_approach4)
  if (!is.null(migration_approach4) && nrow(migration_approach4) > 0) {
    write_csv(migration_approach4,
              file.path(approach4_dir, "results", "approach4_migration_terms.csv"))
    cat("  Found", nrow(migration_approach4), "migration-related terms\n")
    cat("  Top migration term rank:", which(results_approach4$ID == migration_approach4$ID[1]), "\n")
  }
}

cat("\n=================================================================\n")
cat("APPROACH 4 COMPLETE\n")
cat("=================================================================\n")
