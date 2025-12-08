#!/usr/bin/env Rscript

cat("=================================================================\n")
cat("APPROACH 1: Direct Targets Only (Baseline)\n")
cat("=================================================================\n")

# Load workspace and helper functions
source("SRF_Eva_integrated_analysis/focused_enrichment_analysis/scripts/00_helper_functions.R")
load_libraries()
setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top")
load("SRF_Eva_integrated_analysis/focused_enrichment_analysis/data/loaded_data.RData")

# Output directory
base_dir <- "SRF_Eva_integrated_analysis/focused_enrichment_analysis"
approach1_dir <- file.path(base_dir, "approach1_direct_targets")

# Create directories
dir.create(file.path(approach1_dir, "results"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(approach1_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(approach1_dir, "gene_lists"), recursive = TRUE, showWarnings = FALSE)

################################################################################
# ANALYSIS
################################################################################

# TES direct targets
tes_direct_genes <- tes_direct$gene_symbol

cat("TES direct targets:", length(tes_direct_genes), "\n")
cat("Background (expressed genes):", length(all_expressed_genes), "\n")

# Save gene list
write.table(tes_direct_genes,
            file.path(approach1_dir, "gene_lists", "TES_direct_targets.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# GO enrichment
cat("Running GO enrichment...\n")
ego_approach1 <- perform_GO_enrichment(
  gene_list = tes_direct_genes,
  background_genes = all_expressed_genes,
  ont = "BP"
)

results_approach1 <- save_enrichment_results(ego_approach1, "approach1_TES_direct", approach1_dir)

# Extract migration terms
if (!is.null(results_approach1)) {
  migration_approach1 <- extract_migration_terms(results_approach1)
  if (!is.null(migration_approach1) && nrow(migration_approach1) > 0) {
    write_csv(migration_approach1,
              file.path(approach1_dir, "results", "approach1_migration_terms.csv"))
    cat("  Found", nrow(migration_approach1), "migration-related terms\n")
    cat("  Top migration term rank:", which(results_approach1$ID == migration_approach1$ID[1]), "\n")
  }
}

cat("\n=================================================================\n")
cat("APPROACH 1 COMPLETE\n")
cat("=================================================================\n")
