#!/usr/bin/env Rscript

cat("=================================================================\n")
cat("APPROACH 2: Downregulated Direct Targets (Mechanistically Informed)\n")
cat("=================================================================\n")

# Load workspace and helper functions
source("SRF_Eva_integrated_analysis/focused_enrichment_analysis/scripts/00_helper_functions.R")
load_libraries()
setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top")
load("SRF_Eva_integrated_analysis/focused_enrichment_analysis/data/loaded_data.RData")

# Output directory
base_dir <- "SRF_Eva_integrated_analysis/focused_enrichment_analysis"
approach2_dir <- file.path(base_dir, "approach2_downregulated")

# Create directories
dir.create(file.path(approach2_dir, "results"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(approach2_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(approach2_dir, "gene_lists"), recursive = TRUE, showWarnings = FALSE)

################################################################################
# ANALYSIS
################################################################################

# Filter for downregulated genes (TES is a repressor)
tes_direct_down <- tes_direct %>%
  filter(log2FoldChange < 0, padj < 0.05)

tes_direct_down_genes <- tes_direct_down$gene_symbol

cat("TES downregulated direct targets:", length(tes_direct_down_genes), "\n")

# Save gene list
write.table(tes_direct_down_genes,
            file.path(approach2_dir, "gene_lists", "TES_direct_downregulated.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# Save detailed table
write_csv(tes_direct_down %>%
            arrange(padj) %>%
            select(gene_symbol, baseMean, log2FoldChange, padj, tes_bound, tead1_bound),
          file.path(approach2_dir, "gene_lists", "TES_direct_downregulated_detailed.csv"))

# GO enrichment
cat("Running GO enrichment...\n")
ego_approach2 <- perform_GO_enrichment(
  gene_list = tes_direct_down_genes,
  background_genes = all_expressed_genes,
  ont = "BP"
)

results_approach2 <- save_enrichment_results(ego_approach2, "approach2_TES_direct_down", approach2_dir)

# Extract migration terms
if (!is.null(results_approach2)) {
  migration_approach2 <- extract_migration_terms(results_approach2)
  if (!is.null(migration_approach2) && nrow(migration_approach2) > 0) {
    write_csv(migration_approach2,
              file.path(approach2_dir, "results", "approach2_migration_terms.csv"))
    cat("  Found", nrow(migration_approach2), "migration-related terms\n")
    cat("  Top migration term rank:", which(results_approach2$ID == migration_approach2$ID[1]), "\n")
  }
}

cat("\n=================================================================\n")
cat("APPROACH 2 COMPLETE\n")
cat("=================================================================\n")
