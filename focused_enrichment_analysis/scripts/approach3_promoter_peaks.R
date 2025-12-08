#!/usr/bin/env Rscript

cat("=================================================================\n")
cat("APPROACH 3: Promoter Peaks Only\n")
cat("=================================================================\n")

# Load workspace and helper functions
source("SRF_Eva_integrated_analysis/focused_enrichment_analysis/scripts/00_helper_functions.R")
load_libraries()
setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top")
load("SRF_Eva_integrated_analysis/focused_enrichment_analysis/data/loaded_data.RData")

# Output directory
base_dir <- "SRF_Eva_integrated_analysis/focused_enrichment_analysis"
approach3_dir <- file.path(base_dir, "approach3_promoter_peaks")

# Create directories
dir.create(file.path(approach3_dir, "results"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(approach3_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(approach3_dir, "gene_lists"), recursive = TRUE, showWarnings = FALSE)

################################################################################
# ANALYSIS
################################################################################

# Filter for promoter peaks
promoter_keywords <- c("Promoter", "5' UTR", "1st Exon")
tes_promoter_peaks <- tes_peaks_annotated %>%
  filter(grepl(paste(promoter_keywords, collapse = "|"), annotation, ignore.case = TRUE))

cat("Total TES peaks:", nrow(tes_peaks_annotated), "\n")
cat("TES promoter peaks:", nrow(tes_promoter_peaks), "\n")

# Extract genes from promoter peaks
tes_promoter_genes <- extract_genes_from_peaks(tes_promoter_peaks)

cat("Unique genes with promoter peaks:", length(tes_promoter_genes), "\n")

# Intersect with DEGs
tes_promoter_degs <- intersect(tes_promoter_genes, degs$gene_symbol)

cat("Promoter genes that are DEGs:", length(tes_promoter_degs), "\n")

# Save gene lists
write.table(tes_promoter_genes,
            file.path(approach3_dir, "gene_lists", "TES_promoter_peak_genes.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(tes_promoter_degs,
            file.path(approach3_dir, "gene_lists", "TES_promoter_DEGs.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# GO enrichment on promoter DEGs
cat("Running GO enrichment on promoter DEGs...\n")
ego_approach3 <- perform_GO_enrichment(
  gene_list = tes_promoter_degs,
  background_genes = all_expressed_genes,
  ont = "BP"
)

results_approach3 <- save_enrichment_results(ego_approach3, "approach3_TES_promoter_DEGs", approach3_dir)

# Extract migration terms
if (!is.null(results_approach3)) {
  migration_approach3 <- extract_migration_terms(results_approach3)
  if (!is.null(migration_approach3) && nrow(migration_approach3) > 0) {
    write_csv(migration_approach3,
              file.path(approach3_dir, "results", "approach3_migration_terms.csv"))
    cat("  Found", nrow(migration_approach3), "migration-related terms\n")
    cat("  Top migration term rank:", which(results_approach3$ID == migration_approach3$ID[1]), "\n")
  }
}

cat("\n=================================================================\n")
cat("APPROACH 3 COMPLETE\n")
cat("=================================================================\n")
