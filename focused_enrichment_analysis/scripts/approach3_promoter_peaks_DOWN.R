#!/usr/bin/env Rscript

cat("=================================================================\n")
cat("APPROACH 3-DOWN: Promoter Peaks + Downregulated DEGs Only\n")
cat("=================================================================\n")

# Load workspace and helper functions
source("SRF_Eva_integrated_analysis/focused_enrichment_analysis/scripts/00_helper_functions.R")
load_libraries()
setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top")
load("SRF_Eva_integrated_analysis/focused_enrichment_analysis/data/loaded_data.RData")

# Output directory
base_dir <- "SRF_Eva_integrated_analysis/focused_enrichment_analysis"
approach3_down_dir <- file.path(base_dir, "approach3_promoter_peaks_DOWN")

# Create directories
dir.create(file.path(approach3_down_dir, "results"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(approach3_down_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(approach3_down_dir, "gene_lists"), recursive = TRUE, showWarnings = FALSE)

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

# Filter DEGs for downregulated only (TES is a repressor)
degs_down <- degs %>%
  filter(log2FoldChange < 0)

cat("Total DEGs (up+down):", nrow(degs), "\n")
cat("Downregulated DEGs:", nrow(degs_down), "\n")

# Intersect with downregulated DEGs
tes_promoter_degs_down <- intersect(tes_promoter_genes, degs_down$gene_symbol)

cat("Promoter genes that are downregulated DEGs:", length(tes_promoter_degs_down), "\n")

# Save gene lists
write.table(tes_promoter_genes,
            file.path(approach3_down_dir, "gene_lists", "TES_promoter_peak_genes.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(tes_promoter_degs_down,
            file.path(approach3_down_dir, "gene_lists", "TES_promoter_DEGs_DOWN.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# Save detailed table with expression stats
tes_promoter_degs_down_detailed <- degs_down %>%
  filter(gene_symbol %in% tes_promoter_degs_down) %>%
  arrange(padj) %>%
  select(gene_symbol, baseMean, log2FoldChange, padj)

write_csv(tes_promoter_degs_down_detailed,
          file.path(approach3_down_dir, "gene_lists", "TES_promoter_DEGs_DOWN_detailed.csv"))

# GO enrichment on promoter downregulated DEGs
cat("Running GO enrichment on promoter downregulated DEGs...\n")
ego_approach3_down <- perform_GO_enrichment(
  gene_list = tes_promoter_degs_down,
  background_genes = all_expressed_genes,
  ont = "BP"
)

results_approach3_down <- save_enrichment_results(ego_approach3_down, "approach3_TES_promoter_DEGs_DOWN", approach3_down_dir)

# Extract migration terms
if (!is.null(results_approach3_down)) {
  migration_approach3_down <- extract_migration_terms(results_approach3_down)
  if (!is.null(migration_approach3_down) && nrow(migration_approach3_down) > 0) {
    write_csv(migration_approach3_down,
              file.path(approach3_down_dir, "results", "approach3_DOWN_migration_terms.csv"))
    cat("  Found", nrow(migration_approach3_down), "migration-related terms\n")
    cat("  Top migration term rank:", which(results_approach3_down$ID == migration_approach3_down$ID[1]), "\n")
  }
}

cat("\n=================================================================\n")
cat("APPROACH 3-DOWN COMPLETE\n")
cat("Results saved to:", approach3_down_dir, "\n")
cat("=================================================================\n")
