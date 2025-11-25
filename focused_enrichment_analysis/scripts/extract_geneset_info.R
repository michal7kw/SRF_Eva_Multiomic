#!/usr/bin/env Rscript

cat("=================================================================\n")
cat("Extracting Gene Set Information for Approach 6\n")
cat("=================================================================\n\n")

# Load required libraries
suppressPackageStartupMessages({
  library(msigdbr)
  library(dplyr)
  library(readr)
})

# Set working directory
setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top")

cat("Loading MSigDB gene sets...\n")

# Get GO gene sets from msigdbr
msigdb_go <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP")
msigdb_hallmark <- msigdbr(species = "Homo sapiens", collection = "H")

# Gene set 1: Neuron migration (GO:0001764)
cat("\n=== GENE SET 1: NEURON MIGRATION ===\n")
neuron_migration_genes <- msigdb_go %>%
  filter(grepl("NEURON.*MIGRATION|GO:0001764", gs_name, ignore.case = TRUE)) %>%
  pull(gene_symbol) %>%
  unique()

cat("Total neuron migration genes:", length(neuron_migration_genes), "\n")
cat("Gene symbols:\n")
cat(paste(neuron_migration_genes, collapse = ", "), "\n")

# Save to file
write.table(neuron_migration_genes,
  file = "SRF_Eva_integrated_analysis/focused_enrichment_analysis/approach6_migration_focused/gene_lists/neuron_migration_geneset.txt",
  row.names = FALSE, col.names = FALSE, quote = FALSE
)

# Gene set 2: Cell proliferation (GO:0008283)
cat("\n=== GENE SET 2: CELL PROLIFERATION ===\n")
cell_prolif_genes <- msigdb_go %>%
  filter(grepl("CELL_PROLIFERATION|GO:0008283", gs_name, ignore.case = TRUE)) %>%
  pull(gene_symbol) %>%
  unique()

cat("Total cell proliferation genes:", length(cell_prolif_genes), "\n")
cat("Gene symbols:\n")
cat(paste(cell_prolif_genes, collapse = ", "), "\n")

# Save to file
write.table(cell_prolif_genes,
  file = "SRF_Eva_integrated_analysis/focused_enrichment_analysis/approach6_migration_focused/gene_lists/cell_proliferation_geneset.txt",
  row.names = FALSE, col.names = FALSE, quote = FALSE
)

# Gene set 3: Apoptosis from Hallmark
cat("\n=== GENE SET 3: APOPTOSIS (HALLMARK) ===\n")
apoptosis_genes <- msigdb_hallmark %>%
  filter(gs_name == "HALLMARK_APOPTOSIS") %>%
  pull(gene_symbol)

cat("Total apoptosis genes:", length(apoptosis_genes), "\n")
cat("Gene symbols:\n")
cat(paste(apoptosis_genes, collapse = ", "), "\n")

# Save to file
write.table(apoptosis_genes,
  file = "SRF_Eva_integrated_analysis/focused_enrichment_analysis/approach6_migration_focused/gene_lists/apoptosis_geneset.txt",
  row.names = FALSE, col.names = FALSE, quote = FALSE
)

# Create a summary table
cat("\n=== SUMMARY TABLE ===\n")
summary_table <- data.frame(
  GeneSet = c("Neuron Migration (GO:0001764)", "Cell Proliferation (GO:0008283)", "Apoptosis (Hallmark)"),
  GeneCount = c(length(neuron_migration_genes), length(cell_prolif_genes), length(apoptosis_genes)),
  Description = c("Genes involved in neuron migration processes", "Genes involved in cell proliferation", "Genes involved in apoptosis pathways")
)

print(summary_table)

# Save summary table
write_csv(
  summary_table,
  "SRF_Eva_integrated_analysis/focused_enrichment_analysis/approach6_migration_focused/gene_lists/geneset_summary.csv"
)

# Also load TES direct targets if available to show overlap
tes_direct_file <- "SRF_Eva_integrated_analysis/output/results/10_direct_targets/TES_direct_targets_all_genes.csv"
if (file.exists(tes_direct_file)) {
  cat("\n=== TES DIRECT TARGETS OVERLAP ===\n")
  tes_direct <- read_csv(tes_direct_file)
  tes_direct_genes <- tes_direct$gene_symbol

  cat("Total TES direct targets:", length(tes_direct_genes), "\n")

  # Calculate overlaps
  neuron_overlap <- intersect(tes_direct_genes, neuron_migration_genes)
  prolif_overlap <- intersect(tes_direct_genes, cell_prolif_genes)
  apoptosis_overlap <- intersect(tes_direct_genes, apoptosis_genes)

  cat("\nOverlap with TES direct targets:\n")
  cat(
    "  Neuron migration:", length(neuron_overlap), "genes (",
    round(length(neuron_overlap) / length(neuron_migration_genes) * 100, 1), "% of gene set, ",
    round(length(neuron_overlap) / length(tes_direct_genes) * 100, 1), "% of TES targets)\n"
  )
  cat(
    "  Cell proliferation:", length(prolif_overlap), "genes (",
    round(length(prolif_overlap) / length(cell_prolif_genes) * 100, 1), "% of gene set, ",
    round(length(prolif_overlap) / length(tes_direct_genes) * 100, 1), "% of TES targets)\n"
  )
  cat(
    "  Apoptosis:", length(apoptosis_overlap), "genes (",
    round(length(apoptosis_overlap) / length(apoptosis_genes) * 100, 1), "% of gene set, ",
    round(length(apoptosis_overlap) / length(tes_direct_genes) * 100, 1), "% of TES targets)\n"
  )

  # Save overlap genes
  write.table(neuron_overlap,
    file = "SRF_Eva_integrated_analysis/focused_enrichment_analysis/approach6_migration_focused/gene_lists/TES_neuron_migration_overlap.txt",
    row.names = FALSE, col.names = FALSE, quote = FALSE
  )

  write.table(prolif_overlap,
    file = "SRF_Eva_integrated_analysis/focused_enrichment_analysis/approach6_migration_focused/gene_lists/TES_cell_proliferation_overlap.txt",
    row.names = FALSE, col.names = FALSE, quote = FALSE
  )

  write.table(apoptosis_overlap,
    file = "SRF_Eva_integrated_analysis/focused_enrichment_analysis/approach6_migration_focused/gene_lists/TES_apoptosis_overlap.txt",
    row.names = FALSE, col.names = FALSE, quote = FALSE
  )

  # Create overlap summary table
  overlap_table <- data.frame(
    GeneSet = c("Neuron Migration", "Cell Proliferation", "Apoptosis"),
    GeneSetSize = c(length(neuron_migration_genes), length(cell_prolif_genes), length(apoptosis_genes)),
    TES_Overlap = c(length(neuron_overlap), length(prolif_overlap), length(apoptosis_overlap)),
    PercentOfGeneSet = c(
      round(length(neuron_overlap) / length(neuron_migration_genes) * 100, 1),
      round(length(prolif_overlap) / length(cell_prolif_genes) * 100, 1),
      round(length(apoptosis_overlap) / length(apoptosis_genes) * 100, 1)
    ),
    PercentOfTES = c(
      round(length(neuron_overlap) / length(tes_direct_genes) * 100, 1),
      round(length(prolif_overlap) / length(tes_direct_genes) * 100, 1),
      round(length(apoptosis_overlap) / length(tes_direct_genes) * 100, 1)
    )
  )

  cat("\nOverlap Summary Table:\n")
  print(overlap_table)

  write_csv(
    overlap_table,
    "SRF_Eva_integrated_analysis/focused_enrichment_analysis/approach6_migration_focused/gene_lists/overlap_summary.csv"
  )
}

cat("\n=================================================================\n")
cat("Gene set extraction complete!\n")
cat("Files saved in: SRF_Eva_integrated_analysis/focused_enrichment_analysis/approach6_migration_focused/gene_lists/\n")
cat("=================================================================\n")
