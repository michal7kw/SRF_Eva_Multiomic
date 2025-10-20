#!/usr/bin/env Rscript

cat("=================================================================\n")
cat("APPROACH 6: Migration Gene-Focused (Hypothesis-Driven)\n")
cat("=================================================================\n")

# Load workspace and helper functions
source("SRF_Eva_integrated_analysis/focused_enrichment_analysis/scripts/00_helper_functions.R")
load_libraries()
setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top")
load("SRF_Eva_integrated_analysis/focused_enrichment_analysis/data/loaded_data.RData")

# Output directory
base_dir <- "SRF_Eva_integrated_analysis/focused_enrichment_analysis"
approach6_dir <- file.path(base_dir, "approach6_migration_focused")

# Create directories
dir.create(file.path(approach6_dir, "results"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(approach6_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(approach6_dir, "gene_lists"), recursive = TRUE, showWarnings = FALSE)

################################################################################
# LOAD GENE SETS
################################################################################

# Get GO gene sets from msigdbr
msigdb_go <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP")

# Gene set 1: Neuron migration (GO:0001764)
neuron_migration_genes <- msigdb_go %>%
  filter(grepl("NEURON.*MIGRATION|GO:0001764", gs_name, ignore.case = TRUE)) %>%
  pull(gene_symbol) %>%
  unique()

cat("Neuron migration gene set size:", length(neuron_migration_genes), "\n")

# Gene set 2: Cell proliferation (GO:0008283)
cell_prolif_genes <- msigdb_go %>%
  filter(grepl("CELL_PROLIFERATION|GO:0008283", gs_name, ignore.case = TRUE)) %>%
  pull(gene_symbol) %>%
  unique()

cat("Cell proliferation gene set size:", length(cell_prolif_genes), "\n")

# Gene set 3: Apoptosis from Hallmark
apoptosis_genes <- msigdb_hallmark %>%
  filter(gs_name == "HALLMARK_APOPTOSIS") %>%
  pull(gene_symbol)

cat("Apoptosis gene set size:", length(apoptosis_genes), "\n")

################################################################################
# HYPERGEOMETRIC TESTS
################################################################################
cat("\n=== Hypergeometric enrichment tests ===\n")

# TES direct targets
tes_direct_genes <- tes_direct$gene_symbol

# Test 1: Neuron migration
if (length(neuron_migration_genes) > 0) {
  tes_direct_in_neuron_mig <- intersect(tes_direct_genes, neuron_migration_genes)
  cat("\nNeuron migration:\n")
  cat("  TES direct targets in gene set:", length(tes_direct_in_neuron_mig), "\n")

  if (length(tes_direct_in_neuron_mig) > 0) {
    universe_size <- length(all_expressed_genes)
    geneset_in_universe <- length(intersect(neuron_migration_genes, all_expressed_genes))
    tes_direct_size <- length(tes_direct_genes)
    overlap <- length(tes_direct_in_neuron_mig)

    phyper_pval <- phyper(q = overlap - 1,
                          m = geneset_in_universe,
                          n = universe_size - geneset_in_universe,
                          k = tes_direct_size,
                          lower.tail = FALSE)
    cat("  Hypergeometric p-value:", phyper_pval, "\n")

    write.table(tes_direct_in_neuron_mig,
                file.path(approach6_dir, "gene_lists", "TES_direct_neuron_migration.txt"),
                row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
}

# Test 2: Cell proliferation
if (length(cell_prolif_genes) > 0) {
  tes_direct_in_prolif <- intersect(tes_direct_genes, cell_prolif_genes)
  cat("\nCell proliferation:\n")
  cat("  TES direct targets in gene set:", length(tes_direct_in_prolif), "\n")

  if (length(tes_direct_in_prolif) > 0) {
    universe_size <- length(all_expressed_genes)
    geneset_in_universe <- length(intersect(cell_prolif_genes, all_expressed_genes))
    tes_direct_size <- length(tes_direct_genes)
    overlap <- length(tes_direct_in_prolif)

    phyper_pval <- phyper(q = overlap - 1,
                          m = geneset_in_universe,
                          n = universe_size - geneset_in_universe,
                          k = tes_direct_size,
                          lower.tail = FALSE)
    cat("  Hypergeometric p-value:", phyper_pval, "\n")

    write.table(tes_direct_in_prolif,
                file.path(approach6_dir, "gene_lists", "TES_direct_proliferation.txt"),
                row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
}

# Test 3: Apoptosis
if (length(apoptosis_genes) > 0) {
  tes_direct_in_apoptosis <- intersect(tes_direct_genes, apoptosis_genes)
  cat("\nApoptosis:\n")
  cat("  TES direct targets in gene set:", length(tes_direct_in_apoptosis), "\n")

  if (length(tes_direct_in_apoptosis) > 0) {
    universe_size <- length(all_expressed_genes)
    geneset_in_universe <- length(intersect(apoptosis_genes, all_expressed_genes))
    tes_direct_size <- length(tes_direct_genes)
    overlap <- length(tes_direct_in_apoptosis)

    phyper_pval <- phyper(q = overlap - 1,
                          m = geneset_in_universe,
                          n = universe_size - geneset_in_universe,
                          k = tes_direct_size,
                          lower.tail = FALSE)
    cat("  Hypergeometric p-value:", phyper_pval, "\n")

    write.table(tes_direct_in_apoptosis,
                file.path(approach6_dir, "gene_lists", "TES_direct_apoptosis.txt"),
                row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
}

################################################################################
# DETAILED TABLES
################################################################################

# Create detailed tables for each gene set
if (exists("tes_direct_in_neuron_mig") && length(tes_direct_in_neuron_mig) > 0) {
  neuron_mig_detailed <- tes_direct %>%
    filter(gene_symbol %in% tes_direct_in_neuron_mig) %>%
    arrange(padj) %>%
    select(gene_symbol, baseMean, log2FoldChange, padj, tes_bound, tead1_bound)

  write_csv(neuron_mig_detailed,
            file.path(approach6_dir, "gene_lists", "TES_direct_neuron_migration_detailed.csv"))
}

if (exists("tes_direct_in_prolif") && length(tes_direct_in_prolif) > 0) {
  prolif_detailed <- tes_direct %>%
    filter(gene_symbol %in% tes_direct_in_prolif) %>%
    arrange(padj) %>%
    select(gene_symbol, baseMean, log2FoldChange, padj, tes_bound, tead1_bound)

  write_csv(prolif_detailed,
            file.path(approach6_dir, "gene_lists", "TES_direct_proliferation_detailed.csv"))
}

if (exists("tes_direct_in_apoptosis") && length(tes_direct_in_apoptosis) > 0) {
  apoptosis_detailed <- tes_direct %>%
    filter(gene_symbol %in% tes_direct_in_apoptosis) %>%
    arrange(padj) %>%
    select(gene_symbol, baseMean, log2FoldChange, padj, tes_bound, tead1_bound)

  write_csv(apoptosis_detailed,
            file.path(approach6_dir, "gene_lists", "TES_direct_apoptosis_detailed.csv"))
}

################################################################################
# GSEA ANALYSIS
################################################################################
cat("\n=== GSEA analysis ===\n")

# Create ranked list by log2FoldChange * -log10(padj)
degs_ranked <- deseq_results %>%
  filter(!is.na(padj) & !is.na(log2FoldChange) & baseMean > 10) %>%
  mutate(padj_capped = pmax(padj, 1e-300)) %>%
  mutate(rank_metric = log2FoldChange * -log10(padj_capped)) %>%
  filter(is.finite(rank_metric)) %>%
  arrange(desc(rank_metric))

ranked_genes <- setNames(degs_ranked$rank_metric, degs_ranked$gene_symbol)

# Combine the three gene sets of interest into a custom collection
custom_genesets <- bind_rows(
  msigdb_go %>% filter(grepl("NEURON.*MIGRATION|GO:0001764", gs_name, ignore.case = TRUE)) %>%
    mutate(gs_name = "GOBP_NEURON_MIGRATION"),
  msigdb_go %>% filter(grepl("CELL_PROLIFERATION|GO:0008283", gs_name, ignore.case = TRUE)) %>%
    mutate(gs_name = "GOBP_CELL_PROLIFERATION"),
  msigdb_hallmark %>% filter(gs_name == "HALLMARK_APOPTOSIS")
)

cat("Running GSEA on custom gene sets (neuron migration, proliferation, apoptosis)...\n")

if (nrow(custom_genesets) > 0) {
  gsea_custom <- perform_GSEA_msigdb(ranked_genes, custom_genesets, pvalueCutoff = 0.25)

  if (!is.null(gsea_custom) && nrow(gsea_custom) > 0) {
    gsea_results <- as.data.frame(gsea_custom)
    write_csv(gsea_results,
              file.path(approach6_dir, "results", "GSEA_custom_genesets.csv"))

    cat("\nGSEA results:\n")
    print(gsea_results[, c("ID", "NES", "pvalue", "p.adjust")])

    # Plot GSEA results
    pdf(file.path(approach6_dir, "plots", "GSEA_custom_dotplot.pdf"), width = 10, height = 6)
    print(dotplot(gsea_custom, showCategory = 10, title = "GSEA: Migration, Proliferation, Apoptosis"))
    dev.off()

    # Individual GSEA plots for significant gene sets
    for (geneset_id in gsea_results$ID) {
      safe_name <- gsub("[^A-Za-z0-9_]", "_", geneset_id)
      pdf(file.path(approach6_dir, "plots", paste0("GSEA_", safe_name, ".pdf")),
          width = 10, height = 6)
      print(gseaplot2(gsea_custom,
                      geneSetID = geneset_id,
                      title = geneset_id))
      dev.off()
    }
  } else {
    cat("No significant GSEA results found (p < 0.25)\n")
  }
} else {
  cat("Could not create custom gene sets for GSEA\n")
}

################################################################################
# VISUALIZATIONS
################################################################################
cat("\n=== Creating visualizations ===\n")

# Prepare data for enrichment barplot
enrichment_results <- data.frame(
  GeneSet = c("Neuron\nMigration", "Cell\nProliferation", "Apoptosis"),
  Overlap = c(
    if(exists("tes_direct_in_neuron_mig")) length(tes_direct_in_neuron_mig) else 0,
    if(exists("tes_direct_in_prolif")) length(tes_direct_in_prolif) else 0,
    if(exists("tes_direct_in_apoptosis")) length(tes_direct_in_apoptosis) else 0
  ),
  GeneSetSize = c(
    length(neuron_migration_genes),
    length(cell_prolif_genes),
    length(apoptosis_genes)
  ),
  PValue = c(
    if(exists("tes_direct_in_neuron_mig") && length(tes_direct_in_neuron_mig) > 0) {
      phyper(length(tes_direct_in_neuron_mig) - 1,
             length(intersect(neuron_migration_genes, all_expressed_genes)),
             length(all_expressed_genes) - length(intersect(neuron_migration_genes, all_expressed_genes)),
             length(tes_direct_genes), lower.tail = FALSE)
    } else NA,
    if(exists("tes_direct_in_prolif") && length(tes_direct_in_prolif) > 0) {
      phyper(length(tes_direct_in_prolif) - 1,
             length(intersect(cell_prolif_genes, all_expressed_genes)),
             length(all_expressed_genes) - length(intersect(cell_prolif_genes, all_expressed_genes)),
             length(tes_direct_genes), lower.tail = FALSE)
    } else NA,
    if(exists("tes_direct_in_apoptosis") && length(tes_direct_in_apoptosis) > 0) {
      phyper(length(tes_direct_in_apoptosis) - 1,
             length(intersect(apoptosis_genes, all_expressed_genes)),
             length(all_expressed_genes) - length(intersect(apoptosis_genes, all_expressed_genes)),
             length(tes_direct_genes), lower.tail = FALSE)
    } else NA
  )
) %>%
  filter(!is.na(PValue)) %>%
  mutate(
    NegLog10P = -log10(PValue),
    PercentOverlap = (Overlap / GeneSetSize) * 100,
    Significance = case_when(
      PValue < 0.001 ~ "p < 0.001",
      PValue < 0.01 ~ "p < 0.01",
      PValue < 0.05 ~ "p < 0.05",
      TRUE ~ "n.s."
    )
  )

# Plot 1: Enrichment barplot with -log10(p-value)
pdf(file.path(approach6_dir, "plots", "hypergeometric_enrichment.pdf"), width = 8, height = 6)
p1 <- ggplot(enrichment_results, aes(x = reorder(GeneSet, -NegLog10P), y = NegLog10P, fill = GeneSet)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "darkred") +
  geom_text(aes(label = sprintf("p = %.2e", PValue)), vjust = -0.5, size = 3.5) +
  labs(title = "TES Direct Target Enrichment in Biological Gene Sets",
       subtitle = "Hypergeometric test",
       x = "Gene Set",
       y = "-log10(p-value)") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 11),
        plot.title = element_text(size = 14, face = "bold")) +
  scale_fill_brewer(palette = "Set2")
print(p1)
dev.off()

# Plot 2: Overlap counts
pdf(file.path(approach6_dir, "plots", "geneset_overlap_counts.pdf"), width = 8, height = 6)
p2 <- ggplot(enrichment_results, aes(x = reorder(GeneSet, -Overlap), y = Overlap, fill = GeneSet)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = sprintf("%d genes\n(%.1f%%)", Overlap, PercentOverlap)),
            vjust = -0.5, size = 3.5) +
  labs(title = "TES Direct Targets in Biological Gene Sets",
       subtitle = sprintf("Out of %d total TES direct targets", length(tes_direct_genes)),
       x = "Gene Set",
       y = "Number of TES Direct Targets") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 11),
        plot.title = element_text(size = 14, face = "bold")) +
  scale_fill_brewer(palette = "Set2")
print(p2)
dev.off()

# Plot 3: Expression heatmap for top genes in each category
if (exists("tes_direct_in_neuron_mig") && length(tes_direct_in_neuron_mig) > 0) {
  top_neuron_genes <- tes_direct %>%
    filter(gene_symbol %in% tes_direct_in_neuron_mig) %>%
    arrange(padj) %>%
    head(20) %>%
    pull(gene_symbol)

  heatmap_data <- tes_direct %>%
    filter(gene_symbol %in% top_neuron_genes) %>%
    select(gene_symbol, log2FoldChange) %>%
    arrange(log2FoldChange)

  pdf(file.path(approach6_dir, "plots", "neuron_migration_top_genes_heatmap.pdf"),
      width = 6, height = 8)
  p3 <- ggplot(heatmap_data, aes(x = 1, y = reorder(gene_symbol, log2FoldChange),
                                 fill = log2FoldChange)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", log2FoldChange)), size = 3) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                         name = "log2FC") +
    labs(title = "Top 20 Neuron Migration Genes",
         subtitle = "TES Direct Targets",
         x = "", y = "") +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.title = element_text(size = 14, face = "bold"))
  print(p3)
  dev.off()
}

cat("Plots created successfully!\n")

cat("\n=================================================================\n")
cat("APPROACH 6 COMPLETE\n")
cat("=================================================================\n")
