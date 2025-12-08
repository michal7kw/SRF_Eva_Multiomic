#!/usr/bin/env Rscript

cat("=================================================================\n")
cat("APPROACH 6 CORRECTED: Migration Gene-Focused (Hypothesis-Driven)\n")
cat("CORRECTED: Test if migration genes are enriched in TES targets\n")
cat("=================================================================\n")

# Load workspace and helper functions
source("SRF_Eva_integrated_analysis/focused_enrichment_analysis/scripts/00_helper_functions.R")
load_libraries()
setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top")
load("SRF_Eva_integrated_analysis/focused_enrichment_analysis/data/loaded_data.RData")

# Output directory
base_dir <- "SRF_Eva_integrated_analysis/focused_enrichment_analysis"
approach6_dir <- file.path(base_dir, "approach6_migration_focused_CORRECTED")

# Create directories
dir.create(file.path(approach6_dir, "results"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(approach6_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(approach6_dir, "gene_lists"), recursive = TRUE, showWarnings = FALSE)

################################################################################
# LOAD GENE SETS
################################################################################

# Get GO gene sets from msigdbr
msigdb_go <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP")

# Gene set 1: Cell migration (GO:0016477)
cell_migration_genes <- msigdb_go %>%
  filter(grepl("CELL_MIGRATION|GO:0016477", gs_name, ignore.case = TRUE)) %>%
  pull(gene_symbol) %>%
  unique()

cat("Cell migration gene set size:", length(cell_migration_genes), "\n")

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
# CORRECTED HYPERGEOMETRIC TESTS
################################################################################

cat("\n=== CORRECTED Hypergeometric enrichment tests ===\n")
cat("Testing: Are migration/proliferation/apoptosis genes enriched\n")
cat("         among TES direct targets (vs all genes)?\n\n")

# TES direct targets
tes_direct_genes <- tes_direct$gene_symbol

# Universe: all genes with reasonable expression
universe_genes <- all_expressed_genes
universe_size <- length(universe_genes)

cat("Universe (expressed genes):", universe_size, "\n")
cat("TES direct targets:", length(tes_direct_genes), "\n\n")

# CORRECTED Test 1: Cell migration
if (length(cell_migration_genes) > 0) {
  cat("--- Cell Migration ---\n")

  # Genes from migration set that are in universe
  migration_in_universe <- intersect(cell_migration_genes, universe_genes)
  m <- length(migration_in_universe)

  # Genes from migration set that are TES direct targets
  migration_in_tes_targets <- intersect(cell_migration_genes, tes_direct_genes)
  q <- length(migration_in_tes_targets)

  # TES target size
  k <- length(tes_direct_genes)

  # Expected by chance
  expected <- (k * m) / universe_size

  # Fold enrichment
  fold_enrich <- q / expected

  # CORRECTED Hypergeometric test
  # q = number of migration genes in TES targets
  # m = number of migration genes in universe
  # n = number of non-migration genes in universe
  # k = number of TES targets
  phyper_pval <- phyper(q = q - 1,
                        m = m,
                        n = universe_size - m,
                        k = k,
                        lower.tail = FALSE)

  cat("  Migration genes in universe:", m, "\n")
  cat("  Migration genes in TES targets:", q, "\n")
  cat("  Expected by chance:", round(expected, 2), "\n")
  cat("  Fold enrichment:", round(fold_enrich, 3), "\n")
  cat("  Hypergeometric p-value:", format(phyper_pval, scientific = TRUE), "\n\n")

  write.table(migration_in_tes_targets,
              file.path(approach6_dir, "gene_lists", "cell_migration_in_TES_targets.txt"),
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}

# CORRECTED Test 2: Cell proliferation
if (length(cell_prolif_genes) > 0) {
  cat("--- Cell Proliferation ---\n")

  prolif_in_universe <- intersect(cell_prolif_genes, universe_genes)
  m <- length(prolif_in_universe)

  prolif_in_tes_targets <- intersect(cell_prolif_genes, tes_direct_genes)
  q <- length(prolif_in_tes_targets)

  k <- length(tes_direct_genes)
  expected <- (k * m) / universe_size
  fold_enrich <- q / expected

  phyper_pval <- phyper(q = q - 1,
                        m = m,
                        n = universe_size - m,
                        k = k,
                        lower.tail = FALSE)

  cat("  Proliferation genes in universe:", m, "\n")
  cat("  Proliferation genes in TES targets:", q, "\n")
  cat("  Expected by chance:", round(expected, 2), "\n")
  cat("  Fold enrichment:", round(fold_enrich, 3), "\n")
  cat("  Hypergeometric p-value:", format(phyper_pval, scientific = TRUE), "\n\n")

  write.table(prolif_in_tes_targets,
              file.path(approach6_dir, "gene_lists", "cell_proliferation_in_TES_targets.txt"),
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}

# CORRECTED Test 3: Apoptosis
if (length(apoptosis_genes) > 0) {
  cat("--- Apoptosis ---\n")

  apoptosis_in_universe <- intersect(apoptosis_genes, universe_genes)
  m <- length(apoptosis_in_universe)

  apoptosis_in_tes_targets <- intersect(apoptosis_genes, tes_direct_genes)
  q <- length(apoptosis_in_tes_targets)

  k <- length(tes_direct_genes)
  expected <- (k * m) / universe_size
  fold_enrich <- q / expected

  phyper_pval <- phyper(q = q - 1,
                        m = m,
                        n = universe_size - m,
                        k = k,
                        lower.tail = FALSE)

  cat("  Apoptosis genes in universe:", m, "\n")
  cat("  Apoptosis genes in TES targets:", q, "\n")
  cat("  Expected by chance:", round(expected, 2), "\n")
  cat("  Fold enrichment:", round(fold_enrich, 3), "\n")
  cat("  Hypergeometric p-value:", format(phyper_pval, scientific = TRUE), "\n\n")

  write.table(apoptosis_in_tes_targets,
              file.path(approach6_dir, "gene_lists", "apoptosis_in_TES_targets.txt"),
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}

################################################################################
# DETAILED TABLES
################################################################################

cat("=== Creating detailed tables ===\n")

# Cell migration detailed
if (exists("migration_in_tes_targets") && length(migration_in_tes_targets) > 0) {
  migration_detailed <- tes_direct %>%
    filter(gene_symbol %in% migration_in_tes_targets) %>%
    arrange(padj) %>%
    select(gene_symbol, baseMean, log2FoldChange, padj, tes_bound, tead1_bound)

  write_csv(migration_detailed,
            file.path(approach6_dir, "gene_lists", "cell_migration_in_TES_targets_detailed.csv"))
}

# Cell proliferation detailed
if (exists("prolif_in_tes_targets") && length(prolif_in_tes_targets) > 0) {
  prolif_detailed <- tes_direct %>%
    filter(gene_symbol %in% prolif_in_tes_targets) %>%
    arrange(padj) %>%
    select(gene_symbol, baseMean, log2FoldChange, padj, tes_bound, tead1_bound)

  write_csv(prolif_detailed,
            file.path(approach6_dir, "gene_lists", "cell_proliferation_in_TES_targets_detailed.csv"))
}

# Apoptosis detailed
if (exists("apoptosis_in_tes_targets") && length(apoptosis_in_tes_targets) > 0) {
  apoptosis_detailed <- tes_direct %>%
    filter(gene_symbol %in% apoptosis_in_tes_targets) %>%
    arrange(padj) %>%
    select(gene_symbol, baseMean, log2FoldChange, padj, tes_bound, tead1_bound)

  write_csv(apoptosis_detailed,
            file.path(approach6_dir, "gene_lists", "apoptosis_in_TES_targets_detailed.csv"))
}

################################################################################
# GSEA ANALYSIS (SAME AS BEFORE)
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
  msigdb_go %>% filter(grepl("CELL_MIGRATION|GO:0016477", gs_name, ignore.case = TRUE)) %>%
    mutate(gs_name = "GOBP_CELL_MIGRATION"),
  msigdb_go %>% filter(grepl("CELL_PROLIFERATION|GO:0008283", gs_name, ignore.case = TRUE)) %>%
    mutate(gs_name = "GOBP_CELL_PROLIFERATION"),
  msigdb_hallmark %>% filter(gs_name == "HALLMARK_APOPTOSIS")
)

cat("Running GSEA on custom gene sets (migration, proliferation, apoptosis)...\n")

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

    # Individual GSEA plots
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
}

################################################################################
# VISUALIZATIONS
################################################################################

cat("\n=== Creating visualizations ===\n")

# Prepare data for enrichment barplot
enrichment_results <- data.frame(
  GeneSet = c("Cell\nMigration", "Cell\nProliferation", "Apoptosis"),
  Overlap = c(
    if(exists("migration_in_tes_targets")) length(migration_in_tes_targets) else 0,
    if(exists("prolif_in_tes_targets")) length(prolif_in_tes_targets) else 0,
    if(exists("apoptosis_in_tes_targets")) length(apoptosis_in_tes_targets) else 0
  ),
  GeneSetSize = c(
    length(cell_migration_genes),
    length(cell_prolif_genes),
    length(apoptosis_genes)
  ),
  Expected = c(
    (length(tes_direct_genes) * length(intersect(cell_migration_genes, universe_genes))) / universe_size,
    (length(tes_direct_genes) * length(intersect(cell_prolif_genes, universe_genes))) / universe_size,
    (length(tes_direct_genes) * length(intersect(apoptosis_genes, universe_genes))) / universe_size
  ),
  PValue = c(
    if(exists("migration_in_tes_targets") && length(migration_in_tes_targets) > 0) {
      m <- length(intersect(cell_migration_genes, universe_genes))
      q <- length(migration_in_tes_targets)
      k <- length(tes_direct_genes)
      phyper(q - 1, m, universe_size - m, k, lower.tail = FALSE)
    } else NA,
    if(exists("prolif_in_tes_targets") && length(prolif_in_tes_targets) > 0) {
      m <- length(intersect(cell_prolif_genes, universe_genes))
      q <- length(prolif_in_tes_targets)
      k <- length(tes_direct_genes)
      phyper(q - 1, m, universe_size - m, k, lower.tail = FALSE)
    } else NA,
    if(exists("apoptosis_in_tes_targets") && length(apoptosis_in_tes_targets) > 0) {
      m <- length(intersect(apoptosis_genes, universe_genes))
      q <- length(apoptosis_in_tes_targets)
      k <- length(tes_direct_genes)
      phyper(q - 1, m, universe_size - m, k, lower.tail = FALSE)
    } else NA
  )
) %>%
  filter(!is.na(PValue)) %>%
  mutate(
    FoldEnrichment = Overlap / Expected,
    NegLog10P = -log10(PValue),
    PercentOverlap = (Overlap / GeneSetSize) * 100,
    Significance = case_when(
      PValue < 0.001 ~ "***",
      PValue < 0.01 ~ "**",
      PValue < 0.05 ~ "*",
      TRUE ~ "n.s."
    )
  )

# Save enrichment results
write_csv(enrichment_results, file.path(approach6_dir, "results", "hypergeometric_enrichment_results.csv"))

# Plot 1: Fold enrichment barplot
pdf(file.path(approach6_dir, "plots", "hypergeometric_fold_enrichment.pdf"), width = 10, height = 7)

p1 <- ggplot(enrichment_results, aes(x = reorder(GeneSet, -FoldEnrichment),
                                      y = FoldEnrichment, fill = GeneSet)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  geom_text(aes(label = sprintf("%.2fx\n(%d genes)\np=%.2e %s",
                                FoldEnrichment, Overlap, PValue, Significance)),
            vjust = -0.5, size = 3.5) +
  labs(title = "CORRECTED: Enrichment of Migration/Cancer Genes in TES Direct Targets",
       subtitle = "Testing: Are these gene sets enriched among TES targets?",
       x = "Gene Set",
       y = "Fold Enrichment over Expected") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold")) +
  scale_fill_brewer(palette = "Set2") +
  ylim(0, max(enrichment_results$FoldEnrichment) * 1.3)

print(p1)
dev.off()

# Plot 2: -log10(p-value) barplot
pdf(file.path(approach6_dir, "plots", "hypergeometric_pvalue.pdf"), width = 10, height = 7)

p2 <- ggplot(enrichment_results, aes(x = reorder(GeneSet, -NegLog10P),
                                      y = NegLog10P, fill = GeneSet)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "darkred") +
  geom_hline(yintercept = -log10(0.001), linetype = "dashed", color = "darkred", lwd = 1.2) +
  geom_text(aes(label = sprintf("p = %.2e %s", PValue, Significance)),
            vjust = -0.5, size = 3.5) +
  labs(title = "Statistical Significance of Gene Set Enrichment",
       subtitle = "Hypergeometric test (CORRECTED)",
       x = "Gene Set",
       y = "-log10(p-value)") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold")) +
  scale_fill_brewer(palette = "Set2")

print(p2)
dev.off()

################################################################################
# SUMMARY REPORT
################################################################################

summary_text <- paste0(
  "=================================================================\n",
  "APPROACH 6 CORRECTED: MIGRATION GENE-FOCUSED ANALYSIS\n",
  "=================================================================\n\n",
  "IMPORTANT: This is the CORRECTED version.\n",
  "Original Approach 6 tested: Are TES targets enriched in gene sets?\n",
  "CORRECTED version tests: Are gene sets enriched in TES targets?\n\n",
  "Total TES direct targets:", length(tes_direct_genes), "\n",
  "Universe (expressed genes):", universe_size, "\n\n",
  "--- Enrichment Results ---\n"
)

for (i in 1:nrow(enrichment_results)) {
  summary_text <- paste0(
    summary_text,
    sprintf("\n%s:\n", enrichment_results$GeneSet[i]),
    sprintf("  Gene set size: %d\n", enrichment_results$GeneSetSize[i]),
    sprintf("  Overlap with TES targets: %d (%.1f%% of gene set)\n",
            enrichment_results$Overlap[i], enrichment_results$PercentOverlap[i]),
    sprintf("  Expected by chance: %.2f\n", enrichment_results$Expected[i]),
    sprintf("  Fold enrichment: %.2fx\n", enrichment_results$FoldEnrichment[i]),
    sprintf("  P-value: %s %s\n\n", format(enrichment_results$PValue[i], scientific = TRUE),
            enrichment_results$Significance[i])
  )
}

summary_text <- paste0(
  summary_text,
  "=================================================================\n",
  "Analysis complete: ", date(), "\n",
  "=================================================================\n"
)

cat(summary_text)
writeLines(summary_text, file.path(approach6_dir, "results", "analysis_summary.txt"))

cat("\n=================================================================\n")
cat("APPROACH 6 CORRECTED COMPLETE\n")
cat("Results saved to:", approach6_dir, "\n")
cat("=================================================================\n")
