#!/usr/bin/env Rscript

cat("=================================================================\n")
cat("APPROACH 7: Glioblastoma & Cancer Marker Analysis\n")
cat("Focus: Cancer-specific genes and GBM invasion markers\n")
cat("=================================================================\n")

# Load workspace and helper functions
source("SRF_Eva_integrated_analysis/focused_enrichment_analysis/scripts/00_helper_functions.R")
load_libraries()
setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top")
load("SRF_Eva_integrated_analysis/focused_enrichment_analysis/data/loaded_data.RData")

# Output directory
base_dir <- "SRF_Eva_integrated_analysis/focused_enrichment_analysis"
approach7_dir <- file.path(base_dir, "approach7_glioblastoma")

# Create directories
dir.create(file.path(approach7_dir, "results"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(approach7_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(approach7_dir, "gene_lists"), recursive = TRUE, showWarnings = FALSE)

################################################################################
# LOAD CANCER/GLIOBLASTOMA GENE SETS
################################################################################

cat("\n=== Loading cancer and glioblastoma gene sets ===\n")

# 1. MSigDB Hallmark cancer-related gene sets
cat("Loading MSigDB Hallmark cancer gene sets...\n")
hallmark_cancer <- msigdb_hallmark %>%
  filter(grepl("EPITHELIAL_MESENCHYMAL_TRANSITION|ANGIOGENESIS|METASTASIS|HYPOXIA",
               gs_name, ignore.case = TRUE))

emt_genes <- msigdb_hallmark %>%
  filter(gs_name == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION") %>%
  pull(gene_symbol) %>% unique()

angiogenesis_genes <- msigdb_hallmark %>%
  filter(gs_name == "HALLMARK_ANGIOGENESIS") %>%
  pull(gene_symbol) %>% unique()

hypoxia_genes <- msigdb_hallmark %>%
  filter(gs_name == "HALLMARK_HYPOXIA") %>%
  pull(gene_symbol) %>% unique()

cat("  EMT genes:", length(emt_genes), "\n")
cat("  Angiogenesis genes:", length(angiogenesis_genes), "\n")
cat("  Hypoxia genes:", length(hypoxia_genes), "\n")

# 2. KEGG and Reactome glioma pathways
cat("\nLoading KEGG/Reactome glioma pathways...\n")
glioma_kegg <- msigdb_c2 %>%
  filter(grepl("GLIOMA", gs_name, ignore.case = TRUE)) %>%
  pull(gene_symbol) %>% unique()

cancer_pathways <- msigdb_c2 %>%
  filter(grepl("CANCER|METASTASIS|INVASION", gs_name, ignore.case = TRUE)) %>%
  pull(gene_symbol) %>% unique()

cat("  KEGG Glioma genes:", length(glioma_kegg), "\n")
cat("  Cancer pathway genes:", length(cancer_pathways), "\n")

# 3. Curated GBM invasion/migration markers from literature
cat("\nCurating GBM invasion markers from literature...\n")
gbm_invasion_markers <- c(
  # Matrix metalloproteinases
  "MMP2", "MMP9", "MMP14", "MMP1", "MMP7",
  # Receptor tyrosine kinases
  "AXL", "EGFR", "MET", "PDGFRA", "PDGFRB",
  # VEGF signaling
  "VEGFA", "VEGFB", "VEGFC", "KDR", "FLT1",
  # EMT transcription factors
  "TWIST1", "TWIST2", "SNAI1", "SNAI2", "ZEB1", "ZEB2",
  # Mesenchymal markers
  "VIM", "FN1", "CDH2", "ACTA2",
  # Invasion/migration genes
  "ITGAV", "ITGB3", "ITGB1", "ITGA6", "CD44",
  # Integrins and adhesion
  "VCAM1", "ICAM1", "PECAM1",
  # GBM stem cell markers
  "SOX2", "OLIG2", "NES", "PROM1", "BMI1",
  # YAP/TAZ targets known in GBM
  "CTGF", "CYR61", "ANKRD1", "BIRC5", "CCN1", "CCN2",
  # Other invasion genes
  "CXCR4", "CXCL12", "SPP1", "SERPINE1", "PLAU", "PLAUR",
  # Proteases
  "ADAM10", "ADAM17", "ADAM9",
  # FAK signaling
  "PTK2", "PTK2B", "SRC", "PXN",
  # TGF-beta pathway
  "TGFB1", "TGFB2", "TGFBR1", "TGFBR2", "SMAD2", "SMAD3", "SMAD4",
  # Wnt signaling
  "WNT5A", "WNT5B", "FZD2", "DVL2",
  # Notch signaling
  "NOTCH1", "JAG1", "HES1", "HEY1"
)

cat("  Curated GBM invasion markers:", length(gbm_invasion_markers), "\n")

# 4. Mesenchymal GBM subtype markers (from TCGA)
mesenchymal_gbm_markers <- c(
  "CHI3L1", "CD44", "SERPINE1", "TIMP1", "VIM", "FN1",
  "TGFBI", "MET", "AXL", "RELB", "STAT3", "FOSL2",
  "TRADD", "TNFRSF1A", "TAB2", "RELB", "CASP1", "CASP8",
  "C/EBPβ", "TAZ", "BHLHE40"
)

cat("  Mesenchymal GBM markers:", length(mesenchymal_gbm_markers), "\n")

################################################################################
# TES DIRECT TARGETS
################################################################################

tes_direct_genes <- tes_direct$gene_symbol
cat("\nTES direct targets:", length(tes_direct_genes), "\n")

################################################################################
# ENRICHMENT TESTING: HYPERGEOMETRIC TESTS
################################################################################

cat("\n=== Hypergeometric enrichment tests ===\n")

# Universe: all expressed genes
universe_size <- length(all_expressed_genes)

# Function to perform hypergeometric test
test_enrichment <- function(gene_set, gene_set_name) {
  # Genes in gene set that are in universe
  geneset_in_universe <- intersect(gene_set, all_expressed_genes)
  m <- length(geneset_in_universe)

  # Overlap with TES direct targets
  overlap <- intersect(tes_direct_genes, geneset_in_universe)
  k <- length(tes_direct_genes)
  q <- length(overlap)

  if (q == 0) {
    cat(sprintf("  %s: No overlap\n", gene_set_name))
    return(NULL)
  }

  # Hypergeometric test
  pval <- phyper(q = q - 1,
                 m = m,
                 n = universe_size - m,
                 k = k,
                 lower.tail = FALSE)

  # Expected by chance
  expected <- (k * m) / universe_size

  # Fold enrichment
  fold_enrich <- q / expected

  cat(sprintf("  %s:\n", gene_set_name))
  cat(sprintf("    Gene set size (in universe): %d\n", m))
  cat(sprintf("    Overlap with TES targets: %d\n", q))
  cat(sprintf("    Expected by chance: %.2f\n", expected))
  cat(sprintf("    Fold enrichment: %.2f\n", fold_enrich))
  cat(sprintf("    p-value: %s\n", format(pval, scientific = TRUE)))

  # Save gene list
  if (q > 0) {
    write.table(overlap,
                file.path(approach7_dir, "gene_lists",
                         paste0("TES_direct_", gsub("[^A-Za-z0-9_]", "_", gene_set_name), ".txt")),
                row.names = FALSE, col.names = FALSE, quote = FALSE)
  }

  return(data.frame(
    GeneSet = gene_set_name,
    GeneSetSize = m,
    Overlap = q,
    Expected = expected,
    FoldEnrichment = fold_enrich,
    PValue = pval,
    stringsAsFactors = FALSE
  ))
}

# Test each gene set
enrichment_results <- bind_rows(
  test_enrichment(emt_genes, "Hallmark_EMT"),
  test_enrichment(angiogenesis_genes, "Hallmark_Angiogenesis"),
  test_enrichment(hypoxia_genes, "Hallmark_Hypoxia"),
  test_enrichment(glioma_kegg, "KEGG_Glioma"),
  test_enrichment(cancer_pathways, "Cancer_Pathways_Reactome"),
  test_enrichment(gbm_invasion_markers, "GBM_Invasion_Markers_Curated"),
  test_enrichment(mesenchymal_gbm_markers, "Mesenchymal_GBM_Subtype")
) %>%
  mutate(
    NegLog10P = -log10(PValue),
    BH_adjusted = p.adjust(PValue, method = "BH"),
    Significance = case_when(
      PValue < 0.001 ~ "***",
      PValue < 0.01 ~ "**",
      PValue < 0.05 ~ "*",
      TRUE ~ "n.s."
    )
  ) %>%
  arrange(PValue)

# Save enrichment results
write_csv(enrichment_results,
          file.path(approach7_dir, "results", "cancer_gbm_enrichment_results.csv"))

cat("\n=== Enrichment Summary ===\n")
print(enrichment_results)

################################################################################
# GSEA ANALYSIS WITH CANCER GENE SETS
################################################################################

cat("\n=== GSEA analysis on cancer gene sets ===\n")

# Create ranked list by log2FoldChange * -log10(padj)
degs_ranked <- deseq_results %>%
  filter(!is.na(padj) & !is.na(log2FoldChange) & baseMean > 10) %>%
  mutate(padj_capped = pmax(padj, 1e-300)) %>%
  mutate(rank_metric = log2FoldChange * -log10(padj_capped)) %>%
  filter(is.finite(rank_metric)) %>%
  arrange(desc(rank_metric))

ranked_genes <- setNames(degs_ranked$rank_metric, degs_ranked$gene_symbol)

# Combine cancer gene sets
cancer_genesets <- bind_rows(
  msigdb_hallmark %>% filter(grepl("EPITHELIAL_MESENCHYMAL_TRANSITION|ANGIOGENESIS|HYPOXIA", gs_name)),
  msigdb_c2 %>% filter(grepl("GLIOMA|GLIOBLASTOMA", gs_name, ignore.case = TRUE))
)

if (nrow(cancer_genesets) > 0) {
  cat("Running GSEA on", length(unique(cancer_genesets$gs_name)), "cancer gene sets...\n")

  gsea_cancer <- perform_GSEA_msigdb(ranked_genes, cancer_genesets, pvalueCutoff = 0.25)

  if (!is.null(gsea_cancer) && nrow(gsea_cancer) > 0) {
    gsea_results <- as.data.frame(gsea_cancer)
    write_csv(gsea_results,
              file.path(approach7_dir, "results", "GSEA_cancer_genesets.csv"))

    cat("\nGSEA results:\n")
    print(gsea_results[, c("ID", "NES", "pvalue", "p.adjust")])

    # Plot GSEA results
    pdf(file.path(approach7_dir, "plots", "GSEA_cancer_dotplot.pdf"), width = 12, height = 8)
    print(dotplot(gsea_cancer, showCategory = 20, title = "GSEA: Cancer & Glioblastoma Gene Sets"))
    dev.off()

    # Individual GSEA plots for significant gene sets
    for (geneset_id in head(gsea_results$ID, 10)) {
      safe_name <- gsub("[^A-Za-z0-9_]", "_", geneset_id)
      pdf(file.path(approach7_dir, "plots", paste0("GSEA_", safe_name, ".pdf")),
          width = 10, height = 6)
      print(gseaplot2(gsea_cancer,
                      geneSetID = geneset_id,
                      title = geneset_id))
      dev.off()
    }
  } else {
    cat("No significant GSEA results found (p < 0.25)\n")
  }
} else {
  cat("Could not load cancer gene sets for GSEA\n")
}

################################################################################
# GO ENRICHMENT FOR GBM INVASION MARKERS
################################################################################

cat("\n=== GO enrichment for TES targets overlapping with GBM markers ===\n")

gbm_overlap <- intersect(tes_direct_genes, gbm_invasion_markers)

if (length(gbm_overlap) >= 10) {
  cat("Running GO enrichment on", length(gbm_overlap), "GBM marker genes targeted by TES...\n")

  ego_gbm_overlap <- perform_GO_enrichment(
    gene_list = gbm_overlap,
    background_genes = all_expressed_genes,
    ont = "BP",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05
  )

  if (!is.null(ego_gbm_overlap) && nrow(ego_gbm_overlap) > 0) {
    save_enrichment_results(ego_gbm_overlap, "GBM_markers_targeted_by_TES", approach7_dir)
    cat("Found", nrow(ego_gbm_overlap), "enriched GO terms\n")

    # Extract migration terms
    migration_terms <- extract_migration_terms(as.data.frame(ego_gbm_overlap))
    if (!is.null(migration_terms) && nrow(migration_terms) > 0) {
      write_csv(migration_terms,
                file.path(approach7_dir, "results", "GBM_markers_migration_terms.csv"))
      cat("  Migration-related terms:", nrow(migration_terms), "\n")
    }
  }
} else {
  cat("Too few overlapping genes for GO enrichment (", length(gbm_overlap), ")\n")
}

################################################################################
# DETAILED ANNOTATION FOR KEY CANCER/GBM GENES
################################################################################

cat("\n=== Creating detailed tables for key cancer/GBM genes ===\n")

# For each significant gene set, create detailed table
for (i in 1:nrow(enrichment_results)) {
  geneset_name <- enrichment_results$GeneSet[i]

  # Get gene list
  if (grepl("EMT", geneset_name)) {
    genes_in_set <- emt_genes
  } else if (grepl("Angiogenesis", geneset_name)) {
    genes_in_set <- angiogenesis_genes
  } else if (grepl("Hypoxia", geneset_name)) {
    genes_in_set <- hypoxia_genes
  } else if (grepl("Glioma", geneset_name)) {
    genes_in_set <- glioma_kegg
  } else if (grepl("Cancer_Path", geneset_name)) {
    genes_in_set <- cancer_pathways
  } else if (grepl("Invasion", geneset_name)) {
    genes_in_set <- gbm_invasion_markers
  } else if (grepl("Mesenchymal", geneset_name)) {
    genes_in_set <- mesenchymal_gbm_markers
  } else {
    next
  }

  # Find overlap
  overlap_genes <- intersect(tes_direct_genes, genes_in_set)

  if (length(overlap_genes) > 0) {
    detailed_table <- tes_direct %>%
      filter(gene_symbol %in% overlap_genes) %>%
      arrange(padj) %>%
      select(gene_symbol, baseMean, log2FoldChange, padj, tes_bound, tead1_bound)

    write_csv(detailed_table,
              file.path(approach7_dir, "gene_lists",
                       paste0(gsub("[^A-Za-z0-9_]", "_", geneset_name), "_detailed.csv")))
  }
}

################################################################################
# VISUALIZATIONS
################################################################################

cat("\n=== Creating visualizations ===\n")

# Plot 1: Enrichment comparison barplot
pdf(file.path(approach7_dir, "plots", "cancer_gbm_enrichment_barplot.pdf"), width = 12, height = 8)

p1 <- ggplot(enrichment_results, aes(x = reorder(GeneSet, -FoldEnrichment),
                                      y = FoldEnrichment, fill = GeneSet)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  geom_text(aes(label = sprintf("%d genes\n(%.1fx, p=%.2e)",
                                Overlap, FoldEnrichment, PValue)),
            vjust = -0.5, size = 3) +
  labs(title = "TES Direct Target Enrichment in Cancer/Glioblastoma Gene Sets",
       subtitle = sprintf("Out of %d total TES direct targets", length(tes_direct_genes)),
       x = "Gene Set",
       y = "Fold Enrichment over Expected") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        plot.title = element_text(size = 14, face = "bold")) +
  scale_fill_brewer(palette = "Set2") +
  ylim(0, max(enrichment_results$FoldEnrichment) * 1.2)

print(p1)
dev.off()

# Plot 2: -log10(p-value) barplot
pdf(file.path(approach7_dir, "plots", "cancer_gbm_pvalue_barplot.pdf"), width = 12, height = 8)

p2 <- ggplot(enrichment_results, aes(x = reorder(GeneSet, -NegLog10P),
                                      y = NegLog10P, fill = GeneSet)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "darkred") +
  geom_hline(yintercept = -log10(0.001), linetype = "dashed", color = "darkred", lwd = 1.2) +
  geom_text(aes(label = Significance), vjust = -0.5, size = 5) +
  labs(title = "Statistical Significance of Cancer/GBM Gene Set Enrichment",
       subtitle = "Hypergeometric test",
       x = "Gene Set",
       y = "-log10(p-value)") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        plot.title = element_text(size = 14, face = "bold")) +
  scale_fill_brewer(palette = "Set2")

print(p2)
dev.off()

# Plot 3: Bubble plot combining fold enrichment and p-value
pdf(file.path(approach7_dir, "plots", "cancer_gbm_bubble_plot.pdf"), width = 12, height = 8)

p3 <- ggplot(enrichment_results,
             aes(x = FoldEnrichment, y = NegLog10P, size = Overlap, color = GeneSet)) +
  geom_point(alpha = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", alpha = 0.5) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50", alpha = 0.5) +
  geom_text_repel(aes(label = GeneSet), size = 3.5, max.overlaps = 20) +
  scale_size_continuous(name = "Overlap\n(# genes)", range = c(3, 15)) +
  scale_color_brewer(palette = "Set2") +
  labs(title = "Cancer/Glioblastoma Gene Set Enrichment in TES Direct Targets",
       subtitle = "Bubble size = number of overlapping genes",
       x = "Fold Enrichment",
       y = "-log10(p-value)") +
  theme_bw() +
  theme(legend.position = "right",
        plot.title = element_text(size = 14, face = "bold"))

print(p3)
dev.off()

# Plot 4: Heatmap of top GBM invasion markers
gbm_overlap_top <- tes_direct %>%
  filter(gene_symbol %in% gbm_invasion_markers) %>%
  arrange(padj) %>%
  head(30)

if (nrow(gbm_overlap_top) > 0) {
  pdf(file.path(approach7_dir, "plots", "top_gbm_invasion_markers_heatmap.pdf"),
      width = 8, height = 10)

  p4 <- ggplot(gbm_overlap_top,
               aes(x = 1, y = reorder(gene_symbol, -log2FoldChange), fill = log2FoldChange)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", log2FoldChange)), size = 3) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                         name = "log2FC\n(TES vs GFP)") +
    labs(title = "Top 30 GBM Invasion Markers Targeted by TES",
         subtitle = "Curated from literature, ranked by p-value",
         x = "", y = "Gene Symbol") +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.title = element_text(size = 14, face = "bold"))

  print(p4)
  dev.off()
}

################################################################################
# SUMMARY REPORT
################################################################################

cat("\n=== Generating summary report ===\n")

summary_text <- paste0(
  "=================================================================\n",
  "APPROACH 7: GLIOBLASTOMA & CANCER MARKER ANALYSIS SUMMARY\n",
  "=================================================================\n\n",
  "Total TES direct targets analyzed: ", length(tes_direct_genes), "\n",
  "Background universe (expressed genes): ", universe_size, "\n\n",
  "--- Gene Set Enrichment Results ---\n"
)

for (i in 1:nrow(enrichment_results)) {
  summary_text <- paste0(
    summary_text,
    sprintf("\n%d. %s:\n", i, enrichment_results$GeneSet[i]),
    sprintf("   Overlap: %d / %d genes (%.1f%%)\n",
            enrichment_results$Overlap[i],
            enrichment_results$GeneSetSize[i],
            enrichment_results$Overlap[i]/enrichment_results$GeneSetSize[i]*100),
    sprintf("   Fold enrichment: %.2fx\n", enrichment_results$FoldEnrichment[i]),
    sprintf("   P-value: %s %s\n",
            format(enrichment_results$PValue[i], scientific = TRUE),
            enrichment_results$Significance[i]),
    sprintf("   Adj. p-value (BH): %s\n", format(enrichment_results$BH_adjusted[i], scientific = TRUE))
  )
}

summary_text <- paste0(
  summary_text,
  "\n=================================================================\n",
  "KEY FINDINGS:\n",
  "=================================================================\n",
  "1. TES directly targets ", enrichment_results$Overlap[enrichment_results$GeneSet == "GBM_Invasion_Markers_Curated"],
  " of ", enrichment_results$GeneSetSize[enrichment_results$GeneSet == "GBM_Invasion_Markers_Curated"],
  " curated GBM invasion markers\n",
  "2. Strongest enrichment: ", enrichment_results$GeneSet[1],
  " (", sprintf("%.2fx", enrichment_results$FoldEnrichment[1]), ")\n",
  "3. All ", sum(enrichment_results$PValue < 0.05), " gene sets show significant enrichment (p < 0.05)\n\n",
  "Analysis complete: ", date(), "\n",
  "=================================================================\n"
)

cat(summary_text)
writeLines(summary_text, file.path(approach7_dir, "results", "analysis_summary.txt"))

cat("\n=================================================================\n")
cat("APPROACH 7 COMPLETE\n")
cat("Results saved to:", approach7_dir, "\n")
cat("=================================================================\n")
