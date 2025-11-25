#!/usr/bin/env Rscript
#
# DIRECTIONAL GO ENRICHMENT: Separate analysis for UP vs DOWN regulated genes
# Addresses Todo #3: Re-analyze RNAseq and do GO enrichment by dividing pathways
# that are upregulated or downregulated

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(dplyr)
  library(ggplot2)
  library(enrichplot)
  library(readr)
})


setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_1")

RNA_data_dir <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_RNA"

cat("=== DIRECTIONAL GO ENRICHMENT ANALYSIS ===\n")
cat("Separate enrichment for UP vs DOWN regulated genes\n")
cat("Analysis started:", as.character(Sys.time()), "\n\n")

# Create output directory
output_dir <- "output/02_directional_go_enrichment"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# PHASE 1: LOAD AND SPLIT DEG DATA
# =============================================================================

cat("=== PHASE 1: Loading and Splitting DEG Data ===\n")

# Load DESeq2 results
deseq_results <- read.delim(
  file.path(RNA_data_dir, "results/05_deseq2/deseq2_results_TES_vs_GFP.txt"),
  stringsAsFactors = FALSE
)

cat(sprintf("✓ Loaded %d genes from DESeq2\n", nrow(deseq_results)))

# Clean Ensembl IDs
deseq_results$ensembl_id <- gsub("\\..*", "", deseq_results$gene_id)

# Identify significant DEGs
deseq_results$is_significant <- !is.na(deseq_results$padj) &
  deseq_results$padj < 0.05

sig_degs <- deseq_results[deseq_results$is_significant, ]

cat(sprintf("✓ Significant DEGs (padj < 0.05): %d\n", nrow(sig_degs)))

# Split by direction
upregulated <- sig_degs[sig_degs$log2FoldChange > 0, ]
downregulated <- sig_degs[sig_degs$log2FoldChange < 0, ]

cat(sprintf("  - Upregulated genes: %d\n", nrow(upregulated)))
cat(sprintf("  - Downregulated genes: %d\n", nrow(downregulated)))
cat(sprintf("  - UP/DOWN ratio: %.2f\n\n", nrow(upregulated) / nrow(downregulated)))

# Export gene lists
write.csv(upregulated, file.path(output_dir, "upregulated_genes.csv"), row.names = FALSE)
write.csv(downregulated, file.path(output_dir, "downregulated_genes.csv"), row.names = FALSE)

# =============================================================================
# PHASE 2: GO ENRICHMENT FOR UPREGULATED GENES
# =============================================================================

cat("=== PHASE 2: GO Enrichment for Upregulated Genes ===\n")

# Background: all genes with padj values (tested genes)
background_genes <- deseq_results$ensembl_id[!is.na(deseq_results$padj)]
background_genes <- background_genes[!is.na(background_genes)]

cat(sprintf("Background: %d genes tested in DESeq2\n\n", length(background_genes)))

# Convert Ensembl to Entrez (GO enrichment uses Entrez)
upregulated$entrez_id <- mapIds(org.Hs.eg.db,
  keys = upregulated$ensembl_id,
  column = "ENTREZID",
  keytype = "ENSEMBL",
  multiVals = "first"
)

downregulated$entrez_id <- mapIds(org.Hs.eg.db,
  keys = downregulated$ensembl_id,
  column = "ENTREZID",
  keytype = "ENSEMBL",
  multiVals = "first"
)

background_entrez <- mapIds(org.Hs.eg.db,
  keys = background_genes,
  column = "ENTREZID",
  keytype = "ENSEMBL",
  multiVals = "first"
)
background_entrez <- background_entrez[!is.na(background_entrez)]

# Get gene lists (remove NAs)
up_genes <- upregulated$entrez_id[!is.na(upregulated$entrez_id)]
down_genes <- downregulated$entrez_id[!is.na(downregulated$entrez_id)]

cat(sprintf("✓ Upregulated: %d genes with Entrez IDs\n", length(up_genes)))
cat(sprintf("✓ Downregulated: %d genes with Entrez IDs\n", length(down_genes)))
cat(sprintf("✓ Background: %d genes with Entrez IDs\n\n", length(background_entrez)))

# Run GO enrichment for upregulated genes
cat("Running GO enrichment for UPREGULATED genes...\n")
up_go_bp <- enrichGO(
  gene = up_genes,
  universe = background_entrez,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE
)

up_go_mf <- enrichGO(
  gene = up_genes,
  universe = background_entrez,
  OrgDb = org.Hs.eg.db,
  ont = "MF",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE
)

up_go_cc <- enrichGO(
  gene = up_genes,
  universe = background_entrez,
  OrgDb = org.Hs.eg.db,
  ont = "CC",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE
)

cat(sprintf("✓ Upregulated GO enrichment:\n"))
cat(sprintf("  - Biological Process: %d pathways\n", nrow(up_go_bp)))
cat(sprintf("  - Molecular Function: %d pathways\n", nrow(up_go_mf)))
cat(sprintf("  - Cellular Component: %d pathways\n\n", nrow(up_go_cc)))

# =============================================================================
# PHASE 3: GO ENRICHMENT FOR DOWNREGULATED GENES
# =============================================================================

cat("=== PHASE 3: GO Enrichment for Downregulated Genes ===\n")

cat("Running GO enrichment for DOWNREGULATED genes...\n")
down_go_bp <- enrichGO(
  gene = down_genes,
  universe = background_entrez,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE
)

down_go_mf <- enrichGO(
  gene = down_genes,
  universe = background_entrez,
  OrgDb = org.Hs.eg.db,
  ont = "MF",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE
)

down_go_cc <- enrichGO(
  gene = down_genes,
  universe = background_entrez,
  OrgDb = org.Hs.eg.db,
  ont = "CC",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE
)

cat(sprintf("✓ Downregulated GO enrichment:\n"))
cat(sprintf("  - Biological Process: %d pathways\n", nrow(down_go_bp)))
cat(sprintf("  - Molecular Function: %d pathways\n", nrow(down_go_mf)))
cat(sprintf("  - Cellular Component: %d pathways\n\n", nrow(down_go_cc)))

# =============================================================================
# PHASE 4: EXPORT RESULTS
# =============================================================================

cat("=== PHASE 4: Exporting Results ===\n")

# Export GO results
if (nrow(up_go_bp) > 0) {
  write.csv(up_go_bp@result, file.path(output_dir, "upregulated_GO_BP.csv"), row.names = FALSE)
}
if (nrow(up_go_mf) > 0) {
  write.csv(up_go_mf@result, file.path(output_dir, "upregulated_GO_MF.csv"), row.names = FALSE)
}
if (nrow(up_go_cc) > 0) {
  write.csv(up_go_cc@result, file.path(output_dir, "upregulated_GO_CC.csv"), row.names = FALSE)
}

if (nrow(down_go_bp) > 0) {
  write.csv(down_go_bp@result, file.path(output_dir, "downregulated_GO_BP.csv"), row.names = FALSE)
}
if (nrow(down_go_mf) > 0) {
  write.csv(down_go_mf@result, file.path(output_dir, "downregulated_GO_MF.csv"), row.names = FALSE)
}
if (nrow(down_go_cc) > 0) {
  write.csv(down_go_cc@result, file.path(output_dir, "downregulated_GO_CC.csv"), row.names = FALSE)
}

cat("✓ GO enrichment results exported\n\n")

# =============================================================================
# PHASE 5: VISUALIZATIONS
# =============================================================================

cat("=== PHASE 5: Creating Visualizations ===\n")

# Plot 1: Dotplot comparison (BP only for clarity)
if (nrow(up_go_bp) > 0 & nrow(down_go_bp) > 0) {
  cat("Creating comparative dotplots...\n")

  pdf(file.path(output_dir, "01_BP_comparison_dotplot.pdf"), width = 16, height = 12)

  p1_up <- dotplot(up_go_bp, showCategory = 20, title = "Upregulated Genes - GO BP")
  print(p1_up)

  p1_down <- dotplot(down_go_bp, showCategory = 20, title = "Downregulated Genes - GO BP")
  print(p1_down)

  dev.off()
}

# Plot 2: Bar plots
if (nrow(up_go_bp) > 0) {
  pdf(file.path(output_dir, "02_upregulated_barplot.pdf"), width = 12, height = 10)
  p2 <- barplot(up_go_bp, showCategory = 30) +
    labs(title = "Top 30 GO BP - Upregulated Genes")
  print(p2)
  dev.off()
}

if (nrow(down_go_bp) > 0) {
  pdf(file.path(output_dir, "03_downregulated_barplot.pdf"), width = 12, height = 10)
  p3 <- barplot(down_go_bp, showCategory = 30) +
    labs(title = "Top 30 GO BP - Downregulated Genes")
  print(p3)
  dev.off()
}

# Plot 3: Comparison plot combining both
if (nrow(up_go_bp) > 10 & nrow(down_go_bp) > 10) {
  cat("Creating combined comparison plot...\n")

  # Extract top 20 from each
  up_top <- up_go_bp@result %>%
    arrange(p.adjust) %>%
    head(20) %>%
    mutate(
      direction = "Upregulated",
      GeneRatio_numeric = sapply(
        strsplit(GeneRatio, "/"),
        function(x) as.numeric(x[1]) / as.numeric(x[2])
      )
    )

  down_top <- down_go_bp@result %>%
    arrange(p.adjust) %>%
    head(20) %>%
    mutate(
      direction = "Downregulated",
      GeneRatio_numeric = sapply(
        strsplit(GeneRatio, "/"),
        function(x) as.numeric(x[1]) / as.numeric(x[2])
      )
    )

  combined <- rbind(up_top, down_top)

  pdf(file.path(output_dir, "04_combined_comparison.pdf"), width = 14, height = 16)
  p4 <- ggplot(combined, aes(
    x = GeneRatio_numeric, y = reorder(Description, GeneRatio_numeric),
    color = p.adjust, size = Count
  )) +
    geom_point() +
    facet_wrap(~direction, ncol = 2, scales = "free_y") +
    scale_color_gradient(low = "red", high = "blue") +
    labs(
      title = "Top 20 GO BP Pathways: Upregulated vs Downregulated",
      x = "Gene Ratio",
      y = "GO Term",
      color = "Adjusted p-value",
      size = "Gene Count"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      axis.text.y = element_text(size = 8)
    )
  print(p4)
  dev.off()
}

# Plot 4: Enrichment map (network view) if possible
if (nrow(up_go_bp) > 20) {
  cat("Creating enrichment map for upregulated pathways...\n")
  pdf(file.path(output_dir, "05_upregulated_enrichment_map.pdf"), width = 12, height = 10)
  tryCatch(
    {
      up_go_simple <- simplify(up_go_bp, cutoff = 0.7, by = "p.adjust")
      p5 <- emapplot(pairwise_termsim(up_go_simple), showCategory = 30) +
        labs(title = "Enrichment Map - Upregulated Genes")
      print(p5)
    },
    error = function(e) {
      cat("Note: Enrichment map could not be generated\n")
    }
  )
  dev.off()
}

if (nrow(down_go_bp) > 20) {
  cat("Creating enrichment map for downregulated pathways...\n")
  pdf(file.path(output_dir, "06_downregulated_enrichment_map.pdf"), width = 12, height = 10)
  tryCatch(
    {
      down_go_simple <- simplify(down_go_bp, cutoff = 0.7, by = "p.adjust")
      p6 <- emapplot(pairwise_termsim(down_go_simple), showCategory = 30) +
        labs(title = "Enrichment Map - Downregulated Genes")
      print(p6)
    },
    error = function(e) {
      cat("Note: Enrichment map could not be generated\n")
    }
  )
  dev.off()
}

cat("✓ All visualizations created\n\n")

# =============================================================================
# PHASE 6: SUMMARY REPORT
# =============================================================================

cat("=== PHASE 6: Generating Summary Report ===\n")

summary_file <- file.path(output_dir, "DIRECTIONAL_GO_SUMMARY.txt")
cat("DIRECTIONAL GO ENRICHMENT ANALYSIS SUMMARY\n", file = summary_file)
cat("==========================================\n\n", file = summary_file, append = TRUE)
cat(paste("Generated:", Sys.time(), "\n\n"), file = summary_file, append = TRUE)

cat("GENE COUNTS:\n", file = summary_file, append = TRUE)
cat(sprintf("  Total significant DEGs: %d\n", nrow(sig_degs)), file = summary_file, append = TRUE)
cat(
  sprintf(
    "    - Upregulated (log2FC > 0): %d (%.1f%%)\n",
    nrow(upregulated),
    100 * nrow(upregulated) / nrow(sig_degs)
  ),
  file = summary_file, append = TRUE
)
cat(
  sprintf(
    "    - Downregulated (log2FC < 0): %d (%.1f%%)\n\n",
    nrow(downregulated),
    100 * nrow(downregulated) / nrow(sig_degs)
  ),
  file = summary_file, append = TRUE
)

cat("GO ENRICHMENT RESULTS:\n", file = summary_file, append = TRUE)
cat("\nUpregulated Genes:\n", file = summary_file, append = TRUE)
cat(sprintf("  - Biological Process: %d pathways\n", nrow(up_go_bp)), file = summary_file, append = TRUE)
cat(sprintf("  - Molecular Function: %d pathways\n", nrow(up_go_mf)), file = summary_file, append = TRUE)
cat(sprintf("  - Cellular Component: %d pathways\n", nrow(up_go_cc)), file = summary_file, append = TRUE)

cat("\nDownregulated Genes:\n", file = summary_file, append = TRUE)
cat(sprintf("  - Biological Process: %d pathways\n", nrow(down_go_bp)), file = summary_file, append = TRUE)
cat(sprintf("  - Molecular Function: %d pathways\n", nrow(down_go_mf)), file = summary_file, append = TRUE)
cat(sprintf("  - Cellular Component: %d pathways\n", nrow(down_go_cc)), file = summary_file, append = TRUE)

cat("\n\nTOP 10 UPREGULATED PATHWAYS (BP):\n", file = summary_file, append = TRUE)
if (nrow(up_go_bp) > 0) {
  top_up <- up_go_bp@result %>%
    arrange(p.adjust) %>%
    head(10)
  for (i in 1:min(nrow(top_up), 10)) {
    cat(
      sprintf(
        "  %d. %s (FDR=%.2e, Count=%d)\n",
        i, top_up$Description[i], top_up$p.adjust[i], top_up$Count[i]
      ),
      file = summary_file, append = TRUE
    )
  }
}

cat("\n\nTOP 10 DOWNREGULATED PATHWAYS (BP):\n", file = summary_file, append = TRUE)
if (nrow(down_go_bp) > 0) {
  top_down <- down_go_bp@result %>%
    arrange(p.adjust) %>%
    head(10)
  for (i in 1:min(nrow(top_down), 10)) {
    cat(
      sprintf(
        "  %d. %s (FDR=%.2e, Count=%d)\n",
        i, top_down$Description[i], top_down$p.adjust[i], top_down$Count[i]
      ),
      file = summary_file, append = TRUE
    )
  }
}

cat("\n✓ Summary report saved\n\n")

cat("========================================\n")
cat("DIRECTIONAL GO ENRICHMENT COMPLETE\n")
cat("========================================\n")
cat("Completed:", as.character(Sys.time()), "\n")
cat(sprintf("Output directory: %s\n\n", output_dir))
cat("Key files:\n")
cat("  - upregulated_genes.csv\n")
cat("  - downregulated_genes.csv\n")
cat("  - upregulated_GO_BP/MF/CC.csv\n")
cat("  - downregulated_GO_BP/MF/CC.csv\n")
cat("  - Multiple visualization PDFs\n")
cat("  - DIRECTIONAL_GO_SUMMARY.txt\n")
