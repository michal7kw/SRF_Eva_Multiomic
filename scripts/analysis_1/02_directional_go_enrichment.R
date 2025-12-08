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
# PHASE 4b: FILTER FOR CANCER-RELEVANT PATHWAYS
# =============================================================================

cat("=== PHASE 4b: Identifying Cancer-Relevant Pathways ===\n")

# Define cancer-relevant keywords (same as GSEA script for consistency)
cancer_keywords <- c(
  # Cell death pathways
  "apoptosis", "apoptotic", "cell death", "programmed cell death",
  "necrosis", "necrotic", "ferroptosis", "pyroptosis", "anoikis",
  "autophagy", "autophagic", "survival", "viability", "senescence",

  # Cell cycle and proliferation
  "proliferation", "proliferative", "cell cycle", "mitosis", "mitotic",
  "cell division", "growth", "G1/S", "G2/M", "S phase", "M phase",
  "DNA replication", "chromosome", "spindle", "cytokinesis",
  "cyclin", "checkpoint", "DNA repair", "DNA damage",

  # Migration and invasion
  "migration", "migratory", "motility", "invasion", "invasive",
  "chemotaxis", "chemotactic", "cell movement", "locomotion",
  "adhesion", "cell adhesion", "focal adhesion", "cell junction",
  "cytoskeleton", "actin", "tubulin", "microtubule",

  # Angiogenesis and vasculature
  "angiogenesis", "angiogenic", "blood vessel", "vasculature",
  "endothelial", "VEGF", "vascular",

  # Signaling pathways
  "signaling", "signal transduction", "kinase", "phosphorylation",
  "growth factor", "receptor", "activation", "cascade",

  # Metabolism
  "glycolysis", "metabolism", "metabolic", "glucose", "ATP",
  "oxidative", "respiration", "biosynthesis", "catabolic",

  # Transcription and chromatin
  "transcription", "gene expression", "chromatin", "histone",
  "RNA processing", "splicing", "translation",

  # Cancer-specific terms
  "tumor", "cancer", "oncogenic", "transformation",
  "EMT", "epithelial", "mesenchymal", "stemness",
  "Hippo", "YAP", "TEAD", "Wnt", "Notch", "Hedgehog",

  # Stress response
  "stress", "oxidative stress", "hypoxia", "ER stress",
  "unfolded protein", "heat shock", "inflammatory"
)

cancer_pattern <- paste(cancer_keywords, collapse = "|")

# Function to filter for cancer-relevant pathways
filter_cancer_pathways <- function(go_result) {
  if (is.null(go_result) || nrow(go_result) == 0) {
    return(NULL)
  }
  result_df <- go_result@result
  cancer_df <- result_df[grepl(cancer_pattern, result_df$Description, ignore.case = TRUE), ]
  return(cancer_df)
}

# Filter cancer-relevant pathways
up_cancer_bp <- filter_cancer_pathways(up_go_bp)
down_cancer_bp <- filter_cancer_pathways(down_go_bp)

cat(sprintf("✓ Cancer-relevant pathways identified:\n"))
cat(sprintf("  - Upregulated BP: %d (of %d total)\n",
  ifelse(is.null(up_cancer_bp), 0, nrow(up_cancer_bp)), nrow(up_go_bp)))
cat(sprintf("  - Downregulated BP: %d (of %d total)\n\n",
  ifelse(is.null(down_cancer_bp), 0, nrow(down_cancer_bp)), nrow(down_go_bp)))

# Export cancer-relevant pathways
if (!is.null(up_cancer_bp) && nrow(up_cancer_bp) > 0) {
  write.csv(up_cancer_bp, file.path(output_dir, "upregulated_GO_BP_cancer_relevant.csv"), row.names = FALSE)
}
if (!is.null(down_cancer_bp) && nrow(down_cancer_bp) > 0) {
  write.csv(down_cancer_bp, file.path(output_dir, "downregulated_GO_BP_cancer_relevant.csv"), row.names = FALSE)
}

# =============================================================================
# PHASE 5: VISUALIZATIONS
# =============================================================================

cat("=== PHASE 5: Creating Visualizations ===\n")

# =============================================================================
# Plot 1a: TOP PATHWAYS REGARDLESS OF SIGNIFICANCE (for comprehensive view)
# =============================================================================
cat("Creating Plot 1: Top pathways (regardless of significance)...\n")

# Function to create a manual dotplot from enrichResult @result slot
create_top_dotplot <- function(enrich_result, n_show = 20, title = "GO Enrichment") {
  if (is.null(enrich_result) || nrow(enrich_result@result) == 0) {
    return(NULL)
  }

  df <- enrich_result@result %>%
    arrange(p.adjust) %>%
    head(n_show) %>%
    mutate(
      GeneRatio_numeric = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])),
      Description = factor(Description, levels = rev(Description))
    )

  p <- ggplot(df, aes(x = GeneRatio_numeric, y = Description, color = p.adjust, size = Count)) +
    geom_point() +
    scale_color_gradient(low = "red", high = "blue", name = "Adjusted\np-value") +
    scale_size_continuous(name = "Gene\nCount", range = c(3, 10)) +
    labs(
      title = title,
      subtitle = sprintf("Top %d pathways by adjusted p-value", n_show),
      x = "Gene Ratio",
      y = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5),
      axis.text.y = element_text(size = 9)
    )

  return(p)
}

# Create dotplots showing TOP pathways (regardless of q-value cutoff)
if (nrow(up_go_bp@result) > 0 || nrow(down_go_bp@result) > 0) {
  pdf(file.path(output_dir, "01_BP_comparison_dotplot.pdf"), width = 14, height = 12)

  if (nrow(up_go_bp@result) > 0) {
    p_up <- create_top_dotplot(up_go_bp, n_show = 25,
      title = "Upregulated Genes - GO Biological Process")
    if (!is.null(p_up)) print(p_up)
  }

  if (nrow(down_go_bp@result) > 0) {
    p_down <- create_top_dotplot(down_go_bp, n_show = 25,
      title = "Downregulated Genes - GO Biological Process")
    if (!is.null(p_down)) print(p_down)
  }

  dev.off()
  cat("  ✓ Saved: 01_BP_comparison_dotplot.pdf\n")
}

# =============================================================================
# Plot 1b: SIGNIFICANT PATHWAYS ONLY (using enrichplot dotplot)
# =============================================================================
cat("Creating Plot 1b: Significant pathways only (q < 0.05)...\n")

# Count significant results
n_sig_up <- sum(up_go_bp@result$p.adjust < 0.05)
n_sig_down <- sum(down_go_bp@result$p.adjust < 0.05)

cat(sprintf("  Significant GO BP terms (q < 0.05): UP=%d, DOWN=%d\n", n_sig_up, n_sig_down))

if (n_sig_up > 0 || n_sig_down > 0) {
  pdf(file.path(output_dir, "01b_BP_significant_dotplot.pdf"), width = 14, height = 12)

  if (n_sig_up > 0) {
    p1_up <- dotplot(up_go_bp, showCategory = min(25, n_sig_up),
      title = sprintf("Upregulated Genes - GO BP (n=%d significant)", n_sig_up))
    print(p1_up)
  }

  if (n_sig_down > 0) {
    p1_down <- dotplot(down_go_bp, showCategory = min(25, n_sig_down),
      title = sprintf("Downregulated Genes - GO BP (n=%d significant)", n_sig_down))
    print(p1_down)
  }

  dev.off()
  cat("  ✓ Saved: 01b_BP_significant_dotplot.pdf\n")
} else {
  cat("  ⚠ No significant pathways at q < 0.05, skipping 01b plot\n")
}

# Plot 2: Bar plots (Top pathways regardless of significance)
cat("Creating barplots...\n")

# Function to create manual barplot from enrichResult
create_top_barplot <- function(enrich_result, n_show = 30, title = "GO Enrichment") {
  if (is.null(enrich_result) || nrow(enrich_result@result) == 0) {
    return(NULL)
  }

  df <- enrich_result@result %>%
    arrange(p.adjust) %>%
    head(n_show) %>%
    mutate(
      Description = factor(Description, levels = rev(Description)),
      neg_log_padj = -log10(p.adjust),
      is_significant = p.adjust < 0.05
    )

  p <- ggplot(df, aes(x = neg_log_padj, y = Description, fill = is_significant)) +
    geom_bar(stat = "identity") +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "red", linewidth = 0.8) +
    scale_fill_manual(values = c("TRUE" = "#E31A1C", "FALSE" = "#1F78B4"),
      labels = c("Not significant", "Significant (q<0.05)"),
      name = "Significance") +
    labs(
      title = title,
      subtitle = sprintf("Top %d pathways (red line = q=0.05)", n_show),
      x = "-log10(Adjusted p-value)",
      y = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5),
      axis.text.y = element_text(size = 8)
    )

  return(p)
}

if (nrow(up_go_bp@result) > 0) {
  pdf(file.path(output_dir, "02_upregulated_barplot.pdf"), width = 14, height = 12)
  p2 <- create_top_barplot(up_go_bp, n_show = 30, title = "Top 30 GO BP - Upregulated Genes")
  if (!is.null(p2)) print(p2)
  dev.off()
  cat("  ✓ Saved: 02_upregulated_barplot.pdf\n")
}

if (nrow(down_go_bp@result) > 0) {
  pdf(file.path(output_dir, "03_downregulated_barplot.pdf"), width = 14, height = 12)
  p3 <- create_top_barplot(down_go_bp, n_show = 30, title = "Top 30 GO BP - Downregulated Genes")
  if (!is.null(p3)) print(p3)
  dev.off()
  cat("  ✓ Saved: 03_downregulated_barplot.pdf\n")
}

# Plot 3: Comparison plot combining both (use @result to get all pathways)
if (nrow(up_go_bp@result) > 0 & nrow(down_go_bp@result) > 0) {
  cat("Creating combined comparison plot...\n")

  # Extract top 20 from each (regardless of significance)
  up_top <- up_go_bp@result %>%
    arrange(p.adjust) %>%
    head(20) %>%
    mutate(
      direction = "Upregulated",
      GeneRatio_numeric = sapply(
        strsplit(GeneRatio, "/"),
        function(x) as.numeric(x[1]) / as.numeric(x[2])
      ),
      is_significant = p.adjust < 0.05,
      # Truncate long descriptions
      Description_short = ifelse(nchar(Description) > 45,
        paste0(substr(Description, 1, 42), "..."),
        Description)
    )

  down_top <- down_go_bp@result %>%
    arrange(p.adjust) %>%
    head(20) %>%
    mutate(
      direction = "Downregulated",
      GeneRatio_numeric = sapply(
        strsplit(GeneRatio, "/"),
        function(x) as.numeric(x[1]) / as.numeric(x[2])
      ),
      is_significant = p.adjust < 0.05,
      Description_short = ifelse(nchar(Description) > 45,
        paste0(substr(Description, 1, 42), "..."),
        Description)
    )

  combined <- rbind(up_top, down_top)

  pdf(file.path(output_dir, "04_combined_comparison.pdf"), width = 16, height = 14)
  p4 <- ggplot(combined, aes(
    x = GeneRatio_numeric, y = reorder(Description_short, GeneRatio_numeric),
    color = p.adjust, size = Count, shape = is_significant
  )) +
    geom_point() +
    facet_wrap(~direction, ncol = 2, scales = "free_y") +
    scale_color_gradient(low = "red", high = "blue", name = "Adjusted\np-value") +
    scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
      labels = c("q >= 0.05", "q < 0.05"),
      name = "Significance") +
    scale_size_continuous(name = "Gene\nCount", range = c(3, 8)) +
    labs(
      title = "Top 20 GO BP Pathways: Upregulated vs Downregulated",
      subtitle = "Filled circles = significant (q < 0.05), open circles = not significant",
      x = "Gene Ratio",
      y = "GO Term"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5),
      axis.text.y = element_text(size = 8),
      strip.text = element_text(face = "bold", size = 12)
    )
  print(p4)
  dev.off()
  cat("  ✓ Saved: 04_combined_comparison.pdf\n")
}

# Plot 4: Enrichment map (network view) if possible
# Need at least some significant pathways for this
if (n_sig_up > 5) {
  cat("Creating enrichment map for upregulated pathways...\n")
  pdf(file.path(output_dir, "05_upregulated_enrichment_map.pdf"), width = 12, height = 10)
  tryCatch(
    {
      up_go_simple <- simplify(up_go_bp, cutoff = 0.7, by = "p.adjust")
      p5 <- emapplot(pairwise_termsim(up_go_simple), showCategory = 30) +
        labs(title = "Enrichment Map - Upregulated Genes")
      print(p5)
      cat("  ✓ Saved: 05_upregulated_enrichment_map.pdf\n")
    },
    error = function(e) {
      cat("  ⚠ Enrichment map could not be generated for upregulated genes\n")
    }
  )
  dev.off()
} else {
  cat("  ⚠ Not enough significant pathways for upregulated enrichment map\n")
}

if (n_sig_down > 5) {
  cat("Creating enrichment map for downregulated pathways...\n")
  pdf(file.path(output_dir, "06_downregulated_enrichment_map.pdf"), width = 12, height = 10)
  tryCatch(
    {
      down_go_simple <- simplify(down_go_bp, cutoff = 0.7, by = "p.adjust")
      p6 <- emapplot(pairwise_termsim(down_go_simple), showCategory = 30) +
        labs(title = "Enrichment Map - Downregulated Genes")
      print(p6)
      cat("  ✓ Saved: 06_downregulated_enrichment_map.pdf\n")
    },
    error = function(e) {
      cat("  ⚠ Enrichment map could not be generated for downregulated genes\n")
    }
  )
  dev.off()
} else {
  cat("  ⚠ Not enough significant pathways for downregulated enrichment map\n")
}

cat("\n✓ All visualizations created\n\n")

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
cat("  Gene lists:\n")
cat("  - upregulated_genes.csv\n")
cat("  - downregulated_genes.csv\n")
cat("\n  GO Enrichment Results:\n")
cat("  - upregulated_GO_BP/MF/CC.csv\n")
cat("  - downregulated_GO_BP/MF/CC.csv\n")
cat("  - *_cancer_relevant.csv (filtered by keywords)\n")
cat("\n  Visualizations:\n")
cat("  - 01_BP_comparison_dotplot.pdf (TOP 25 pathways - all)\n")
cat("  - 01b_BP_significant_dotplot.pdf (significant q<0.05 only)\n")
cat("  - 02_upregulated_barplot.pdf (TOP 30 with significance line)\n")
cat("  - 03_downregulated_barplot.pdf (TOP 30 with significance line)\n")
cat("  - 04_combined_comparison.pdf (side-by-side UP vs DOWN)\n")
cat("  - 05/06_enrichment_map.pdf (network view, if enough significant)\n")
cat("\n  Summary:\n")
cat("  - DIRECTIONAL_GO_SUMMARY.txt\n")
cat(sprintf("\nNote: Upregulated genes had %d significant GO BP terms (q<0.05)\n", n_sig_up))
cat(sprintf("      Downregulated genes had %d significant GO BP terms (q<0.05)\n", n_sig_down))
