#!/usr/bin/env Rscript

# 05_functional_enrichment.R
# Module 5: Functional Enrichment
# Part of the SRF_Eva_integrated_analysis pipeline
#
# FIXED: Absolute paths

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(dplyr)
})

# ===================== Configuration =====================

# Base directories - use absolute paths throughout
BASE_DIR <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top"
SCRIPT_DIR <- file.path(BASE_DIR, "SRF_Eva_integrated_analysis/scripts/analysis_2")

# Input: Annotated Peaks from Module 1
ANNOTATION_FILE <- file.path(SCRIPT_DIR,
                             "results/01_peak_classification/Master_Peak_Annotation.csv")

# Output
OUTPUT_DIR <- file.path(SCRIPT_DIR, "results/05_functional_enrichment")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ===================== Main Analysis =====================

cat("=== Module 5: Functional Enrichment ===\n\n")

# 1. Load Data
cat("Loading annotated peaks...\n")
if (!file.exists(ANNOTATION_FILE)) {
  stop(paste("Annotation file not found:", ANNOTATION_FILE,
             "\nPlease run Module 1 (01_peak_classification.R) first."))
}
anno_df <- read.csv(ANNOTATION_FILE, stringsAsFactors = FALSE)

cat(sprintf("  Loaded %d annotated peaks\n", nrow(anno_df)))

# 2. Perform Enrichment by Category
cat("\nPerforming enrichment analysis...\n")

categories <- unique(anno_df$category)
cat(sprintf("  Categories: %s\n", paste(categories, collapse = ", ")))

enrich_results <- list()

for (cat_name in categories) {
  cat(sprintf("\n  Processing %s...\n", cat_name))

  # Get genes for this category (Entrez IDs)
  genes <- unique(anno_df$geneId[anno_df$category == cat_name])
  genes <- genes[!is.na(genes)]
  genes <- as.character(genes)

  cat(sprintf("    Unique genes: %d\n", length(genes)))

  if (length(genes) < 10) {
    cat("    Skipping (too few genes for enrichment)\n")
    next
  }

  # GO Enrichment (Biological Process)
  tryCatch({
    ego <- enrichGO(gene = genes,
                    OrgDb = org.Hs.eg.db,
                    keyType = "ENTREZID",
                    ont = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05,
                    readable = TRUE)

    if (!is.null(ego) && nrow(ego) > 0) {
      cat(sprintf("    Enriched GO terms: %d\n", nrow(ego)))

      # Save results
      write.csv(as.data.frame(ego),
                file.path(OUTPUT_DIR, paste0("GO_BP_", cat_name, ".csv")),
                row.names = FALSE)

      # Store for comparison
      enrich_results[[cat_name]] <- ego
    } else {
      cat("    No significant GO terms found\n")
    }
  }, error = function(e) {
    cat(sprintf("    Error in enrichment: %s\n", e$message))
  })
}

# 3. Visualization
cat("\n\nGenerating enrichment plots...\n")

# Dotplot for each category
for (cat_name in names(enrich_results)) {
  cat(sprintf("  Creating dotplot for %s...\n", cat_name))

  pdf(file.path(OUTPUT_DIR, paste0("Dotplot_", cat_name, ".pdf")),
      width = 10, height = 8)
  p <- dotplot(enrich_results[[cat_name]],
               showCategory = 20,
               title = paste("GO:BP Enrichment -", cat_name))
  print(p)
  dev.off()
}

# Compare Clusters (if multiple categories have results)
if (length(enrich_results) > 1) {
  cat("\n  Comparing categories...\n")

  # Create a list of gene vectors (Entrez IDs)
  gene_list <- lapply(categories, function(x) {
    genes <- unique(anno_df$geneId[anno_df$category == x])
    genes <- genes[!is.na(genes)]
    as.character(genes)
  })
  names(gene_list) <- categories

  # Remove empty categories
  gene_list <- gene_list[sapply(gene_list, length) >= 10]

  if (length(gene_list) > 1) {
    tryCatch({
      # CompareCluster
      ck <- compareCluster(geneCluster = gene_list,
                           fun = "enrichGO",
                           OrgDb = org.Hs.eg.db,
                           keyType = "ENTREZID",
                           ont = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05)

      if (!is.null(ck) && nrow(ck) > 0) {
        pdf(file.path(OUTPUT_DIR, "Comparison_Dotplot.pdf"),
            width = 12, height = 10)
        p <- dotplot(ck,
                     showCategory = 10,
                     title = "Functional Comparison of Peak Categories")
        print(p)
        dev.off()

        # Save comparison results
        write.csv(as.data.frame(ck),
                  file.path(OUTPUT_DIR, "Comparison_All_Categories.csv"),
                  row.names = FALSE)

        cat("    Comparison plot generated\n")
      }
    }, error = function(e) {
      cat(sprintf("    Error in comparison: %s\n", e$message))
    })
  }
}

# 4. Summary
cat("\n=== Summary ===\n")
for (cat_name in names(enrich_results)) {
  n_terms <- nrow(enrich_results[[cat_name]])
  cat(sprintf("  %s: %d enriched GO:BP terms\n", cat_name, n_terms))
}

cat("\nAnalysis complete. Results saved to", OUTPUT_DIR, "\n")
