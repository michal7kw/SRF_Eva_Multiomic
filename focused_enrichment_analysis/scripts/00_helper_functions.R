#!/usr/bin/env Rscript
# Helper functions for enrichment analysis

################################################################################
# LIBRARY LOADING
################################################################################
load_libraries <- function() {
  suppressPackageStartupMessages({
    library(readr)
    library(ggplot2)
    library(GenomicRanges)
    library(ChIPseeker)
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    library(enrichplot)
    library(ggpubr)
    library(ComplexHeatmap)
    library(circlize)
    library(VennDiagram)
    library(gridExtra)
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(msigdbr)
    library(dplyr)  # Load dplyr LAST to avoid namespace conflicts
  })
}

################################################################################
# GO ENRICHMENT FUNCTION
################################################################################
perform_GO_enrichment <- function(gene_list, background_genes, ont = "BP",
                                   pvalueCutoff = 0.05, qvalueCutoff = 0.05,
                                   minGSSize = 10, maxGSSize = 500) {

  # Convert gene symbols to Entrez IDs
  gene_entrez <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID",
                      OrgDb = org.Hs.eg.db, drop = TRUE)

  background_entrez <- bitr(background_genes, fromType = "SYMBOL", toType = "ENTREZID",
                            OrgDb = org.Hs.eg.db, drop = TRUE)

  # Run enrichGO
  ego <- enrichGO(gene = gene_entrez$ENTREZID,
                  universe = background_entrez$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = ont,
                  pAdjustMethod = "BH",
                  pvalueCutoff = pvalueCutoff,
                  qvalueCutoff = qvalueCutoff,
                  minGSSize = minGSSize,
                  maxGSSize = maxGSSize,
                  readable = TRUE)

  return(ego)
}

################################################################################
# GSEA FUNCTION
################################################################################
perform_GSEA_msigdb <- function(ranked_genes, gene_sets, pvalueCutoff = 0.05) {

  # Convert to Entrez IDs
  gene_mapping <- bitr(names(ranked_genes), fromType = "SYMBOL", toType = "ENTREZID",
                       OrgDb = org.Hs.eg.db, drop = TRUE)

  # Create ranked list with Entrez IDs
  ranked_entrez <- ranked_genes[gene_mapping$SYMBOL]
  names(ranked_entrez) <- as.character(gene_mapping$ENTREZID)
  ranked_entrez <- sort(ranked_entrez, decreasing = TRUE)

  # Prepare gene sets
  gene_sets_df <- as.data.frame(gene_sets)

  # Check which column name is available
  if ("ncbi_gene" %in% colnames(gene_sets_df)) {
    msigdb_t2g <- gene_sets_df %>% dplyr::select(gs_name, gene = ncbi_gene)
  } else if ("entrez_gene" %in% colnames(gene_sets_df)) {
    msigdb_t2g <- gene_sets_df %>% dplyr::select(gs_name, gene = entrez_gene)
  } else if ("human_entrez_gene" %in% colnames(gene_sets_df)) {
    msigdb_t2g <- gene_sets_df %>% dplyr::select(gs_name, gene = human_entrez_gene)
  } else if ("entrez_id" %in% colnames(gene_sets_df)) {
    msigdb_t2g <- gene_sets_df %>% dplyr::select(gs_name, gene = entrez_id)
  } else {
    cat("Available columns in gene_sets:", paste(colnames(gene_sets_df), collapse=", "), "\n")
    stop("Cannot find Entrez gene column in msigdbr results")
  }

  # Ensure gene column is character
  msigdb_t2g$gene <- as.character(msigdb_t2g$gene)

  # Debug: check gene ID formats match
  cat("DEBUG - Sample ranked gene IDs:", paste(head(names(ranked_entrez), 10), collapse=", "), "\n")
  cat("DEBUG - Sample TERM2GENE IDs:", paste(head(unique(msigdb_t2g$gene), 10), collapse=", "), "\n")
  cat("DEBUG - Overlap check:", sum(names(ranked_entrez) %in% msigdb_t2g$gene), "genes overlap out of", length(ranked_entrez), "\n")

  # Run GSEA
  gsea_result <- GSEA(geneList = ranked_entrez,
                      TERM2GENE = msigdb_t2g,
                      pvalueCutoff = pvalueCutoff,
                      pAdjustMethod = "BH",
                      minGSSize = 15,
                      maxGSSize = 500)

  return(gsea_result)
}

################################################################################
# SAVE ENRICHMENT RESULTS
################################################################################
save_enrichment_results <- function(enrich_obj, prefix, output_dir) {

  if (is.null(enrich_obj) || nrow(enrich_obj) == 0) {
    cat("  No significant enrichment found for", prefix, "\n")
    return(NULL)
  }

  # Save full results
  results_df <- as.data.frame(enrich_obj)
  write_csv(results_df, file.path(output_dir, "results", paste0(prefix, "_GO_enrichment.csv")))

  # Save top 20 results
  if (nrow(results_df) > 0) {
    top20 <- head(results_df, 20)
    write_csv(top20, file.path(output_dir, "results", paste0(prefix, "_top20_terms.csv")))
  }

  # Generate plots
  if (nrow(enrich_obj) >= 10) {
    cat("  Generating plots for", prefix, "with", nrow(enrich_obj), "terms...\n")

    # Determine number of categories to show
    n_categories <- min(20, nrow(enrich_obj))

    # Calculate dynamic height based on number of terms
    # Base height + 0.5 inches per term (ensures adequate spacing)
    plot_height <- max(14, 6 + (n_categories * 0.5))

    cat("  Plot dimensions: width=15, height=", plot_height, "\n")

    # Dot plot
    p1 <- dotplot(enrich_obj, showCategory = n_categories) +
      ggtitle(paste(prefix, "- Top", n_categories, "GO Terms")) +
      theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
            axis.text.y = element_text(size = 8, lineheight = 1.2, margin = margin(r = 5)),
            axis.text.x = element_text(size = 10),
            axis.title = element_text(size = 11),
            strip.text.y = element_text(size = 10),
            plot.margin = margin(5, 20, 5, 5)) +
      scale_y_discrete(labels = function(x) {
        # Wrap long labels at 55 characters
        sapply(x, function(label) {
          if (nchar(label) > 55) {
            paste(strwrap(label, width = 55), collapse = "\n")
          } else {
            label
          }
        })
      })

    # Save as PDF
    pdf_file <- file.path(output_dir, "plots", paste0(prefix, "_dotplot.pdf"))
    cat("  Saving PDF:", pdf_file, "\n")
    ggsave(pdf_file, p1, width = 15, height = plot_height, device = "pdf")

    # Save as PNG
    png_file <- file.path(output_dir, "plots", paste0(prefix, "_dotplot.png"))
    cat("  Saving PNG:", png_file, "\n")
    ggsave(png_file, p1, width = 15, height = plot_height, device = "png", dpi = 300)

    # Bar plot
    p2 <- barplot(enrich_obj, showCategory = n_categories) +
      ggtitle(paste(prefix, "- Top", n_categories, "GO Terms")) +
      theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
            axis.text.y = element_text(size = 8, lineheight = 1.2, margin = margin(r = 5)),
            axis.text.x = element_text(size = 10),
            axis.title = element_text(size = 11),
            strip.text.y = element_text(size = 10),
            plot.margin = margin(5, 20, 5, 5)) +
      scale_y_discrete(labels = function(x) {
        # Wrap long labels at 55 characters
        sapply(x, function(label) {
          if (nchar(label) > 55) {
            paste(strwrap(label, width = 55), collapse = "\n")
          } else {
            label
          }
        })
      })

    # Save as PDF
    pdf_file <- file.path(output_dir, "plots", paste0(prefix, "_barplot.pdf"))
    cat("  Saving PDF:", pdf_file, "\n")
    ggsave(pdf_file, p2, width = 15, height = plot_height, device = "pdf")

    # Save as PNG
    png_file <- file.path(output_dir, "plots", paste0(prefix, "_barplot.png"))
    cat("  Saving PNG:", png_file, "\n")
    ggsave(png_file, p2, width = 15, height = plot_height, device = "png", dpi = 300)

    cat("  Plots saved successfully!\n")
  }

  return(results_df)
}

################################################################################
# EXTRACT MIGRATION TERMS
################################################################################
extract_migration_terms <- function(enrich_results) {
  if (is.null(enrich_results) || nrow(enrich_results) == 0) return(NULL)

  migration_keywords <- c("migration", "motility", "adhesion", "invasion",
                          "metastasis", "locomotion", "chemotaxis")

  pattern <- paste(migration_keywords, collapse = "|")
  migration_terms <- enrich_results %>%
    filter(grepl(pattern, Description, ignore.case = TRUE)) %>%
    arrange(qvalue)

  return(migration_terms)
}

################################################################################
# EXTRACT GENES FROM PEAKS
################################################################################
extract_genes_from_peaks <- function(peaks_df) {
  # Check if SYMBOL column exists
  if ("SYMBOL" %in% colnames(peaks_df)) {
    genes <- unique(peaks_df$SYMBOL[!is.na(peaks_df$SYMBOL)])
  } else if ("geneId" %in% colnames(peaks_df)) {
    # Convert Entrez IDs to gene symbols
    entrez_ids <- unique(peaks_df$geneId[!is.na(peaks_df$geneId)])
    # Remove non-numeric IDs
    entrez_ids <- entrez_ids[grepl("^[0-9]+$", entrez_ids)]

    if (length(entrez_ids) > 0) {
      gene_mapping <- bitr(entrez_ids, fromType = "ENTREZID", toType = "SYMBOL",
                          OrgDb = org.Hs.eg.db, drop = TRUE)
      genes <- unique(gene_mapping$SYMBOL)
    } else {
      genes <- character(0)
    }
  } else {
    cat("WARNING: No gene identifier column found in peaks\n")
    genes <- character(0)
  }

  return(genes)
}
