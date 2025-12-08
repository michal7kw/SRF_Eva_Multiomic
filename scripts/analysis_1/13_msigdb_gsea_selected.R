#!/usr/bin/env Rscript
#
# COMPREHENSIVE MSIGDB GSEA ANALYSIS
# Uses fgsea package with MSigDB gene sets from .grp format
# Analyzes ALL genes ranked by fold change
#
# Author: Michal Kubacki
# Date: 2025-10-13

suppressPackageStartupMessages({
  library(fgsea)
  library(dplyr)
  library(ggplot2)
  library(org.Hs.eg.db)
  library(readr)
  library(stringr)
  library(tidyr)
  library(pheatmap)
  library(RColorBrewer)
})

setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_1")

cat("=============================================================================\n")
cat("  COMPREHENSIVE MSIGDB GSEA ANALYSIS\n")
cat("  Using fgsea with Molecular Signatures Database gene sets\n")
cat("=============================================================================\n")
cat("Analysis started:", as.character(Sys.time()), "\n\n")

# =============================================================================
# CONFIGURATION
# =============================================================================

# Paths
GENE_SET_DIR <- "../../GENE_SETS_selected"
RNA_SEQ_FILE <- "../../../SRF_Eva_RNA/results/05_deseq2/deseq2_results_TES_vs_GFP.txt"
OUTPUT_DIR <- "output/13_msigdb_gsea_selected"

# Analysis parameters
MIN_GENE_SET_SIZE <- 10
MAX_GENE_SET_SIZE <- 500
N_PERMUTATIONS <- 10000
FDR_THRESHOLD <- 0.05
N_PROC <- 8 # Parallel processing cores

# Create output directory
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# FUNCTION: Load .grp gene set files
# =============================================================================

load_grp_file <- function(filepath) {
  # Read .grp file (first line is pathway name, rest are genes)
  lines <- readLines(filepath)

  # Skip comment lines and empty lines
  lines <- lines[!grepl("^#", lines) & nchar(trimws(lines)) > 0]

  if (length(lines) < 2) {
    return(NULL)
  }

  pathway_name <- lines[1]
  genes <- lines[-1] # All lines except first

  return(list(name = pathway_name, genes = genes))
}


# =============================================================================
# PHASE 1: LOAD AND PREPARE GENE SETS
# =============================================================================

cat("\n=== PHASE 1: Loading MSigDB Gene Sets ===\n")

# Find all .grp files
grp_files <- list.files(GENE_SET_DIR, pattern = "\\.grp$", full.names = TRUE)

if (length(grp_files) == 0) {
  stop("ERROR: No .grp files found in ", GENE_SET_DIR)
}

cat(sprintf("Found %d .grp gene set files\n", length(grp_files)))

# Load all gene sets
msigdb_gene_sets <- list()

for (grp_file in grp_files) {
  gene_set <- load_grp_file(grp_file)

  if (!is.null(gene_set)) {
    cat(sprintf(
      "  Loaded: %s (%d genes)\n",
      gene_set$name, length(gene_set$genes)
    ))
    msigdb_gene_sets[[gene_set$name]] <- gene_set$genes
  }
}

# Filter by gene set size
initial_count <- length(msigdb_gene_sets)
msigdb_gene_sets <- msigdb_gene_sets[
  sapply(msigdb_gene_sets, length) >= MIN_GENE_SET_SIZE &
    sapply(msigdb_gene_sets, length) <= MAX_GENE_SET_SIZE
]

cat(sprintf("\n✓ Loaded %d gene sets\n", initial_count))
cat(sprintf(
  "✓ After size filtering (%d-%d genes): %d gene sets\n",
  MIN_GENE_SET_SIZE, MAX_GENE_SET_SIZE, length(msigdb_gene_sets)
))

if (length(msigdb_gene_sets) == 0) {
  stop("ERROR: No valid gene sets remain after filtering")
}

# =============================================================================
# PHASE 2: PREPARE RANKED GENE LIST FROM RNA-SEQ
# =============================================================================

cat("\n=== PHASE 2: Creating Ranked Gene List ===\n")

# Load RNA-seq results (ALL genes)
rna_all <- read.delim(RNA_SEQ_FILE, stringsAsFactors = FALSE)
cat(sprintf("✓ Loaded %d genes from RNA-seq\n", nrow(rna_all)))

# Clean and prepare for ranking
rna_all$ensembl_id <- gsub("\\..*", "", rna_all$gene_id)

# Get gene symbols
rna_all$gene_symbol <- mapIds(org.Hs.eg.db,
  keys = rna_all$ensembl_id,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

# Remove genes without symbols or log2FC
rna_ranked <- rna_all %>%
  dplyr::filter(!is.na(gene_symbol) & !is.na(log2FoldChange)) %>%
  arrange(desc(log2FoldChange))

cat(sprintf("✓ Ranked %d genes with valid symbols and fold changes\n", nrow(rna_ranked)))

# Create ranked list (gene symbol → log2FC)
gene_ranks <- setNames(rna_ranked$log2FoldChange, rna_ranked$gene_symbol)

# Remove duplicates (keep first = highest ranking)
gene_ranks <- gene_ranks[!duplicated(names(gene_ranks))]

cat(sprintf("✓ Final ranked list: %d unique genes\n", length(gene_ranks)))
cat(sprintf("  Range: %.2f to %.2f log2FC\n", min(gene_ranks), max(gene_ranks)))
cat(sprintf("  Most upregulated: %s (%.2f)\n", names(gene_ranks)[1], gene_ranks[1]))
cat(sprintf(
  "  Most downregulated: %s (%.2f)\n",
  names(gene_ranks)[length(gene_ranks)],
  gene_ranks[length(gene_ranks)]
))

# =============================================================================
# PHASE 3: RUN FGSEA
# =============================================================================

cat("\n=== PHASE 3: Running Gene Set Enrichment Analysis ===\n")
cat(sprintf("Parameters:\n"))
cat(sprintf("  Gene sets: %d\n", length(msigdb_gene_sets)))
cat(sprintf("  Ranked genes: %d\n", length(gene_ranks)))
cat(sprintf("  Permutations: %d\n", N_PERMUTATIONS))
cat(sprintf("  Parallel cores: %d\n", N_PROC))
cat("\nThis may take several minutes...\n")

set.seed(42) # For reproducibility

# Run fgsea
fgsea_results <- fgsea(
  pathways = msigdb_gene_sets,
  stats = gene_ranks,
  minSize = MIN_GENE_SET_SIZE,
  maxSize = MAX_GENE_SET_SIZE,
  nproc = N_PROC,
  nPermSimple = N_PERMUTATIONS
)

# Add collection information from pathway names
fgsea_results <- fgsea_results %>%
  mutate(
    collection = case_when(
      grepl("^HALLMARK_", pathway) ~ "Hallmark",
      grepl("^GOBP_", pathway) ~ "GO Biological Process",
      grepl("^GOCC_", pathway) ~ "GO Cellular Component",
      grepl("^GOMF_", pathway) ~ "GO Molecular Function",
      grepl("^KEGG_", pathway) ~ "KEGG",
      grepl("^REACTOME_", pathway) ~ "Reactome",
      grepl("^WP_", pathway) ~ "WikiPathways",
      grepl("^BIOCARTA_", pathway) ~ "BioCarta",
      grepl("^PID_", pathway) ~ "PID",
      # C6 Oncogenic signatures (pattern: NAME_UP/DN.V#_UP/DN)
      grepl("_(UP|DN)\\.V[0-9]_(UP|DN)$", pathway) ~ "C6 Oncogenic",
      TRUE ~ "Other"
    )
  )

# Filter significant results
fgsea_sig <- fgsea_results %>%
  dplyr::filter(padj < FDR_THRESHOLD) %>%
  arrange(padj)

cat(sprintf("\n✓ GSEA complete!\n"))
cat(sprintf("  Total gene sets tested: %d\n", nrow(fgsea_results)))
cat(sprintf("  Significant gene sets (FDR < %.2f): %d\n", FDR_THRESHOLD, nrow(fgsea_sig)))
cat(sprintf("  Positively enriched (NES > 0): %d\n", sum(fgsea_sig$NES > 0)))
cat(sprintf("  Negatively enriched (NES < 0): %d\n", sum(fgsea_sig$NES < 0)))

# Print breakdown by collection
cat("\nBreakdown by collection:\n")
collection_summary <- fgsea_sig %>%
  group_by(collection) %>%
  summarise(
    total = n(),
    up = sum(NES > 0),
    down = sum(NES < 0),
    .groups = "drop"
  ) %>%
  arrange(desc(total))

for (i in 1:nrow(collection_summary)) {
  cat(sprintf(
    "  %s: %d (%d up, %d down)\n",
    collection_summary$collection[i],
    collection_summary$total[i],
    collection_summary$up[i],
    collection_summary$down[i]
  ))
}

# =============================================================================
# PHASE 4: EXPORT RESULTS
# =============================================================================

cat("\n=== PHASE 4: Exporting Results ===\n")

# Convert list columns to character strings for CSV export
fgsea_results_export <- fgsea_results %>%
  mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";"))

fgsea_sig_export <- fgsea_sig %>%
  mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";"))

# Export all results
write.csv(fgsea_results_export,
  file.path(OUTPUT_DIR, "01_msigdb_all_pathways.csv"),
  row.names = FALSE
)

write.csv(fgsea_sig_export,
  file.path(OUTPUT_DIR, "02_msigdb_significant_pathways.csv"),
  row.names = FALSE
)

# Export by collection
for (coll in unique(fgsea_sig$collection)) {
  coll_data <- fgsea_sig %>%
    dplyr::filter(collection == coll) %>%
    mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";"))

  if (nrow(coll_data) > 0) {
    filename <- paste0("03_", gsub(" ", "_", tolower(coll)), "_pathways.csv")
    write.csv(coll_data, file.path(OUTPUT_DIR, filename), row.names = FALSE)
  }
}

cat(sprintf("✓ Exported %d result files\n", 2 + length(unique(fgsea_sig$collection))))

# =============================================================================
# PHASE 5: VISUALIZATIONS
# =============================================================================

cat("\n=== PHASE 5: Creating Visualizations ===\n")

# Plot 1: Top pathways barplot
if (nrow(fgsea_sig) > 0) {
  cat("Creating top pathways barplot...\n")

  top_n <- 25
  top_pathways <- fgsea_sig %>%
    arrange(desc(abs(NES))) %>%
    head(top_n) %>%
    mutate(
      pathway_short = str_trunc(gsub("_", " ", pathway), width = 60),
      direction = ifelse(NES > 0, "Upregulated", "Downregulated")
    )

  pdf(file.path(OUTPUT_DIR, "plot_01_top_pathways_barplot.pdf"),
    width = 14, height = 10
  )
  p1 <- ggplot(top_pathways, aes(x = reorder(pathway_short, NES), y = NES, fill = direction)) +
    geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
    coord_flip() +
    labs(
      title = sprintf("Top %d MSigDB Pathways (by |NES|)", top_n),
      subtitle = paste0("TES vs GFP | FDR < ", FDR_THRESHOLD),
      x = "Pathway",
      y = "Normalized Enrichment Score (NES)"
    ) +
    scale_fill_manual(values = c("Upregulated" = "#E31A1C", "Downregulated" = "#1F78B4")) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      axis.text.y = element_text(size = 9),
      legend.position = "bottom"
    )
  print(p1)
  dev.off()
}

# Plot 2: Collection overview
if (nrow(collection_summary) > 0) {
  cat("Creating collection overview plot...\n")

  collection_long <- collection_summary %>%
    pivot_longer(cols = c(up, down), names_to = "direction", values_to = "count") %>%
    mutate(direction = ifelse(direction == "up", "Upregulated", "Downregulated"))

  pdf(file.path(OUTPUT_DIR, "plot_02_collection_overview.pdf"),
    width = 12, height = 8
  )
  p2 <- ggplot(collection_long, aes(x = reorder(collection, total), y = count, fill = direction)) +
    geom_bar(stat = "identity", position = "dodge", color = "black", linewidth = 0.3) +
    coord_flip() +
    labs(
      title = "MSigDB Collection Overview",
      subtitle = paste0("Significant pathways (FDR < ", FDR_THRESHOLD, ")"),
      x = "Collection",
      y = "Number of Pathways"
    ) +
    scale_fill_manual(values = c("Upregulated" = "#E31A1C", "Downregulated" = "#1F78B4")) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      legend.position = "bottom"
    )
  print(p2)
  dev.off()
}

# Plot 3: NES vs -log10(padj) scatter plot
if (nrow(fgsea_sig) > 0) {
  cat("Creating volcano-style scatter plot...\n")

  fgsea_sig_plot <- fgsea_sig %>%
    mutate(
      neg_log_padj = -log10(padj),
      label = ifelse(rank(padj) <= 10, gsub("_", " ", pathway), "")
    )

  pdf(file.path(OUTPUT_DIR, "plot_03_enrichment_scatter.pdf"),
    width = 14, height = 10
  )
  p3 <- ggplot(fgsea_sig_plot, aes(x = NES, y = neg_log_padj, color = collection, size = size)) +
    geom_point(alpha = 0.6) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
    ggrepel::geom_text_repel(aes(label = label), size = 3, max.overlaps = 15) +
    labs(
      title = "MSigDB Pathway Enrichment",
      subtitle = "Bubble size = gene set size",
      x = "Normalized Enrichment Score (NES)",
      y = "-log10(Adjusted P-value)",
      color = "Collection",
      size = "Gene Set Size"
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      legend.position = "right"
    ) +
    scale_color_brewer(palette = "Set1")
  print(p3)
  dev.off()
}

# Plot 4: Detailed enrichment plots for top pathways
if (nrow(fgsea_sig) > 0) {
  cat("Creating detailed enrichment plots...\n")

  top_for_enrichment <- min(12, nrow(fgsea_sig))
  top_pathway_names <- fgsea_sig %>%
    arrange(padj) %>%
    head(top_for_enrichment) %>%
    pull(pathway)

  pdf(file.path(OUTPUT_DIR, "plot_04_detailed_enrichment_curves.pdf"),
    width = 12, height = 8
  )
  for (pathway_name in top_pathway_names) {
    p <- plotEnrichment(msigdb_gene_sets[[pathway_name]], gene_ranks) +
      labs(title = gsub("_", " ", pathway_name)) +
      theme_bw(base_size = 14) +
      theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
    print(p)
  }
  dev.off()
}

# Plot 5: Heatmap of top pathways (if many significant)
if (nrow(fgsea_sig) >= 20) {
  cat("Creating pathway heatmap...\n")

  # Get top pathways
  top_pathways_heat <- fgsea_sig %>%
    arrange(padj) %>%
    head(50) %>%
    mutate(pathway_short = str_trunc(gsub("_", " ", pathway), width = 60))

  # Create matrix for heatmap
  heat_matrix <- matrix(top_pathways_heat$NES, ncol = 1)
  rownames(heat_matrix) <- top_pathways_heat$pathway_short
  colnames(heat_matrix) <- "TES vs GFP"

  # Add significance annotation
  annotation_row <- data.frame(
    Collection = top_pathways_heat$collection,
    row.names = top_pathways_heat$pathway_short
  )

  pdf(file.path(OUTPUT_DIR, "plot_05_pathway_heatmap.pdf"),
    width = 10, height = 14
  )
  pheatmap(heat_matrix,
    cluster_cols = FALSE,
    cluster_rows = TRUE,
    annotation_row = annotation_row,
    color = colorRampPalette(c("#1F78B4", "white", "#E31A1C"))(100),
    breaks = seq(-max(abs(top_pathways_heat$NES)),
      max(abs(top_pathways_heat$NES)),
      length.out = 101
    ),
    main = "Top 50 MSigDB Pathways",
    fontsize = 8,
    fontsize_row = 7
  )
  dev.off()
}

cat("✓ All visualizations created\n")

# =============================================================================
# PHASE 6: SUMMARY REPORT
# =============================================================================

cat("\n=== PHASE 6: Generating Summary Report ===\n")

summary_file <- file.path(OUTPUT_DIR, "MSIGDB_GSEA_SUMMARY.txt")

cat("", file = summary_file) # Clear file
cat("=============================================================================\n", file = summary_file, append = TRUE)
cat("  COMPREHENSIVE MSIGDB GSEA ANALYSIS - SUMMARY REPORT\n", file = summary_file, append = TRUE)
cat("=============================================================================\n\n", file = summary_file, append = TRUE)
cat(paste("Generated:", Sys.time(), "\n\n"), file = summary_file, append = TRUE)

cat("INPUT DATA:\n", file = summary_file, append = TRUE)
cat(sprintf("  RNA-seq file: %s\n", RNA_SEQ_FILE), file = summary_file, append = TRUE)
cat(sprintf("  Gene sets directory: %s\n", GENE_SET_DIR), file = summary_file, append = TRUE)
cat(sprintf("  Number of .grp files: %d\n\n", length(grp_files)), file = summary_file, append = TRUE)

cat("PARAMETERS:\n", file = summary_file, append = TRUE)
cat(sprintf("  Gene set size range: %d - %d genes\n", MIN_GENE_SET_SIZE, MAX_GENE_SET_SIZE), file = summary_file, append = TRUE)
cat(sprintf("  Permutations: %d\n", N_PERMUTATIONS), file = summary_file, append = TRUE)
cat(sprintf("  FDR threshold: %.2f\n", FDR_THRESHOLD), file = summary_file, append = TRUE)
cat(sprintf("  Ranking metric: log2FoldChange (TES vs GFP)\n\n"), file = summary_file, append = TRUE)

cat("GENE SET STATISTICS:\n", file = summary_file, append = TRUE)
cat(sprintf("  Total gene sets tested: %d\n", nrow(fgsea_results)), file = summary_file, append = TRUE)
cat(sprintf("  Genes in ranked list: %d\n\n", length(gene_ranks)), file = summary_file, append = TRUE)

cat("ENRICHMENT RESULTS:\n", file = summary_file, append = TRUE)
cat(sprintf("  Significant pathways (FDR < %.2f): %d\n", FDR_THRESHOLD, nrow(fgsea_sig)), file = summary_file, append = TRUE)
cat(sprintf("  Positively enriched (NES > 0): %d\n", sum(fgsea_sig$NES > 0)), file = summary_file, append = TRUE)
cat(sprintf("  Negatively enriched (NES < 0): %d\n\n", sum(fgsea_sig$NES < 0)), file = summary_file, append = TRUE)

cat("COLLECTION BREAKDOWN:\n", file = summary_file, append = TRUE)
for (i in 1:nrow(collection_summary)) {
  cat(
    sprintf(
      "  %s: %d pathways (%d up, %d down)\n",
      collection_summary$collection[i],
      collection_summary$total[i],
      collection_summary$up[i],
      collection_summary$down[i]
    ),
    file = summary_file, append = TRUE
  )
}
cat("\n", file = summary_file, append = TRUE)

# Top 20 upregulated
cat("TOP 20 UPREGULATED PATHWAYS:\n", file = summary_file, append = TRUE)
if (nrow(fgsea_sig) > 0) {
  top_up <- fgsea_sig %>%
    dplyr::filter(NES > 0) %>%
    arrange(desc(NES)) %>%
    head(20)

  for (i in 1:nrow(top_up)) {
    cat(
      sprintf(
        "  %2d. [%s] %s\n      NES=%.2f, FDR=%.2e, size=%d\n",
        i, top_up$collection[i], top_up$pathway[i],
        top_up$NES[i], top_up$padj[i], top_up$size[i]
      ),
      file = summary_file, append = TRUE
    )
  }
}
cat("\n", file = summary_file, append = TRUE)

# Top 20 downregulated
cat("TOP 20 DOWNREGULATED PATHWAYS:\n", file = summary_file, append = TRUE)
if (nrow(fgsea_sig) > 0) {
  top_down <- fgsea_sig %>%
    dplyr::filter(NES < 0) %>%
    arrange(NES) %>%
    head(20)

  for (i in 1:nrow(top_down)) {
    cat(
      sprintf(
        "  %2d. [%s] %s\n      NES=%.2f, FDR=%.2e, size=%d\n",
        i, top_down$collection[i], top_down$pathway[i],
        top_down$NES[i], top_down$padj[i], top_down$size[i]
      ),
      file = summary_file, append = TRUE
    )
  }
}
cat("\n", file = summary_file, append = TRUE)

cat("OUTPUT FILES:\n", file = summary_file, append = TRUE)
cat("  CSV files:\n", file = summary_file, append = TRUE)
cat("    - 01_msigdb_all_pathways.csv (all tested pathways)\n", file = summary_file, append = TRUE)
cat("    - 02_msigdb_significant_pathways.csv (FDR < 0.05)\n", file = summary_file, append = TRUE)
cat("    - 03_*_pathways.csv (by collection)\n", file = summary_file, append = TRUE)
cat("  Plots:\n", file = summary_file, append = TRUE)
cat("    - plot_01_top_pathways_barplot.pdf\n", file = summary_file, append = TRUE)
cat("    - plot_02_collection_overview.pdf\n", file = summary_file, append = TRUE)
cat("    - plot_03_enrichment_scatter.pdf\n", file = summary_file, append = TRUE)
cat("    - plot_04_detailed_enrichment_curves.pdf\n", file = summary_file, append = TRUE)
if (nrow(fgsea_sig) >= 20) {
  cat("    - plot_05_pathway_heatmap.pdf\n", file = summary_file, append = TRUE)
}
cat("\n", file = summary_file, append = TRUE)

cat("=============================================================================\n", file = summary_file, append = TRUE)

cat("✓ Summary report saved\n")

# =============================================================================
# COMPLETION
# =============================================================================

cat("\n=============================================================================\n")
cat("  COMPREHENSIVE MSIGDB GSEA ANALYSIS COMPLETE\n")
cat("=============================================================================\n")
cat("Completed:", as.character(Sys.time()), "\n")
cat(sprintf("Output directory: %s\n\n", OUTPUT_DIR))

cat("Key results:\n")
cat(sprintf("  • Tested %d MSigDB gene sets\n", nrow(fgsea_results)))
cat(sprintf("  • Found %d significant pathways (FDR < %.2f)\n", nrow(fgsea_sig), FDR_THRESHOLD))
cat(sprintf(
  "  • %d upregulated, %d downregulated\n",
  sum(fgsea_sig$NES > 0), sum(fgsea_sig$NES < 0)
))
cat("\nThis analysis uses true rank-based GSEA with the fgsea package.\n")
cat("All genes contribute to enrichment score calculation.\n\n")

# Print session info for reproducibility
cat("Session Info:\n")
sessionInfo()
