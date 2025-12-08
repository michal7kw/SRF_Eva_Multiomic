#!/usr/bin/env Rscript
#
# MSIGDB GSEA ANALYSIS - ORGANIZED BY COLLECTION
# Runs GSEA and saves results in separate folders for each MSigDB collection
# Uses fgsea package with MSigDB gene sets from .grp format
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
cat("  MSIGDB GSEA ANALYSIS - ORGANIZED BY COLLECTION\n")
cat("  Results will be saved in separate folders per MSigDB collection\n")
cat("=============================================================================\n")
cat("Analysis started:", as.character(Sys.time()), "\n\n")

# =============================================================================
# CONFIGURATION
# =============================================================================

# Paths
GENE_SET_DIR <- "../../GENE_SETS"
RNA_SEQ_FILE <- "../../../SRF_Eva_RNA/results/05_deseq2/deseq2_results_TES_vs_GFP.txt"
BASE_OUTPUT_DIR <- "output/12_msigdb_by_collection"

# Analysis parameters
MIN_GENE_SET_SIZE <- 10
MAX_GENE_SET_SIZE <- 500
N_PERMUTATIONS <- 10000
FDR_THRESHOLD <- 0.05
N_PROC <- 8

# Create base output directory
dir.create(BASE_OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# FUNCTION: Load .grp gene set files
# =============================================================================

load_grp_file <- function(filepath) {
  lines <- readLines(filepath)
  lines <- lines[!grepl("^#", lines) & nchar(trimws(lines)) > 0]

  if (length(lines) < 2) {
    return(NULL)
  }

  pathway_name <- lines[1]
  genes <- lines[-1]

  return(list(name = pathway_name, genes = genes))
}

# =============================================================================
# FUNCTION: Determine collection from pathway name
# =============================================================================

get_collection_name <- function(pathway_name) {
  collection <- case_when(
    grepl("^HALLMARK_", pathway_name) ~ "01_Hallmark",
    grepl("^GOBP_", pathway_name) ~ "02_GO_Biological_Process",
    grepl("^GOCC_", pathway_name) ~ "03_GO_Cellular_Component",
    grepl("^GOMF_", pathway_name) ~ "04_GO_Molecular_Function",
    grepl("^KEGG_", pathway_name) ~ "05_KEGG_Pathways",
    grepl("^REACTOME_", pathway_name) ~ "06_Reactome_Pathways",
    grepl("^WP_", pathway_name) ~ "07_WikiPathways",
    grepl("^BIOCARTA_", pathway_name) ~ "08_BioCarta",
    grepl("^PID_", pathway_name) ~ "09_PID",
    # C6 Oncogenic signatures
    grepl("_(UP|DN)\\.V[0-9]_(UP|DN)$", pathway_name) ~ "10_C6_Oncogenic",
    # C3 Cell Atlas (CAR = curated cancer cell atlas)
    grepl("^(CAHOY_|CORDENONSI_|CAR_)", pathway_name) ~ "11_C3_Cancer_Atlas",
    # C4 Computational (MODULE = CM cancer modules, GCM = cancer gene neighborhoods)
    grepl("^(CGN_|CM_|MODULE_|GCM_)", pathway_name) ~ "12_C4_Computational",
    # C2 Curated gene sets (MORF, GNF2, etc.)
    grepl("^(MORF_|GNF2_)", pathway_name) ~ "14_C2_Curated",
    TRUE ~ "13_Other"
  )
  return(collection)
}

# =============================================================================
# PHASE 1: LOAD AND ORGANIZE GENE SETS BY COLLECTION
# =============================================================================

cat("\n=== PHASE 1: Loading and Organizing Gene Sets ===\n")

grp_files <- list.files(GENE_SET_DIR, pattern = "\\.grp$", full.names = TRUE)

if (length(grp_files) == 0) {
  stop("ERROR: No .grp files found in ", GENE_SET_DIR)
}

cat(sprintf("Found %d .grp gene set files\n", length(grp_files)))

# Load all gene sets and organize by collection
all_gene_sets <- list()
collection_map <- list()

for (grp_file in grp_files) {
  gene_set <- load_grp_file(grp_file)

  if (!is.null(gene_set)) {
    # Determine collection
    collection <- get_collection_name(gene_set$name)

    # Store gene set
    all_gene_sets[[gene_set$name]] <- gene_set$genes

    # Track collection membership
    if (!collection %in% names(collection_map)) {
      collection_map[[collection]] <- c()
    }
    collection_map[[collection]] <- c(collection_map[[collection]], gene_set$name)
  }
}

# Filter by gene set size
initial_count <- length(all_gene_sets)
all_gene_sets <- all_gene_sets[
  sapply(all_gene_sets, length) >= MIN_GENE_SET_SIZE &
    sapply(all_gene_sets, length) <= MAX_GENE_SET_SIZE
]

cat(sprintf("\n✓ Loaded %d gene sets\n", initial_count))
cat(sprintf(
  "✓ After size filtering (%d-%d genes): %d gene sets\n",
  MIN_GENE_SET_SIZE, MAX_GENE_SET_SIZE, length(all_gene_sets)
))

# Print collection summary
cat("\nGene sets by collection:\n")
for (coll in sort(names(collection_map))) {
  valid_sets <- intersect(collection_map[[coll]], names(all_gene_sets))
  cat(sprintf("  %s: %d gene sets\n", coll, length(valid_sets)))
}

# =============================================================================
# PHASE 2: PREPARE RANKED GENE LIST FROM RNA-SEQ
# =============================================================================

cat("\n=== PHASE 2: Creating Ranked Gene List ===\n")

rna_all <- read.delim(RNA_SEQ_FILE, stringsAsFactors = FALSE)
cat(sprintf("✓ Loaded %d genes from RNA-seq\n", nrow(rna_all)))

rna_all$ensembl_id <- gsub("\\..*", "", rna_all$gene_id)

rna_all$gene_symbol <- mapIds(org.Hs.eg.db,
  keys = rna_all$ensembl_id,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

rna_ranked <- rna_all %>%
  filter(!is.na(gene_symbol) & !is.na(log2FoldChange)) %>%
  arrange(desc(log2FoldChange))

cat(sprintf("✓ Ranked %d genes with valid symbols and fold changes\n", nrow(rna_ranked)))

gene_ranks <- setNames(rna_ranked$log2FoldChange, rna_ranked$gene_symbol)
gene_ranks <- gene_ranks[!duplicated(names(gene_ranks))]

cat(sprintf("✓ Final ranked list: %d unique genes\n", length(gene_ranks)))
cat(sprintf("  Range: %.2f to %.2f log2FC\n", min(gene_ranks), max(gene_ranks)))

# =============================================================================
# PHASE 3: RUN FGSEA ON ALL GENE SETS
# =============================================================================

cat("\n=== PHASE 3: Running GSEA on All Gene Sets ===\n")
cat(sprintf(
  "Testing %d gene sets with %d permutations...\n",
  length(all_gene_sets), N_PERMUTATIONS
))

set.seed(42)

fgsea_all_results <- fgsea(
  pathways = all_gene_sets,
  stats = gene_ranks,
  minSize = MIN_GENE_SET_SIZE,
  maxSize = MAX_GENE_SET_SIZE,
  nproc = N_PROC,
  nPermSimple = N_PERMUTATIONS
)

# Add collection information
fgsea_all_results$collection <- sapply(fgsea_all_results$pathway, get_collection_name)

cat(sprintf("\n✓ GSEA complete for all gene sets!\n"))
cat(sprintf("  Total gene sets tested: %d\n", nrow(fgsea_all_results)))
cat(sprintf(
  "  Significant (FDR < %.2f): %d\n",
  FDR_THRESHOLD, sum(fgsea_all_results$padj < FDR_THRESHOLD, na.rm = TRUE)
))

# =============================================================================
# PHASE 4: PROCESS AND SAVE RESULTS BY COLLECTION
# =============================================================================

cat("\n=== PHASE 4: Processing Results by Collection ===\n")

# Get unique collections from results
collections_in_results <- unique(fgsea_all_results$collection)

for (collection in sort(collections_in_results)) {
  cat(sprintf("\nProcessing collection: %s\n", collection))

  # Filter results for this collection
  collection_results <- fgsea_all_results %>%
    dplyr::filter(collection == !!collection)

  collection_sig <- collection_results %>%
    dplyr::filter(padj < FDR_THRESHOLD) %>%
    arrange(padj)

  cat(sprintf("  Gene sets in collection: %d\n", nrow(collection_results)))
  cat(sprintf(
    "  Significant: %d (%d up, %d down)\n",
    nrow(collection_sig),
    sum(collection_sig$NES > 0),
    sum(collection_sig$NES < 0)
  ))

  # Create output directory for this collection
  output_dir <- file.path(BASE_OUTPUT_DIR, collection)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # Export CSV files
  collection_results_export <- collection_results %>%
    mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";"))

  collection_sig_export <- collection_sig %>%
    mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";"))

  write.csv(collection_results_export,
    file.path(output_dir, "all_pathways.csv"),
    row.names = FALSE
  )

  write.csv(collection_sig_export,
    file.path(output_dir, "significant_pathways.csv"),
    row.names = FALSE
  )

  # Create visualizations if there are significant results
  if (nrow(collection_sig) > 0) {
    # Plot 1: Top pathways barplot
    top_n <- min(25, nrow(collection_sig))
    top_pathways <- collection_sig %>%
      arrange(desc(abs(NES))) %>%
      head(top_n) %>%
      mutate(
        pathway_short = str_trunc(gsub("_", " ", pathway), width = 60),
        direction = ifelse(NES > 0, "Upregulated", "Downregulated")
      )

    pdf(file.path(output_dir, "plot_01_top_pathways.pdf"),
      width = 14, height = 10
    )
    p1 <- ggplot(top_pathways, aes(x = reorder(pathway_short, NES), y = NES, fill = direction)) +
      geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
      coord_flip() +
      labs(
        title = paste0(collection, " - Top Pathways"),
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

    # Plot 2: Volcano-style plot
    collection_sig_plot <- collection_sig %>%
      mutate(
        neg_log_padj = -log10(padj),
        label = ifelse(rank(padj) <= 10, str_trunc(gsub("_", " ", pathway), 40), "")
      )

    pdf(file.path(output_dir, "plot_02_enrichment_scatter.pdf"),
      width = 14, height = 10
    )
    p2 <- ggplot(collection_sig_plot, aes(x = NES, y = neg_log_padj, size = size)) +
      geom_point(alpha = 0.6, color = "#E31A1C") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
      ggrepel::geom_text_repel(aes(label = label), size = 3, max.overlaps = 15) +
      labs(
        title = paste0(collection, " - Pathway Enrichment"),
        subtitle = "Bubble size = gene set size",
        x = "Normalized Enrichment Score (NES)",
        y = "-log10(Adjusted P-value)",
        size = "Gene Set Size"
      ) +
      theme_bw(base_size = 14) +
      theme(
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5),
        legend.position = "right"
      )
    print(p2)
    dev.off()

    # Plot 3: Enrichment curves for top pathways
    top_for_curves <- min(12, nrow(collection_sig))
    top_pathway_names <- collection_sig %>%
      arrange(padj) %>%
      head(top_for_curves) %>%
      pull(pathway)

    pdf(file.path(output_dir, "plot_03_enrichment_curves.pdf"),
      width = 12, height = 8
    )
    for (pathway_name in top_pathway_names) {
      if (pathway_name %in% names(all_gene_sets)) {
        p <- plotEnrichment(all_gene_sets[[pathway_name]], gene_ranks) +
          labs(title = gsub("_", " ", pathway_name)) +
          theme_bw(base_size = 14) +
          theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
        print(p)
      }
    }
    dev.off()
  }

  # Create summary report
  summary_file <- file.path(output_dir, "SUMMARY.txt")
  cat("", file = summary_file)
  cat(paste0(collection, " - GSEA SUMMARY\n"), file = summary_file, append = TRUE)
  cat(paste0(strrep("=", nchar(collection) + 17), "\n\n"), file = summary_file, append = TRUE)
  cat(paste("Generated:", Sys.time(), "\n\n"), file = summary_file, append = TRUE)

  cat("RESULTS:\n", file = summary_file, append = TRUE)
  cat(sprintf("  Total gene sets tested: %d\n", nrow(collection_results)), file = summary_file, append = TRUE)
  cat(sprintf("  Significant (FDR < %.2f): %d\n", FDR_THRESHOLD, nrow(collection_sig)), file = summary_file, append = TRUE)
  cat(sprintf("  Upregulated (NES > 0): %d\n", sum(collection_sig$NES > 0)), file = summary_file, append = TRUE)
  cat(sprintf("  Downregulated (NES < 0): %d\n\n", sum(collection_sig$NES < 0)), file = summary_file, append = TRUE)

  if (nrow(collection_sig) > 0) {
    cat("TOP 10 UPREGULATED:\n", file = summary_file, append = TRUE)
    top_up <- collection_sig %>%
      filter(NES > 0) %>%
      arrange(desc(NES)) %>%
      head(10)
    for (i in 1:nrow(top_up)) {
      cat(
        sprintf(
          "  %2d. %s (NES=%.2f, FDR=%.2e)\n",
          i, top_up$pathway[i], top_up$NES[i], top_up$padj[i]
        ),
        file = summary_file, append = TRUE
      )
    }

    cat("\nTOP 10 DOWNREGULATED:\n", file = summary_file, append = TRUE)
    top_down <- collection_sig %>%
      filter(NES < 0) %>%
      arrange(NES) %>%
      head(10)
    for (i in 1:nrow(top_down)) {
      cat(
        sprintf(
          "  %2d. %s (NES=%.2f, FDR=%.2e)\n",
          i, top_down$pathway[i], top_down$NES[i], top_down$padj[i]
        ),
        file = summary_file, append = TRUE
      )
    }
  }

  cat(sprintf("\n  ✓ Saved results to: %s\n", output_dir))
}

# =============================================================================
# PHASE 5: CREATE MASTER SUMMARY
# =============================================================================

cat("\n=== PHASE 5: Creating Master Summary ===\n")

master_summary_file <- file.path(BASE_OUTPUT_DIR, "MASTER_SUMMARY.txt")

cat("", file = master_summary_file)
cat("=============================================================================\n", file = master_summary_file, append = TRUE)
cat("  MSIGDB GSEA ANALYSIS - MASTER SUMMARY\n", file = master_summary_file, append = TRUE)
cat("=============================================================================\n\n", file = master_summary_file, append = TRUE)
cat(paste("Generated:", Sys.time(), "\n\n"), file = master_summary_file, append = TRUE)

cat("OVERALL STATISTICS:\n", file = master_summary_file, append = TRUE)
cat(sprintf("  Total gene sets tested: %d\n", nrow(fgsea_all_results)), file = master_summary_file, append = TRUE)
cat(
  sprintf(
    "  Significant (FDR < %.2f): %d\n",
    FDR_THRESHOLD, sum(fgsea_all_results$padj < FDR_THRESHOLD, na.rm = TRUE)
  ),
  file = master_summary_file, append = TRUE
)
cat(sprintf("  Collections analyzed: %d\n\n", length(collections_in_results)),
  file = master_summary_file, append = TRUE
)

cat("RESULTS BY COLLECTION:\n", file = master_summary_file, append = TRUE)
for (coll_name in sort(collections_in_results)) {
  coll_results <- fgsea_all_results %>% dplyr::filter(collection == coll_name)
  coll_sig <- coll_results %>% dplyr::filter(padj < FDR_THRESHOLD)

  cat(sprintf("  %s:\n", coll_name), file = master_summary_file, append = TRUE)
  cat(
    sprintf(
      "    Tested: %d | Significant: %d | Up: %d | Down: %d\n",
      nrow(coll_results), nrow(coll_sig),
      sum(coll_sig$NES > 0), sum(coll_sig$NES < 0)
    ),
    file = master_summary_file, append = TRUE
  )
  cat(sprintf("    Output: %s\n\n", file.path(BASE_OUTPUT_DIR, coll_name)),
    file = master_summary_file, append = TRUE
  )
}

cat("\n✓ Master summary saved\n")

# Export combined results
fgsea_all_results_export <- fgsea_all_results %>%
  mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";"))

write.csv(fgsea_all_results_export,
  file.path(BASE_OUTPUT_DIR, "ALL_COLLECTIONS_combined.csv"),
  row.names = FALSE
)

fgsea_all_sig_export <- fgsea_all_results %>%
  filter(padj < FDR_THRESHOLD) %>%
  mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";"))

write.csv(fgsea_all_sig_export,
  file.path(BASE_OUTPUT_DIR, "ALL_COLLECTIONS_significant.csv"),
  row.names = FALSE
)

# =============================================================================
# COMPLETION
# =============================================================================

cat("\n=============================================================================\n")
cat("  ANALYSIS COMPLETE\n")
cat("=============================================================================\n")
cat("Completed:", as.character(Sys.time()), "\n")
cat(sprintf("Base output directory: %s\n\n", BASE_OUTPUT_DIR))

cat("Results organized by collection:\n")
for (collection in sort(collections_in_results)) {
  cat(sprintf("  - %s/\n", file.path(BASE_OUTPUT_DIR, collection)))
}

cat("\nMaster files:\n")
cat(sprintf("  - %s\n", file.path(BASE_OUTPUT_DIR, "MASTER_SUMMARY.txt")))
cat(sprintf("  - %s\n", file.path(BASE_OUTPUT_DIR, "ALL_COLLECTIONS_combined.csv")))
cat(sprintf("  - %s\n\n", file.path(BASE_OUTPUT_DIR, "ALL_COLLECTIONS_significant.csv")))

cat("Session Info:\n")
sessionInfo()
