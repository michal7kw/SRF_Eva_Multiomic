#!/usr/bin/env Rscript

cat("=================================================================\n")
cat("Starting Comprehensive Enrichment Analysis\n")
cat("=================================================================\n\n")

# Load required libraries
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
  library(dplyr) # Load dplyr LAST to avoid namespace conflicts
})

# Set working directory
setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top")
base_dir <- "SRF_Eva_integrated_analysis/focused_enrichment_analysis"

# Output directories
dirs <- list(
  approach1 = file.path(base_dir, "approach1_direct_targets"),
  approach2 = file.path(base_dir, "approach2_downregulated"),
  approach3 = file.path(base_dir, "approach3_promoter_peaks"),
  approach4 = file.path(base_dir, "approach4_high_confidence"),
  approach5 = file.path(base_dir, "approach5_diffbind"),
  # approach6 = file.path(base_dir, "approach6_migration_focused"),
  # tier1 = file.path(base_dir, "tier1_progression"),
  # tier2 = file.path(base_dir, "tier2_validation")
)

# Create all directories
for (dir_path in dirs) {
  dir.create(file.path(dir_path, "results"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dir_path, "plots"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dir_path, "gene_lists"), recursive = TRUE, showWarnings = FALSE)
}

################################################################################
# LOAD DATA
################################################################################
cat("\n=== Loading Data ===\n")

# RNA-seq differential expression results
cat("Loading RNA-seq DEG data...\n")
deseq_results <- read.table(
  "SRF_Eva_RNA/results/05_deseq2/deseq2_results_TES_vs_GFP.txt",
  header = TRUE,
  sep = "\t", # Tab-delimited, not comma
  stringsAsFactors = FALSE
)

# Integrative analysis results
cat("Loading integrative analysis results...\n")
tes_direct <- read_csv("SRF_Eva_integrated_analysis/output/results/10_direct_targets/TES_direct_targets_all_genes.csv")
tead1_direct <- read_csv("SRF_Eva_integrated_analysis/output/results/10_direct_targets/TEAD1_direct_targets_all_genes.csv")
tes_specific <- read_csv("SRF_Eva_integrated_analysis/output/results/10_direct_targets/TES_specific_targets_all_genes.csv")

# Cut&Tag annotated peaks
cat("Loading Cut&Tag peak annotations...\n")
tes_peaks_annotated <- read_csv("SRF_Eva_CUTandTAG/results/07_analysis_narrow/TES_peaks_annotated.csv")
tead1_peaks_annotated <- read_csv("SRF_Eva_CUTandTAG/results/07_analysis_narrow/TEAD1_peaks_annotated.csv")

# DiffBind results
cat("Loading DiffBind results...\n")
diffbind_tes_vs_tesmut <- read_csv("SRF_Eva_CUTandTAG/results/07_analysis_narrow/DiffBind_TES_vs_TESmut.csv")

# MSigDB Hallmark gene sets
cat("Loading MSigDB Hallmark gene sets...\n")
msigdb_hallmark <- msigdbr(species = "Homo sapiens", collection = "H")
msigdb_c2 <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME")

cat("Data loading complete!\n")
cat("  DEGs: ", nrow(deseq_results), "\n")
cat("  TES direct targets: ", nrow(tes_direct), "\n")
cat("  TEAD1 direct targets: ", nrow(tead1_direct), "\n")
cat("  TES peaks annotated: ", nrow(tes_peaks_annotated), "\n")

################################################################################
# HELPER FUNCTIONS
################################################################################

# Function to perform GO enrichment with proper background
perform_GO_enrichment <- function(gene_list, background_genes, ont = "BP",
                                  pvalueCutoff = 0.05, qvalueCutoff = 0.05,
                                  minGSSize = 10, maxGSSize = 500) {
  # Convert gene symbols to Entrez IDs
  gene_entrez <- bitr(gene_list,
    fromType = "SYMBOL", toType = "ENTREZID",
    OrgDb = org.Hs.eg.db, drop = TRUE
  )

  background_entrez <- bitr(background_genes,
    fromType = "SYMBOL", toType = "ENTREZID",
    OrgDb = org.Hs.eg.db, drop = TRUE
  )

  # Run enrichGO
  ego <- enrichGO(
    gene = gene_entrez$ENTREZID,
    universe = background_entrez$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = ont,
    pAdjustMethod = "BH",
    pvalueCutoff = pvalueCutoff,
    qvalueCutoff = qvalueCutoff,
    minGSSize = minGSSize,
    maxGSSize = maxGSSize,
    readable = TRUE
  )

  return(ego)
}

# Function to perform GSEA using MSigDB
perform_GSEA_msigdb <- function(ranked_genes, gene_sets, pvalueCutoff = 0.05) {
  # Convert to Entrez IDs
  gene_mapping <- bitr(names(ranked_genes),
    fromType = "SYMBOL", toType = "ENTREZID",
    OrgDb = org.Hs.eg.db, drop = TRUE
  )

  # Create ranked list with Entrez IDs
  ranked_entrez <- ranked_genes[gene_mapping$SYMBOL]
  names(ranked_entrez) <- as.character(gene_mapping$ENTREZID) # ENSURE CHARACTER
  ranked_entrez <- sort(ranked_entrez, decreasing = TRUE)

  # Prepare gene sets
  # msigdbr columns: gs_name, entrez_gene OR human_entrez_gene OR gene_symbol
  gene_sets_df <- as.data.frame(gene_sets)

  # Check which column name is available
  # Modern msigdbr uses: ncbi_gene (not entrez_gene)
  if ("ncbi_gene" %in% colnames(gene_sets_df)) {
    msigdb_t2g <- gene_sets_df %>% dplyr::select(gs_name, gene = ncbi_gene)
  } else if ("entrez_gene" %in% colnames(gene_sets_df)) {
    msigdb_t2g <- gene_sets_df %>% dplyr::select(gs_name, gene = entrez_gene)
  } else if ("human_entrez_gene" %in% colnames(gene_sets_df)) {
    msigdb_t2g <- gene_sets_df %>% dplyr::select(gs_name, gene = human_entrez_gene)
  } else if ("entrez_id" %in% colnames(gene_sets_df)) {
    msigdb_t2g <- gene_sets_df %>% dplyr::select(gs_name, gene = entrez_id)
  } else {
    # Print available columns for debugging
    cat("Available columns in gene_sets:", paste(colnames(gene_sets_df), collapse = ", "), "\n")
    stop("Cannot find Entrez gene column in msigdbr results")
  }

  # Ensure gene column is character
  msigdb_t2g$gene <- as.character(msigdb_t2g$gene)

  # Debug: check gene ID formats match
  cat("DEBUG - Sample ranked gene IDs:", paste(head(names(ranked_entrez), 10), collapse = ", "), "\n")
  cat("DEBUG - Sample TERM2GENE IDs:", paste(head(unique(msigdb_t2g$gene), 10), collapse = ", "), "\n")
  cat("DEBUG - Overlap check:", sum(names(ranked_entrez) %in% msigdb_t2g$gene), "genes overlap out of", length(ranked_entrez), "\n")

  # Run GSEA
  gsea_result <- GSEA(
    geneList = ranked_entrez,
    TERM2GENE = msigdb_t2g,
    pvalueCutoff = pvalueCutoff,
    pAdjustMethod = "BH",
    minGSSize = 15,
    maxGSSize = 500
  )

  return(gsea_result)
}

# Function to save enrichment results
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
    # Dot plot
    p1 <- dotplot(enrich_obj, showCategory = 20) +
      ggtitle(paste(prefix, "- Top 20 GO Terms")) +
      theme(plot.title = element_text(hjust = 0.5))
    ggsave(file.path(output_dir, "plots", paste0(prefix, "_dotplot.pdf")),
      p1,
      width = 10, height = 8
    )

    # Bar plot
    p2 <- barplot(enrich_obj, showCategory = 20) +
      ggtitle(paste(prefix, "- Top 20 GO Terms"))
    ggsave(file.path(output_dir, "plots", paste0(prefix, "_barplot.pdf")),
      p2,
      width = 10, height = 8
    )
  }

  return(results_df)
}

# Function to extract migration-related terms
extract_migration_terms <- function(enrich_results) {
  if (is.null(enrich_results) || nrow(enrich_results) == 0) {
    return(NULL)
  }

  migration_keywords <- c(
    "migration", "motility", "adhesion", "invasion",
    "metastasis", "locomotion", "chemotaxis"
  )

  pattern <- paste(migration_keywords, collapse = "|")
  migration_terms <- enrich_results %>%
    filter(grepl(pattern, Description, ignore.case = TRUE)) %>%
    arrange(qvalue)

  return(migration_terms)
}

# Function to extract gene symbols from annotated peaks
# ChIPseeker annotated peaks have geneId (Entrez) but not SYMBOL
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
      gene_mapping <- bitr(entrez_ids,
        fromType = "ENTREZID", toType = "SYMBOL",
        OrgDb = org.Hs.eg.db, drop = TRUE
      )
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

################################################################################
# APPROACH 1: DIRECT TARGETS ONLY
################################################################################
cat("\n=================================================================\n")
cat("APPROACH 1: Direct Targets Only (Baseline)\n")
cat("=================================================================\n")

approach1_dir <- dirs$approach1

# TES direct targets
tes_direct_genes <- tes_direct$gene_symbol

# Background: all genes with reasonable expression
all_expressed_genes <- deseq_results$gene_symbol[!is.na(deseq_results$baseMean) &
  deseq_results$baseMean > 10 &
  !is.na(deseq_results$gene_symbol)]

cat("TES direct targets:", length(tes_direct_genes), "\n")
cat("Background (expressed genes):", length(all_expressed_genes), "\n")

# Save gene list
write.table(tes_direct_genes,
  file.path(approach1_dir, "gene_lists", "TES_direct_targets.txt"),
  row.names = FALSE, col.names = FALSE, quote = FALSE
)

# GO enrichment
cat("Running GO enrichment...\n")
ego_approach1 <- perform_GO_enrichment(
  gene_list = tes_direct_genes,
  background_genes = all_expressed_genes,
  ont = "BP"
)

results_approach1 <- save_enrichment_results(ego_approach1, "approach1_TES_direct", approach1_dir)

# Extract migration terms
if (!is.null(results_approach1)) {
  migration_approach1 <- extract_migration_terms(results_approach1)
  if (!is.null(migration_approach1) && nrow(migration_approach1) > 0) {
    write_csv(
      migration_approach1,
      file.path(approach1_dir, "results", "approach1_migration_terms.csv")
    )
    cat("  Found", nrow(migration_approach1), "migration-related terms\n")
    cat("  Top migration term rank:", which(results_approach1$ID == migration_approach1$ID[1]), "\n")
  }
}

################################################################################
# APPROACH 2: DOWNREGULATED DIRECT TARGETS
################################################################################
cat("\n=================================================================\n")
cat("APPROACH 2: Downregulated Direct Targets (Mechanistically Informed)\n")
cat("=================================================================\n")

approach2_dir <- dirs$approach2

# Filter for downregulated genes (TES is a repressor)
tes_direct_down <- tes_direct %>%
  filter(log2FoldChange < 0, padj < 0.05)

tes_direct_down_genes <- tes_direct_down$gene_symbol

cat("TES downregulated direct targets:", length(tes_direct_down_genes), "\n")

# Save gene list
write.table(tes_direct_down_genes,
  file.path(approach2_dir, "gene_lists", "TES_direct_downregulated.txt"),
  row.names = FALSE, col.names = FALSE, quote = FALSE
)

# Save detailed table
write_csv(
  tes_direct_down %>%
    arrange(padj) %>%
    select(gene_symbol, baseMean, log2FoldChange, padj, tes_bound, tead1_bound),
  file.path(approach2_dir, "gene_lists", "TES_direct_downregulated_detailed.csv")
)

# GO enrichment
cat("Running GO enrichment...\n")
ego_approach2 <- perform_GO_enrichment(
  gene_list = tes_direct_down_genes,
  background_genes = all_expressed_genes,
  ont = "BP"
)

results_approach2 <- save_enrichment_results(ego_approach2, "approach2_TES_direct_down", approach2_dir)

# Extract migration terms
if (!is.null(results_approach2)) {
  migration_approach2 <- extract_migration_terms(results_approach2)
  if (!is.null(migration_approach2) && nrow(migration_approach2) > 0) {
    write_csv(
      migration_approach2,
      file.path(approach2_dir, "results", "approach2_migration_terms.csv")
    )
    cat("  Found", nrow(migration_approach2), "migration-related terms\n")
    cat("  Top migration term rank:", which(results_approach2$ID == migration_approach2$ID[1]), "\n")
  }
}

################################################################################
# APPROACH 3: PROMOTER PEAKS ONLY
################################################################################
cat("\n=================================================================\n")
cat("APPROACH 3: Promoter Peaks Only\n")
cat("=================================================================\n")

approach3_dir <- dirs$approach3

# Filter for promoter peaks
promoter_keywords <- c("Promoter", "5' UTR", "1st Exon")
tes_promoter_peaks <- tes_peaks_annotated %>%
  filter(grepl(paste(promoter_keywords, collapse = "|"), annotation, ignore.case = TRUE))

cat("Total TES peaks:", nrow(tes_peaks_annotated), "\n")
cat("TES promoter peaks:", nrow(tes_promoter_peaks), "\n")

# Extract genes from promoter peaks
tes_promoter_genes <- extract_genes_from_peaks(tes_promoter_peaks)

cat("Unique genes with promoter peaks:", length(tes_promoter_genes), "\n")

# Intersect with DEGs
degs <- deseq_results %>%
  filter(!is.na(padj) & padj < 0.05 & !is.na(gene_symbol))
tes_promoter_degs <- intersect(tes_promoter_genes, degs$gene_symbol)

cat("Promoter genes that are DEGs:", length(tes_promoter_degs), "\n")

# Save gene lists
write.table(tes_promoter_genes,
  file.path(approach3_dir, "gene_lists", "TES_promoter_peak_genes.txt"),
  row.names = FALSE, col.names = FALSE, quote = FALSE
)

write.table(tes_promoter_degs,
  file.path(approach3_dir, "gene_lists", "TES_promoter_DEGs.txt"),
  row.names = FALSE, col.names = FALSE, quote = FALSE
)

# GO enrichment on promoter DEGs
cat("Running GO enrichment on promoter DEGs...\n")
ego_approach3 <- perform_GO_enrichment(
  gene_list = tes_promoter_degs,
  background_genes = all_expressed_genes,
  ont = "BP"
)

results_approach3 <- save_enrichment_results(ego_approach3, "approach3_TES_promoter_DEGs", approach3_dir)

# Extract migration terms
if (!is.null(results_approach3)) {
  migration_approach3 <- extract_migration_terms(results_approach3)
  if (!is.null(migration_approach3) && nrow(migration_approach3) > 0) {
    write_csv(
      migration_approach3,
      file.path(approach3_dir, "results", "approach3_migration_terms.csv")
    )
    cat("  Found", nrow(migration_approach3), "migration-related terms\n")
    cat("  Top migration term rank:", which(results_approach3$ID == migration_approach3$ID[1]), "\n")
  }
}

################################################################################
# APPROACH 4: HIGH-CONFIDENCE PEAKS
################################################################################
cat("\n=================================================================\n")
cat("APPROACH 4: High-Confidence Peaks (Top 50% by signal)\n")
cat("=================================================================\n")

approach4_dir <- dirs$approach4

# Filter for high-confidence peaks (top 50% by qValue or fold enrichment)
# Note: In narrowPeak format from ChIPseeker:
#  V7 = signalValue (fold enrichment)
#  V8 = pValue (-log10)
#  V9 = qValue (-log10)

# Check which columns are available
if ("V9" %in% colnames(tes_peaks_annotated)) {
  # Use qValue (-log10, so higher is better)
  cat("Using V9 (qValue -log10) for filtering...\n")
  median_qval <- median(tes_peaks_annotated$V9, na.rm = TRUE)
  tes_highconf_peaks <- tes_peaks_annotated %>%
    filter(V9 >= median_qval)
} else if ("V7" %in% colnames(tes_peaks_annotated)) {
  # Use signalValue (fold enrichment)
  cat("Using V7 (signalValue) for filtering...\n")
  median_signal <- median(tes_peaks_annotated$V7, na.rm = TRUE)
  tes_highconf_peaks <- tes_peaks_annotated %>%
    filter(V7 >= median_signal)
} else if ("V8" %in% colnames(tes_peaks_annotated)) {
  # Use pValue (-log10, so higher is better)
  cat("Using V8 (pValue -log10) for filtering...\n")
  median_pval <- median(tes_peaks_annotated$V8, na.rm = TRUE)
  tes_highconf_peaks <- tes_peaks_annotated %>%
    filter(V8 >= median_pval)
} else {
  cat("WARNING: No suitable quality columns found, using all peaks\n")
  tes_highconf_peaks <- tes_peaks_annotated
}

cat("Total TES peaks:", nrow(tes_peaks_annotated), "\n")
cat("High-confidence peaks:", nrow(tes_highconf_peaks), "\n")

# Extract genes
tes_highconf_genes <- extract_genes_from_peaks(tes_highconf_peaks)
cat("Unique genes with high-conf peaks:", length(tes_highconf_genes), "\n")

# Intersect with DEGs
tes_highconf_degs <- intersect(tes_highconf_genes, degs$gene_symbol)
cat("High-conf genes that are DEGs:", length(tes_highconf_degs), "\n")

# Save gene lists
write.table(tes_highconf_genes,
  file.path(approach4_dir, "gene_lists", "TES_highconf_peak_genes.txt"),
  row.names = FALSE, col.names = FALSE, quote = FALSE
)

write.table(tes_highconf_degs,
  file.path(approach4_dir, "gene_lists", "TES_highconf_DEGs.txt"),
  row.names = FALSE, col.names = FALSE, quote = FALSE
)

# GO enrichment
cat("Running GO enrichment on high-confidence DEGs...\n")
ego_approach4 <- perform_GO_enrichment(
  gene_list = tes_highconf_degs,
  background_genes = all_expressed_genes,
  ont = "BP"
)

results_approach4 <- save_enrichment_results(ego_approach4, "approach4_TES_highconf_DEGs", approach4_dir)

# Extract migration terms
if (!is.null(results_approach4)) {
  migration_approach4 <- extract_migration_terms(results_approach4)
  if (!is.null(migration_approach4) && nrow(migration_approach4) > 0) {
    write_csv(
      migration_approach4,
      file.path(approach4_dir, "results", "approach4_migration_terms.csv")
    )
    cat("  Found", nrow(migration_approach4), "migration-related terms\n")
    cat("  Top migration term rank:", which(results_approach4$ID == migration_approach4$ID[1]), "\n")
  }
}

################################################################################
# APPROACH 5: DIFFERENTIAL BINDING + DIFFERENTIAL EXPRESSION
################################################################################
cat("\n=================================================================\n")
cat("APPROACH 5: Differential Binding + Differential Expression\n")
cat("=================================================================\n")

approach5_dir <- dirs$approach5

# Filter for TES-specific differential peaks (TES > TESmut)
if ("FDR" %in% colnames(diffbind_tes_vs_tesmut)) {
  tes_specific_diffbind <- diffbind_tes_vs_tesmut %>%
    filter(FDR < 0.05, Fold > 1.5)

  cat("Total differential sites:", nrow(diffbind_tes_vs_tesmut), "\n")
  cat("TES-specific sites (vs TESmut):", nrow(tes_specific_diffbind), "\n")

  # DiffBind results don't have gene annotations
  # Overlap with annotated peaks to get genes
  cat("Overlapping differential peaks with annotated peaks...\n")

  # Create GRanges for differential peaks
  diffbind_gr <- GRanges(
    seqnames = tes_specific_diffbind$seqnames,
    ranges = IRanges(
      start = tes_specific_diffbind$start,
      end = tes_specific_diffbind$end
    )
  )

  # Create GRanges for annotated peaks
  annotated_gr <- GRanges(
    seqnames = tes_peaks_annotated$seqnames,
    ranges = IRanges(
      start = tes_peaks_annotated$start,
      end = tes_peaks_annotated$end
    )
  )

  # Find overlaps
  overlaps <- findOverlaps(diffbind_gr, annotated_gr)
  overlapping_annotated_idx <- subjectHits(overlaps)

  # Get genes from overlapping annotated peaks
  overlapping_peaks <- tes_peaks_annotated[overlapping_annotated_idx, ]
  tes_diffbind_genes <- extract_genes_from_peaks(overlapping_peaks)

  if (length(tes_diffbind_genes) > 0) {
    cat("Unique genes at differential sites:", length(tes_diffbind_genes), "\n")

    # Intersect with downregulated DEGs (since TES is repressor)
    tes_diffbind_down <- intersect(tes_diffbind_genes, tes_direct_down_genes)
    cat("Differential binding + downregulated:", length(tes_diffbind_down), "\n")

    # DEBUG: Check what we have
    cat("DEBUG - TES diffbind genes:", length(tes_diffbind_genes), "\n")
    cat("DEBUG - TES direct down genes:", length(tes_direct_down_genes), "\n")
    cat("DEBUG - Overlap (diffbind + down):", length(tes_diffbind_down), "\n")

    # Save gene lists
    if (length(tes_diffbind_genes) > 0) {
      write.table(tes_diffbind_genes,
        file.path(approach5_dir, "gene_lists", "TES_diffbind_genes.txt"),
        row.names = FALSE, col.names = FALSE, quote = FALSE
      )
    }

    if (length(tes_diffbind_down) > 0) {
      write.table(tes_diffbind_down,
        file.path(approach5_dir, "gene_lists", "TES_diffbind_down.txt"),
        row.names = FALSE, col.names = FALSE, quote = FALSE
      )
    }

    # GO enrichment - try with all diffbind genes if down is too few
    genes_for_enrichment <- if (length(tes_diffbind_down) >= 10) {
      cat("Using differential binding + downregulated genes for enrichment\n")
      tes_diffbind_down
    } else if (length(tes_diffbind_genes) >= 10) {
      cat("Too few downregulated genes, using all differential binding genes instead\n")
      tes_diffbind_genes
    } else {
      cat("Too few genes for enrichment analysis\n")
      NULL
    }

    if (!is.null(genes_for_enrichment) && length(genes_for_enrichment) >= 10) {
      cat("Running GO enrichment on", length(genes_for_enrichment), "genes...\n")
      ego_approach5 <- perform_GO_enrichment(
        gene_list = genes_for_enrichment,
        background_genes = all_expressed_genes,
        ont = "BP"
      )

      results_approach5 <- save_enrichment_results(ego_approach5, "approach5_diffbind", approach5_dir)

      # Extract migration terms
      if (!is.null(results_approach5)) {
        migration_approach5 <- extract_migration_terms(results_approach5)
        if (!is.null(migration_approach5) && nrow(migration_approach5) > 0) {
          write_csv(
            migration_approach5,
            file.path(approach5_dir, "results", "approach5_migration_terms.csv")
          )
          cat("  Found", nrow(migration_approach5), "migration-related terms\n")
          cat("  Top migration term rank:", which(results_approach5$ID == migration_approach5$ID[1]), "\n")
        }
      }
    } else {
      results_approach5 <- NULL
    }
  }
} else {
  cat("DiffBind results not in expected format, skipping Approach 5\n")
  results_approach5 <- NULL
}

################################################################################
# APPROACH 6: MIGRATION GENE-FOCUSED (Hypothesis-Driven)
################################################################################
# cat("\n=================================================================\n")
# cat("APPROACH 6: Migration Gene-Focused (Hypothesis-Driven)\n")
# cat("=================================================================\n")

# approach6_dir <- dirs$approach6

# # Load specific gene sets of interest:
# # 1. GOBP_NEURON_MIGRATION (from GO)
# # 2. CELL_PROLIFERATION_GO_0008283 (from GO)
# # 3. HALLMARK_APOPTOSIS (from MSigDB Hallmark)

# # Get GO gene sets from msigdbr
# msigdb_go <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP")

# # Gene set 1: Neuron migration (GO:0001764)
# neuron_migration_genes <- msigdb_go %>%
#   filter(grepl("NEURON.*MIGRATION|GO:0001764", gs_name, ignore.case = TRUE)) %>%
#   pull(gene_symbol) %>%
#   unique()

# cat("Neuron migration gene set size:", length(neuron_migration_genes), "\n")

# # Gene set 2: Cell proliferation (GO:0008283)
# cell_prolif_genes <- msigdb_go %>%
#   filter(grepl("CELL_PROLIFERATION|GO:0008283", gs_name, ignore.case = TRUE)) %>%
#   pull(gene_symbol) %>%
#   unique()

# cat("Cell proliferation gene set size:", length(cell_prolif_genes), "\n")

# # Gene set 3: Apoptosis from Hallmark
# apoptosis_genes <- msigdb_hallmark %>%
#   filter(gs_name == "HALLMARK_APOPTOSIS") %>%
#   pull(gene_symbol)

# cat("Apoptosis gene set size:", length(apoptosis_genes), "\n")

# # Hypergeometric tests for each gene set
# cat("\n=== Hypergeometric enrichment tests ===\n")

# # Test 1: Neuron migration
# if (length(neuron_migration_genes) > 0) {
#   tes_direct_in_neuron_mig <- intersect(tes_direct_genes, neuron_migration_genes)
#   cat("\nNeuron migration:\n")
#   cat("  TES direct targets in gene set:", length(tes_direct_in_neuron_mig), "\n")

#   if (length(tes_direct_in_neuron_mig) > 0) {
#     universe_size <- length(all_expressed_genes)
#     geneset_in_universe <- length(intersect(neuron_migration_genes, all_expressed_genes))
#     tes_direct_size <- length(tes_direct_genes)
#     overlap <- length(tes_direct_in_neuron_mig)

#     phyper_pval <- phyper(q = overlap - 1,
#                           m = geneset_in_universe,
#                           n = universe_size - geneset_in_universe,
#                           k = tes_direct_size,
#                           lower.tail = FALSE)
#     cat("  Hypergeometric p-value:", phyper_pval, "\n")

#     write.table(tes_direct_in_neuron_mig,
#                 file.path(approach6_dir, "gene_lists", "TES_direct_neuron_migration.txt"),
#                 row.names = FALSE, col.names = FALSE, quote = FALSE)
#   }
# }

# # Test 2: Cell proliferation
# if (length(cell_prolif_genes) > 0) {
#   tes_direct_in_prolif <- intersect(tes_direct_genes, cell_prolif_genes)
#   cat("\nCell proliferation:\n")
#   cat("  TES direct targets in gene set:", length(tes_direct_in_prolif), "\n")

#   if (length(tes_direct_in_prolif) > 0) {
#     universe_size <- length(all_expressed_genes)
#     geneset_in_universe <- length(intersect(cell_prolif_genes, all_expressed_genes))
#     tes_direct_size <- length(tes_direct_genes)
#     overlap <- length(tes_direct_in_prolif)

#     phyper_pval <- phyper(q = overlap - 1,
#                           m = geneset_in_universe,
#                           n = universe_size - geneset_in_universe,
#                           k = tes_direct_size,
#                           lower.tail = FALSE)
#     cat("  Hypergeometric p-value:", phyper_pval, "\n")

#     write.table(tes_direct_in_prolif,
#                 file.path(approach6_dir, "gene_lists", "TES_direct_proliferation.txt"),
#                 row.names = FALSE, col.names = FALSE, quote = FALSE)
#   }
# }

# # Test 3: Apoptosis
# if (length(apoptosis_genes) > 0) {
#   tes_direct_in_apoptosis <- intersect(tes_direct_genes, apoptosis_genes)
#   cat("\nApoptosis:\n")
#   cat("  TES direct targets in gene set:", length(tes_direct_in_apoptosis), "\n")

#   if (length(tes_direct_in_apoptosis) > 0) {
#     universe_size <- length(all_expressed_genes)
#     geneset_in_universe <- length(intersect(apoptosis_genes, all_expressed_genes))
#     tes_direct_size <- length(tes_direct_genes)
#     overlap <- length(tes_direct_in_apoptosis)

#     phyper_pval <- phyper(q = overlap - 1,
#                           m = geneset_in_universe,
#                           n = universe_size - geneset_in_universe,
#                           k = tes_direct_size,
#                           lower.tail = FALSE)
#     cat("  Hypergeometric p-value:", phyper_pval, "\n")

#     write.table(tes_direct_in_apoptosis,
#                 file.path(approach6_dir, "gene_lists", "TES_direct_apoptosis.txt"),
#                 row.names = FALSE, col.names = FALSE, quote = FALSE)
#   }
# }

# # Create detailed tables for each gene set
# if (exists("tes_direct_in_neuron_mig") && length(tes_direct_in_neuron_mig) > 0) {
#   neuron_mig_detailed <- tes_direct %>%
#     filter(gene_symbol %in% tes_direct_in_neuron_mig) %>%
#     arrange(padj) %>%
#     select(gene_symbol, baseMean, log2FoldChange, padj, tes_bound, tead1_bound)

#   write_csv(neuron_mig_detailed,
#             file.path(approach6_dir, "gene_lists", "TES_direct_neuron_migration_detailed.csv"))
# }

# if (exists("tes_direct_in_prolif") && length(tes_direct_in_prolif) > 0) {
#   prolif_detailed <- tes_direct %>%
#     filter(gene_symbol %in% tes_direct_in_prolif) %>%
#     arrange(padj) %>%
#     select(gene_symbol, baseMean, log2FoldChange, padj, tes_bound, tead1_bound)

#   write_csv(prolif_detailed,
#             file.path(approach6_dir, "gene_lists", "TES_direct_proliferation_detailed.csv"))
# }

# if (exists("tes_direct_in_apoptosis") && length(tes_direct_in_apoptosis) > 0) {
#   apoptosis_detailed <- tes_direct %>%
#     filter(gene_symbol %in% tes_direct_in_apoptosis) %>%
#     arrange(padj) %>%
#     select(gene_symbol, baseMean, log2FoldChange, padj, tes_bound, tead1_bound)

#   write_csv(apoptosis_detailed,
#             file.path(approach6_dir, "gene_lists", "TES_direct_apoptosis_detailed.csv"))
# }

# # GSEA analysis using ranked genes
# cat("\n=== GSEA analysis ===\n")

# # Create ranked list by log2FoldChange * -log10(padj)
# # Cap padj at a minimum value to avoid Inf values
# degs_ranked <- deseq_results %>%
#   filter(!is.na(padj) & !is.na(log2FoldChange) & baseMean > 10) %>%
#   mutate(padj_capped = pmax(padj, 1e-300)) %>%  # Prevent -log10 from going to Inf
#   mutate(rank_metric = log2FoldChange * -log10(padj_capped)) %>%
#   filter(is.finite(rank_metric)) %>%  # Remove any remaining non-finite values
#   arrange(desc(rank_metric))

# ranked_genes <- setNames(degs_ranked$rank_metric, degs_ranked$gene_symbol)

# # Combine the three gene sets of interest into a custom collection
# custom_genesets <- bind_rows(
#   msigdb_go %>% filter(grepl("NEURON.*MIGRATION|GO:0001764", gs_name, ignore.case = TRUE)) %>%
#     mutate(gs_name = "GOBP_NEURON_MIGRATION"),
#   msigdb_go %>% filter(grepl("CELL_PROLIFERATION|GO:0008283", gs_name, ignore.case = TRUE)) %>%
#     mutate(gs_name = "GOBP_CELL_PROLIFERATION"),
#   msigdb_hallmark %>% filter(gs_name == "HALLMARK_APOPTOSIS")
# )

# cat("Running GSEA on custom gene sets (neuron migration, proliferation, apoptosis)...\n")

# if (nrow(custom_genesets) > 0) {
#   gsea_custom <- perform_GSEA_msigdb(ranked_genes, custom_genesets, pvalueCutoff = 0.25)

#   if (!is.null(gsea_custom) && nrow(gsea_custom) > 0) {
#     gsea_results <- as.data.frame(gsea_custom)
#     write_csv(gsea_results,
#               file.path(approach6_dir, "results", "GSEA_custom_genesets.csv"))

#     cat("\nGSEA results:\n")
#     print(gsea_results[, c("ID", "NES", "pvalue", "p.adjust")])

#     # Plot GSEA results
#     pdf(file.path(approach6_dir, "plots", "GSEA_custom_dotplot.pdf"), width = 10, height = 6)
#     print(dotplot(gsea_custom, showCategory = 10, title = "GSEA: Migration, Proliferation, Apoptosis"))
#     dev.off()

#     # Individual GSEA plots for significant gene sets
#     for (geneset_id in gsea_results$ID) {
#       safe_name <- gsub("[^A-Za-z0-9_]", "_", geneset_id)
#       pdf(file.path(approach6_dir, "plots", paste0("GSEA_", safe_name, ".pdf")),
#           width = 10, height = 6)
#       print(gseaplot2(gsea_custom,
#                       geneSetID = geneset_id,
#                       title = geneset_id))
#       dev.off()
#     }
#   } else {
#     cat("No significant GSEA results found (p < 0.25)\n")
#   }
# } else {
#   cat("Could not create custom gene sets for GSEA\n")
# }

# # Create visualizations for the hypergeometric enrichment results
# cat("\n=== Creating visualizations ===\n")

# # Prepare data for enrichment barplot
# enrichment_results <- data.frame(
#   GeneSet = c("Neuron\nMigration", "Cell\nProliferation", "Apoptosis"),
#   Overlap = c(
#     if(exists("tes_direct_in_neuron_mig")) length(tes_direct_in_neuron_mig) else 0,
#     if(exists("tes_direct_in_prolif")) length(tes_direct_in_prolif) else 0,
#     if(exists("tes_direct_in_apoptosis")) length(tes_direct_in_apoptosis) else 0
#   ),
#   GeneSetSize = c(
#     length(neuron_migration_genes),
#     length(cell_prolif_genes),
#     length(apoptosis_genes)
#   ),
#   PValue = c(
#     if(exists("tes_direct_in_neuron_mig") && length(tes_direct_in_neuron_mig) > 0) {
#       phyper(length(tes_direct_in_neuron_mig) - 1,
#              length(intersect(neuron_migration_genes, all_expressed_genes)),
#              length(all_expressed_genes) - length(intersect(neuron_migration_genes, all_expressed_genes)),
#              length(tes_direct_genes), lower.tail = FALSE)
#     } else NA,
#     if(exists("tes_direct_in_prolif") && length(tes_direct_in_prolif) > 0) {
#       phyper(length(tes_direct_in_prolif) - 1,
#              length(intersect(cell_prolif_genes, all_expressed_genes)),
#              length(all_expressed_genes) - length(intersect(cell_prolif_genes, all_expressed_genes)),
#              length(tes_direct_genes), lower.tail = FALSE)
#     } else NA,
#     if(exists("tes_direct_in_apoptosis") && length(tes_direct_in_apoptosis) > 0) {
#       phyper(length(tes_direct_in_apoptosis) - 1,
#              length(intersect(apoptosis_genes, all_expressed_genes)),
#              length(all_expressed_genes) - length(intersect(apoptosis_genes, all_expressed_genes)),
#              length(tes_direct_genes), lower.tail = FALSE)
#     } else NA
#   )
# ) %>%
#   filter(!is.na(PValue)) %>%
#   mutate(
#     NegLog10P = -log10(PValue),
#     PercentOverlap = (Overlap / GeneSetSize) * 100,
#     Significance = case_when(
#       PValue < 0.001 ~ "p < 0.001",
#       PValue < 0.01 ~ "p < 0.01",
#       PValue < 0.05 ~ "p < 0.05",
#       TRUE ~ "n.s."
#     )
#   )

# # Plot 1: Enrichment barplot with -log10(p-value)
# pdf(file.path(approach6_dir, "plots", "hypergeometric_enrichment.pdf"), width = 8, height = 6)
# p1 <- ggplot(enrichment_results, aes(x = reorder(GeneSet, -NegLog10P), y = NegLog10P, fill = GeneSet)) +
#   geom_bar(stat = "identity") +
#   geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
#   geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "darkred") +
#   geom_text(aes(label = sprintf("p = %.2e", PValue)), vjust = -0.5, size = 3.5) +
#   labs(title = "TES Direct Target Enrichment in Biological Gene Sets",
#        subtitle = "Hypergeometric test",
#        x = "Gene Set",
#        y = "-log10(p-value)") +
#   theme_bw() +
#   theme(legend.position = "none",
#         axis.text.x = element_text(size = 11),
#         plot.title = element_text(size = 14, face = "bold")) +
#   scale_fill_brewer(palette = "Set2")
# print(p1)
# dev.off()

# # Plot 2: Overlap counts
# pdf(file.path(approach6_dir, "plots", "geneset_overlap_counts.pdf"), width = 8, height = 6)
# p2 <- ggplot(enrichment_results, aes(x = reorder(GeneSet, -Overlap), y = Overlap, fill = GeneSet)) +
#   geom_bar(stat = "identity") +
#   geom_text(aes(label = sprintf("%d genes\n(%.1f%%)", Overlap, PercentOverlap)),
#             vjust = -0.5, size = 3.5) +
#   labs(title = "TES Direct Targets in Biological Gene Sets",
#        subtitle = sprintf("Out of %d total TES direct targets", length(tes_direct_genes)),
#        x = "Gene Set",
#        y = "Number of TES Direct Targets") +
#   theme_bw() +
#   theme(legend.position = "none",
#         axis.text.x = element_text(size = 11),
#         plot.title = element_text(size = 14, face = "bold")) +
#   scale_fill_brewer(palette = "Set2")
# print(p2)
# dev.off()

# # Plot 3: Expression heatmap for top genes in each category
# if (exists("tes_direct_in_neuron_mig") && length(tes_direct_in_neuron_mig) > 0) {
#   top_neuron_genes <- tes_direct %>%
#     filter(gene_symbol %in% tes_direct_in_neuron_mig) %>%
#     arrange(padj) %>%
#     head(20) %>%
#     pull(gene_symbol)

#   heatmap_data <- tes_direct %>%
#     filter(gene_symbol %in% top_neuron_genes) %>%
#     select(gene_symbol, log2FoldChange) %>%
#     arrange(log2FoldChange)

#   pdf(file.path(approach6_dir, "plots", "neuron_migration_top_genes_heatmap.pdf"),
#       width = 6, height = 8)
#   p3 <- ggplot(heatmap_data, aes(x = 1, y = reorder(gene_symbol, log2FoldChange),
#                                  fill = log2FoldChange)) +
#     geom_tile(color = "white") +
#     geom_text(aes(label = sprintf("%.2f", log2FoldChange)), size = 3) +
#     scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
#                          name = "log2FC") +
#     labs(title = "Top 20 Neuron Migration Genes",
#          subtitle = "TES Direct Targets",
#          x = "", y = "") +
#     theme_minimal() +
#     theme(axis.text.x = element_blank(),
#           axis.ticks.x = element_blank(),
#           plot.title = element_text(size = 14, face = "bold"))
#   print(p3)
#   dev.off()
# }

# cat("Plots created successfully!\n")

################################################################################
# TIER 1: PROGRESSION OF SPECIFICITY
################################################################################
# cat("\n=================================================================\n")
# cat("TIER 1: Show Progression of Specificity\n")
# cat("=================================================================\n")

# tier1_dir <- dirs$tier1

# # Compile all approaches
# progression_data <- data.frame(
#   Approach = c("All TES peaks",
#                "Direct targets (any)",
#                "Direct downregulated",
#                "Promoter DEGs",
#                "High-conf DEGs",
#                "DiffBind + Down"),
#   N_genes = c(
#     length(unique(tes_peaks_annotated$SYMBOL[!is.na(tes_peaks_annotated$SYMBOL)])),
#     length(tes_direct_genes),
#     length(tes_direct_down_genes),
#     length(tes_promoter_degs),
#     length(tes_highconf_degs),
#     ifelse(exists("tes_diffbind_down"), length(tes_diffbind_down), NA)
#   ),
#   N_sig_terms = c(
#     NA,  # Would need to load original analysis
#     ifelse(!is.null(results_approach1), sum(results_approach1$qvalue < 0.05), 0),
#     ifelse(!is.null(results_approach2), sum(results_approach2$qvalue < 0.05), 0),
#     ifelse(!is.null(results_approach3), sum(results_approach3$qvalue < 0.05), 0),
#     ifelse(!is.null(results_approach4), sum(results_approach4$qvalue < 0.05), 0),
#     ifelse(exists("results_approach5") && !is.null(results_approach5),
#            sum(results_approach5$qvalue < 0.05), NA)
#   ),
#   Top_migration_rank = c(
#     NA,
#     ifelse(!is.null(migration_approach1) && nrow(migration_approach1) > 0,
#            which(results_approach1$ID == migration_approach1$ID[1])[1], NA),
#     ifelse(!is.null(migration_approach2) && nrow(migration_approach2) > 0,
#            which(results_approach2$ID == migration_approach2$ID[1])[1], NA),
#     ifelse(!is.null(migration_approach3) && nrow(migration_approach3) > 0,
#            which(results_approach3$ID == migration_approach3$ID[1])[1], NA),
#     ifelse(!is.null(migration_approach4) && nrow(migration_approach4) > 0,
#            which(results_approach4$ID == migration_approach4$ID[1])[1], NA),
#     ifelse(exists("migration_approach5") && !is.null(migration_approach5) && nrow(migration_approach5) > 0,
#            which(results_approach5$ID == migration_approach5$ID[1])[1], NA)
#   )
# )

# write_csv(progression_data,
#           file.path(tier1_dir, "results", "progression_summary.csv"))

# cat("\n=== PROGRESSION SUMMARY ===\n")
# print(progression_data)

# # Create comparison plot
# p_progression <- ggplot(progression_data %>% filter(!is.na(Top_migration_rank)),
#                         aes(x = reorder(Approach, -Top_migration_rank),
#                             y = Top_migration_rank,
#                             fill = N_genes)) +
#   geom_bar(stat = "identity") +
#   geom_text(aes(label = N_genes), hjust = -0.2) +
#   coord_flip() +
#   scale_fill_gradient(low = "lightblue", high = "darkblue") +
#   labs(title = "Progression of Specificity: Migration Term Ranking",
#        subtitle = "Lower rank = stronger enrichment",
#        x = "Approach",
#        y = "Rank of Top Migration Term",
#        fill = "N genes") +
#   theme_minimal() +
#   theme(plot.title = element_text(hjust = 0.5, face = "bold"),
#         plot.subtitle = element_text(hjust = 0.5))

# ggsave(file.path(tier1_dir, "plots", "progression_comparison.pdf"),
#        p_progression, width = 10, height = 6)

################################################################################
# TIER 2: VALIDATION WITH ORTHOGONAL APPROACH
################################################################################
# cat("\n=================================================================\n")
# cat("TIER 2: Validation with Orthogonal Approach\n")
# cat("=================================================================\n")

# tier2_dir <- dirs$tier2

# # Compare: Direct targets vs DiffBind-derived targets
# # This validates that our findings are robust across different methods

# if (exists("tes_diffbind_down") && length(tes_diffbind_down) > 0) {

#   # Venn diagram: Direct Down vs DiffBind Down
#   venn_data <- list(
#     "Direct\nDownregulated" = tes_direct_down_genes,
#     "DiffBind\n+ Down" = tes_diffbind_down
#   )

#   pdf(file.path(tier2_dir, "plots", "validation_venn.pdf"), width = 8, height = 8)
#   venn.plot <- venn.diagram(
#     x = venn_data,
#     filename = NULL,
#     fill = c("lightblue", "pink"),
#     alpha = 0.5,
#     category.names = names(venn_data),
#     main = "Validation: Direct vs DiffBind Approaches"
#   )
#   grid.draw(venn.plot)
#   dev.off()

#   # High-confidence genes (in both)
#   high_conf_validated <- intersect(tes_direct_down_genes, tes_diffbind_down)
#   cat("High-confidence validated genes:", length(high_conf_validated), "\n")

#   write.table(high_conf_validated,
#               file.path(tier2_dir, "gene_lists", "validated_highconf_genes.txt"),
#               row.names = FALSE, col.names = FALSE, quote = FALSE)

#   # Detailed table
#   validated_detailed <- tes_direct_down %>%
#     filter(gene_symbol %in% high_conf_validated) %>%
#     arrange(padj)

#   write_csv(validated_detailed,
#             file.path(tier2_dir, "gene_lists", "validated_highconf_detailed.csv"))

#   # Enrichment on validated genes
#   if (length(high_conf_validated) >= 10) {
#     cat("Running GO enrichment on validated high-confidence genes...\n")
#     ego_tier2 <- perform_GO_enrichment(
#       gene_list = high_conf_validated,
#       background_genes = all_expressed_genes,
#       ont = "BP"
#     )

#     results_tier2 <- save_enrichment_results(ego_tier2, "tier2_validated", tier2_dir)

#     # Extract migration terms
#     if (!is.null(results_tier2)) {
#       migration_tier2 <- extract_migration_terms(results_tier2)
#       if (!is.null(migration_tier2) && nrow(migration_tier2) > 0) {
#         write_csv(migration_tier2,
#                   file.path(tier2_dir, "results", "tier2_migration_terms.csv"))
#         cat("  Found", nrow(migration_tier2), "migration-related terms\n")
#         cat("  Top migration term rank:", which(results_tier2$ID == migration_tier2$ID[1]), "\n")
#       }
#     }
#   }
# }

################################################################################
# FINAL SUMMARY
################################################################################
cat("\n=================================================================\n")
cat("ANALYSIS COMPLETE - SUMMARY\n")
cat("=================================================================\n\n")

cat("All approaches have been completed!\n\n")
cat("Results are saved in:\n")
cat("  ", base_dir, "\n\n")

cat("Key findings:\n")
cat("  1. Direct targets: ", length(tes_direct_genes), " genes\n")
cat("  2. Downregulated direct: ", length(tes_direct_down_genes), " genes\n")
cat("  3. Promoter DEGs: ", length(tes_promoter_degs), " genes\n")
cat("  4. High-confidence DEGs: ", length(tes_highconf_degs), " genes\n")
# if (exists("tes_diffbind_down")) {
#   cat("  5. DiffBind + Down: ", length(tes_diffbind_down), " genes\n")
# }
# if (exists("high_conf_validated")) {
#   cat("  6. Validated high-conf: ", length(high_conf_validated), " genes\n")
# }

cat("\n=================================================================\n")
cat("Check individual approach folders for detailed results!\n")
cat("=================================================================\n")
