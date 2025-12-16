#!/usr/bin/env Rscript
#
# GENOMICPLOT METAGENE PROFILES
# Creates gene body profiles showing:
#   - Upstream (-2kb) → 5'UTR → CDS → 3'UTR → Downstream (+1kb)
# For both binding (TES/TEAD1) and methylation (TES vs GFP)
#

suppressPackageStartupMessages({
    library(GenomicPlot)
    library(GenomicFeatures)
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    library(org.Hs.eg.db)
    library(GenomicRanges)
    library(rtracklayer)
    library(dplyr)
    library(ggplot2)
})

setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_1")

cat("==========================================================\n")
cat("GENOMICPLOT METAGENE PROFILES\n")
cat("==========================================================\n")
cat("Analysis started:", as.character(Sys.time()), "\n\n")

# =============================================================================
# PATH CONFIGURATION
# =============================================================================

# BigWig files - Binding (Cut&Tag)
CUTANDTAG_DIR <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_CUTandTAG/results/06_bigwig"
TES_BINDING_BW <- file.path(CUTANDTAG_DIR, "TES_comb.bw")
TEAD1_BINDING_BW <- file.path(CUTANDTAG_DIR, "TEAD1_comb.bw")

# BigWig files - Methylation (meDIP)
MEDIP_DIR <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP/results/05_bigwig"
TES_METH_BW <- file.path(MEDIP_DIR, "TES_combined_RPKM.bw")
GFP_METH_BW <- file.path(MEDIP_DIR, "GFP_combined_RPKM.bw")

# DESeq2 results for gene filtering
DESEQ2_FILE <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_RNA/results/05_deseq2/deseq2_results_TES_vs_GFP.txt"

# Output directory
OUTPUT_BASE <- "output/23_genomicplot_metagene"
dir.create(OUTPUT_BASE, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# PHASE 1: PREPARE GENE ANNOTATIONS
# =============================================================================

cat("=== PHASE 1: Preparing Gene Annotations ===\n")

# Load TxDb
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
cat("  Loaded TxDb for hg38\n")

# Get 5-part gene features using GenomicPlot's helper
# This creates features with: promoter, 5UTR, CDS, 3UTR, downstream
cat("  Extracting 5-part gene features...\n")

# Get genes with proper structure
genes <- genes(txdb)
transcripts <- transcriptsBy(txdb, by = "gene")
cds <- cdsBy(txdb, by = "tx", use.names = TRUE)
utr5 <- fiveUTRsByTranscript(txdb, use.names = TRUE)
utr3 <- threeUTRsByTranscript(txdb, use.names = TRUE)

# Create the 5-part gene feature object using GenomicPlot
# First, let's prepare the data for plot_5parts_metagene
gf5 <- prepare_5parts_genomic_features(
    txdb = txdb,
    meta = TRUE,  # Use metagene (without introns)
    nbins = 100,
    fiveP = -2000,  # 2kb upstream
    threeP = 1000,  # 1kb downstream
    longest = TRUE  # Use longest transcript per gene
)

cat(sprintf("  Created 5-part features for %d genes\n", length(gf5)))

cat("\n")

# =============================================================================
# PHASE 2: LOAD DESEQ2 AND DEFINE GENE SETS
# =============================================================================

cat("=== PHASE 2: Defining Gene Sets ===\n")

# Load DESeq2 results
deseq2 <- read.delim(DESEQ2_FILE, stringsAsFactors = FALSE)
deseq2$gene_id_clean <- gsub("\\..*", "", deseq2$gene_id)

# Convert Ensembl to Entrez (GenomicPlot uses Entrez IDs)
ensembl_to_entrez <- mapIds(org.Hs.eg.db,
                             keys = deseq2$gene_id_clean,
                             column = "ENTREZID",
                             keytype = "ENSEMBL",
                             multiVals = "first")
deseq2$entrez_id <- ensembl_to_entrez[deseq2$gene_id_clean]

# All expressed genes (non-NA padj)
all_expressed <- deseq2 %>%
    filter(!is.na(padj), !is.na(entrez_id))
cat(sprintf("  All expressed genes with Entrez IDs: %d\n", nrow(all_expressed)))

# DEGs DOWN
degs_down <- deseq2 %>%
    filter(!is.na(padj), padj < 0.05, log2FoldChange < 0, !is.na(entrez_id))
cat(sprintf("  DEGs DOWN with Entrez IDs: %d\n", nrow(degs_down)))

# Filter gene features for each set
all_genes_ids <- all_expressed$entrez_id
degs_down_ids <- degs_down$entrez_id

# Subset the 5-part features
gf5_all <- gf5[names(gf5) %in% all_genes_ids]
gf5_down <- gf5[names(gf5) %in% degs_down_ids]

cat(sprintf("  Genes in 5-part features (all): %d\n", length(gf5_all)))
cat(sprintf("  Genes in 5-part features (DEGs DOWN): %d\n", length(gf5_down)))

cat("\n")

# =============================================================================
# PHASE 3: SET UP BIGWIG IMPORT PARAMETERS
# =============================================================================

cat("=== PHASE 3: Setting Up Import Parameters ===\n")

# Import parameters for BigWig files
bwImportParams <- setImportParams(
    offset = 0,
    fix_width = 0,
    fix_point = "start",
    norm = TRUE,
    useScore = FALSE,
    outRle = TRUE,
    useSizeFactor = FALSE,
    genome = "hg38"
)

cat("  BigWig import parameters configured\n\n")

# =============================================================================
# PHASE 4: CREATE BINDING PROFILES
# =============================================================================

cat("=== PHASE 4: Creating Binding Profiles ===\n")

# Define binding BigWig files
binding_files <- c(TES_BINDING_BW, TEAD1_BINDING_BW)
names(binding_files) <- c("TES", "TEAD1")

# Check files exist
for (f in binding_files) {
    if (!file.exists(f)) {
        stop(sprintf("BigWig file not found: %s", f))
    }
}
cat("  All binding BigWig files found\n")

# Plot for ALL GENES - Binding
cat("  Generating binding profile for all genes...\n")
tryCatch({
    plot_5parts_metagene(
        queryFiles = binding_files,
        gFeatures_list = list("All_Genes" = gf5_all),
        inputFiles = NULL,
        scale = FALSE,
        verbose = FALSE,
        transform = NA,
        smooth = TRUE,
        stranded = FALSE,
        outPrefix = file.path(OUTPUT_BASE, "binding_all_genes"),
        importParams = bwImportParams,
        heatmap = FALSE,
        rmOutlier = 0.01,
        nc = 4
    )
    cat("  Created: binding_all_genes.pdf\n")
}, error = function(e) {
    cat(sprintf("  Warning: Could not create binding_all_genes plot: %s\n", e$message))
})

# Plot for DEGs DOWN - Binding
cat("  Generating binding profile for DEGs DOWN...\n")
tryCatch({
    plot_5parts_metagene(
        queryFiles = binding_files,
        gFeatures_list = list("DEGs_DOWN" = gf5_down),
        inputFiles = NULL,
        scale = FALSE,
        verbose = FALSE,
        transform = NA,
        smooth = TRUE,
        stranded = FALSE,
        outPrefix = file.path(OUTPUT_BASE, "binding_DEGs_DOWN"),
        importParams = bwImportParams,
        heatmap = FALSE,
        rmOutlier = 0.01,
        nc = 4
    )
    cat("  Created: binding_DEGs_DOWN.pdf\n")
}, error = function(e) {
    cat(sprintf("  Warning: Could not create binding_DEGs_DOWN plot: %s\n", e$message))
})

cat("\n")

# =============================================================================
# PHASE 5: CREATE METHYLATION PROFILES
# =============================================================================

cat("=== PHASE 5: Creating Methylation Profiles ===\n")

# Define methylation BigWig files
meth_files <- c(TES_METH_BW, GFP_METH_BW)
names(meth_files) <- c("TES_meDIP", "GFP_meDIP")

# Check files exist
for (f in meth_files) {
    if (!file.exists(f)) {
        stop(sprintf("BigWig file not found: %s", f))
    }
}
cat("  All methylation BigWig files found\n")

# Plot for ALL GENES - Methylation
cat("  Generating methylation profile for all genes...\n")
tryCatch({
    plot_5parts_metagene(
        queryFiles = meth_files,
        gFeatures_list = list("All_Genes" = gf5_all),
        inputFiles = NULL,
        scale = FALSE,
        verbose = FALSE,
        transform = NA,
        smooth = TRUE,
        stranded = FALSE,
        outPrefix = file.path(OUTPUT_BASE, "methylation_all_genes"),
        importParams = bwImportParams,
        heatmap = FALSE,
        rmOutlier = 0.01,
        nc = 4
    )
    cat("  Created: methylation_all_genes.pdf\n")
}, error = function(e) {
    cat(sprintf("  Warning: Could not create methylation_all_genes plot: %s\n", e$message))
})

# Plot for DEGs DOWN - Methylation
cat("  Generating methylation profile for DEGs DOWN...\n")
tryCatch({
    plot_5parts_metagene(
        queryFiles = meth_files,
        gFeatures_list = list("DEGs_DOWN" = gf5_down),
        inputFiles = NULL,
        scale = FALSE,
        verbose = FALSE,
        transform = NA,
        smooth = TRUE,
        stranded = FALSE,
        outPrefix = file.path(OUTPUT_BASE, "methylation_DEGs_DOWN"),
        importParams = bwImportParams,
        heatmap = FALSE,
        rmOutlier = 0.01,
        nc = 4
    )
    cat("  Created: methylation_DEGs_DOWN.pdf\n")
}, error = function(e) {
    cat(sprintf("  Warning: Could not create methylation_DEGs_DOWN plot: %s\n", e$message))
})

cat("\n")

# =============================================================================
# PHASE 6: CREATE COMBINED COMPARISON PLOTS
# =============================================================================

cat("=== PHASE 6: Creating Combined Comparison Plots ===\n")

# Compare ALL vs DEGs DOWN for binding
cat("  Generating binding comparison (All vs DEGs DOWN)...\n")
tryCatch({
    plot_5parts_metagene(
        queryFiles = binding_files,
        gFeatures_list = list("All_Genes" = gf5_all, "DEGs_DOWN" = gf5_down),
        inputFiles = NULL,
        scale = FALSE,
        verbose = FALSE,
        transform = NA,
        smooth = TRUE,
        stranded = FALSE,
        outPrefix = file.path(OUTPUT_BASE, "binding_comparison"),
        importParams = bwImportParams,
        heatmap = FALSE,
        rmOutlier = 0.01,
        nc = 4
    )
    cat("  Created: binding_comparison.pdf\n")
}, error = function(e) {
    cat(sprintf("  Warning: Could not create binding_comparison plot: %s\n", e$message))
})

# Compare ALL vs DEGs DOWN for methylation
cat("  Generating methylation comparison (All vs DEGs DOWN)...\n")
tryCatch({
    plot_5parts_metagene(
        queryFiles = meth_files,
        gFeatures_list = list("All_Genes" = gf5_all, "DEGs_DOWN" = gf5_down),
        inputFiles = NULL,
        scale = FALSE,
        verbose = FALSE,
        transform = NA,
        smooth = TRUE,
        stranded = FALSE,
        outPrefix = file.path(OUTPUT_BASE, "methylation_comparison"),
        importParams = bwImportParams,
        heatmap = FALSE,
        rmOutlier = 0.01,
        nc = 4
    )
    cat("  Created: methylation_comparison.pdf\n")
}, error = function(e) {
    cat(sprintf("  Warning: Could not create methylation_comparison plot: %s\n", e$message))
})

cat("\n")

# =============================================================================
# COMPLETION
# =============================================================================

cat("==========================================================\n")
cat("GENOMICPLOT METAGENE PROFILES COMPLETE\n")
cat("==========================================================\n")
cat("Completed:", as.character(Sys.time()), "\n")
cat(sprintf("Output directory: %s\n\n", OUTPUT_BASE))
cat("Generated files:\n")
cat("  Binding profiles:\n")
cat("    - binding_all_genes.pdf\n")
cat("    - binding_DEGs_DOWN.pdf\n")
cat("    - binding_comparison.pdf\n")
cat("  Methylation profiles:\n")
cat("    - methylation_all_genes.pdf\n")
cat("    - methylation_DEGs_DOWN.pdf\n")
cat("    - methylation_comparison.pdf\n")
