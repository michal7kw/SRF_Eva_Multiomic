#!/usr/bin/env Rscript
#
# PEAKED GENES METAGENE PROFILES - R PREPARATION SCRIPT
# Creates BED files with gene body coordinates for deepTools scale-regions
# Only includes genes that have Cut&Tag peaks
#
# Version 1: All genes with peaks (no RNA-seq filtering)
# Version 2: DEGs with peaks (UP/DOWN separated)
#

suppressPackageStartupMessages({
    library(GenomicFeatures)
    library(org.Hs.eg.db)
    library(GenomicRanges)
    library(dplyr)
    library(rtracklayer)
})

setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_1")

cat("==========================================================\n")
cat("PEAKED GENES METAGENE - BED FILE PREPARATION\n")
cat("==========================================================\n")
cat("Analysis started:", as.character(Sys.time()), "\n\n")

# =============================================================================
# PATH CONFIGURATION
# =============================================================================

# Peak annotation files
TES_PEAKS_FILE <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_CUTandTAG/results/07_analysis_narrow/TES_peaks_annotated.csv"
TEAD1_PEAKS_FILE <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_CUTandTAG/results/07_analysis_narrow/TEAD1_peaks_annotated.csv"

# DESeq2 results (for Version 2)
DESEQ2_FILE <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_RNA/results/05_deseq2/deseq2_results_TES_vs_GFP.txt"

# GTF annotation
GTF_FILE <- "/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/annotation/gencode.v44.annotation.gtf"

# Output directory
OUTPUT_BASE <- "output/30_peaked_genes_metagene"
dir.create(OUTPUT_BASE, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# PHASE 1: LOAD GENE ANNOTATIONS FROM GTF
# =============================================================================

cat("=== PHASE 1: Loading Gene Annotations ===\n")

# Load GTF and extract gene information
cat("  Loading GTF annotation...\n")
gtf <- rtracklayer::import(GTF_FILE)

# Get gene-level information
genes_gtf <- gtf[gtf$type == "gene"]
genes_df <- as.data.frame(genes_gtf)

# Filter for protein-coding genes on standard chromosomes
genes_df <- genes_df %>%
    filter(gene_type == "protein_coding") %>%
    filter(seqnames %in% paste0("chr", c(1:22, "X", "Y")))

cat(sprintf("  Protein-coding genes on standard chromosomes: %d\n", nrow(genes_df)))

# Clean gene IDs (remove version)
genes_df$gene_id_clean <- gsub("\\..*", "", genes_df$gene_id)

# Create lookup table: gene_name -> gene coordinates
genes_lookup <- genes_df %>%
    select(gene_name, seqnames, start, end, strand, gene_id_clean) %>%
    distinct(gene_name, .keep_all = TRUE)

cat(sprintf("  Unique gene symbols with coordinates: %d\n", nrow(genes_lookup)))

cat("\n")

# =============================================================================
# PHASE 2: LOAD PEAK ANNOTATIONS AND EXTRACT PEAKED GENES
# =============================================================================

cat("=== PHASE 2: Loading Peak Annotations ===\n")

# Load TES peaks
cat("  Loading TES peaks...\n")
tes_peaks <- read.csv(TES_PEAKS_FILE, stringsAsFactors = FALSE)
cat(sprintf("    Total TES peaks: %d\n", nrow(tes_peaks)))

# Load TEAD1 peaks
cat("  Loading TEAD1 peaks...\n")
tead1_peaks <- read.csv(TEAD1_PEAKS_FILE, stringsAsFactors = FALSE)
cat(sprintf("    Total TEAD1 peaks: %d\n", nrow(tead1_peaks)))

# Extract unique gene IDs from peak annotations (Entrez IDs)
tes_entrez_ids <- unique(tes_peaks$geneId[!is.na(tes_peaks$geneId)])
tead1_entrez_ids <- unique(tead1_peaks$geneId[!is.na(tead1_peaks$geneId)])

cat(sprintf("  Unique Entrez IDs with TES peaks: %d\n", length(tes_entrez_ids)))
cat(sprintf("  Unique Entrez IDs with TEAD1 peaks: %d\n", length(tead1_entrez_ids)))

# Convert Entrez IDs to gene symbols
cat("  Converting Entrez IDs to gene symbols...\n")
tes_symbols <- mapIds(org.Hs.eg.db,
                      keys = as.character(tes_entrez_ids),
                      column = "SYMBOL",
                      keytype = "ENTREZID",
                      multiVals = "first")
tes_symbols <- unique(tes_symbols[!is.na(tes_symbols)])

tead1_symbols <- mapIds(org.Hs.eg.db,
                        keys = as.character(tead1_entrez_ids),
                        column = "SYMBOL",
                        keytype = "ENTREZID",
                        multiVals = "first")
tead1_symbols <- unique(tead1_symbols[!is.na(tead1_symbols)])

cat(sprintf("  Gene symbols with TES peaks: %d\n", length(tes_symbols)))
cat(sprintf("  Gene symbols with TEAD1 peaks: %d\n", length(tead1_symbols)))

cat("\n")

# =============================================================================
# PHASE 3: LOAD DESEQ2 FOR DEG FILTERING (VERSION 2)
# =============================================================================

cat("=== PHASE 3: Loading DESeq2 Results ===\n")

deseq2 <- read.delim(DESEQ2_FILE, stringsAsFactors = FALSE)
deseq2$gene_id_clean <- gsub("\\..*", "", deseq2$gene_id)

# DEGs UP
degs_up <- deseq2 %>%
    filter(!is.na(padj), padj < 0.05, log2FoldChange > 0)
cat(sprintf("  DEGs UP: %d\n", nrow(degs_up)))

# DEGs DOWN
degs_down <- deseq2 %>%
    filter(!is.na(padj), padj < 0.05, log2FoldChange < 0)
cat(sprintf("  DEGs DOWN: %d\n", nrow(degs_down)))

# Get gene symbols for DEGs
degs_up_symbols <- unique(degs_up$gene_symbol[!is.na(degs_up$gene_symbol)])
degs_down_symbols <- unique(degs_down$gene_symbol[!is.na(degs_down$gene_symbol)])

cat("\n")

# =============================================================================
# PHASE 4: CREATE GENE SETS
# =============================================================================

cat("=== PHASE 4: Creating Gene Sets ===\n")

# VERSION 1: All genes with peaks (no DEG filter)
tes_peaked_genes <- intersect(tes_symbols, genes_lookup$gene_name)
tead1_peaked_genes <- intersect(tead1_symbols, genes_lookup$gene_name)

cat(sprintf("  VERSION 1 (All peaked genes):\n"))
cat(sprintf("    TES-peaked genes with coordinates: %d\n", length(tes_peaked_genes)))
cat(sprintf("    TEAD1-peaked genes with coordinates: %d\n", length(tead1_peaked_genes)))

# VERSION 2: DEGs with peaks
tes_degs_up_peaked <- intersect(tes_symbols, degs_up_symbols)
tes_degs_up_peaked <- intersect(tes_degs_up_peaked, genes_lookup$gene_name)

tes_degs_down_peaked <- intersect(tes_symbols, degs_down_symbols)
tes_degs_down_peaked <- intersect(tes_degs_down_peaked, genes_lookup$gene_name)

tead1_degs_up_peaked <- intersect(tead1_symbols, degs_up_symbols)
tead1_degs_up_peaked <- intersect(tead1_degs_up_peaked, genes_lookup$gene_name)

tead1_degs_down_peaked <- intersect(tead1_symbols, degs_down_symbols)
tead1_degs_down_peaked <- intersect(tead1_degs_down_peaked, genes_lookup$gene_name)

cat(sprintf("  VERSION 2 (DEGs with peaks):\n"))
cat(sprintf("    TES DEGs UP with peaks: %d\n", length(tes_degs_up_peaked)))
cat(sprintf("    TES DEGs DOWN with peaks: %d\n", length(tes_degs_down_peaked)))
cat(sprintf("    TEAD1 DEGs UP with peaks: %d\n", length(tead1_degs_up_peaked)))
cat(sprintf("    TEAD1 DEGs DOWN with peaks: %d\n", length(tead1_degs_down_peaked)))

cat("\n")

# =============================================================================
# PHASE 5: CREATE BED FILES FOR DEEPTOOLS
# =============================================================================

cat("=== PHASE 5: Creating BED Files ===\n")

# Function to create BED6 file for gene bodies
create_gene_body_bed <- function(gene_symbols, gene_lookup, output_file) {
    # Get coordinates for these genes
    gene_data <- gene_lookup %>%
        filter(gene_name %in% gene_symbols)

    if (nrow(gene_data) == 0) {
        cat(sprintf("  WARNING: No genes found for %s\n", basename(output_file)))
        return(0)
    }

    bed <- data.frame(
        chr = gene_data$seqnames,
        start = gene_data$start - 1,  # BED is 0-based
        end = gene_data$end,
        name = gene_data$gene_name,
        score = 0,
        strand = gene_data$strand,
        stringsAsFactors = FALSE
    )

    # Filter valid entries
    bed <- bed[!is.na(bed$start) & !is.na(bed$end), ]
    bed <- bed[bed$start >= 0, ]

    # Sort by chromosome and position
    bed <- bed[order(bed$chr, bed$start), ]

    # Remove duplicates
    bed <- bed[!duplicated(paste(bed$chr, bed$start, bed$end)), ]

    write.table(bed, output_file, quote = FALSE, sep = "\t",
                row.names = FALSE, col.names = FALSE)

    cat(sprintf("  Created: %s (%d genes)\n", basename(output_file), nrow(bed)))
    return(nrow(bed))
}

# VERSION 1: All peaked genes
cat("\n  VERSION 1 - All peaked genes:\n")
n_tes <- create_gene_body_bed(tes_peaked_genes, genes_lookup,
                               file.path(OUTPUT_BASE, "TES_peaked_genes.bed"))
n_tead1 <- create_gene_body_bed(tead1_peaked_genes, genes_lookup,
                                 file.path(OUTPUT_BASE, "TEAD1_peaked_genes.bed"))

# VERSION 2: DEGs with peaks
cat("\n  VERSION 2 - DEGs with peaks:\n")
n_tes_up <- create_gene_body_bed(tes_degs_up_peaked, genes_lookup,
                                  file.path(OUTPUT_BASE, "TES_DEGs_UP_peaked.bed"))
n_tes_down <- create_gene_body_bed(tes_degs_down_peaked, genes_lookup,
                                    file.path(OUTPUT_BASE, "TES_DEGs_DOWN_peaked.bed"))
n_tead1_up <- create_gene_body_bed(tead1_degs_up_peaked, genes_lookup,
                                    file.path(OUTPUT_BASE, "TEAD1_DEGs_UP_peaked.bed"))
n_tead1_down <- create_gene_body_bed(tead1_degs_down_peaked, genes_lookup,
                                      file.path(OUTPUT_BASE, "TEAD1_DEGs_DOWN_peaked.bed"))

cat("\n")

# =============================================================================
# SUMMARY
# =============================================================================

cat("==========================================================\n")
cat("BED FILE PREPARATION COMPLETE\n")
cat("==========================================================\n")
cat("Completed:", as.character(Sys.time()), "\n")
cat(sprintf("Output directory: %s\n\n", OUTPUT_BASE))

cat("VERSION 1 - All peaked genes:\n")
cat(sprintf("  - TES_peaked_genes.bed (%d genes)\n", n_tes))
cat(sprintf("  - TEAD1_peaked_genes.bed (%d genes)\n", n_tead1))

cat("\nVERSION 2 - DEGs with peaks:\n")
cat(sprintf("  - TES_DEGs_UP_peaked.bed (%d genes)\n", n_tes_up))
cat(sprintf("  - TES_DEGs_DOWN_peaked.bed (%d genes)\n", n_tes_down))
cat(sprintf("  - TEAD1_DEGs_UP_peaked.bed (%d genes)\n", n_tead1_up))
cat(sprintf("  - TEAD1_DEGs_DOWN_peaked.bed (%d genes)\n", n_tead1_down))

cat("\nNext: Run deepTools to generate metagene profiles\n")
