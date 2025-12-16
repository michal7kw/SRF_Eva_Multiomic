#!/usr/bin/env Rscript
#
# GENE BODY METAGENE PROFILES - R PREPARATION SCRIPT
# Creates BED files with gene body coordinates for deepTools scale-regions
# Includes: TSS, gene body (scaled), TES
#

suppressPackageStartupMessages({
    library(GenomicFeatures)
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    library(org.Hs.eg.db)
    library(GenomicRanges)
    library(dplyr)
    library(rtracklayer)
})

setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_1")

cat("==========================================================\n")
cat("GENE BODY METAGENE - BED FILE PREPARATION\n")
cat("==========================================================\n")
cat("Analysis started:", as.character(Sys.time()), "\n\n")

# =============================================================================
# PATH CONFIGURATION
# =============================================================================

# DESeq2 results
DESEQ2_FILE <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_RNA/results/05_deseq2/deseq2_results_TES_vs_GFP.txt"

# GTF annotation
GTF_FILE <- "/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/annotation/gencode.v44.annotation.gtf"

# Output directory
OUTPUT_BASE <- "output/23_genebody_metagene"
dir.create(OUTPUT_BASE, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# PHASE 1: LOAD GENE ANNOTATIONS
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

# Clean gene IDs
genes_df$gene_id_clean <- gsub("\\..*", "", genes_df$gene_id)

cat("\n")

# =============================================================================
# PHASE 2: LOAD DESEQ2 AND DEFINE GENE SETS
# =============================================================================

cat("=== PHASE 2: Defining Gene Sets ===\n")

# Load DESeq2 results
deseq2 <- read.delim(DESEQ2_FILE, stringsAsFactors = FALSE)
deseq2$gene_id_clean <- gsub("\\..*", "", deseq2$gene_id)

# All expressed genes
all_expressed <- deseq2 %>%
    filter(!is.na(padj))
cat(sprintf("  All expressed genes: %d\n", nrow(all_expressed)))

# DEGs DOWN
degs_down <- deseq2 %>%
    filter(!is.na(padj), padj < 0.05, log2FoldChange < 0)
cat(sprintf("  DEGs DOWN: %d\n", nrow(degs_down)))

# Merge with gene coordinates
all_genes_coords <- genes_df %>%
    filter(gene_id_clean %in% all_expressed$gene_id_clean)
cat(sprintf("  All expressed with coordinates: %d\n", nrow(all_genes_coords)))

degs_down_coords <- genes_df %>%
    filter(gene_id_clean %in% degs_down$gene_id_clean)
cat(sprintf("  DEGs DOWN with coordinates: %d\n", nrow(degs_down_coords)))

cat("\n")

# =============================================================================
# PHASE 3: CREATE BED FILES FOR DEEPTOOLS
# =============================================================================

cat("=== PHASE 3: Creating BED Files ===\n")

# Function to create BED6 file for gene bodies
create_gene_body_bed <- function(gene_data, output_file) {
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

# Create BED files
n_all <- create_gene_body_bed(all_genes_coords,
                               file.path(OUTPUT_BASE, "all_expressed_genes.bed"))
n_down <- create_gene_body_bed(degs_down_coords,
                                file.path(OUTPUT_BASE, "DEGs_DOWN_genes.bed"))

cat("\n")

# =============================================================================
# SUMMARY
# =============================================================================

cat("==========================================================\n")
cat("BED FILE PREPARATION COMPLETE\n")
cat("==========================================================\n")
cat("Completed:", as.character(Sys.time()), "\n")
cat(sprintf("Output directory: %s\n\n", OUTPUT_BASE))
cat("Created files:\n")
cat(sprintf("  - all_expressed_genes.bed (%d genes)\n", n_all))
cat(sprintf("  - DEGs_DOWN_genes.bed (%d genes)\n", n_down))
cat("\nNext: Run deepTools to generate metagene profiles\n")
