#!/usr/bin/env Rscript
#
# METAGENE PROFILES: TES/TEAD1 Binding + Methylation
# Creates BED files for:
#   1. All expressed genes
#   2. Downregulated DEGs (DEGs DOWN)
#
# For use with deepTools to generate binding + methylation metagene profiles

suppressPackageStartupMessages({
    library(dplyr)
    library(GenomicRanges)
    library(rtracklayer)
})

setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_1")

cat("=== METAGENE BINDING + METHYLATION PROFILES ===\n")
cat("Creating BED files for deepTools visualization\n")
cat("Analysis started:", as.character(Sys.time()), "\n\n")

# =============================================================================
# PATH CONFIGURATION
# =============================================================================

# Input: DESeq2 results
DESEQ2_FILE <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_RNA/results/05_deseq2/deseq2_results_TES_vs_GFP.txt"

# Reference annotation
GTF_FILE <- "/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/annotation/gencode.v44.annotation.gtf"

# Output
OUTPUT_BASE <- "output/21_metagene_binding_methylation"
dir.create(OUTPUT_BASE, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# PHASE 1: LOAD GENE ANNOTATIONS
# =============================================================================

cat("=== PHASE 1: Loading Gene Annotations ===\n")

cat(sprintf("Loading GTF: %s\n", GTF_FILE))
gtf <- rtracklayer::import(GTF_FILE)

# Filter for genes only and get coordinates
genes_df <- as.data.frame(gtf) %>%
    filter(type == "gene") %>%
    select(seqnames, start, end, strand, gene_id, gene_name)

# Clean gene IDs (remove version)
genes_df$gene_id_clean <- gsub("\\..*", "", genes_df$gene_id)

cat(sprintf("✓ Loaded %d gene annotations\n\n", nrow(genes_df)))

# =============================================================================
# PHASE 2: LOAD DESEQ2 RESULTS
# =============================================================================

cat("=== PHASE 2: Loading DESeq2 Results ===\n")

deseq2 <- read.delim(DESEQ2_FILE, header = TRUE, stringsAsFactors = FALSE)
cat(sprintf("✓ Loaded %d genes from DESeq2\n", nrow(deseq2)))

# Clean gene IDs
deseq2$gene_id_clean <- gsub("\\..*", "", deseq2$gene_id)

# Define gene sets
# 1. All expressed genes (genes with non-NA padj, meaning they passed filtering)
all_expressed <- deseq2 %>%
    filter(!is.na(padj)) %>%
    select(gene_id_clean, log2FoldChange, padj, gene_symbol)

cat(sprintf("✓ All expressed genes: %d\n", nrow(all_expressed)))

# 2. DEGs DOWN (significantly downregulated: padj < 0.05 & log2FC < 0)
degs_down <- deseq2 %>%
    filter(!is.na(padj), padj < 0.05, log2FoldChange < 0) %>%
    select(gene_id_clean, log2FoldChange, padj, gene_symbol)

cat(sprintf("✓ DEGs DOWN (padj<0.05, log2FC<0): %d\n\n", nrow(degs_down)))

# =============================================================================
# PHASE 3: CREATE BED FILES
# =============================================================================

cat("=== PHASE 3: Creating BED Files ===\n")

# Function to create TSS-centered BED file
create_tss_bed <- function(gene_data, genes_annotation, output_file, description) {
    # Merge with annotation
    merged <- merge(gene_data, genes_annotation,
        by.x = "gene_id_clean", by.y = "gene_id_clean",
        all.x = TRUE
    )

    # Remove genes without coordinates
    merged <- merged[!is.na(merged$start), ]

    if (nrow(merged) == 0) {
        cat(sprintf("Warning: No genes matched for %s\n", basename(output_file)))
        return(0)
    }

    # Calculate TSS position
    merged$tss <- ifelse(merged$strand == "+", merged$start, merged$end)

    # Create BED format for TSS (single bp for reference-point mode)
    # We'll use a small region around TSS for the BED file
    bed <- data.frame(
        chr = merged$seqnames,
        start = merged$tss - 1,  # BED is 0-based
        end = merged$tss,
        name = ifelse(is.na(merged$gene_symbol) | merged$gene_symbol == "",
                      merged$gene_id_clean, merged$gene_symbol),
        score = 0,
        strand = merged$strand
    )

    # Filter for standard chromosomes
    valid_chr <- paste0("chr", c(1:22, "X", "Y"))
    bed <- bed[bed$chr %in% valid_chr, ]

    # Remove duplicates (in case of multiple gene entries)
    bed <- bed[!duplicated(paste(bed$chr, bed$start, bed$end, bed$strand)), ]

    # Sort by chromosome and position
    bed <- bed[order(bed$chr, bed$start), ]

    # Write BED file
    write.table(bed, output_file,
        quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE
    )

    cat(sprintf("✓ Created %s: %d TSS regions (%s)\n",
                basename(output_file), nrow(bed), description))
    return(nrow(bed))
}

# Create BED files
n_all <- create_tss_bed(all_expressed, genes_df,
                         file.path(OUTPUT_BASE, "all_expressed_genes_TSS.bed"),
                         "all expressed genes")

n_down <- create_tss_bed(degs_down, genes_df,
                          file.path(OUTPUT_BASE, "DEGs_DOWN_TSS.bed"),
                          "downregulated DEGs")

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n========================================\n")
cat("BED FILE CREATION COMPLETE\n")
cat("========================================\n")
cat("Completed:", as.character(Sys.time()), "\n")
cat(sprintf("Output directory: %s\n\n", OUTPUT_BASE))
cat("Created files:\n")
cat(sprintf("  - all_expressed_genes_TSS.bed (%d regions)\n", n_all))
cat(sprintf("  - DEGs_DOWN_TSS.bed (%d regions)\n", n_down))
cat("\nNext step: Run the shell script to compute matrices and generate plots\n")
