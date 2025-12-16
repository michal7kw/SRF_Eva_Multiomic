#!/usr/bin/env Rscript
#
# GENE BODY METHYLATION AT DMR-ASSOCIATED GENES
# Creates BED files for genes stratified by DMR status:
#   - Genes with hypermethylated DMRs in gene body
#   - Genes with hypomethylated DMRs in gene body
#   - Genes without DMRs (control)
#

suppressPackageStartupMessages({
    library(GenomicRanges)
    library(rtracklayer)
    library(dplyr)
})

setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_1")

cat("==========================================================\n")
cat("GENE BODY METHYLATION - DMR-ASSOCIATED GENES\n")
cat("==========================================================\n")
cat("Analysis started:", as.character(Sys.time()), "\n\n")

# =============================================================================
# PATH CONFIGURATION
# =============================================================================

# DESeq2 results
DESEQ2_FILE <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_RNA/results/05_deseq2/deseq2_results_TES_vs_GFP.txt"

# GTF annotation
GTF_FILE <- "/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/annotation/gencode.v44.annotation.gtf"

# DMR files
DMR_FILE <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP/results/07_differential_MEDIPS/TES_vs_GFP_DMRs_FDR05_FC2.csv"

# Output directory
OUTPUT_BASE <- "output/25_genebody_methylation_DMR_genes"
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

# Create GRanges for genes
genes_gr <- GRanges(
    seqnames = genes_df$seqnames,
    ranges = IRanges(start = genes_df$start, end = genes_df$end),
    strand = genes_df$strand,
    gene_name = genes_df$gene_name,
    gene_id = genes_df$gene_id_clean
)

cat("\n")

# =============================================================================
# PHASE 2: LOAD DMRs
# =============================================================================

cat("=== PHASE 2: Loading DMRs ===\n")

dmr_data <- read.csv(DMR_FILE, stringsAsFactors = FALSE)
cat(sprintf("  Total DMRs: %d\n", nrow(dmr_data)))

# Separate hyper and hypo methylated
# Note: Column is 'logFC' and 'stop' instead of 'end'
dmr_hyper <- dmr_data %>% filter(logFC > 0)  # TES > GFP (hypermethylated)
dmr_hypo <- dmr_data %>% filter(logFC < 0)   # TES < GFP (hypomethylated)

cat(sprintf("  Hypermethylated DMRs: %d\n", nrow(dmr_hyper)))
cat(sprintf("  Hypomethylated DMRs: %d\n", nrow(dmr_hypo)))

# Create GRanges for DMRs
# Note: DMR file uses "1" format, genes use "chr1" format - add prefix
create_dmr_gr <- function(dmr_df) {
    # Add "chr" prefix if not present
    chr_names <- dmr_df$chr
    if (!grepl("^chr", chr_names[1])) {
        chr_names <- paste0("chr", chr_names)
    }

    GRanges(
        seqnames = chr_names,
        ranges = IRanges(start = dmr_df$start, end = dmr_df$stop),
        logFC = dmr_df$logFC,
        FDR = dmr_df$FDR
    )
}

dmr_hyper_gr <- create_dmr_gr(dmr_hyper)
dmr_hypo_gr <- create_dmr_gr(dmr_hypo)

cat("\n")

# =============================================================================
# PHASE 3: FIND GENES WITH DMRs IN GENE BODY
# =============================================================================

cat("=== PHASE 3: Finding Genes with DMRs in Gene Body ===\n")

# Find overlaps between genes and DMRs
# Using gene body (not just promoter)
hyper_overlaps <- findOverlaps(genes_gr, dmr_hyper_gr)
hypo_overlaps <- findOverlaps(genes_gr, dmr_hypo_gr)

# Get unique genes with each type of DMR
genes_with_hyper_dmr <- unique(queryHits(hyper_overlaps))
genes_with_hypo_dmr <- unique(queryHits(hypo_overlaps))

cat(sprintf("  Genes with hypermethylated DMRs in gene body: %d\n", length(genes_with_hyper_dmr)))
cat(sprintf("  Genes with hypomethylated DMRs in gene body: %d\n", length(genes_with_hypo_dmr)))

# Genes with only hypermethylated DMRs (not hypomethylated)
genes_hyper_only <- setdiff(genes_with_hyper_dmr, genes_with_hypo_dmr)
cat(sprintf("  Genes with ONLY hypermethylated DMRs: %d\n", length(genes_hyper_only)))

# Genes with only hypomethylated DMRs (not hypermethylated)
genes_hypo_only <- setdiff(genes_with_hypo_dmr, genes_with_hyper_dmr)
cat(sprintf("  Genes with ONLY hypomethylated DMRs: %d\n", length(genes_hypo_only)))

# Genes with no DMRs at all
all_dmr_overlaps <- findOverlaps(genes_gr, c(dmr_hyper_gr, dmr_hypo_gr))
genes_with_any_dmr <- unique(queryHits(all_dmr_overlaps))
genes_no_dmr <- setdiff(1:length(genes_gr), genes_with_any_dmr)
cat(sprintf("  Genes without any DMRs: %d\n", length(genes_no_dmr)))

cat("\n")

# =============================================================================
# PHASE 4: LOAD DESEQ2 AND FILTER FOR EXPRESSED GENES
# =============================================================================

cat("=== PHASE 4: Filtering for Expressed Genes ===\n")

deseq2 <- read.delim(DESEQ2_FILE, stringsAsFactors = FALSE)
deseq2$gene_id_clean <- gsub("\\..*", "", deseq2$gene_id)

# Get expressed genes
expressed_genes <- deseq2 %>%
    filter(!is.na(padj)) %>%
    pull(gene_id_clean)

cat(sprintf("  Expressed genes (non-NA padj): %d\n", length(expressed_genes)))

# Filter gene sets to only expressed genes
gene_ids <- genes_gr$gene_id

# Hypermethylated genes (expressed)
hyper_gene_ids <- gene_ids[genes_hyper_only]
hyper_expressed <- hyper_gene_ids[hyper_gene_ids %in% expressed_genes]
cat(sprintf("  Hypermethylated genes (expressed): %d\n", length(hyper_expressed)))

# Hypomethylated genes (expressed)
hypo_gene_ids <- gene_ids[genes_hypo_only]
hypo_expressed <- hypo_gene_ids[hypo_gene_ids %in% expressed_genes]
cat(sprintf("  Hypomethylated genes (expressed): %d\n", length(hypo_expressed)))

# No DMR genes (expressed) - sample to match size
no_dmr_gene_ids <- gene_ids[genes_no_dmr]
no_dmr_expressed <- no_dmr_gene_ids[no_dmr_gene_ids %in% expressed_genes]
# Sample to reasonable size for comparison
set.seed(42)
if (length(no_dmr_expressed) > 5000) {
    no_dmr_expressed <- sample(no_dmr_expressed, 5000)
}
cat(sprintf("  No-DMR genes (expressed, sampled): %d\n", length(no_dmr_expressed)))

cat("\n")

# =============================================================================
# PHASE 5: CREATE BED FILES FOR DEEPTOOLS
# =============================================================================

cat("=== PHASE 5: Creating BED Files ===\n")

# Function to create BED6 file for gene bodies
create_gene_body_bed <- function(gene_ids_to_include, genes_df, output_file) {
    selected <- genes_df[genes_df$gene_id_clean %in% gene_ids_to_include, ]

    bed <- data.frame(
        chr = selected$seqnames,
        start = selected$start - 1,  # BED is 0-based
        end = selected$end,
        name = selected$gene_name,
        score = 0,
        strand = selected$strand,
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

# Create BED files for each gene set
n_hyper <- create_gene_body_bed(hyper_expressed, genes_df,
                                 file.path(OUTPUT_BASE, "genes_with_hypermethylated_DMRs.bed"))
n_hypo <- create_gene_body_bed(hypo_expressed, genes_df,
                                file.path(OUTPUT_BASE, "genes_with_hypomethylated_DMRs.bed"))
n_no_dmr <- create_gene_body_bed(no_dmr_expressed, genes_df,
                                  file.path(OUTPUT_BASE, "genes_without_DMRs.bed"))

# Also create DEGs DOWN with hypermethylated DMRs
degs_down <- deseq2 %>%
    filter(!is.na(padj), padj < 0.05, log2FoldChange < 0) %>%
    pull(gene_id_clean)

degs_down_hyper <- intersect(hyper_expressed, degs_down)
cat(sprintf("\n  DEGs DOWN with hypermethylated DMRs: %d\n", length(degs_down_hyper)))

if (length(degs_down_hyper) >= 50) {
    n_degs_hyper <- create_gene_body_bed(degs_down_hyper, genes_df,
                                          file.path(OUTPUT_BASE, "DEGs_DOWN_hypermethylated.bed"))
}

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
cat(sprintf("  - genes_with_hypermethylated_DMRs.bed (%d genes)\n", n_hyper))
cat(sprintf("  - genes_with_hypomethylated_DMRs.bed (%d genes)\n", n_hypo))
cat(sprintf("  - genes_without_DMRs.bed (%d genes)\n", n_no_dmr))
if (exists("n_degs_hyper")) {
    cat(sprintf("  - DEGs_DOWN_hypermethylated.bed (%d genes)\n", n_degs_hyper))
}
cat("\nNext: Run deepTools to generate gene body methylation profiles\n")
