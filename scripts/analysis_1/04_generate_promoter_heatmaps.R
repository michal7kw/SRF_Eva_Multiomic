#!/usr/bin/env Rscript
#
# GENERATE PROMOTER HEATMAPS: Create BED files for deepTools visualization
# Creates BED files for TES and TEAD1 UP/DOWN regulated gene promoters

suppressPackageStartupMessages({
    library(dplyr)
    library(GenomicRanges)
    library(rtracklayer)
})

setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_1")

cat("=== GENERATE PROMOTER HEATMAPS ===\n")
cat("Creating BED files for deepTools visualization\n")
cat("Analysis started:", as.character(Sys.time()), "\n\n")

# =============================================================================
# PATH CONFIGURATION
# =============================================================================
# Input paths from integrative analysis (script 11 - all genes version)
INPUT_BASE <- "output/11_final_integrative_analysis_all_genes"
INPUT_DIRECT_TARGETS <- file.path(INPUT_BASE, "direct_targets")

# Output paths
OUTPUT_BASE <- "output/04_generate_promoter_heatmaps"

# Reference annotation
GTF_FILE <- "/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/annotation/gencode.v44.annotation.gtf"

# Create output directory
dir.create(OUTPUT_BASE, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# PHASE 1: LOAD GENE ANNOTATIONS
# =============================================================================

cat("=== PHASE 1: Loading Gene Annotations ===\n")

# Load GTF annotation
cat(sprintf("Loading GTF: %s\n", GTF_FILE))
gtf <- rtracklayer::import(GTF_FILE)

# Convert to data frame and filter for genes only
genes_df <- as.data.frame(gtf) %>%
    filter(type == "gene") %>%
    select(seqnames, start, end, strand, gene_id, gene_name)

# Clean gene IDs (remove version)
genes_df$gene_id_clean <- gsub("\\..*", "", genes_df$gene_id)

cat(sprintf("✓ Loaded %d gene annotations\n\n", nrow(genes_df)))

# =============================================================================
# PHASE 2: LOAD INTEGRATIVE ANALYSIS RESULTS
# =============================================================================

cat("=== PHASE 2: Loading Integrative Analysis Results ===\n")
cat(sprintf("Input directory: %s\n", INPUT_DIRECT_TARGETS))

# Load direct target gene lists
tes_direct <- read.csv(file.path(INPUT_DIRECT_TARGETS, "TES_direct_targets_all_genes.csv"))
tead1_direct <- read.csv(file.path(INPUT_DIRECT_TARGETS, "TEAD1_direct_targets_all_genes.csv"))

cat(sprintf("✓ Loaded %d TES direct targets\n", nrow(tes_direct)))
cat(sprintf("✓ Loaded %d TEAD1 direct targets\n\n", nrow(tead1_direct)))

# Split by direction
tes_up <- tes_direct[tes_direct$log2FoldChange > 0, ]
tes_down <- tes_direct[tes_direct$log2FoldChange < 0, ]
tead1_up <- tead1_direct[tead1_direct$log2FoldChange > 0, ]
tead1_down <- tead1_direct[tead1_direct$log2FoldChange < 0, ]

cat(sprintf("TES UP: %d genes, TES DOWN: %d genes\n", nrow(tes_up), nrow(tes_down)))
cat(sprintf("TEAD1 UP: %d genes, TEAD1 DOWN: %d genes\n\n", nrow(tead1_up), nrow(tead1_down)))

# =============================================================================
# PHASE 3: CREATE BED FILES
# =============================================================================

cat("=== PHASE 3: Creating BED Files ===\n")

# Function to create BED file
create_bed <- function(gene_data, genes_annotation, output_file) {
    # Merge with annotation
    merged <- merge(gene_data, genes_annotation,
        by.x = "ensembl_id", by.y = "gene_id_clean",
        all.x = TRUE
    )

    # Remove genes without coordinates
    merged <- merged[!is.na(merged$start), ]

    if (nrow(merged) == 0) {
        cat(sprintf("Warning: No genes matched for %s\n", basename(output_file)))
        return(0)
    }

    # Define promoter regions (TSS +/- 2kb)
    merged$promoter_start <- ifelse(merged$strand == "+",
        merged$start - 2000,
        merged$end - 2000
    )
    merged$promoter_end <- ifelse(merged$strand == "+",
        merged$start + 2000,
        merged$end + 2000
    )

    # Ensure valid coordinates
    merged$promoter_start[merged$promoter_start < 0] <- 0

    # Create BED format (chr, start, end, name, score, strand)
    bed <- data.frame(
        chr = merged$seqnames,
        start = merged$promoter_start,
        end = merged$promoter_end,
        name = merged$gene_name,
        score = abs(merged$log2FoldChange) * 100, # Scale for visibility
        strand = merged$strand
    )

    # Sort by chromosome and start
    bed <- bed[order(bed$chr, bed$start), ]

    # Write BED file
    write.table(bed, output_file,
        quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE
    )

    cat(sprintf("✓ Created %s: %d regions\n", basename(output_file), nrow(bed)))
    return(nrow(bed))
}

# Create BED files for each category
n_tes_up <- create_bed(tes_up, genes_df, file.path(OUTPUT_BASE, "TES_UP_promoters.bed"))
n_tes_down <- create_bed(tes_down, genes_df, file.path(OUTPUT_BASE, "TES_DOWN_promoters.bed"))
n_tead1_up <- create_bed(tead1_up, genes_df, file.path(OUTPUT_BASE, "TEAD1_UP_promoters.bed"))
n_tead1_down <- create_bed(tead1_down, genes_df, file.path(OUTPUT_BASE, "TEAD1_DOWN_promoters.bed"))

cat("\n========================================\n")
cat("BED FILE CREATION COMPLETE\n")
cat("========================================\n")
cat("Completed:", as.character(Sys.time()), "\n")
cat(sprintf("Input directory: %s\n", INPUT_DIRECT_TARGETS))
cat(sprintf("Output directory: %s\n\n", OUTPUT_BASE))
cat("Created files:\n")
cat("  - TES_UP_promoters.bed\n")
cat("  - TES_DOWN_promoters.bed\n")
cat("  - TEAD1_UP_promoters.bed\n")
cat("  - TEAD1_DOWN_promoters.bed\n")
