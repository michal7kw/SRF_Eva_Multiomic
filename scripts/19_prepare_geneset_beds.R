#!/usr/bin/env Rscript

# ============================================================================
# Prepare Gene Set BED Files for Metagene Analysis
# ============================================================================
# Purpose: Convert gene lists to BED format for deepTools metagene analysis
# Creates BED files for:
#   1. Downregulated DEGs (from meDIP/results/12_gene_sets/)
#   2. TES_degs.txt custom gene list (gene symbols → genomic coordinates)
#
# Output: Gene body coordinates (TSS to TES) for scale-regions mode
#         Promoter coordinates (TSS ± 2kb) for reference-point mode
# ============================================================================

suppressPackageStartupMessages({
    library(rtracklayer)
    library(GenomicRanges)
    library(GenomicFeatures)
    library(org.Hs.eg.db)
    library(dplyr)
    library(stringr)
})

# Set working directory
setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis")

# Create output directory
dir.create("results/19_metagene_beds", recursive = TRUE, showWarnings = FALSE)

cat("============================================\n")
cat("Prepare Gene Set BED Files\n")
cat("Started:", format(Sys.time()), "\n")
cat("============================================\n\n")

# ============================================================================
# Load annotation resources
# ============================================================================

cat("Loading annotation resources...\n")

# Load GTF for gene coordinates
gtf_file <- "/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/annotation/gencode.v44.annotation.gtf"
if (!file.exists(gtf_file)) {
    stop("GTF file not found: ", gtf_file)
}

# Create TxDb from GTF
txdb <- makeTxDbFromGFF(gtf_file, format = "gtf")
genes_gr <- genes(txdb)

# Get gene symbol mappings
# IMPORTANT: Strip version numbers from Ensembl IDs (ENSG00000000003.15 -> ENSG00000000003)
# org.Hs.eg.db requires unversioned Ensembl IDs
ensembl_ids_clean <- gsub("\\..*", "", names(genes_gr))

gene_symbols <- mapIds(org.Hs.eg.db,
    keys = ensembl_ids_clean,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
)

# Add gene symbols to GRanges
# Note: gene_symbols is indexed by clean IDs, so we need to use clean IDs for lookup
genes_gr$gene_symbol <- gene_symbols[ensembl_ids_clean]
genes_gr$ensembl_id <- names(genes_gr)
genes_gr$ensembl_id_clean <- ensembl_ids_clean

cat(sprintf("  Loaded %d genes from GTF\n", length(genes_gr)))
cat(sprintf("  Mapped %d genes to symbols\n", sum(!is.na(genes_gr$gene_symbol))))

# ============================================================================
# Function to create BED files from gene list
# ============================================================================

create_gene_beds <- function(gene_list, gene_type = "symbol", name_prefix, genes_ref = genes_gr) {
    cat(sprintf("\nProcessing gene list: %s (%d genes, type: %s)\n",
                name_prefix, length(gene_list), gene_type))

    # NOTE: BigWig files use Ensembl naming (1, 2, 3...) not UCSC (chr1, chr2, chr3...)
    # We need to convert chromosome names to Ensembl format for deepTools compatibility

    # Match genes
    if (gene_type == "symbol") {
        # Clean gene symbols (remove whitespace, uppercase)
        gene_list <- toupper(trimws(gene_list))
        gene_list <- gene_list[gene_list != ""]

        # Match by symbol
        matched_idx <- which(toupper(genes_ref$gene_symbol) %in% gene_list)
        matched_genes <- genes_ref[matched_idx]

        # Report matching
        matched_symbols <- unique(genes_ref$gene_symbol[matched_idx])
        missing <- setdiff(gene_list, toupper(matched_symbols))
        cat(sprintf("  Matched: %d/%d genes\n", length(matched_symbols), length(gene_list)))
        if (length(missing) > 0 && length(missing) <= 20) {
            cat(sprintf("  Missing: %s\n", paste(missing, collapse = ", ")))
        } else if (length(missing) > 20) {
            cat(sprintf("  Missing: %d genes (showing first 10: %s...)\n",
                        length(missing), paste(head(missing, 10), collapse = ", ")))
        }

    } else if (gene_type == "ensembl") {
        # Match by Ensembl ID (remove version if present)
        gene_list_clean <- gsub("\\..*", "", gene_list)

        # Use pre-computed clean IDs if available, otherwise compute them
        if ("ensembl_id_clean" %in% names(mcols(genes_ref))) {
            genes_ref_clean <- genes_ref$ensembl_id_clean
        } else {
            genes_ref_clean <- gsub("\\..*", "", names(genes_ref))
        }

        matched_idx <- which(genes_ref_clean %in% gene_list_clean)
        matched_genes <- genes_ref[matched_idx]

        cat(sprintf("  Matched: %d/%d genes\n", length(matched_idx), length(gene_list)))
    } else {
        stop("Unknown gene_type: ", gene_type)
    }

    if (length(matched_genes) < 5) {
        cat("  WARNING: Too few genes matched. Skipping.\n")
        return(NULL)
    }

    # Create gene body BED (TSS to TES)
    # IMPORTANT: Convert UCSC chromosome names (chr1, chr2) to Ensembl format (1, 2)
    # This is required because BigWig files use Ensembl naming
    chr_names <- as.character(seqnames(matched_genes))
    chr_names <- gsub("^chr", "", chr_names)  # Remove "chr" prefix for Ensembl compatibility

    gene_body <- data.frame(
        chr = chr_names,
        start = start(matched_genes) - 1,  # BED is 0-based
        end = end(matched_genes),
        name = ifelse(is.na(matched_genes$gene_symbol),
                      matched_genes$ensembl_id, matched_genes$gene_symbol),
        score = 0,
        strand = as.character(strand(matched_genes))
    )

    # Create promoter BED (TSS ± 2kb)
    tss <- ifelse(gene_body$strand == "+", gene_body$start, gene_body$end)
    promoter <- data.frame(
        chr = gene_body$chr,  # Already converted to Ensembl naming
        start = pmax(0, tss - 2000),
        end = tss + 2000,
        name = gene_body$name,
        score = 0,
        strand = gene_body$strand
    )

    # Write BED files
    gene_body_file <- sprintf("results/19_metagene_beds/%s_gene_body.bed", name_prefix)
    promoter_file <- sprintf("results/19_metagene_beds/%s_promoter.bed", name_prefix)

    write.table(gene_body, gene_body_file, sep = "\t", quote = FALSE,
                row.names = FALSE, col.names = FALSE)
    write.table(promoter, promoter_file, sep = "\t", quote = FALSE,
                row.names = FALSE, col.names = FALSE)

    cat(sprintf("  Created: %s (%d regions)\n", gene_body_file, nrow(gene_body)))
    cat(sprintf("  Created: %s (%d regions)\n", promoter_file, nrow(promoter)))

    return(list(
        gene_body = gene_body_file,
        promoter = promoter_file,
        n_genes = nrow(gene_body),
        matched_genes = matched_genes
    ))
}

# ============================================================================
# Process TES_degs.txt gene list
# ============================================================================

cat("\n=== Processing TES_degs.txt ===\n")

tes_degs_file <- "data/TES_degs.txt"
if (file.exists(tes_degs_file)) {
    tes_degs <- readLines(tes_degs_file)
    tes_degs <- tes_degs[tes_degs != ""]  # Remove empty lines

    cat(sprintf("Read %d genes from TES_degs.txt\n", length(tes_degs)))

    result_tes <- create_gene_beds(
        gene_list = tes_degs,
        gene_type = "symbol",
        name_prefix = "TES_degs"
    )
} else {
    cat("WARNING: TES_degs.txt not found at", tes_degs_file, "\n")
}

# ============================================================================
# Process downregulated DEGs from meDIP gene sets
# ============================================================================

cat("\n=== Processing Downregulated DEGs ===\n")

downreg_promoter_file <- "../meDIP/results/12_gene_sets/downregulated_promoters.bed"
if (file.exists(downreg_promoter_file)) {
    # Read existing BED file
    downreg_bed <- read.table(downreg_promoter_file, sep = "\t", header = FALSE,
                               stringsAsFactors = FALSE)

    # Extract gene IDs (column 4 contains gene names/IDs)
    downreg_genes <- downreg_bed[, 4]

    # Check if these are Ensembl IDs or symbols
    if (grepl("^ENSG", downreg_genes[1])) {
        gene_type <- "ensembl"
    } else {
        gene_type <- "symbol"
    }

    cat(sprintf("Read %d downregulated genes from meDIP BED file\n", length(downreg_genes)))

    result_downreg <- create_gene_beds(
        gene_list = downreg_genes,
        gene_type = gene_type,
        name_prefix = "downregulated_degs"
    )
} else {
    cat("WARNING: Downregulated promoters file not found at", downreg_promoter_file, "\n")

    # Alternative: try to get from DESeq2 results directly
    deseq_file <- "../SRF_Eva_RNA/results/05_deseq2/deseq2_results_TES_vs_GFP.txt"
    if (file.exists(deseq_file)) {
        cat("  Using DESeq2 results as alternative source...\n")
        deseq <- read.table(deseq_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

        # Filter for downregulated genes (padj < 0.05, log2FC < -1)
        downreg <- deseq[!is.na(deseq$padj) & deseq$padj < 0.05 &
                         !is.na(deseq$log2FoldChange) & deseq$log2FoldChange < -1, ]

        cat(sprintf("  Found %d downregulated genes (padj<0.05, log2FC<-1)\n", nrow(downreg)))

        if (nrow(downreg) > 0) {
            # Extract Ensembl IDs
            downreg_genes <- downreg$gene_id
            if (is.null(downreg_genes)) downreg_genes <- rownames(downreg)

            result_downreg <- create_gene_beds(
                gene_list = downreg_genes,
                gene_type = "ensembl",
                name_prefix = "downregulated_degs"
            )
        }
    } else {
        cat("  WARNING: DESeq2 results also not found.\n")
    }
}

# ============================================================================
# Process upregulated DEGs
# ============================================================================

cat("\n=== Processing Upregulated DEGs ===\n")

upreg_promoter_file <- "../meDIP/results/12_gene_sets/upregulated_promoters.bed"
if (file.exists(upreg_promoter_file)) {
    upreg_bed <- read.table(upreg_promoter_file, sep = "\t", header = FALSE,
                             stringsAsFactors = FALSE)
    upreg_genes <- upreg_bed[, 4]

    if (grepl("^ENSG", upreg_genes[1])) {
        gene_type <- "ensembl"
    } else {
        gene_type <- "symbol"
    }

    cat(sprintf("Read %d upregulated genes from meDIP BED file\n", length(upreg_genes)))

    result_upreg <- create_gene_beds(
        gene_list = upreg_genes,
        gene_type = gene_type,
        name_prefix = "upregulated_degs"
    )
} else {
    cat("INFO: Upregulated promoters file not found at", upreg_promoter_file, "\n")
}

# ============================================================================
# Summary
# ============================================================================

cat("\n============================================\n")
cat("Gene Set BED Files Created\n")
cat("============================================\n\n")

cat("Output directory: results/19_metagene_beds/\n\n")

bed_files <- list.files("results/19_metagene_beds", pattern = "\\.bed$", full.names = TRUE)
if (length(bed_files) > 0) {
    cat("Created files:\n")
    for (f in bed_files) {
        n_regions <- length(readLines(f))
        cat(sprintf("  - %s (%d regions)\n", basename(f), n_regions))
    }
}

cat("\nBED file types:\n")
cat("  - *_gene_body.bed: Full gene coordinates (for scale-regions mode)\n")
cat("  - *_promoter.bed: Promoter regions TSS ± 2kb (for reference-point mode)\n")

cat("\nFinished:", format(Sys.time()), "\n")
