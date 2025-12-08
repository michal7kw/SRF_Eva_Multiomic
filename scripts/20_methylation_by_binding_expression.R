#!/usr/bin/env Rscript

# ============================================================================
# Methylation Analysis by Binding and Expression Status
# ============================================================================
# Purpose: Create comprehensive gene classifications combining:
#   - Cut&Tag binding status (TES-bound, TEAD1-bound, both, neither)
#   - Expression status (upregulated, downregulated, unchanged)
#
# Output: BED files for each gene category for methylation visualization
#   - Promoter regions (TSS ± 2kb)
#   - Gene body regions (TSS to TES)
# ============================================================================

suppressPackageStartupMessages({
    library(GenomicRanges)
    library(GenomicFeatures)
    library(rtracklayer)
    library(org.Hs.eg.db)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
})

# Set working directory
setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis")

# Create output directory
dir.create("results/20_methylation_binding_expression", recursive = TRUE, showWarnings = FALSE)
outdir <- "results/20_methylation_binding_expression"

cat("============================================\n")
cat("Methylation by Binding and Expression Analysis\n")
cat("Started:", format(Sys.time()), "\n")
cat("============================================\n\n")

# ============================================================================
# Load data
# ============================================================================

cat("=== Loading data ===\n\n")

# Load RNA-seq differential expression results
deseq_file <- "../SRF_Eva_RNA/results/05_deseq2/deseq2_results_TES_vs_GFP.txt"
if (!file.exists(deseq_file)) {
    stop("DESeq2 results not found: ", deseq_file)
}
deseq <- read.table(deseq_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
cat(sprintf("  Loaded DESeq2 results: %d genes\n", nrow(deseq)))

# Get gene ID column
if ("gene_id" %in% colnames(deseq)) {
    deseq$ensembl_id <- gsub("\\..*", "", deseq$gene_id)
} else {
    deseq$ensembl_id <- gsub("\\..*", "", rownames(deseq))
    deseq$gene_id <- rownames(deseq)
}

# Load Cut&Tag peak annotations
tes_peaks_file <- "../SRF_Eva_CUTandTAG/results/07_analysis_narrow/TES_peaks_annotated.csv"
tead1_peaks_file <- "../SRF_Eva_CUTandTAG/results/07_analysis_narrow/TEAD1_peaks_annotated.csv"

# Check if files exist
if (!file.exists(tes_peaks_file)) {
    warning("TES peaks not found: ", tes_peaks_file)
    tes_peaks <- NULL
} else {
    tes_peaks <- read.csv(tes_peaks_file, stringsAsFactors = FALSE)
    cat(sprintf("  Loaded TES peaks: %d peaks\n", nrow(tes_peaks)))
}

if (!file.exists(tead1_peaks_file)) {
    warning("TEAD1 peaks not found: ", tead1_peaks_file)
    tead1_peaks <- NULL
} else {
    tead1_peaks <- read.csv(tead1_peaks_file, stringsAsFactors = FALSE)
    cat(sprintf("  Loaded TEAD1 peaks: %d peaks\n", nrow(tead1_peaks)))
}

# Load gene annotation from GTF
gtf_file <- "/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/annotation/gencode.v44.annotation.gtf"
txdb <- makeTxDbFromGFF(gtf_file, format = "gtf")
genes_gr <- genes(txdb)

# Clean Ensembl IDs (remove version numbers) before mapping
clean_ensembl_ids <- gsub("\\..*", "", names(genes_gr))

# Map to gene symbols using cleaned IDs
gene_symbols <- mapIds(org.Hs.eg.db,
    keys = clean_ensembl_ids,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
)
genes_gr$gene_symbol <- gene_symbols
genes_gr$ensembl_id <- clean_ensembl_ids

cat(sprintf("  Loaded gene annotations: %d genes\n", length(genes_gr)))

# ============================================================================
# Extract genes bound by TES/TEAD1
# ============================================================================

cat("\n=== Identifying bound genes ===\n")

# Function to extract bound genes from peak annotation
get_bound_genes <- function(peaks_df, tf_name) {
    if (is.null(peaks_df) || nrow(peaks_df) == 0) {
        cat(sprintf("  WARNING: No peaks available for %s\n", tf_name))
        return(character(0))
    }

    # Get promoter peaks
    promoter_peaks <- peaks_df[grepl("Promoter", peaks_df$annotation, ignore.case = TRUE), ]

    # Get gene IDs
    if ("geneId" %in% colnames(promoter_peaks)) {
        gene_ids <- unique(promoter_peaks$geneId)
    } else if ("GENEID" %in% colnames(promoter_peaks)) {
        gene_ids <- unique(promoter_peaks$GENEID)
    } else {
        # Try to find any gene-related column
        gene_cols <- grep("gene", colnames(promoter_peaks), ignore.case = TRUE, value = TRUE)
        if (length(gene_cols) > 0) {
            gene_ids <- unique(promoter_peaks[, gene_cols[1]])
        } else {
            cat(sprintf("  WARNING: No gene ID column found for %s\n", tf_name))
            return(character(0))
        }
    }

    # Convert Entrez to Ensembl if needed
    if (all(grepl("^[0-9]+$", head(na.omit(gene_ids))))) {
        # These are Entrez IDs, convert to Ensembl
        ensembl_ids <- mapIds(org.Hs.eg.db,
            keys = as.character(gene_ids),
            column = "ENSEMBL",
            keytype = "ENTREZID",
            multiVals = "first"
        )
        gene_ids <- unique(na.omit(ensembl_ids))
    }

    # Clean Ensembl IDs (remove version)
    gene_ids <- gsub("\\..*", "", gene_ids)
    gene_ids <- unique(na.omit(gene_ids))

    cat(sprintf("  %s: %d genes with promoter binding\n", tf_name, length(gene_ids)))
    return(gene_ids)
}

# Get bound genes
tes_bound_genes <- get_bound_genes(tes_peaks, "TES")
tead1_bound_genes <- get_bound_genes(tead1_peaks, "TEAD1")

# ============================================================================
# Classify genes by binding and expression
# ============================================================================

cat("\n=== Classifying genes ===\n")

# Expression thresholds
padj_cutoff <- 0.05
log2fc_up <- 1      # |log2FC| > 1 for up/down
log2fc_down <- -1

# Classify expression
deseq <- deseq %>%
    mutate(
        expression_status = case_when(
            is.na(padj) ~ "not_tested",
            padj >= padj_cutoff ~ "unchanged",
            log2FoldChange > log2fc_up ~ "upregulated",
            log2FoldChange < log2fc_down ~ "downregulated",
            TRUE ~ "unchanged"
        )
    )

# Classify binding
deseq <- deseq %>%
    mutate(
        tes_bound = ensembl_id %in% tes_bound_genes,
        tead1_bound = ensembl_id %in% tead1_bound_genes,
        binding_status = case_when(
            tes_bound & tead1_bound ~ "TES_TEAD1_bound",
            tes_bound & !tead1_bound ~ "TES_only_bound",
            !tes_bound & tead1_bound ~ "TEAD1_only_bound",
            TRUE ~ "Neither_bound"
        )
    )

# Combined classification
deseq <- deseq %>%
    mutate(
        combined_class = paste(binding_status, expression_status, sep = "_")
    )

# Summary statistics
cat("\nBinding status summary:\n")
print(table(deseq$binding_status))

cat("\nExpression status summary:\n")
print(table(deseq$expression_status))

cat("\nCombined classification summary:\n")
print(table(deseq$binding_status, deseq$expression_status))

# Save classification table
classification_file <- file.path(outdir, "gene_classification_summary.csv")
write.csv(deseq, classification_file, row.names = FALSE)
cat(sprintf("\nSaved classification table: %s\n", classification_file))

# ============================================================================
# Create BED files for each category
# ============================================================================

cat("\n=== Creating BED files ===\n")

# Function to create BED file for a gene set
create_bed_for_genes <- function(gene_ids, name, genes_ref = genes_gr) {
    # Clean IDs
    gene_ids <- gsub("\\..*", "", gene_ids)

    # Match to reference
    matched_idx <- which(gsub("\\..*", "", names(genes_ref)) %in% gene_ids)

    if (length(matched_idx) < 5) {
        cat(sprintf("  %s: Too few genes (%d) - skipping\n", name, length(matched_idx)))
        return(NULL)
    }

    matched_genes <- genes_ref[matched_idx]

    # IMPORTANT: Convert UCSC chromosome names (chr1, chr2) to Ensembl format (1, 2)
    # This is required because BigWig files use Ensembl naming
    chr_names <- as.character(seqnames(matched_genes))
    chr_names <- gsub("^chr", "", chr_names)  # Remove "chr" prefix for Ensembl compatibility

    # Gene body BED
    gene_body <- data.frame(
        chr = chr_names,
        start = start(matched_genes) - 1,
        end = end(matched_genes),
        name = ifelse(is.na(matched_genes$gene_symbol),
                      gsub("\\..*", "", names(matched_genes)),
                      matched_genes$gene_symbol),
        score = 0,
        strand = as.character(strand(matched_genes))
    )

    # Promoter BED (TSS ± 2kb)
    tss <- ifelse(gene_body$strand == "+", gene_body$start + 1, gene_body$end)
    promoter <- data.frame(
        chr = gene_body$chr,  # Already converted to Ensembl naming
        start = pmax(0, tss - 2000),
        end = tss + 2000,
        name = gene_body$name,
        score = 0,
        strand = gene_body$strand
    )

    # Write files
    body_file <- file.path(outdir, "beds", sprintf("%s_gene_body.bed", name))
    prom_file <- file.path(outdir, "beds", sprintf("%s_promoter.bed", name))

    dir.create(file.path(outdir, "beds"), showWarnings = FALSE)

    write.table(gene_body, body_file, sep = "\t", quote = FALSE,
                row.names = FALSE, col.names = FALSE)
    write.table(promoter, prom_file, sep = "\t", quote = FALSE,
                row.names = FALSE, col.names = FALSE)

    cat(sprintf("  %s: %d genes\n", name, nrow(gene_body)))
    cat(sprintf("    -> %s\n", body_file))
    cat(sprintf("    -> %s\n", prom_file))

    return(list(body = body_file, promoter = prom_file, n = nrow(gene_body)))
}

# Create BED files for binding categories
cat("\nBinding categories:\n")
for (status in c("TES_only_bound", "TEAD1_only_bound", "TES_TEAD1_bound", "Neither_bound")) {
    genes <- deseq$ensembl_id[deseq$binding_status == status]
    create_bed_for_genes(genes, status)
}

# Create BED files for expression categories
cat("\nExpression categories:\n")
for (status in c("upregulated", "downregulated", "unchanged")) {
    genes <- deseq$ensembl_id[deseq$expression_status == status]
    create_bed_for_genes(genes, status)
}

# Create BED files for combined binding + expression
cat("\nCombined binding + expression categories:\n")

# TES-bound up/down
create_bed_for_genes(
    deseq$ensembl_id[deseq$tes_bound & deseq$expression_status == "upregulated"],
    "TES_bound_upregulated"
)
create_bed_for_genes(
    deseq$ensembl_id[deseq$tes_bound & deseq$expression_status == "downregulated"],
    "TES_bound_downregulated"
)

# TEAD1-bound up/down
create_bed_for_genes(
    deseq$ensembl_id[deseq$tead1_bound & deseq$expression_status == "upregulated"],
    "TEAD1_bound_upregulated"
)
create_bed_for_genes(
    deseq$ensembl_id[deseq$tead1_bound & deseq$expression_status == "downregulated"],
    "TEAD1_bound_downregulated"
)

# Both bound up/down
create_bed_for_genes(
    deseq$ensembl_id[deseq$binding_status == "TES_TEAD1_bound" & deseq$expression_status == "upregulated"],
    "Both_bound_upregulated"
)
create_bed_for_genes(
    deseq$ensembl_id[deseq$binding_status == "TES_TEAD1_bound" & deseq$expression_status == "downregulated"],
    "Both_bound_downregulated"
)

# Unbound DEGs (indirect effects)
create_bed_for_genes(
    deseq$ensembl_id[deseq$binding_status == "Neither_bound" & deseq$expression_status == "upregulated"],
    "Indirect_upregulated"
)
create_bed_for_genes(
    deseq$ensembl_id[deseq$binding_status == "Neither_bound" & deseq$expression_status == "downregulated"],
    "Indirect_downregulated"
)

# ============================================================================
# Create summary plot
# ============================================================================

cat("\n=== Creating summary visualization ===\n")

# Create summary data for plotting
summary_df <- deseq %>%
    filter(expression_status %in% c("upregulated", "downregulated")) %>%
    group_by(binding_status, expression_status) %>%
    summarise(n = n(), .groups = "drop")

if (nrow(summary_df) > 0) {
    # Bar plot
    p <- ggplot(summary_df, aes(x = binding_status, y = n, fill = expression_status)) +
        geom_bar(stat = "identity", position = "dodge") +
        scale_fill_manual(values = c("downregulated" = "steelblue", "upregulated" = "firebrick")) +
        labs(
            title = "DEGs by Binding Status",
            x = "Binding Status",
            y = "Number of Genes",
            fill = "Expression"
        ) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

    ggsave(file.path(outdir, "binding_expression_summary.pdf"), p, width = 10, height = 6)
    ggsave(file.path(outdir, "binding_expression_summary.png"), p, width = 10, height = 6, dpi = 300)
    cat("  Created summary plots\n")
}

# ============================================================================
# Summary
# ============================================================================

cat("\n============================================\n")
cat("Gene Classification Complete\n")
cat("============================================\n\n")

cat(sprintf("Output directory: %s/\n\n", outdir))

cat("Files created:\n")
cat("  - gene_classification_summary.csv: Full classification table\n")
cat("  - beds/*_gene_body.bed: Gene body coordinates\n")
cat("  - beds/*_promoter.bed: Promoter coordinates (TSS ± 2kb)\n")
cat("  - binding_expression_summary.pdf: Summary visualization\n")

cat("\nBED file categories:\n")
bed_files <- list.files(file.path(outdir, "beds"), pattern = "\\.bed$")
for (f in bed_files) {
    n <- length(readLines(file.path(outdir, "beds", f)))
    cat(sprintf("  - %s: %d regions\n", f, n))
}

cat("\nFinished:", format(Sys.time()), "\n")
