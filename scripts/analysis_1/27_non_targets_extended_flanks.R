#!/usr/bin/env Rscript
#
# NON-TARGET GENES ANALYSIS
# =============================================================================
# Identify genes that are NOT bound by TES or TEAD1 (non-targets)
# and stratify by DMR status to see if methylation changes occur
# independently of direct TF binding.
#
# This tests whether methylation changes at non-target genes are:
# 1. Similar to target genes (suggesting indirect/global effects)
# 2. Different (suggesting methylation requires binding)
# =============================================================================

suppressPackageStartupMessages({
    library(GenomicRanges)
    library(rtracklayer)
    library(dplyr)
})

setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_1")

cat("==========================================================\n")
cat("NON-TARGET GENES - DMR STRATIFIED ANALYSIS\n")
cat("==========================================================\n")
cat("Analysis started:", as.character(Sys.time()), "\n\n")

# =============================================================================
# PATH CONFIGURATION
# =============================================================================

# Peak files (consensus peaks from merged replicates)
TES_PEAKS <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_CUTandTAG/results/05_peaks_narrow/TES_peaks.narrowPeak"
TEAD1_PEAKS <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_CUTandTAG/results/05_peaks_narrow/TEAD1_peaks.narrowPeak"

# DESeq2 results
DESEQ2_FILE <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_RNA/results/05_deseq2/deseq2_results_TES_vs_GFP.txt"

# GTF annotation
GTF_FILE <- "/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/annotation/gencode.v44.annotation.gtf"

# DMR files
DMR_FILE <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP/results/07_differential_MEDIPS/TES_vs_GFP_DMRs_FDR05_FC2.csv"

# Output directory
OUTPUT_BASE <- "output/27_non_targets_extended_flanks"
dir.create(OUTPUT_BASE, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# PHASE 1: LOAD GENE ANNOTATIONS
# =============================================================================

cat("=== PHASE 1: Loading Gene Annotations ===\n")

gtf <- rtracklayer::import(GTF_FILE)
genes_gtf <- gtf[gtf$type == "gene"]
genes_df <- as.data.frame(genes_gtf)

# Filter for protein-coding genes on standard chromosomes
genes_df <- genes_df %>%
    filter(gene_type == "protein_coding") %>%
    filter(seqnames %in% paste0("chr", c(1:22, "X", "Y")))

cat(sprintf("  Protein-coding genes: %d\n", nrow(genes_df)))

genes_df$gene_id_clean <- gsub("\\..*", "", genes_df$gene_id)

# Create GRanges for genes
genes_gr <- GRanges(
    seqnames = genes_df$seqnames,
    ranges = IRanges(start = genes_df$start, end = genes_df$end),
    strand = genes_df$strand,
    gene_name = genes_df$gene_name,
    gene_id = genes_df$gene_id_clean
)

# Create promoter regions (+/- 2kb from TSS)
promoters_gr <- promoters(genes_gr, upstream = 2000, downstream = 2000)

cat("\n")

# =============================================================================
# PHASE 2: LOAD TES/TEAD1 PEAKS
# =============================================================================

cat("=== PHASE 2: Loading TES/TEAD1 Binding Peaks ===\n")

# Load narrowPeak files
load_peaks <- function(peak_file, name) {
    peaks <- read.table(peak_file, header = FALSE, stringsAsFactors = FALSE,
                        col.names = c("chr", "start", "end", "name", "score",
                                      "strand", "signalValue", "pValue", "qValue", "peak"))
    cat(sprintf("  %s peaks: %d\n", name, nrow(peaks)))

    # Add "chr" prefix if not present (peaks use "1", genes use "chr1")
    chr_names <- peaks$chr
    if (!grepl("^chr", chr_names[1])) {
        chr_names <- paste0("chr", chr_names)
    }

    GRanges(
        seqnames = chr_names,
        ranges = IRanges(start = peaks$start + 1, end = peaks$end),  # Convert to 1-based
        score = peaks$score,
        signalValue = peaks$signalValue
    )
}

tes_peaks_gr <- load_peaks(TES_PEAKS, "TES")
tead1_peaks_gr <- load_peaks(TEAD1_PEAKS, "TEAD1")

# Combine all TES/TEAD1 peaks
all_binding_peaks <- c(tes_peaks_gr, tead1_peaks_gr)
cat(sprintf("  Combined TES+TEAD1 peaks: %d\n", length(all_binding_peaks)))

cat("\n")

# =============================================================================
# PHASE 3: IDENTIFY NON-TARGET GENES
# =============================================================================

cat("=== PHASE 3: Identifying Non-Target Genes ===\n")

# Find genes with peaks overlapping their promoter region
promoter_overlaps <- findOverlaps(promoters_gr, all_binding_peaks)
genes_with_binding <- unique(queryHits(promoter_overlaps))

cat(sprintf("  Genes with TES/TEAD1 binding at promoter (+/-2kb): %d\n", length(genes_with_binding)))

# Non-target genes = genes WITHOUT any TES/TEAD1 binding at promoter
non_target_indices <- setdiff(1:length(genes_gr), genes_with_binding)
cat(sprintf("  Non-target genes (no binding): %d\n", length(non_target_indices)))

cat("\n")

# =============================================================================
# PHASE 4: LOAD DMRs
# =============================================================================

cat("=== PHASE 4: Loading DMRs ===\n")

dmr_data <- read.csv(DMR_FILE, stringsAsFactors = FALSE)
cat(sprintf("  Total DMRs: %d\n", nrow(dmr_data)))

dmr_hyper <- dmr_data %>% filter(logFC > 0)
dmr_hypo <- dmr_data %>% filter(logFC < 0)

cat(sprintf("  Hypermethylated DMRs: %d\n", nrow(dmr_hyper)))
cat(sprintf("  Hypomethylated DMRs: %d\n", nrow(dmr_hypo)))

# Create GRanges for DMRs (add chr prefix if needed)
create_dmr_gr <- function(dmr_df) {
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
# PHASE 5: STRATIFY NON-TARGETS BY DMR STATUS
# =============================================================================

cat("=== PHASE 5: Stratifying Non-Targets by DMR Status ===\n")

# Get non-target genes GRanges
non_target_gr <- genes_gr[non_target_indices]

# Find overlaps with DMRs
hyper_overlaps <- findOverlaps(non_target_gr, dmr_hyper_gr)
hypo_overlaps <- findOverlaps(non_target_gr, dmr_hypo_gr)

# Get indices within non_target_gr
non_target_with_hyper <- unique(queryHits(hyper_overlaps))
non_target_with_hypo <- unique(queryHits(hypo_overlaps))

# Non-targets with ONLY hyper DMRs
non_target_hyper_only <- setdiff(non_target_with_hyper, non_target_with_hypo)
# Non-targets with ONLY hypo DMRs
non_target_hypo_only <- setdiff(non_target_with_hypo, non_target_with_hyper)
# Non-targets with no DMRs
all_dmr_overlaps <- findOverlaps(non_target_gr, c(dmr_hyper_gr, dmr_hypo_gr))
non_target_with_any_dmr <- unique(queryHits(all_dmr_overlaps))
non_target_no_dmr <- setdiff(1:length(non_target_gr), non_target_with_any_dmr)

cat(sprintf("  Non-targets with hypermethylated DMRs only: %d\n", length(non_target_hyper_only)))
cat(sprintf("  Non-targets with hypomethylated DMRs only: %d\n", length(non_target_hypo_only)))
cat(sprintf("  Non-targets without any DMRs: %d\n", length(non_target_no_dmr)))

cat("\n")

# =============================================================================
# PHASE 6: FILTER FOR EXPRESSED GENES
# =============================================================================

cat("=== PHASE 6: Filtering for Expressed Genes ===\n")

deseq2 <- read.delim(DESEQ2_FILE, stringsAsFactors = FALSE)
deseq2$gene_id_clean <- gsub("\\..*", "", deseq2$gene_id)

expressed_genes <- deseq2 %>%
    filter(!is.na(padj)) %>%
    pull(gene_id_clean)

cat(sprintf("  Expressed genes (non-NA padj): %d\n", length(expressed_genes)))

# Filter non-target gene sets to expressed genes
non_target_gene_ids <- non_target_gr$gene_id

# Hyper non-targets (expressed)
hyper_ids <- non_target_gene_ids[non_target_hyper_only]
hyper_expressed <- hyper_ids[hyper_ids %in% expressed_genes]
cat(sprintf("  Non-target + Hypermethylated (expressed): %d\n", length(hyper_expressed)))

# Hypo non-targets (expressed)
hypo_ids <- non_target_gene_ids[non_target_hypo_only]
hypo_expressed <- hypo_ids[hypo_ids %in% expressed_genes]
cat(sprintf("  Non-target + Hypomethylated (expressed): %d\n", length(hypo_expressed)))

# No-DMR non-targets (expressed) - sample if too many
no_dmr_ids <- non_target_gene_ids[non_target_no_dmr]
no_dmr_expressed <- no_dmr_ids[no_dmr_ids %in% expressed_genes]
set.seed(42)
if (length(no_dmr_expressed) > 5000) {
    no_dmr_expressed <- sample(no_dmr_expressed, 5000)
}
cat(sprintf("  Non-target + No DMR (expressed, sampled): %d\n", length(no_dmr_expressed)))

cat("\n")

# =============================================================================
# PHASE 7: CREATE BED FILES
# =============================================================================

cat("=== PHASE 7: Creating BED Files ===\n")

create_gene_body_bed <- function(gene_ids_to_include, genes_df, output_file) {
    selected <- genes_df[genes_df$gene_id_clean %in% gene_ids_to_include, ]

    bed <- data.frame(
        chr = selected$seqnames,
        start = selected$start - 1,
        end = selected$end,
        name = selected$gene_name,
        score = 0,
        strand = selected$strand,
        stringsAsFactors = FALSE
    )

    bed <- bed[!is.na(bed$start) & !is.na(bed$end), ]
    bed <- bed[bed$start >= 0, ]
    bed <- bed[order(bed$chr, bed$start), ]
    bed <- bed[!duplicated(paste(bed$chr, bed$start, bed$end)), ]

    write.table(bed, output_file, quote = FALSE, sep = "\t",
                row.names = FALSE, col.names = FALSE)

    cat(sprintf("  Created: %s (%d genes)\n", basename(output_file), nrow(bed)))
    return(nrow(bed))
}

# Create BED files
n_hyper <- create_gene_body_bed(hyper_expressed, genes_df,
                                 file.path(OUTPUT_BASE, "non_targets_hypermethylated.bed"))
n_hypo <- create_gene_body_bed(hypo_expressed, genes_df,
                                file.path(OUTPUT_BASE, "non_targets_hypomethylated.bed"))
n_no_dmr <- create_gene_body_bed(no_dmr_expressed, genes_df,
                                  file.path(OUTPUT_BASE, "non_targets_no_dmr.bed"))

# Also save gene counts for later use
counts <- data.frame(
    category = c("non_target_hyper", "non_target_hypo", "non_target_no_dmr"),
    n_genes = c(n_hyper, n_hypo, n_no_dmr)
)
write.csv(counts, file.path(OUTPUT_BASE, "gene_counts.csv"), row.names = FALSE)

cat("\n")

# =============================================================================
# SUMMARY
# =============================================================================

cat("==========================================================\n")
cat("NON-TARGET GENE SETS CREATED\n")
cat("==========================================================\n")
cat("Completed:", as.character(Sys.time()), "\n")
cat(sprintf("Output directory: %s\n\n", OUTPUT_BASE))
cat("Created files:\n")
cat(sprintf("  - non_targets_hypermethylated.bed (%d genes)\n", n_hyper))
cat(sprintf("  - non_targets_hypomethylated.bed (%d genes)\n", n_hypo))
cat(sprintf("  - non_targets_no_dmr.bed (%d genes)\n", n_no_dmr))
cat("\n")
cat("These are genes WITHOUT TES/TEAD1 binding at promoters,\n")
cat("stratified by whether they have DMRs in the gene body.\n")
cat("\n")
