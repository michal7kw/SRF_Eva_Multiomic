#!/usr/bin/env Rscript
#
# DEGs DOWN STRATIFIED BY TES/TEAD1 BINDING - METAGENE ANALYSIS
# =============================================================================
#
# Creates 3 gene sets for metagene analysis:
#   1. DEGs DOWN WITH TES/TEAD1 binding at promoter (±2kb from TSS)
#   2. DEGs DOWN WITHOUT TES/TEAD1 binding at promoter
#   3. Random expressed genes (same N as group 1) WITHOUT binding (control)
#
# Purpose: Test whether TES-induced methylation requires direct binding
#   - Group 1: Expect differential MeDIP (TES > GFP)
#   - Groups 2 & 3: Expect similar MeDIP (TES ≈ GFP) if methylation requires binding
#
# =============================================================================

suppressPackageStartupMessages({
    library(GenomicRanges)
    library(rtracklayer)
    library(dplyr)
})

setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_1")

cat("==========================================================\n")
cat("DEGs DOWN STRATIFIED BY TES/TEAD1 BINDING\n")
cat("==========================================================\n")
cat("Analysis started:", as.character(Sys.time()), "\n\n")

# =============================================================================
# PATH CONFIGURATION
# =============================================================================

# Peak files (narrow peaks from merged replicates)
TES_PEAKS <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_CUTandTAG/results/05_peaks_narrow/TES_peaks.narrowPeak"
TEAD1_PEAKS <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_CUTandTAG/results/05_peaks_narrow/TEAD1_peaks.narrowPeak"

# DESeq2 results
DESEQ2_FILE <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_RNA/results/05_deseq2/deseq2_results_TES_vs_GFP.txt"

# GTF annotation
GTF_FILE <- "/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/annotation/gencode.v44.annotation.gtf"

# Output directory
OUTPUT_BASE <- "output/29_degs_down_binding_methylation"
dir.create(OUTPUT_BASE, showWarnings = FALSE, recursive = TRUE)

# Parameters
PROMOTER_UPSTREAM <- 2000   # bp upstream of TSS
PROMOTER_DOWNSTREAM <- 2000 # bp downstream of TSS

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

cat(sprintf("  Protein-coding genes on standard chromosomes: %d\n", nrow(genes_df)))

# Clean gene IDs (remove version)
genes_df$gene_id_clean <- gsub("\\..*", "", genes_df$gene_id)

# Create GRanges for genes
genes_gr <- GRanges(
    seqnames = genes_df$seqnames,
    ranges = IRanges(start = genes_df$start, end = genes_df$end),
    strand = genes_df$strand,
    gene_name = genes_df$gene_name,
    gene_id = genes_df$gene_id_clean
)

# Create promoter regions (±2kb from TSS)
promoters_gr <- promoters(genes_gr, upstream = PROMOTER_UPSTREAM, downstream = PROMOTER_DOWNSTREAM)

cat(sprintf("  Created promoter regions: ±%d bp from TSS\n", PROMOTER_UPSTREAM))
cat("\n")

# =============================================================================
# PHASE 2: LOAD TES/TEAD1 PEAKS
# =============================================================================

cat("=== PHASE 2: Loading TES/TEAD1 Binding Peaks ===\n")

# Function to load narrowPeak files
load_peaks <- function(peak_file, name) {
    peaks <- read.table(peak_file, header = FALSE, stringsAsFactors = FALSE,
                        col.names = c("chr", "start", "end", "name", "score",
                                      "strand", "signalValue", "pValue", "qValue", "peak"))
    cat(sprintf("  %s peaks loaded: %d\n", name, nrow(peaks)))

    # Add "chr" prefix if not present
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
# PHASE 3: IDENTIFY GENES WITH/WITHOUT PROMOTER BINDING
# =============================================================================

cat("=== PHASE 3: Identifying Genes with Promoter Binding ===\n")

# Find overlaps between promoters and peaks
promoter_overlaps <- findOverlaps(promoters_gr, all_binding_peaks)
genes_with_binding_idx <- unique(queryHits(promoter_overlaps))
genes_without_binding_idx <- setdiff(1:length(genes_gr), genes_with_binding_idx)

cat(sprintf("  Genes with TES/TEAD1 binding at promoter (±2kb): %d\n", length(genes_with_binding_idx)))
cat(sprintf("  Genes without TES/TEAD1 binding at promoter: %d\n", length(genes_without_binding_idx)))

# Get gene IDs for each group
bound_gene_ids <- genes_gr$gene_id[genes_with_binding_idx]
unbound_gene_ids <- genes_gr$gene_id[genes_without_binding_idx]

cat("\n")

# =============================================================================
# PHASE 4: LOAD DESEQ2 AND DEFINE GENE SETS
# =============================================================================

cat("=== PHASE 4: Loading DESeq2 Results and Defining Gene Sets ===\n")

deseq2 <- read.delim(DESEQ2_FILE, stringsAsFactors = FALSE)
deseq2$gene_id_clean <- gsub("\\..*", "", deseq2$gene_id)

cat(sprintf("  Total genes in DESeq2: %d\n", nrow(deseq2)))

# All expressed genes (non-NA padj)
expressed_genes <- deseq2 %>%
    filter(!is.na(padj)) %>%
    pull(gene_id_clean)
cat(sprintf("  Expressed genes (non-NA padj): %d\n", length(expressed_genes)))

# DEGs DOWN: significantly downregulated (padj < 0.05 & log2FC < 0)
degs_down <- deseq2 %>%
    filter(!is.na(padj), padj < 0.05, log2FoldChange < 0) %>%
    pull(gene_id_clean)
cat(sprintf("  DEGs DOWN (padj<0.05, log2FC<0): %d\n", length(degs_down)))

# =============================================================================
# PHASE 5: CREATE 3 GENE SETS
# =============================================================================

cat("\n=== PHASE 5: Creating 3 Gene Sets ===\n")

# GROUP 1: DEGs DOWN WITH TES/TEAD1 binding at promoter
degs_down_with_binding <- intersect(degs_down, bound_gene_ids)
cat(sprintf("  GROUP 1 - DEGs DOWN with binding: %d genes\n", length(degs_down_with_binding)))

# GROUP 2: DEGs DOWN WITHOUT TES/TEAD1 binding at promoter
degs_down_without_binding <- intersect(degs_down, unbound_gene_ids)
cat(sprintf("  GROUP 2 - DEGs DOWN without binding: %d genes\n", length(degs_down_without_binding)))

# GROUP 3: Random expressed genes (same N as group 1) WITHOUT binding
# Exclude DEGs to get truly "random" non-DE genes
all_de_genes <- deseq2 %>%
    filter(!is.na(padj), padj < 0.05) %>%
    pull(gene_id_clean)

non_de_expressed_unbound <- setdiff(
    intersect(expressed_genes, unbound_gene_ids),
    all_de_genes
)
cat(sprintf("  Non-DE expressed genes without binding (pool for random): %d\n", length(non_de_expressed_unbound)))

# Sample same number as group 1
set.seed(42)  # For reproducibility
n_to_sample <- length(degs_down_with_binding)
if (length(non_de_expressed_unbound) >= n_to_sample) {
    random_control <- sample(non_de_expressed_unbound, n_to_sample)
} else {
    cat("  WARNING: Not enough non-DE unbound genes. Using all available.\n")
    random_control <- non_de_expressed_unbound
}
cat(sprintf("  GROUP 3 - Random control (matched N, no binding): %d genes\n", length(random_control)))

cat("\n")

# =============================================================================
# PHASE 6: CREATE BED FILES
# =============================================================================

cat("=== PHASE 6: Creating BED Files ===\n")

# Function to create gene body BED file
create_gene_body_bed <- function(gene_ids_to_include, genes_df, output_file, description) {
    selected <- genes_df[genes_df$gene_id_clean %in% gene_ids_to_include, ]

    if (nrow(selected) == 0) {
        cat(sprintf("  WARNING: No genes matched for %s\n", description))
        return(0)
    }

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

    cat(sprintf("  Created: %s (%d genes) - %s\n",
                basename(output_file), nrow(bed), description))
    return(nrow(bed))
}

# Create BED files for each group
n_group1 <- create_gene_body_bed(
    degs_down_with_binding, genes_df,
    file.path(OUTPUT_BASE, "DEGs_DOWN_with_binding.bed"),
    "DEGs DOWN + TES/TEAD1 binding"
)

n_group2 <- create_gene_body_bed(
    degs_down_without_binding, genes_df,
    file.path(OUTPUT_BASE, "DEGs_DOWN_without_binding.bed"),
    "DEGs DOWN + NO binding"
)

n_group3 <- create_gene_body_bed(
    random_control, genes_df,
    file.path(OUTPUT_BASE, "random_control_no_binding.bed"),
    "Random control (no binding)"
)

# Save gene lists as CSV for reference
gene_list_df <- data.frame(
    gene_id = c(degs_down_with_binding, degs_down_without_binding, random_control),
    group = c(
        rep("DEGs_DOWN_with_binding", length(degs_down_with_binding)),
        rep("DEGs_DOWN_without_binding", length(degs_down_without_binding)),
        rep("random_control_no_binding", length(random_control))
    )
)
write.csv(gene_list_df, file.path(OUTPUT_BASE, "gene_lists.csv"), row.names = FALSE)

# Save summary counts
counts_df <- data.frame(
    group = c("DEGs_DOWN_with_binding", "DEGs_DOWN_without_binding", "random_control_no_binding"),
    n_genes = c(n_group1, n_group2, n_group3),
    description = c(
        "DEGs DOWN with TES/TEAD1 binding at promoter (±2kb)",
        "DEGs DOWN without TES/TEAD1 binding at promoter",
        "Random expressed non-DE genes without binding (matched N)"
    )
)
write.csv(counts_df, file.path(OUTPUT_BASE, "gene_counts.csv"), row.names = FALSE)

cat("\n")

# =============================================================================
# PHASE 7: ADDITIONAL STATISTICS
# =============================================================================

cat("=== PHASE 7: Additional Statistics ===\n")

# Get log2FC values for each group
deseq2_subset <- deseq2[deseq2$gene_id_clean %in% c(degs_down_with_binding, degs_down_without_binding), ]
deseq2_subset$group <- ifelse(
    deseq2_subset$gene_id_clean %in% degs_down_with_binding,
    "with_binding",
    "without_binding"
)

# Summary statistics
group_stats <- deseq2_subset %>%
    group_by(group) %>%
    summarise(
        n = n(),
        mean_log2FC = mean(log2FoldChange, na.rm = TRUE),
        median_log2FC = median(log2FoldChange, na.rm = TRUE),
        sd_log2FC = sd(log2FoldChange, na.rm = TRUE),
        .groups = "drop"
    )

cat("\n  Log2FC statistics for DEGs DOWN:\n")
print(as.data.frame(group_stats), row.names = FALSE)

# Statistical test comparing log2FC between groups
wilcox_result <- wilcox.test(
    log2FoldChange ~ group,
    data = deseq2_subset
)
cat(sprintf("\n  Wilcoxon test (log2FC: with_binding vs without_binding):\n"))
cat(sprintf("    p-value: %.4e\n", wilcox_result$p.value))

# Save statistics
write.csv(group_stats, file.path(OUTPUT_BASE, "log2FC_statistics.csv"), row.names = FALSE)

cat("\n")

# =============================================================================
# SUMMARY
# =============================================================================

cat("==========================================================\n")
cat("GENE SET CREATION COMPLETE\n")
cat("==========================================================\n")
cat("Completed:", as.character(Sys.time()), "\n")
cat(sprintf("Output directory: %s\n\n", OUTPUT_BASE))

cat("Created files:\n")
cat(sprintf("  1. DEGs_DOWN_with_binding.bed (%d genes)\n", n_group1))
cat(sprintf("     → DEGs DOWN with TES/TEAD1 binding at promoter\n"))
cat(sprintf("     → Expected: DIFFERENTIAL MeDIP (TES > GFP)\n\n"))

cat(sprintf("  2. DEGs_DOWN_without_binding.bed (%d genes)\n", n_group2))
cat(sprintf("     → DEGs DOWN without TES/TEAD1 binding\n"))
cat(sprintf("     → Expected: SIMILAR MeDIP if methylation requires binding\n\n"))

cat(sprintf("  3. random_control_no_binding.bed (%d genes)\n", n_group3))
cat(sprintf("     → Random expressed genes without binding (matched N)\n"))
cat(sprintf("     → Expected: SIMILAR MeDIP (baseline control)\n\n"))

cat("  4. gene_lists.csv - Full gene ID list by group\n")
cat("  5. gene_counts.csv - Summary counts\n")
cat("  6. log2FC_statistics.csv - Expression statistics\n")

cat("\nNext step: Run the shell script to compute metagene profiles\n")
cat("  sbatch 29_degs_down_binding_methylation.sh\n")
