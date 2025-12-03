#!/usr/bin/env Rscript

################################################################################
# Script: 10_integrative.R
# Purpose: Multi-omics integration of meDIP, Cut&Tag, and RNA-seq
#
# Description:
#   Integrates DNA methylation (meDIP), TF binding (Cut&Tag), and gene
#   expression (RNA-seq) to classify genes by regulatory mechanisms and
#   test the hypothesis that TES recruits DNA methyltransferases.
#
# Integration Strategy:
#   1. Load annotated DMRs from meDIP analysis (step 08 output)
#   2. Load TES/TEAD1 binding peaks from Cut&Tag
#   3. Load differential expression from RNA-seq
#   4. Identify genes with:
#      - Promoter DMRs (±2kb from TSS)
#      - TF binding at promoters
#      - Expression changes
#   5. Classify genes into regulatory categories
#
# Biological Hypothesis:
#   TES recruits DNA methyltransferases to TEAD1 target genes, causing
#   promoter hypermethylation and transcriptional repression.
#
# Expected Pattern (if hypothesis is correct):
#   Direct epigenetic targets show:
#   - TES binding at promoter (Cut&Tag)
#   - Increased promoter methylation (meDIP hypermethylation)
#   - Decreased gene expression (RNA-seq downregulation)
#
# Key Questions This Analysis Answers:
#   1. Do TES-hypermethylated genes overlap with TES binding sites?
#      → Tests if methylation is direct (at binding sites) or indirect
#   2. Are hypermethylated genes downregulated?
#      → Tests functional consequence of methylation
#   3. Is methylation TEAD1-dependent (TES-specific, not TESmut)?
#      → Tests if TEAD1 binding is required for methylation recruitment
#
# Input:
#   - meDIP annotated DMRs: ../results/08_annotation/*_annotated.csv
#   - Cut&Tag peaks: ../../SRF_Eva_CUTandTAG/results/11_combined_replicates_narrow/peaks/
#   - RNA-seq DEGs: ../../SRF_Eva_RNA/results/05_deseq2/deseq2_results_TES_vs_GFP.txt
#
# Output:
#   - Integrated gene list: integrated_genes.csv
#   - Direct epigenetic targets: direct_epigenetic_targets.csv
#   - Overlap Venn diagram: binding_methylation_expression_venn.png
#   - Heatmap: integrated_regulation_heatmap.pdf
#   - Summary statistics: integration_summary.txt
#
################################################################################

cat("====================================================\n")
cat("Multi-Omics Integration: meDIP + Cut&Tag + RNA-seq\n")
cat("====================================================\n")
cat(paste("Start:", Sys.time(), "\n\n"))

suppressPackageStartupMessages({
    library(GenomicRanges)
    library(ChIPseeker)
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    library(org.Hs.eg.db)
    library(rtracklayer)
    library(ggplot2)
    library(pheatmap)
    library(VennDiagram)
    library(dplyr)
})

# Define paths
base_dir <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top"
medip_annotation_dir <- file.path(base_dir, "meDIP/results/08_annotation")
medip_dmr_dir <- file.path(base_dir, "meDIP/results/07_differential_MEDIPS")
cuttag_dir <- file.path(base_dir, "SRF_Eva_CUTandTAG/results")
rnaseq_dir <- file.path(base_dir, "SRF_Eva_RNA/results/05_deseq2")
out_dir <- file.path(base_dir, "SRF_Eva_Multiomic/scripts/results/10_integrative")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
setwd(out_dir)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

cat("Loading datasets...\n\n")

################################################################################
# 1. Load RNA-seq differential expression
################################################################################

cat("========================================\n")
cat("1. Loading RNA-seq Data\n")
cat("========================================\n")

rnaseq_file <- file.path(rnaseq_dir, "deseq2_results_TES_vs_GFP.txt")
if (file.exists(rnaseq_file)) {
    degs <- read.delim(rnaseq_file, stringsAsFactors = FALSE)
    # Use "gene_symbol" not "gene_name"
    degs_sig <- degs[!is.na(degs$padj) & degs$padj < 0.05 & abs(degs$log2FoldChange) > 1, ]

    cat(paste("  Total genes tested:", nrow(degs), "\n"))
    cat(paste("  Significant DEGs (FDR<0.05, |FC|>2):", nrow(degs_sig), "\n"))
    cat(paste("  Upregulated (TES > GFP):", sum(degs_sig$log2FoldChange > 0), "\n"))
    cat(paste("  Downregulated (TES < GFP):", sum(degs_sig$log2FoldChange < 0), "\n\n"))

    # Create gene-level expression summary
    rna_genes <- data.frame(
        gene = degs$gene_symbol,
        log2FC = degs$log2FoldChange,
        padj = degs$padj,
        is_DEG = !is.na(degs$padj) & degs$padj < 0.05 & abs(degs$log2FoldChange) > 1,
        direction = ifelse(degs$log2FoldChange > 0, "Up", "Down"),
        stringsAsFactors = FALSE
    )
} else {
    cat("  ERROR: RNA-seq file not found at:", rnaseq_file, "\n")
    cat("  Please check that RNA-seq analysis has been completed.\n\n")
    quit(status = 1)
}

################################################################################
# 2. Load Cut&Tag TES and TEAD1 binding peaks
################################################################################

cat("========================================\n")
cat("2. Loading Cut&Tag Binding Peaks\n")
cat("========================================\n")

# Load TES peaks
tes_peak_file <- file.path(cuttag_dir, "11_combined_replicates_narrow/peaks/TES_combined_peaks.narrowPeak")
if (file.exists(tes_peak_file)) {
    tes_peaks <- import(tes_peak_file, format = "narrowPeak")
    # CRITICAL FIX: Convert chromosome names from Ensembl (1,2,3) to UCSC (chr1,chr2,chr3)
    seqlevelsStyle(tes_peaks) <- "UCSC"
    tes_anno <- annotatePeak(tes_peaks, tssRegion = c(-2000, 500), TxDb = txdb, annoDb = "org.Hs.eg.db")
    tes_anno_df <- as.data.frame(tes_anno)

    # Extract genes with TES binding at promoters
    tes_promoter_genes <- tes_anno_df %>%
        filter(grepl("Promoter", annotation)) %>%
        pull(SYMBOL) %>%
        unique()

    cat(paste("  TES peaks:", length(tes_peaks), "\n"))
    cat(paste("  TES promoter-bound genes:", length(tes_promoter_genes), "\n"))
} else {
    cat("  WARNING: TES peaks not found at:", tes_peak_file, "\n")
    cat("  Continuing without TES binding data...\n")
    tes_peaks <- NULL
    tes_promoter_genes <- c()
}

# Load TESmut peaks (for comparison)
tesmut_peak_file <- file.path(cuttag_dir, "11_combined_replicates_narrow/peaks/TESmut_combined_peaks.narrowPeak")
if (file.exists(tesmut_peak_file)) {
    tesmut_peaks <- import(tesmut_peak_file, format = "narrowPeak")
    # CRITICAL FIX: Convert chromosome names
    seqlevelsStyle(tesmut_peaks) <- "UCSC"
    tesmut_anno <- annotatePeak(tesmut_peaks, tssRegion = c(-2000, 500), TxDb = txdb, annoDb = "org.Hs.eg.db")
    tesmut_anno_df <- as.data.frame(tesmut_anno)

    tesmut_promoter_genes <- tesmut_anno_df %>%
        filter(grepl("Promoter", annotation)) %>%
        pull(SYMBOL) %>%
        unique()

    cat(paste("  TESmut peaks:", length(tesmut_peaks), "\n"))
    cat(paste("  TESmut promoter-bound genes:", length(tesmut_promoter_genes), "\n"))
} else {
    cat("  WARNING: TESmut peaks not found\n")
    tesmut_promoter_genes <- c()
}

# Load TEAD1 peaks (endogenous TEAD1)
tead1_peak_file <- file.path(cuttag_dir, "11_combined_replicates_narrow/peaks/TEAD1_combined_peaks.narrowPeak")
if (file.exists(tead1_peak_file)) {
    tead1_peaks <- import(tead1_peak_file, format = "narrowPeak")
    # CRITICAL FIX: Convert chromosome names
    seqlevelsStyle(tead1_peaks) <- "UCSC"
    tead1_anno <- annotatePeak(tead1_peaks, tssRegion = c(-2000, 500), TxDb = txdb, annoDb = "org.Hs.eg.db")
    tead1_anno_df <- as.data.frame(tead1_anno)

    tead1_promoter_genes <- tead1_anno_df %>%
        filter(grepl("Promoter", annotation)) %>%
        pull(SYMBOL) %>%
        unique()

    cat(paste("  TEAD1 peaks:", length(tead1_peaks), "\n"))
    cat(paste("  TEAD1 promoter-bound genes:", length(tead1_promoter_genes), "\n\n"))
} else {
    cat("  WARNING: TEAD1 peaks not found\n\n")
    tead1_promoter_genes <- c()
}

################################################################################
# 3. Load meDIP Annotated Promoter DMRs
################################################################################

cat("========================================\n")
cat("3. Loading meDIP Promoter DMRs\n")
cat("========================================\n")

# Try to use annotated DMRs from step 08 (preferred)
tes_promoter_dmr_file <- file.path(medip_annotation_dir, "TES_vs_GFP_promoter_DMRs.csv")

if (file.exists(tes_promoter_dmr_file)) {
    cat("Using annotated promoter DMRs from step 08...\n")
    dmrs_tes <- read.csv(tes_promoter_dmr_file, stringsAsFactors = FALSE)

    # Extract genes with promoter methylation changes
    dmr_genes_all <- unique(dmrs_tes$SYMBOL[!is.na(dmrs_tes$SYMBOL)])
    dmr_genes_hyper <- unique(dmrs_tes$SYMBOL[!is.na(dmrs_tes$SYMBOL) & dmrs_tes$logFC > 1])
    dmr_genes_hypo <- unique(dmrs_tes$SYMBOL[!is.na(dmrs_tes$SYMBOL) & dmrs_tes$logFC < -1])

    cat(paste("  TES vs GFP promoter DMRs:", nrow(dmrs_tes), "\n"))
    cat(paste("  Genes with promoter DMRs:", length(dmr_genes_all), "\n"))
    cat(paste("  Hypermethylated genes:", length(dmr_genes_hyper), "\n"))
    cat(paste("  Hypomethylated genes:", length(dmr_genes_hypo), "\n"))
} else {
    # Fallback: Load raw MEDIPS output and filter for promoter DMRs
    cat("Annotated DMRs not found, loading raw MEDIPS output...\n")
    dmr_file <- file.path(medip_dmr_dir, "TES_vs_GFP_DMRs_FDR05_FC2.csv")

    if (file.exists(dmr_file)) {
        dmrs <- read.csv(dmr_file, stringsAsFactors = FALSE)

        # Convert to GRanges and annotate
        dmrs_gr <- makeGRangesFromDataFrame(dmrs,
            seqnames.field = "chr", # FIXED
            start.field = "start", # FIXED
            end.field = "stop", # FIXED
            keep.extra.columns = TRUE
        )

        dmr_anno <- annotatePeak(dmrs_gr, tssRegion = c(-2000, 500), TxDb = txdb, annoDb = "org.Hs.eg.db")
        dmrs_tes <- as.data.frame(dmr_anno)

        # Filter for promoter DMRs
        dmrs_tes <- dmrs_tes[grepl("Promoter", dmrs_tes$annotation), ]

        dmr_genes_all <- unique(dmrs_tes$SYMBOL[!is.na(dmrs_tes$SYMBOL)])
        dmr_genes_hyper <- unique(dmrs_tes$SYMBOL[!is.na(dmrs_tes$SYMBOL) & dmrs_tes$logFC > 1]) # FIXED
        dmr_genes_hypo <- unique(dmrs_tes$SYMBOL[!is.na(dmrs_tes$SYMBOL) & dmrs_tes$logFC < -1]) # FIXED

        cat(paste("  Promoter DMRs (TES vs GFP):", nrow(dmrs_tes), "\n"))
        cat(paste("  Genes with promoter DMRs:", length(dmr_genes_all), "\n"))
        cat(paste("  Hypermethylated:", length(dmr_genes_hyper), "\n"))
        cat(paste("  Hypomethylated:", length(dmr_genes_hypo), "\n"))
    } else {
        cat("  ERROR: DMR file not found at:", dmr_file, "\n")
        cat("  Please run MEDIPS analysis (step 07) and annotation (step 08) first.\n\n")
        quit(status = 1)
    }
}

# Load TESmut DMRs for comparison
tesmut_promoter_dmr_file <- file.path(medip_annotation_dir, "TESmut_vs_GFP_promoter_DMRs.csv")
if (file.exists(tesmut_promoter_dmr_file)) {
    dmrs_tesmut <- read.csv(tesmut_promoter_dmr_file, stringsAsFactors = FALSE)
    dmr_genes_tesmut <- unique(dmrs_tesmut$SYMBOL[!is.na(dmrs_tesmut$SYMBOL)])
    cat(paste("  TESmut promoter DMR genes:", length(dmr_genes_tesmut), "\n\n"))
} else {
    cat("  WARNING: TESmut annotated DMRs not found\n\n")
    dmr_genes_tesmut <- c()
}

################################################################################
# 4. Integrate Datasets
################################################################################

cat("========================================\n")
cat("4. Integrating Datasets\n")
cat("========================================\n\n")

# Create master gene list
all_genes <- unique(c(rna_genes$gene, tes_promoter_genes, dmr_genes_all, tead1_promoter_genes))
all_genes <- all_genes[!is.na(all_genes) & all_genes != ""]

cat(paste("Total unique genes in analysis:", length(all_genes), "\n\n"))

# Build integrated data frame
integrated <- data.frame(
    gene = all_genes,
    stringsAsFactors = FALSE
)

# Add TF binding information
integrated$has_TES_binding <- integrated$gene %in% tes_promoter_genes
integrated$has_TESmut_binding <- integrated$gene %in% tesmut_promoter_genes
integrated$has_TEAD1_binding <- integrated$gene %in% tead1_promoter_genes

# Add methylation information
integrated$has_promoter_DMR <- integrated$gene %in% dmr_genes_all
integrated$is_hypermethylated <- integrated$gene %in% dmr_genes_hyper
integrated$is_hypomethylated <- integrated$gene %in% dmr_genes_hypo
integrated$has_TESmut_DMR <- integrated$gene %in% dmr_genes_tesmut

# Add expression data
match_idx <- match(integrated$gene, rna_genes$gene) # vector of numbers c(1, 2, 3, NA)
integrated$log2FC_expression <- NA
integrated$padj_expression <- NA
integrated$is_DEG <- FALSE
integrated$expression_direction <- NA

valid_match <- !is.na(match_idx)
integrated$log2FC_expression[valid_match] <- rna_genes$log2FC[match_idx[valid_match]]
integrated$padj_expression[valid_match] <- rna_genes$padj[match_idx[valid_match]]
integrated$is_DEG[valid_match] <- rna_genes$is_DEG[match_idx[valid_match]]
integrated$expression_direction[valid_match] <- rna_genes$direction[match_idx[valid_match]]

# Classify regulatory mechanisms
integrated$regulatory_class <- "Not regulated"

# Direct epigenetic regulation: TES binding + promoter methylation + expression change
integrated$regulatory_class[
    integrated$has_TES_binding &
        integrated$has_promoter_DMR &
        integrated$is_DEG
] <- "Direct epigenetic"

# Direct non-epigenetic: TES binding + expression change (no methylation)
integrated$regulatory_class[
    integrated$has_TES_binding &
        !integrated$has_promoter_DMR &
        integrated$is_DEG
] <- "Direct non-epigenetic"

# Indirect methylation: Methylation + expression change (no TES binding)
integrated$regulatory_class[
    !integrated$has_TES_binding &
        integrated$has_promoter_DMR &
        integrated$is_DEG
] <- "Indirect methylation"

# Indirect regulation: Expression change only
integrated$regulatory_class[
    !integrated$has_TES_binding &
        !integrated$has_promoter_DMR &
        integrated$is_DEG
] <- "Indirect expression"

# Poised regulation: TES binding + methylation (no expression change yet)
integrated$regulatory_class[
    integrated$has_TES_binding &
        integrated$has_promoter_DMR &
        !integrated$is_DEG
] <- "Poised (bound + methylated)"

# Add detailed subclassification for direct epigenetic targets
integrated$epigenetic_subclass <- NA
direct_epi_idx <- integrated$regulatory_class == "Direct epigenetic"

if (sum(direct_epi_idx) > 0) {
    integrated$epigenetic_subclass[direct_epi_idx] <- paste0(
        ifelse(integrated$is_hypermethylated[direct_epi_idx], "Hyper", "Hypo"),
        "methylated + ",
        integrated$expression_direction[direct_epi_idx],
        "regulated"
    )
}

# Save integrated results
write.csv(integrated, "integrated_genes.csv", row.names = FALSE)
cat("Saved: integrated_genes.csv\n\n")

# Extract and save direct epigenetic targets (most important)
direct_targets <- integrated[integrated$regulatory_class == "Direct epigenetic", ]
direct_targets <- direct_targets[order(-abs(direct_targets$log2FC_expression)), ]
write.csv(direct_targets, "direct_epigenetic_targets.csv", row.names = FALSE)
cat(paste("Saved:", nrow(direct_targets), "direct epigenetic targets -> direct_epigenetic_targets.csv\n\n"))

################################################################################
# 5. Summary Statistics
################################################################################

cat("========================================\n")
cat("Integration Summary\n")
cat("========================================\n\n")

cat(paste("Total genes analyzed:", nrow(integrated), "\n\n"))

cat("Regulatory categories:\n")
print(table(integrated$regulatory_class))
cat("\n")

if (nrow(direct_targets) > 0) {
    cat("Direct epigenetic targets breakdown:\n")
    print(table(direct_targets$epigenetic_subclass))
    cat("\n")
}

# Overlap statistics
cat("Overlap statistics:\n")
cat(paste("  TES binding only:", sum(integrated$has_TES_binding & !integrated$has_promoter_DMR & !integrated$is_DEG), "\n"))
cat(paste("  Methylation only:", sum(!integrated$has_TES_binding & integrated$has_promoter_DMR & !integrated$is_DEG), "\n"))
cat(paste("  Expression only:", sum(!integrated$has_TES_binding & !integrated$has_promoter_DMR & integrated$is_DEG), "\n"))
cat(paste("  Binding + Methylation:", sum(integrated$has_TES_binding & integrated$has_promoter_DMR), "\n"))
cat(paste("  Binding + Expression:", sum(integrated$has_TES_binding & integrated$is_DEG), "\n"))
cat(paste("  Methylation + Expression:", sum(integrated$has_promoter_DMR & integrated$is_DEG), "\n"))
cat(paste("  All three (direct epigenetic):", nrow(direct_targets), "\n\n"))

################################################################################
# 6. Visualizations
################################################################################

cat("========================================\n")
cat("Generating Visualizations\n")
cat("========================================\n\n")

# Venn diagram
cat("Creating Venn diagram...\n")
venn_list <- list(
    "TES Binding\n(promoter)" = integrated$gene[integrated$has_TES_binding],
    "Promoter\nMethylation" = integrated$gene[integrated$has_promoter_DMR],
    "Expression\nChange" = integrated$gene[integrated$is_DEG]
)

venn.diagram(
    x = venn_list,
    filename = "binding_methylation_expression_venn.png",
    output = TRUE,
    imagetype = "png",
    height = 3000,
    width = 3000,
    resolution = 300,
    col = c("#E41A1C", "#377EB8", "#4DAF4A"),
    fill = c(alpha("#E41A1C", 0.3), alpha("#377EB8", 0.3), alpha("#4DAF4A", 0.3)),
    cat.col = c("#E41A1C", "#377EB8", "#4DAF4A"),
    cat.cex = 1.5,
    cex = 1.5,
    margin = 0.1,
    main = "Multi-Omics Integration",
    main.cex = 2
)
cat("Saved: binding_methylation_expression_venn.png\n\n")

# Heatmap of top direct epigenetic targets
if (nrow(direct_targets) >= 5) {
    cat("Creating heatmap of top direct epigenetic targets...\n")

    # Select top 50 by expression fold change
    top_targets <- head(direct_targets, 50)

    # Create matrix for heatmap
    heatmap_mat <- data.frame(
        TES_binding = ifelse(top_targets$has_TES_binding, 1, 0),
        Hypermethylation = ifelse(top_targets$is_hypermethylated, 1,
            ifelse(top_targets$is_hypomethylated, -1, 0)
        ),
        Expression_FC = top_targets$log2FC_expression / max(abs(top_targets$log2FC_expression), na.rm = TRUE)
    )
    rownames(heatmap_mat) <- top_targets$gene

    pdf("direct_epigenetic_targets_heatmap.pdf", width = 8, height = 12)
    pheatmap(as.matrix(heatmap_mat),
        cluster_cols = FALSE,
        cluster_rows = TRUE,
        show_rownames = TRUE,
        color = colorRampPalette(c("blue", "white", "red"))(100),
        main = "Direct Epigenetic Targets (Top 50 by |FC|)",
        fontsize_row = 6,
        angle_col = 45
    )
    dev.off()
    cat("Saved: direct_epigenetic_targets_heatmap.pdf\n\n")
}

################################################################################
# 7. Hypothesis Testing
################################################################################

cat("========================================\n")
cat("Hypothesis Testing\n")
cat("========================================\n\n")

cat("Testing Hypothesis: TES recruits methylation to TEAD1 target genes\n\n")

# Test 1: Do TES-bound genes show promoter hypermethylation?
tes_bound <- integrated[integrated$has_TES_binding, ]
tes_hyper_rate <- sum(tes_bound$is_hypermethylated) / nrow(tes_bound)
cat(paste(
    "1. TES-bound genes with hypermethylation:",
    sum(tes_bound$is_hypermethylated), "/", nrow(tes_bound),
    "(", round(100 * tes_hyper_rate, 1), "%)\n"
))

# Test 2: Are hypermethylated genes downregulated?
hyper_genes <- integrated[integrated$is_hypermethylated, ]
hyper_down_rate <- sum(hyper_genes$is_DEG & hyper_genes$log2FC_expression < 0, na.rm = TRUE) /
    sum(hyper_genes$is_DEG, na.rm = TRUE)
cat(paste(
    "2. Hypermethylated DEGs that are downregulated:",
    sum(hyper_genes$is_DEG & hyper_genes$log2FC_expression < 0, na.rm = TRUE), "/",
    sum(hyper_genes$is_DEG, na.rm = TRUE),
    "(", round(100 * hyper_down_rate, 1), "%)\n"
))

# Test 3: Is methylation TEAD1-dependent? (TES-specific, not TESmut)
tes_specific_dmr <- setdiff(dmr_genes_all, dmr_genes_tesmut)
cat(paste(
    "3. TES-specific DMR genes (not in TESmut):",
    length(tes_specific_dmr), "/", length(dmr_genes_all),
    "(", round(100 * length(tes_specific_dmr) / length(dmr_genes_all), 1), "%)\n\n"
))

# Test 4: Direct epigenetic targets match expected pattern?
expected_pattern <- direct_targets %>%
    filter(is_hypermethylated & log2FC_expression < 0)
cat(paste("4. Direct epigenetic targets matching expected pattern:\n"))
cat(paste(
    "   (Hypermethylated + Downregulated):",
    nrow(expected_pattern), "/", nrow(direct_targets),
    "(", round(100 * nrow(expected_pattern) / nrow(direct_targets), 1), "%)\n\n"
))

################################################################################
# 8. Save Summary Report
################################################################################

sink("integration_summary.txt")
cat("====================================================\n")
cat("Multi-Omics Integration Summary Report\n")
cat("====================================================\n")
cat(paste("Generated:", Sys.time(), "\n\n"))

cat("INPUT DATA:\n")
cat(paste("  RNA-seq DEGs:", sum(rna_genes$is_DEG), "\n"))
cat(paste("  TES binding (promoters):", length(tes_promoter_genes), "\n"))
cat(paste("  TESmut binding (promoters):", length(tesmut_promoter_genes), "\n"))
cat(paste("  TEAD1 binding (promoters):", length(tead1_promoter_genes), "\n"))
cat(paste("  Promoter DMRs (TES vs GFP):", length(dmr_genes_all), "\n"))
cat(paste("    - Hypermethylated:", length(dmr_genes_hyper), "\n"))
cat(paste("    - Hypomethylated:", length(dmr_genes_hypo), "\n\n"))

cat("INTEGRATION RESULTS:\n")
cat(paste("  Total genes analyzed:", nrow(integrated), "\n\n"))
cat("  Regulatory categories:\n")
print(table(integrated$regulatory_class))
cat("\n")

cat("HYPOTHESIS TESTING:\n")
cat(paste("  1. TES-bound genes with hypermethylation: ", round(100 * tes_hyper_rate, 1), "%\n", sep = ""))
cat(paste("  2. Hypermethylated DEGs downregulated: ", round(100 * hyper_down_rate, 1), "%\n", sep = ""))
cat(paste("  3. TES-specific DMR genes: ", round(100 * length(tes_specific_dmr) / length(dmr_genes_all), 1), "%\n", sep = ""))
if (nrow(direct_targets) > 0) {
    cat(paste("  4. Direct targets matching pattern: ", round(100 * nrow(expected_pattern) / nrow(direct_targets), 1), "%\n\n", sep = ""))
}

cat("KEY FINDINGS:\n")
if (nrow(expected_pattern) > 0) {
    cat("  ✓ Direct epigenetic regulation detected\n")
    cat("  ✓ Pattern: TES binding → Hypermethylation → Repression\n")
    cat(paste("  ✓", nrow(expected_pattern), "genes follow expected mechanism\n\n"))
} else {
    cat("  ⚠ Few direct epigenetic targets detected\n")
    cat("  ⚠ Methylation may be indirect or require further investigation\n\n")
}

cat("TOP 10 DIRECT EPIGENETIC TARGETS:\n")
if (nrow(direct_targets) > 0) {
    top10 <- head(direct_targets[, c(
        "gene", "log2FC_expression", "padj_expression",
        "is_hypermethylated", "has_TES_binding"
    )], 10)
    print(top10, row.names = FALSE)
}

sink()

cat("Saved: integration_summary.txt\n\n")

################################################################################
# Final Output
################################################################################

cat("====================================================\n")
cat("Integration Complete!\n")
cat("====================================================\n")
cat(paste("End:", Sys.time(), "\n\n"))

cat("Output files:\n")
cat("  1. integrated_genes.csv - Full integrated dataset\n")
cat("  2. direct_epigenetic_targets.csv - Direct TES targets with methylation\n")
cat("  3. binding_methylation_expression_venn.png - Overlap visualization\n")
cat("  4. direct_epigenetic_targets_heatmap.pdf - Heatmap of top targets\n")
cat("  5. integration_summary.txt - Statistical summary\n\n")

cat("Next steps:\n")
cat("  1. Review direct_epigenetic_targets.csv for candidate genes\n")
cat("  2. Validate top targets with bisulfite sequencing\n")
cat("  3. Check if pattern is TEAD1-dependent (compare with TESmut)\n")
cat("  4. Perform pathway enrichment on direct epigenetic targets\n\n")

if (nrow(direct_targets) == 0) {
    cat("WARNING: No direct epigenetic targets found!\n")
    cat("This could mean:\n")
    cat("  - Methylation changes are indirect (not at TES binding sites)\n")
    cat("  - Need more sensitive thresholds (check FDR and FC cutoffs)\n")
    cat("  - TES mechanism is non-epigenetic (direct transcriptional regulation)\n\n")
}
