#!/usr/bin/env Rscript

################################################################################
# Script: integrative_TES_TEAD1.R
# Purpose: Multi-omics integration focusing on TES vs TEAD1 comparison
#
# Description:
#   Integrates DNA methylation (meDIP), TF binding (Cut&Tag), and gene
#   expression (RNA-seq) to classify genes by regulatory mechanisms.
#   This version focuses exclusively on TES vs TEAD1 comparison,
#   excluding TESmut samples from the analysis.
#
# Integration Strategy:
#   1. Load annotated DMRs from meDIP analysis (step 08 output)
#   2. Load TES and TEAD1 binding peaks from Cut&Tag
#   3. Load differential expression from RNA-seq
#   4. Identify genes with:
#      - Promoter DMRs (±2kb from TSS)
#      - TF binding at promoters (TES and/or TEAD1)
#      - Expression changes
#   5. Classify genes into regulatory categories
#
# Biological Hypothesis:
#   TES binding at TEAD1 target sites recruits DNA methyltransferases,
#   causing promoter hypermethylation and transcriptional repression.
#
# Input:
#   - meDIP annotated DMRs: ../results/08_annotation/*_annotated.csv
#   - Cut&Tag peaks: ../../SRF_Eva_CUTandTAG/results/11_combined_replicates_narrow/peaks/
#   - RNA-seq DEGs: ../../SRF_Eva_RNA/results/05_deseq2/deseq2_results_TES_vs_GFP.txt
#
# Output:
#   - Integrated gene list: integrated_genes_TES_TEAD1.csv
#   - Direct epigenetic targets: direct_epigenetic_targets_TES_TEAD1.csv
#   - Overlap Venn diagram: TES_TEAD1_binding_methylation_expression_venn.png
#   - Heatmap: TES_TEAD1_regulation_heatmap.pdf
#   - Summary statistics: TES_TEAD1_integration_summary.txt
#
################################################################################

cat("====================================================\n")
cat("Multi-Omics Integration: TES vs TEAD1 Analysis\n")
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
    library(ComplexHeatmap)
    library(circlize)
    library(VennDiagram)
    library(dplyr)
})

# Define paths
base_dir <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top"
medip_annotation_dir <- file.path(base_dir, "meDIP/results/08_annotation")
medip_dmr_dir <- file.path(base_dir, "meDIP/results/07_differential_MEDIPS")
cuttag_dir <- file.path(base_dir, "SRF_Eva_CUTandTAG/results")
rnaseq_dir <- file.path(base_dir, "SRF_Eva_RNA/results/05_deseq2")
out_dir <- file.path(base_dir, "SRF_Eva_integrated_analysis/scripts/integrative_TES_TEAD1")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
setwd(out_dir)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

cat("Loading datasets (TES and TEAD1 only)...\n\n")

################################################################################
# 1. Load RNA-seq differential expression
################################################################################

cat("========================================\n")
cat("1. Loading RNA-seq Data\n")
cat("========================================\n")

rnaseq_file <- file.path(rnaseq_dir, "deseq2_results_TES_vs_GFP.txt")
if (file.exists(rnaseq_file)) {
    degs <- read.delim(rnaseq_file, stringsAsFactors = FALSE)
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
# 2. Load Cut&Tag TES and TEAD1 binding peaks (excluding TESmut)
################################################################################

cat("========================================\n")
cat("2. Loading Cut&Tag Binding Peaks\n")
cat("========================================\n")

# Load TES peaks
tes_peak_file <- file.path(cuttag_dir, "11_combined_replicates_narrow/peaks/TES_combined_peaks.narrowPeak")
if (file.exists(tes_peak_file)) {
    tes_peaks <- import(tes_peak_file, format = "narrowPeak")
    # Convert chromosome names from Ensembl (1,2,3) to UCSC (chr1,chr2,chr3)
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

# Load TEAD1 peaks (endogenous TEAD1)
tead1_peak_file <- file.path(cuttag_dir, "11_combined_replicates_narrow/peaks/TEAD1_combined_peaks.narrowPeak")
if (file.exists(tead1_peak_file)) {
    tead1_peaks <- import(tead1_peak_file, format = "narrowPeak")
    # Convert chromosome names
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

# Load annotated DMRs from step 08 (preferred)
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
    cat(paste("  Hypomethylated genes:", length(dmr_genes_hypo), "\n\n"))
} else {
    # Fallback: Load raw MEDIPS output and filter for promoter DMRs
    cat("Annotated DMRs not found, loading raw MEDIPS output...\n")
    dmr_file <- file.path(medip_dmr_dir, "TES_vs_GFP_DMRs_FDR05_FC2.csv")

    if (file.exists(dmr_file)) {
        dmrs <- read.csv(dmr_file, stringsAsFactors = FALSE)

        # Convert to GRanges and annotate
        dmrs_gr <- makeGRangesFromDataFrame(dmrs,
            seqnames.field = "chr",
            start.field = "start",
            end.field = "stop",
            keep.extra.columns = TRUE
        )

        dmr_anno <- annotatePeak(dmrs_gr, tssRegion = c(-2000, 500), TxDb = txdb, annoDb = "org.Hs.eg.db")
        dmrs_tes <- as.data.frame(dmr_anno)

        # Filter for promoter DMRs
        dmrs_tes <- dmrs_tes[grepl("Promoter", dmrs_tes$annotation), ]

        dmr_genes_all <- unique(dmrs_tes$SYMBOL[!is.na(dmrs_tes$SYMBOL)])
        dmr_genes_hyper <- unique(dmrs_tes$SYMBOL[!is.na(dmrs_tes$SYMBOL) & dmrs_tes$logFC > 1])
        dmr_genes_hypo <- unique(dmrs_tes$SYMBOL[!is.na(dmrs_tes$SYMBOL) & dmrs_tes$logFC < -1])

        cat(paste("  Promoter DMRs (TES vs GFP):", nrow(dmrs_tes), "\n"))
        cat(paste("  Genes with promoter DMRs:", length(dmr_genes_all), "\n"))
        cat(paste("  Hypermethylated:", length(dmr_genes_hyper), "\n"))
        cat(paste("  Hypomethylated:", length(dmr_genes_hypo), "\n\n"))
    } else {
        cat("  ERROR: DMR file not found at:", dmr_file, "\n")
        cat("  Please run MEDIPS analysis (step 07) and annotation (step 08) first.\n\n")
        quit(status = 1)
    }
}

################################################################################
# 4. Integrate Datasets (TES and TEAD1 only)
################################################################################

cat("========================================\n")
cat("4. Integrating Datasets (TES vs TEAD1)\n")
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

# Add TF binding information (TES and TEAD1 only)
integrated$has_TES_binding <- integrated$gene %in% tes_promoter_genes
integrated$has_TEAD1_binding <- integrated$gene %in% tead1_promoter_genes

# Add methylation information
integrated$has_promoter_DMR <- integrated$gene %in% dmr_genes_all
integrated$is_hypermethylated <- integrated$gene %in% dmr_genes_hyper
integrated$is_hypomethylated <- integrated$gene %in% dmr_genes_hypo

# Add expression data
match_idx <- match(integrated$gene, rna_genes$gene)
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
] <- "Direct epigenetic (TES)"

# Direct non-epigenetic: TES binding + expression change (no methylation)
integrated$regulatory_class[
    integrated$has_TES_binding &
        !integrated$has_promoter_DMR &
        integrated$is_DEG
] <- "Direct non-epigenetic (TES)"

# TEAD1-specific regulation
integrated$regulatory_class[
    !integrated$has_TES_binding &
        integrated$has_TEAD1_binding &
        integrated$is_DEG
] <- "TEAD1-specific"

# Co-bound by TES and TEAD1
integrated$regulatory_class[
    integrated$has_TES_binding &
        integrated$has_TEAD1_binding &
        integrated$is_DEG
] <- "TES+TEAD1 co-regulated"

# Indirect methylation: Methylation + expression change (no binding)
integrated$regulatory_class[
    !integrated$has_TES_binding &
        !integrated$has_TEAD1_binding &
        integrated$has_promoter_DMR &
        integrated$is_DEG
] <- "Indirect methylation"

# Indirect regulation: Expression change only
integrated$regulatory_class[
    !integrated$has_TES_binding &
        !integrated$has_TEAD1_binding &
        !integrated$has_promoter_DMR &
        integrated$is_DEG
] <- "Indirect expression"

# Poised regulation: Binding + methylation (no expression change yet)
integrated$regulatory_class[
    (integrated$has_TES_binding | integrated$has_TEAD1_binding) &
        integrated$has_promoter_DMR &
        !integrated$is_DEG
] <- "Poised (bound + methylated)"

# Add detailed subclassification for direct epigenetic targets
integrated$epigenetic_subclass <- NA
direct_epi_idx <- integrated$regulatory_class == "Direct epigenetic (TES)"

if (sum(direct_epi_idx) > 0) {
    integrated$epigenetic_subclass[direct_epi_idx] <- paste0(
        ifelse(integrated$is_hypermethylated[direct_epi_idx], "Hyper", "Hypo"),
        "methylated + ",
        integrated$expression_direction[direct_epi_idx],
        "regulated"
    )
}

# Save integrated results
write.csv(integrated, "integrated_genes_TES_TEAD1.csv", row.names = FALSE)
cat("Saved: integrated_genes_TES_TEAD1.csv\n\n")

# Extract and save direct epigenetic targets
direct_targets <- integrated[integrated$regulatory_class == "Direct epigenetic (TES)", ]
direct_targets <- direct_targets[order(-abs(direct_targets$log2FC_expression)), ]
write.csv(direct_targets, "direct_epigenetic_targets_TES_TEAD1.csv", row.names = FALSE)
cat(paste("Saved:", nrow(direct_targets), "direct epigenetic targets -> direct_epigenetic_targets_TES_TEAD1.csv\n\n"))

################################################################################
# 5. Summary Statistics
################################################################################

cat("========================================\n")
cat("Integration Summary (TES vs TEAD1)\n")
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

# TES vs TEAD1 comparison
cat("TES vs TEAD1 binding comparison:\n")
cat(paste("  TES-only bound genes:", sum(integrated$has_TES_binding & !integrated$has_TEAD1_binding), "\n"))
cat(paste("  TEAD1-only bound genes:", sum(!integrated$has_TES_binding & integrated$has_TEAD1_binding), "\n"))
cat(paste("  Co-bound by TES and TEAD1:", sum(integrated$has_TES_binding & integrated$has_TEAD1_binding), "\n"))
cat(paste("  Neither TES nor TEAD1:", sum(!integrated$has_TES_binding & !integrated$has_TEAD1_binding), "\n\n"))

# Overlap statistics
cat("Overlap statistics:\n")
cat(paste("  TES binding only:", sum(integrated$has_TES_binding & !integrated$has_promoter_DMR & !integrated$is_DEG), "\n"))
cat(paste("  Methylation only:", sum(!integrated$has_TES_binding & integrated$has_promoter_DMR & !integrated$is_DEG), "\n"))
cat(paste("  Expression only:", sum(!integrated$has_TES_binding & !integrated$has_promoter_DMR & integrated$is_DEG), "\n"))
cat(paste("  TES binding + Methylation:", sum(integrated$has_TES_binding & integrated$has_promoter_DMR), "\n"))
cat(paste("  TES binding + Expression:", sum(integrated$has_TES_binding & integrated$is_DEG), "\n"))
cat(paste("  Methylation + Expression:", sum(integrated$has_promoter_DMR & integrated$is_DEG), "\n"))
cat(paste("  All three (direct epigenetic):", nrow(direct_targets), "\n\n"))

################################################################################
# 6. Visualizations
################################################################################

cat("========================================\n")
cat("Generating Visualizations\n")
cat("========================================\n\n")

# Venn diagram - TES binding, methylation, expression
cat("Creating Venn diagram (TES binding + methylation + expression)...\n")
venn_list <- list(
    "TES Binding\n(promoter)" = integrated$gene[integrated$has_TES_binding],
    "Promoter\nMethylation" = integrated$gene[integrated$has_promoter_DMR],
    "Expression\nChange" = integrated$gene[integrated$is_DEG]
)

venn.diagram(
    x = venn_list,
    filename = "TES_TEAD1_binding_methylation_expression_venn.png",
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
    main = "Multi-Omics Integration (TES vs TEAD1)",
    main.cex = 2
)
cat("Saved: TES_TEAD1_binding_methylation_expression_venn.png\n\n")

# Additional Venn diagram - TES vs TEAD1 binding
cat("Creating TES vs TEAD1 binding comparison Venn diagram...\n")
venn_binding <- list(
    "TES\nBinding" = integrated$gene[integrated$has_TES_binding],
    "TEAD1\nBinding" = integrated$gene[integrated$has_TEAD1_binding]
)

venn.diagram(
    x = venn_binding,
    filename = "TES_vs_TEAD1_binding_venn.png",
    output = TRUE,
    imagetype = "png",
    height = 2000,
    width = 2000,
    resolution = 300,
    col = c("#E41A1C", "#984EA3"),
    fill = c(alpha("#E41A1C", 0.3), alpha("#984EA3", 0.3)),
    cat.col = c("#E41A1C", "#984EA3"),
    cat.cex = 1.5,
    cex = 1.5,
    margin = 0.1,
    main = "TES vs TEAD1 Binding Overlap",
    main.cex = 2
)
cat("Saved: TES_vs_TEAD1_binding_venn.png\n\n")

# Heatmap of top direct epigenetic targets using ComplexHeatmap with independent color scales
if (nrow(direct_targets) >= 5) {
    cat("Creating heatmap of top direct epigenetic targets...\n")

    # Select top 50 by expression fold change
    top_targets <- head(direct_targets, 50)
    n_genes <- nrow(top_targets)

    # Prepare data for each column separately
    tes_binding_vec <- ifelse(top_targets$has_TES_binding, 1, 0)
    tead1_binding_vec <- ifelse(top_targets$has_TEAD1_binding, 1, 0)
    methylation_vec <- ifelse(top_targets$is_hypermethylated, 1,
        ifelse(top_targets$is_hypomethylated, -1, 0)
    )
    expression_fc_vec <- top_targets$log2FC_expression
    rownames_vec <- top_targets$gene

    # Create hierarchical clustering based on expression FC
    if (n_genes > 1) {
        hc <- hclust(dist(expression_fc_vec))
        row_order <- hc$order
    } else {
        row_order <- 1
    }

    # Reorder all vectors
    tes_binding_ordered <- tes_binding_vec[row_order]
    tead1_binding_ordered <- tead1_binding_vec[row_order]
    methylation_ordered <- methylation_vec[row_order]
    expression_ordered <- expression_fc_vec[row_order]
    rownames_ordered <- rownames_vec[row_order]

    # Define independent color scales

    # 1. TES Binding column (binary: 0/1) - white/red scale
    col_tes <- colorRamp2(c(0, 1), c("white", "#E41A1C"))
    ht_tes <- Heatmap(
        matrix(tes_binding_ordered, ncol = 1),
        name = "TES\nBinding",
        col = col_tes,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = FALSE,
        column_names_rot = 45,
        column_names_gp = gpar(fontsize = 10),
        column_labels = "TES_binding",
        width = unit(1.5, "cm"),
        heatmap_legend_param = list(
            title = "TES\nBinding",
            at = c(0, 1),
            labels = c("No", "Yes"),
            legend_height = unit(2, "cm")
        )
    )

    # 2. TEAD1 Binding column (binary: 0/1) - white/purple scale
    col_tead1 <- colorRamp2(c(0, 1), c("white", "#984EA3"))
    ht_tead1 <- Heatmap(
        matrix(tead1_binding_ordered, ncol = 1),
        name = "TEAD1\nBinding",
        col = col_tead1,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = FALSE,
        column_names_rot = 45,
        column_names_gp = gpar(fontsize = 10),
        column_labels = "TEAD1_binding",
        width = unit(1.5, "cm"),
        heatmap_legend_param = list(
            title = "TEAD1\nBinding",
            at = c(0, 1),
            labels = c("No", "Yes"),
            legend_height = unit(2, "cm")
        )
    )

    # 3. Methylation column (ternary: -1/0/1) - blue/white/red scale
    col_meth <- colorRamp2(c(-1, 0, 1), c("#377EB8", "white", "#E41A1C"))
    ht_meth <- Heatmap(
        matrix(methylation_ordered, ncol = 1),
        name = "Methylation",
        col = col_meth,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = FALSE,
        column_names_rot = 45,
        column_names_gp = gpar(fontsize = 10),
        column_labels = "Methylation",
        width = unit(1.5, "cm"),
        heatmap_legend_param = list(
            title = "Methylation",
            at = c(-1, 0, 1),
            labels = c("Hypo", "None", "Hyper"),
            legend_height = unit(2, "cm")
        )
    )

    # 4. Expression FC column (continuous) - scale based on actual data range
    fc_range <- max(abs(expression_ordered), na.rm = TRUE)
    col_expr <- colorRamp2(
        c(-fc_range, 0, fc_range),
        c("#2166AC", "white", "#B2182B")
    )
    ht_expr <- Heatmap(
        matrix(expression_ordered, ncol = 1),
        name = "Expression\nlog2FC",
        col = col_expr,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = TRUE,
        row_names_side = "right",
        row_names_gp = gpar(fontsize = 6),
        column_names_rot = 45,
        column_names_gp = gpar(fontsize = 10),
        column_labels = "Expression_FC",
        width = unit(1.5, "cm"),
        row_labels = rownames_ordered,
        heatmap_legend_param = list(
            title = "Expression\nlog2FC",
            legend_height = unit(3, "cm")
        )
    )

    # Combine heatmaps
    ht_list <- ht_tes + ht_tead1 + ht_meth + ht_expr

    png("TES_TEAD1_direct_epigenetic_targets_heatmap.png", width = 10, height = 12, units = "in", res = 300)
    draw(ht_list,
        column_title = "Direct Epigenetic Targets - TES vs TEAD1 (Top 50)",
        column_title_gp = gpar(fontsize = 14, fontface = "bold"),
        heatmap_legend_side = "right",
        padding = unit(c(2, 2, 2, 10), "mm")
    )
    dev.off()
    cat("Saved: TES_TEAD1_direct_epigenetic_targets_heatmap.png\n\n")
}

# Helper function to create heatmap for any gene set
create_multiomics_heatmap <- function(gene_data, output_file, title, max_genes = 50) {
    if (nrow(gene_data) < 5) {
        cat(paste("  Skipping", output_file, "- fewer than 5 genes\n"))
        return(invisible(NULL))
    }

    # Select top genes by expression fold change
    top_data <- head(gene_data[order(-abs(gene_data$log2FC_expression)), ], max_genes)
    n_genes <- nrow(top_data)

    # Prepare data for each column
    tes_binding_vec <- ifelse(top_data$has_TES_binding, 1, 0)
    tead1_binding_vec <- ifelse(top_data$has_TEAD1_binding, 1, 0)
    methylation_vec <- ifelse(top_data$is_hypermethylated, 1,
        ifelse(top_data$is_hypomethylated, -1, 0)
    )
    expression_fc_vec <- top_data$log2FC_expression
    rownames_vec <- top_data$gene

    # Hierarchical clustering based on expression FC
    if (n_genes > 1) {
        hc <- hclust(dist(expression_fc_vec))
        row_order <- hc$order
    } else {
        row_order <- 1
    }

    # Reorder all vectors
    tes_binding_ordered <- tes_binding_vec[row_order]
    tead1_binding_ordered <- tead1_binding_vec[row_order]
    methylation_ordered <- methylation_vec[row_order]
    expression_ordered <- expression_fc_vec[row_order]
    rownames_ordered <- rownames_vec[row_order]

    # Define color scales
    col_tes <- colorRamp2(c(0, 1), c("white", "#E41A1C"))
    col_tead1 <- colorRamp2(c(0, 1), c("white", "#984EA3"))
    col_meth <- colorRamp2(c(-1, 0, 1), c("#377EB8", "white", "#E41A1C"))
    fc_range <- max(abs(expression_ordered), na.rm = TRUE)
    col_expr <- colorRamp2(c(-fc_range, 0, fc_range), c("#2166AC", "white", "#B2182B"))

    # Create heatmaps
    ht_tes <- Heatmap(
        matrix(tes_binding_ordered, ncol = 1),
        name = "TES\nBinding",
        col = col_tes,
        cluster_rows = FALSE, cluster_columns = FALSE,
        show_row_names = FALSE,
        column_names_rot = 45, column_names_gp = gpar(fontsize = 10),
        column_labels = "TES_binding",
        width = unit(1.5, "cm"),
        heatmap_legend_param = list(title = "TES\nBinding", at = c(0, 1), labels = c("No", "Yes"), legend_height = unit(2, "cm"))
    )

    ht_tead1 <- Heatmap(
        matrix(tead1_binding_ordered, ncol = 1),
        name = "TEAD1\nBinding",
        col = col_tead1,
        cluster_rows = FALSE, cluster_columns = FALSE,
        show_row_names = FALSE,
        column_names_rot = 45, column_names_gp = gpar(fontsize = 10),
        column_labels = "TEAD1_binding",
        width = unit(1.5, "cm"),
        heatmap_legend_param = list(title = "TEAD1\nBinding", at = c(0, 1), labels = c("No", "Yes"), legend_height = unit(2, "cm"))
    )

    ht_meth <- Heatmap(
        matrix(methylation_ordered, ncol = 1),
        name = "Methylation",
        col = col_meth,
        cluster_rows = FALSE, cluster_columns = FALSE,
        show_row_names = FALSE,
        column_names_rot = 45, column_names_gp = gpar(fontsize = 10),
        column_labels = "Methylation",
        width = unit(1.5, "cm"),
        heatmap_legend_param = list(title = "Methylation", at = c(-1, 0, 1), labels = c("Hypo", "None", "Hyper"), legend_height = unit(2, "cm"))
    )

    ht_expr <- Heatmap(
        matrix(expression_ordered, ncol = 1),
        name = "Expression\nlog2FC",
        col = col_expr,
        cluster_rows = FALSE, cluster_columns = FALSE,
        show_row_names = TRUE,
        row_names_side = "right", row_names_gp = gpar(fontsize = 6),
        column_names_rot = 45, column_names_gp = gpar(fontsize = 10),
        column_labels = "Expression_FC",
        width = unit(1.5, "cm"),
        row_labels = rownames_ordered,
        heatmap_legend_param = list(title = "Expression\nlog2FC", legend_height = unit(3, "cm"))
    )

    # Combine and save
    ht_list <- ht_tes + ht_tead1 + ht_meth + ht_expr

    png(output_file, width = 10, height = 12, units = "in", res = 300)
    draw(ht_list,
        column_title = title,
        column_title_gp = gpar(fontsize = 14, fontface = "bold"),
        heatmap_legend_side = "right",
        padding = unit(c(2, 2, 2, 10), "mm")
    )
    dev.off()
    cat(paste("Saved:", output_file, "\n"))
}

# Heatmap 2: Co-bound targets (TES + TEAD1)
cat("Creating heatmap of co-bound targets (TES + TEAD1)...\n")
cobound_targets <- integrated[integrated$has_TES_binding & integrated$has_TEAD1_binding & integrated$is_DEG, ]
cobound_targets <- cobound_targets[order(-abs(cobound_targets$log2FC_expression)), ]
if (nrow(cobound_targets) >= 5) {
    create_multiomics_heatmap(
        cobound_targets,
        "TES_TEAD1_cobound_targets_heatmap.png",
        "Co-bound Targets (TES + TEAD1 Binding + DEG, Top 50)"
    )
} else {
    cat("  Skipping co-bound heatmap - fewer than 5 genes\n")
}
cat("\n")

# Heatmap 3: Top upregulated DEGs
cat("Creating heatmap of top upregulated DEGs...\n")
upregulated_degs <- integrated[integrated$is_DEG & !is.na(integrated$log2FC_expression) & integrated$log2FC_expression > 0, ]
upregulated_degs <- upregulated_degs[order(-upregulated_degs$log2FC_expression), ]
if (nrow(upregulated_degs) >= 5) {
    create_multiomics_heatmap(
        upregulated_degs,
        "TES_TEAD1_top_upregulated_heatmap.png",
        "Top Upregulated DEGs (TES vs GFP, Top 50)"
    )
} else {
    cat("  Skipping upregulated heatmap - fewer than 5 genes\n")
}
cat("\n")

# Heatmap 4: Top downregulated DEGs
cat("Creating heatmap of top downregulated DEGs...\n")
downregulated_degs <- integrated[integrated$is_DEG & !is.na(integrated$log2FC_expression) & integrated$log2FC_expression < 0, ]
downregulated_degs <- downregulated_degs[order(downregulated_degs$log2FC_expression), ]
if (nrow(downregulated_degs) >= 5) {
    create_multiomics_heatmap(
        downregulated_degs,
        "TES_TEAD1_top_downregulated_heatmap.png",
        "Top Downregulated DEGs (TES vs GFP, Top 50)"
    )
} else {
    cat("  Skipping downregulated heatmap - fewer than 5 genes\n")
}
cat("\n")

################################################################################
# 7. Hypothesis Testing (TES vs TEAD1 focus)
################################################################################

cat("========================================\n")
cat("Hypothesis Testing (TES vs TEAD1)\n")
cat("========================================\n\n")

cat("Testing: Does TES recruit methylation to TEAD1-bound genes?\n\n")

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

# Test 3: Do TES and TEAD1 co-localize at promoters?
cobound_genes <- integrated[integrated$has_TES_binding & integrated$has_TEAD1_binding, ]
cobound_rate <- nrow(cobound_genes) / sum(integrated$has_TES_binding | integrated$has_TEAD1_binding)
cat(paste(
    "3. Genes co-bound by TES and TEAD1:",
    nrow(cobound_genes), "/", sum(integrated$has_TES_binding | integrated$has_TEAD1_binding),
    "(", round(100 * cobound_rate, 1), "%)\n"
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

sink("TES_TEAD1_integration_summary.txt")
cat("====================================================\n")
cat("Multi-Omics Integration Summary: TES vs TEAD1\n")
cat("====================================================\n")
cat(paste("Generated:", Sys.time(), "\n\n"))

cat("INPUT DATA:\n")
cat(paste("  RNA-seq DEGs:", sum(rna_genes$is_DEG), "\n"))
cat(paste("  TES binding (promoters):", length(tes_promoter_genes), "\n"))
cat(paste("  TEAD1 binding (promoters):", length(tead1_promoter_genes), "\n"))
cat(paste("  Promoter DMRs (TES vs GFP):", length(dmr_genes_all), "\n"))
cat(paste("    - Hypermethylated:", length(dmr_genes_hyper), "\n"))
cat(paste("    - Hypomethylated:", length(dmr_genes_hypo), "\n\n"))

cat("INTEGRATION RESULTS:\n")
cat(paste("  Total genes analyzed:", nrow(integrated), "\n\n"))
cat("  Regulatory categories:\n")
print(table(integrated$regulatory_class))
cat("\n")

cat("TES vs TEAD1 COMPARISON:\n")
cat(paste("  TES-only bound:", sum(integrated$has_TES_binding & !integrated$has_TEAD1_binding), "\n"))
cat(paste("  TEAD1-only bound:", sum(!integrated$has_TES_binding & integrated$has_TEAD1_binding), "\n"))
cat(paste("  Co-bound by both:", sum(integrated$has_TES_binding & integrated$has_TEAD1_binding), "\n\n"))

cat("HYPOTHESIS TESTING:\n")
cat(paste("  1. TES-bound genes with hypermethylation: ", round(100 * tes_hyper_rate, 1), "%\n", sep = ""))
cat(paste("  2. Hypermethylated DEGs downregulated: ", round(100 * hyper_down_rate, 1), "%\n", sep = ""))
cat(paste("  3. TES-TEAD1 co-binding rate: ", round(100 * cobound_rate, 1), "%\n", sep = ""))
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
        "is_hypermethylated", "has_TES_binding", "has_TEAD1_binding"
    )], 10)
    print(top10, row.names = FALSE)
}

sink()

cat("Saved: TES_TEAD1_integration_summary.txt\n\n")

################################################################################
# Final Output
################################################################################

cat("====================================================\n")
cat("Integration Complete (TES vs TEAD1 Analysis)!\n")
cat("====================================================\n")
cat(paste("End:", Sys.time(), "\n\n"))

cat("Output files:\n")
cat("  1. integrated_genes_TES_TEAD1.csv - Full integrated dataset\n")
cat("  2. direct_epigenetic_targets_TES_TEAD1.csv - Direct TES targets with methylation\n")
cat("  3. TES_TEAD1_binding_methylation_expression_venn.png - Overlap visualization\n")
cat("  4. TES_vs_TEAD1_binding_venn.png - TES vs TEAD1 binding comparison\n")
cat("  5. TES_TEAD1_direct_epigenetic_targets_heatmap.png - Heatmap of direct epigenetic targets\n")
cat("  6. TES_TEAD1_cobound_targets_heatmap.png - Heatmap of TES+TEAD1 co-bound DEGs\n")
cat("  7. TES_TEAD1_top_upregulated_heatmap.png - Heatmap of top upregulated DEGs\n")
cat("  8. TES_TEAD1_top_downregulated_heatmap.png - Heatmap of top downregulated DEGs\n")
cat("  9. TES_TEAD1_integration_summary.txt - Statistical summary\n\n")

cat("Next steps:\n")
cat("  1. Review direct_epigenetic_targets_TES_TEAD1.csv for candidate genes\n")
cat("  2. Compare TES-specific vs TEAD1-specific vs co-bound targets\n")
cat("  3. Validate top targets with bisulfite sequencing\n")
cat("  4. Perform pathway enrichment on regulatory categories\n\n")

if (nrow(direct_targets) == 0) {
    cat("WARNING: No direct epigenetic targets found!\n")
    cat("This could mean:\n")
    cat("  - Methylation changes are indirect (not at TES binding sites)\n")
    cat("  - Need more sensitive thresholds (check FDR and FC cutoffs)\n")
    cat("  - TES mechanism is non-epigenetic (direct transcriptional regulation)\n\n")
}
