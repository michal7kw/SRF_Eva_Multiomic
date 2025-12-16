#!/usr/bin/env Rscript
#
# PEAK-DMR MAPPING ANALYSIS
# Maps TES/TEAD1 binding peaks to DMRs
# Outputs: Positions, distances, gene annotations, overlaps
#

suppressPackageStartupMessages({
    library(GenomicRanges)
    library(rtracklayer)
    library(ChIPseeker)
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    library(org.Hs.eg.db)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(VennDiagram)
    library(RColorBrewer)
    library(pheatmap)
    library(gridExtra)
})

setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_1")

cat("==========================================================\n")
cat("PEAK-DMR MAPPING ANALYSIS\n")
cat("==========================================================\n")
cat("Analysis started:", as.character(Sys.time()), "\n\n")

# =============================================================================
# PATH CONFIGURATION
# =============================================================================

# Input files
TES_PEAKS <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_CUTandTAG/results/11_combined_replicates_narrow/peaks/TES_consensus_peaks.bed"
TEAD1_PEAKS <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_CUTandTAG/results/11_combined_replicates_narrow/peaks/TEAD1_consensus_peaks.bed"
DMR_FILE <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP/results/07_differential_MEDIPS/TES_vs_GFP_DMRs_FDR05_FC2.csv"

# Output directory
OUTPUT_BASE <- "output/22_peak_DMR_mapping"
dir.create(OUTPUT_BASE, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTPUT_BASE, "tables"), showWarnings = FALSE)
dir.create(file.path(OUTPUT_BASE, "plots"), showWarnings = FALSE)

# Parameters
OVERLAP_DISTANCE <- 0  # For direct overlaps
PROXIMITY_DISTANCE <- 5000  # 5kb for nearby features

# =============================================================================
# PHASE 1: LOAD DATA
# =============================================================================

cat("=== PHASE 1: Loading Data ===\n")

# Load TES peaks
tes_peaks_df <- read.table(TES_PEAKS, header = FALSE, stringsAsFactors = FALSE,
                            col.names = c("chr", "start", "end", "replicate_count", "replicates"))
tes_peaks_df$chr <- paste0("chr", tes_peaks_df$chr)
tes_peaks <- makeGRangesFromDataFrame(tes_peaks_df, keep.extra.columns = TRUE)
cat(sprintf("  TES peaks loaded: %d\n", length(tes_peaks)))

# Load TEAD1 peaks
tead1_peaks_df <- read.table(TEAD1_PEAKS, header = FALSE, stringsAsFactors = FALSE,
                              col.names = c("chr", "start", "end", "replicate_count", "replicates"))
tead1_peaks_df$chr <- paste0("chr", tead1_peaks_df$chr)
tead1_peaks <- makeGRangesFromDataFrame(tead1_peaks_df, keep.extra.columns = TRUE)
cat(sprintf("  TEAD1 peaks loaded: %d\n", length(tead1_peaks)))

# Load DMRs
dmr_df <- read.csv(DMR_FILE, stringsAsFactors = FALSE)
dmr_df$chr <- paste0("chr", dmr_df$chr)
dmr_df$direction <- ifelse(dmr_df$logFC > 0, "Hyper", "Hypo")
dmrs <- makeGRangesFromDataFrame(dmr_df,
                                  seqnames.field = "chr",
                                  start.field = "start",
                                  end.field = "stop",
                                  keep.extra.columns = TRUE)
cat(sprintf("  DMRs loaded: %d (Hyper: %d, Hypo: %d)\n",
            length(dmrs), sum(dmr_df$direction == "Hyper"), sum(dmr_df$direction == "Hypo")))

cat("\n")

# =============================================================================
# PHASE 2: FIND OVERLAPS
# =============================================================================

cat("=== PHASE 2: Finding Overlaps ===\n")

# Direct overlaps
tes_dmr_overlaps <- findOverlaps(tes_peaks, dmrs)
tead1_dmr_overlaps <- findOverlaps(tead1_peaks, dmrs)

cat(sprintf("  TES peaks overlapping DMRs: %d (%.1f%%)\n",
            length(unique(queryHits(tes_dmr_overlaps))),
            100 * length(unique(queryHits(tes_dmr_overlaps))) / length(tes_peaks)))
cat(sprintf("  TEAD1 peaks overlapping DMRs: %d (%.1f%%)\n",
            length(unique(queryHits(tead1_dmr_overlaps))),
            100 * length(unique(queryHits(tead1_dmr_overlaps))) / length(tead1_peaks)))
cat(sprintf("  DMRs overlapping TES peaks: %d (%.1f%%)\n",
            length(unique(subjectHits(tes_dmr_overlaps))),
            100 * length(unique(subjectHits(tes_dmr_overlaps))) / length(dmrs)))
cat(sprintf("  DMRs overlapping TEAD1 peaks: %d (%.1f%%)\n",
            length(unique(subjectHits(tead1_dmr_overlaps))),
            100 * length(unique(subjectHits(tead1_dmr_overlaps))) / length(dmrs)))

# Proximity analysis (within 5kb)
tes_peaks_extended <- resize(tes_peaks, width = width(tes_peaks) + 2 * PROXIMITY_DISTANCE, fix = "center")
tead1_peaks_extended <- resize(tead1_peaks, width = width(tead1_peaks) + 2 * PROXIMITY_DISTANCE, fix = "center")

tes_dmr_nearby <- findOverlaps(tes_peaks_extended, dmrs)
tead1_dmr_nearby <- findOverlaps(tead1_peaks_extended, dmrs)

cat(sprintf("\n  TES peaks with DMR within 5kb: %d (%.1f%%)\n",
            length(unique(queryHits(tes_dmr_nearby))),
            100 * length(unique(queryHits(tes_dmr_nearby))) / length(tes_peaks)))
cat(sprintf("  TEAD1 peaks with DMR within 5kb: %d (%.1f%%)\n",
            length(unique(queryHits(tead1_dmr_nearby))),
            100 * length(unique(queryHits(tead1_dmr_nearby))) / length(tead1_peaks)))

cat("\n")

# =============================================================================
# PHASE 3: CALCULATE DISTANCES
# =============================================================================

cat("=== PHASE 3: Calculating Distances ===\n")

# Find nearest DMR for each peak
tes_nearest <- distanceToNearest(tes_peaks, dmrs)
tead1_nearest <- distanceToNearest(tead1_peaks, dmrs)

# Add distance info to peaks
tes_peaks$nearest_dmr_idx <- NA
tes_peaks$nearest_dmr_distance <- NA
tes_peaks$nearest_dmr_idx[queryHits(tes_nearest)] <- subjectHits(tes_nearest)
tes_peaks$nearest_dmr_distance[queryHits(tes_nearest)] <- mcols(tes_nearest)$distance

tead1_peaks$nearest_dmr_idx <- NA
tead1_peaks$nearest_dmr_distance <- NA
tead1_peaks$nearest_dmr_idx[queryHits(tead1_nearest)] <- subjectHits(tead1_nearest)
tead1_peaks$nearest_dmr_distance[queryHits(tead1_nearest)] <- mcols(tead1_nearest)$distance

cat(sprintf("  Median distance from TES peak to nearest DMR: %.0f bp\n",
            median(tes_peaks$nearest_dmr_distance, na.rm = TRUE)))
cat(sprintf("  Median distance from TEAD1 peak to nearest DMR: %.0f bp\n",
            median(tead1_peaks$nearest_dmr_distance, na.rm = TRUE)))

cat("\n")

# =============================================================================
# PHASE 4: GENE ANNOTATION
# =============================================================================

cat("=== PHASE 4: Annotating to Genes ===\n")

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Annotate TES peaks
cat("  Annotating TES peaks...\n")
tes_anno <- annotatePeak(tes_peaks, TxDb = txdb, annoDb = "org.Hs.eg.db", verbose = FALSE)
tes_anno_df <- as.data.frame(tes_anno)

# Annotate TEAD1 peaks
cat("  Annotating TEAD1 peaks...\n")
tead1_anno <- annotatePeak(tead1_peaks, TxDb = txdb, annoDb = "org.Hs.eg.db", verbose = FALSE)
tead1_anno_df <- as.data.frame(tead1_anno)

# Annotate DMRs
cat("  Annotating DMRs...\n")
dmr_anno <- annotatePeak(dmrs, TxDb = txdb, annoDb = "org.Hs.eg.db", verbose = FALSE)
dmr_anno_df <- as.data.frame(dmr_anno)

cat(sprintf("  TES peaks annotated to %d unique genes\n", length(unique(na.omit(tes_anno_df$SYMBOL)))))
cat(sprintf("  TEAD1 peaks annotated to %d unique genes\n", length(unique(na.omit(tead1_anno_df$SYMBOL)))))
cat(sprintf("  DMRs annotated to %d unique genes\n", length(unique(na.omit(dmr_anno_df$SYMBOL)))))

cat("\n")

# =============================================================================
# PHASE 5: CREATE COMPREHENSIVE TABLES
# =============================================================================

cat("=== PHASE 5: Creating Comprehensive Tables ===\n")

# Table 1: TES peaks with DMR information
tes_table <- data.frame(
    peak_id = paste0("TES_peak_", 1:length(tes_peaks)),
    chr = as.character(seqnames(tes_peaks)),
    start = start(tes_peaks),
    end = end(tes_peaks),
    width = width(tes_peaks),
    replicate_count = tes_peaks$replicate_count,
    gene = tes_anno_df$SYMBOL,
    gene_id = tes_anno_df$geneId,
    annotation = tes_anno_df$annotation,
    distance_to_TSS = tes_anno_df$distanceToTSS,
    nearest_DMR_distance = tes_peaks$nearest_dmr_distance,
    has_overlapping_DMR = 1:length(tes_peaks) %in% queryHits(tes_dmr_overlaps),
    has_nearby_DMR_5kb = 1:length(tes_peaks) %in% queryHits(tes_dmr_nearby),
    stringsAsFactors = FALSE
)

# Add DMR info for overlapping peaks
tes_table$overlapping_DMR_logFC <- NA
tes_table$overlapping_DMR_direction <- NA
for (i in which(tes_table$has_overlapping_DMR)) {
    dmr_idx <- subjectHits(tes_dmr_overlaps)[queryHits(tes_dmr_overlaps) == i][1]
    tes_table$overlapping_DMR_logFC[i] <- dmrs$logFC[dmr_idx]
    tes_table$overlapping_DMR_direction[i] <- dmrs$direction[dmr_idx]
}

# Table 2: TEAD1 peaks with DMR information
tead1_table <- data.frame(
    peak_id = paste0("TEAD1_peak_", 1:length(tead1_peaks)),
    chr = as.character(seqnames(tead1_peaks)),
    start = start(tead1_peaks),
    end = end(tead1_peaks),
    width = width(tead1_peaks),
    replicate_count = tead1_peaks$replicate_count,
    gene = tead1_anno_df$SYMBOL,
    gene_id = tead1_anno_df$geneId,
    annotation = tead1_anno_df$annotation,
    distance_to_TSS = tead1_anno_df$distanceToTSS,
    nearest_DMR_distance = tead1_peaks$nearest_dmr_distance,
    has_overlapping_DMR = 1:length(tead1_peaks) %in% queryHits(tead1_dmr_overlaps),
    has_nearby_DMR_5kb = 1:length(tead1_peaks) %in% queryHits(tead1_dmr_nearby),
    stringsAsFactors = FALSE
)

# Add DMR info for overlapping peaks
tead1_table$overlapping_DMR_logFC <- NA
tead1_table$overlapping_DMR_direction <- NA
for (i in which(tead1_table$has_overlapping_DMR)) {
    dmr_idx <- subjectHits(tead1_dmr_overlaps)[queryHits(tead1_dmr_overlaps) == i][1]
    tead1_table$overlapping_DMR_logFC[i] <- dmrs$logFC[dmr_idx]
    tead1_table$overlapping_DMR_direction[i] <- dmrs$direction[dmr_idx]
}

# Table 3: DMRs with peak information
dmr_table <- data.frame(
    dmr_id = paste0("DMR_", 1:length(dmrs)),
    chr = as.character(seqnames(dmrs)),
    start = start(dmrs),
    end = end(dmrs),
    width = width(dmrs),
    logFC = dmrs$logFC,
    FDR = dmrs$FDR,
    direction = dmrs$direction,
    CpG_count = dmrs$CpG_count,
    gene = dmr_anno_df$SYMBOL,
    gene_id = dmr_anno_df$geneId,
    annotation = dmr_anno_df$annotation,
    distance_to_TSS = dmr_anno_df$distanceToTSS,
    has_TES_peak = 1:length(dmrs) %in% subjectHits(tes_dmr_overlaps),
    has_TEAD1_peak = 1:length(dmrs) %in% subjectHits(tead1_dmr_overlaps),
    stringsAsFactors = FALSE
)

# Find nearest peaks for DMRs
dmr_to_tes <- distanceToNearest(dmrs, tes_peaks)
dmr_to_tead1 <- distanceToNearest(dmrs, tead1_peaks)

dmr_table$nearest_TES_peak_distance <- NA
dmr_table$nearest_TEAD1_peak_distance <- NA
dmr_table$nearest_TES_peak_distance[queryHits(dmr_to_tes)] <- mcols(dmr_to_tes)$distance
dmr_table$nearest_TEAD1_peak_distance[queryHits(dmr_to_tead1)] <- mcols(dmr_to_tead1)$distance

# Write tables
write.csv(tes_table, file.path(OUTPUT_BASE, "tables", "TES_peaks_with_DMR_info.csv"), row.names = FALSE)
write.csv(tead1_table, file.path(OUTPUT_BASE, "tables", "TEAD1_peaks_with_DMR_info.csv"), row.names = FALSE)
write.csv(dmr_table, file.path(OUTPUT_BASE, "tables", "DMRs_with_peak_info.csv"), row.names = FALSE)

cat("  Created: TES_peaks_with_DMR_info.csv\n")
cat("  Created: TEAD1_peaks_with_DMR_info.csv\n")
cat("  Created: DMRs_with_peak_info.csv\n")

# Table 4: Gene-level summary
cat("\n  Creating gene-level summary...\n")

# Get unique genes from all sources
all_genes <- unique(c(
    na.omit(tes_anno_df$SYMBOL),
    na.omit(tead1_anno_df$SYMBOL),
    na.omit(dmr_anno_df$SYMBOL)
))

gene_summary <- data.frame(
    gene = all_genes,
    has_TES_peak = all_genes %in% tes_anno_df$SYMBOL,
    has_TEAD1_peak = all_genes %in% tead1_anno_df$SYMBOL,
    has_DMR = all_genes %in% dmr_anno_df$SYMBOL,
    stringsAsFactors = FALSE
)

# Count features per gene
gene_summary$TES_peak_count <- sapply(gene_summary$gene, function(g) sum(tes_anno_df$SYMBOL == g, na.rm = TRUE))
gene_summary$TEAD1_peak_count <- sapply(gene_summary$gene, function(g) sum(tead1_anno_df$SYMBOL == g, na.rm = TRUE))
gene_summary$DMR_count <- sapply(gene_summary$gene, function(g) sum(dmr_anno_df$SYMBOL == g, na.rm = TRUE))

# Classify genes
gene_summary$category <- "Other"
gene_summary$category[gene_summary$has_TES_peak & gene_summary$has_DMR] <- "TES + DMR"
gene_summary$category[gene_summary$has_TEAD1_peak & gene_summary$has_DMR] <- "TEAD1 + DMR"
gene_summary$category[gene_summary$has_TES_peak & gene_summary$has_TEAD1_peak & gene_summary$has_DMR] <- "TES + TEAD1 + DMR"
gene_summary$category[gene_summary$has_TES_peak & !gene_summary$has_DMR] <- "TES only"
gene_summary$category[gene_summary$has_TEAD1_peak & !gene_summary$has_DMR] <- "TEAD1 only"
gene_summary$category[!gene_summary$has_TES_peak & !gene_summary$has_TEAD1_peak & gene_summary$has_DMR] <- "DMR only"

write.csv(gene_summary, file.path(OUTPUT_BASE, "tables", "gene_level_summary.csv"), row.names = FALSE)
cat("  Created: gene_level_summary.csv\n")

# Print category summary
cat("\n  Gene category distribution:\n")
print(table(gene_summary$category))

cat("\n")

# =============================================================================
# PHASE 6: CREATE VISUALIZATIONS
# =============================================================================

cat("=== PHASE 6: Creating Visualizations ===\n")

# Plot 1: Venn Diagram - Peaks and DMRs at gene level
cat("  Creating Venn diagram...\n")

tes_genes <- unique(na.omit(tes_anno_df$SYMBOL))
tead1_genes <- unique(na.omit(tead1_anno_df$SYMBOL))
dmr_genes <- unique(na.omit(dmr_anno_df$SYMBOL))

# Three-way Venn
pdf(file.path(OUTPUT_BASE, "plots", "venn_genes_peaks_DMRs.pdf"), width = 10, height = 8)

venn_data <- list(
    TES = tes_genes,
    TEAD1 = tead1_genes,
    DMR = dmr_genes
)

venn.plot <- venn.diagram(
    x = venn_data,
    category.names = c(paste0("TES peaks\n(", length(tes_genes), " genes)"),
                       paste0("TEAD1 peaks\n(", length(tead1_genes), " genes)"),
                       paste0("DMRs\n(", length(dmr_genes), " genes)")),
    filename = NULL,
    output = TRUE,
    main = "Genes with TES/TEAD1 Peaks and DMRs",
    main.cex = 1.5,
    fill = c("#E41A1C", "#377EB8", "#4DAF4A"),
    alpha = 0.5,
    cat.cex = 1.2,
    cat.fontface = "bold",
    margin = 0.1
)
grid.draw(venn.plot)
dev.off()

cat("  Created: venn_genes_peaks_DMRs.pdf\n")

# Plot 2: Distance distribution
cat("  Creating distance distribution plots...\n")

pdf(file.path(OUTPUT_BASE, "plots", "distance_distributions.pdf"), width = 12, height = 8)

par(mfrow = c(2, 2))

# TES to nearest DMR
hist(log10(tes_peaks$nearest_dmr_distance + 1), breaks = 50,
     main = "Distance: TES Peak to Nearest DMR",
     xlab = "log10(Distance + 1) [bp]", col = "#E41A1C", border = "white")
abline(v = log10(5001), col = "blue", lty = 2, lwd = 2)
legend("topright", "5kb threshold", col = "blue", lty = 2, lwd = 2)

# TEAD1 to nearest DMR
hist(log10(tead1_peaks$nearest_dmr_distance + 1), breaks = 50,
     main = "Distance: TEAD1 Peak to Nearest DMR",
     xlab = "log10(Distance + 1) [bp]", col = "#377EB8", border = "white")
abline(v = log10(5001), col = "blue", lty = 2, lwd = 2)

# DMR to nearest TES peak
hist(log10(dmr_table$nearest_TES_peak_distance + 1), breaks = 50,
     main = "Distance: DMR to Nearest TES Peak",
     xlab = "log10(Distance + 1) [bp]", col = "#4DAF4A", border = "white")
abline(v = log10(5001), col = "blue", lty = 2, lwd = 2)

# DMR to nearest TEAD1 peak
hist(log10(dmr_table$nearest_TEAD1_peak_distance + 1), breaks = 50,
     main = "Distance: DMR to Nearest TEAD1 Peak",
     xlab = "log10(Distance + 1) [bp]", col = "#984EA3", border = "white")
abline(v = log10(5001), col = "blue", lty = 2, lwd = 2)

dev.off()

cat("  Created: distance_distributions.pdf\n")

# Plot 3: Bar chart of overlaps
cat("  Creating overlap summary bar chart...\n")

overlap_summary <- data.frame(
    Category = c("TES peaks with DMR", "TES peaks without DMR",
                 "TEAD1 peaks with DMR", "TEAD1 peaks without DMR",
                 "DMRs with TES peak", "DMRs without TES peak",
                 "DMRs with TEAD1 peak", "DMRs without TEAD1 peak"),
    Count = c(sum(tes_table$has_overlapping_DMR), sum(!tes_table$has_overlapping_DMR),
              sum(tead1_table$has_overlapping_DMR), sum(!tead1_table$has_overlapping_DMR),
              sum(dmr_table$has_TES_peak), sum(!dmr_table$has_TES_peak),
              sum(dmr_table$has_TEAD1_peak), sum(!dmr_table$has_TEAD1_peak)),
    Type = rep(c("With", "Without"), 4),
    Feature = rep(c("TES peaks", "TES peaks", "TEAD1 peaks", "TEAD1 peaks",
                    "DMRs (vs TES)", "DMRs (vs TES)", "DMRs (vs TEAD1)", "DMRs (vs TEAD1)"), each = 1)
)

p_overlap <- ggplot(overlap_summary, aes(x = Feature, y = Count, fill = Type)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("With" = "#4DAF4A", "Without" = "#E41A1C")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Peak-DMR Overlap Summary",
         x = "", y = "Count", fill = "Overlap Status")

ggsave(file.path(OUTPUT_BASE, "plots", "overlap_summary_barplot.pdf"), p_overlap, width = 10, height = 6)
cat("  Created: overlap_summary_barplot.pdf\n")

# Plot 4: DMR direction at peaks
cat("  Creating DMR direction analysis...\n")

dmr_direction_at_peaks <- data.frame(
    Peak_Type = c(rep("TES", 2), rep("TEAD1", 2)),
    DMR_Direction = rep(c("Hyper", "Hypo"), 2),
    Count = c(
        sum(tes_table$overlapping_DMR_direction == "Hyper", na.rm = TRUE),
        sum(tes_table$overlapping_DMR_direction == "Hypo", na.rm = TRUE),
        sum(tead1_table$overlapping_DMR_direction == "Hyper", na.rm = TRUE),
        sum(tead1_table$overlapping_DMR_direction == "Hypo", na.rm = TRUE)
    )
)

p_direction <- ggplot(dmr_direction_at_peaks, aes(x = Peak_Type, y = Count, fill = DMR_Direction)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("Hyper" = "#E41A1C", "Hypo" = "#377EB8")) +
    theme_minimal() +
    labs(title = "DMR Direction at Peak Sites",
         subtitle = "Hyper = increased methylation in TES vs GFP",
         x = "Peak Type", y = "Number of overlapping DMRs", fill = "DMR Direction")

ggsave(file.path(OUTPUT_BASE, "plots", "DMR_direction_at_peaks.pdf"), p_direction, width = 8, height = 6)
cat("  Created: DMR_direction_at_peaks.pdf\n")

# Plot 5: Genomic annotation distribution
cat("  Creating genomic annotation plots...\n")

pdf(file.path(OUTPUT_BASE, "plots", "genomic_annotations.pdf"), width = 14, height = 10)

par(mfrow = c(2, 2))

# Simplify annotations
simplify_anno <- function(anno) {
    anno <- gsub(" \\(.*", "", anno)
    anno <- gsub("Distal Intergenic", "Intergenic", anno)
    return(anno)
}

tes_anno_simple <- simplify_anno(tes_anno_df$annotation)
tead1_anno_simple <- simplify_anno(tead1_anno_df$annotation)
dmr_anno_simple <- simplify_anno(dmr_anno_df$annotation)

barplot(sort(table(tes_anno_simple), decreasing = TRUE),
        las = 2, col = "#E41A1C", main = "TES Peak Annotations", cex.names = 0.8)

barplot(sort(table(tead1_anno_simple), decreasing = TRUE),
        las = 2, col = "#377EB8", main = "TEAD1 Peak Annotations", cex.names = 0.8)

barplot(sort(table(dmr_anno_simple), decreasing = TRUE),
        las = 2, col = "#4DAF4A", main = "DMR Annotations", cex.names = 0.8)

# Combined comparison
anno_comparison <- data.frame(
    Annotation = rep(c("Promoter", "Intron", "Exon", "Intergenic", "3' UTR", "5' UTR"), 3),
    Feature = rep(c("TES", "TEAD1", "DMR"), each = 6),
    Percentage = c(
        100 * sum(grepl("Promoter", tes_anno_simple)) / length(tes_anno_simple),
        100 * sum(grepl("Intron", tes_anno_simple)) / length(tes_anno_simple),
        100 * sum(grepl("Exon", tes_anno_simple)) / length(tes_anno_simple),
        100 * sum(grepl("Intergenic", tes_anno_simple)) / length(tes_anno_simple),
        100 * sum(grepl("3' UTR", tes_anno_simple)) / length(tes_anno_simple),
        100 * sum(grepl("5' UTR", tes_anno_simple)) / length(tes_anno_simple),
        100 * sum(grepl("Promoter", tead1_anno_simple)) / length(tead1_anno_simple),
        100 * sum(grepl("Intron", tead1_anno_simple)) / length(tead1_anno_simple),
        100 * sum(grepl("Exon", tead1_anno_simple)) / length(tead1_anno_simple),
        100 * sum(grepl("Intergenic", tead1_anno_simple)) / length(tead1_anno_simple),
        100 * sum(grepl("3' UTR", tead1_anno_simple)) / length(tead1_anno_simple),
        100 * sum(grepl("5' UTR", tead1_anno_simple)) / length(tead1_anno_simple),
        100 * sum(grepl("Promoter", dmr_anno_simple)) / length(dmr_anno_simple),
        100 * sum(grepl("Intron", dmr_anno_simple)) / length(dmr_anno_simple),
        100 * sum(grepl("Exon", dmr_anno_simple)) / length(dmr_anno_simple),
        100 * sum(grepl("Intergenic", dmr_anno_simple)) / length(dmr_anno_simple),
        100 * sum(grepl("3' UTR", dmr_anno_simple)) / length(dmr_anno_simple),
        100 * sum(grepl("5' UTR", dmr_anno_simple)) / length(dmr_anno_simple)
    )
)

# Create matrix: rows = features (TES, TEAD1, DMR), cols = annotations
anno_matrix <- matrix(anno_comparison$Percentage, nrow = 3, byrow = TRUE)
colnames(anno_matrix) <- c("Promoter", "Intron", "Exon", "Intergenic", "3' UTR", "5' UTR")
rownames(anno_matrix) <- c("TES", "TEAD1", "DMR")

barplot(anno_matrix, beside = TRUE,
        col = c("#E41A1C", "#377EB8", "#4DAF4A"),
        legend.text = rownames(anno_matrix),
        args.legend = list(x = "topright"),
        main = "Genomic Distribution Comparison",
        ylab = "Percentage", las = 2)

dev.off()

cat("  Created: genomic_annotations.pdf\n")

# Plot 6: Scatter plot - Peak position vs nearest DMR distance
cat("  Creating peak-DMR relationship scatter...\n")

pdf(file.path(OUTPUT_BASE, "plots", "peak_DMR_scatter.pdf"), width = 12, height = 6)

par(mfrow = c(1, 2))

# TES peaks colored by whether they overlap a DMR
plot(1:nrow(tes_table), log10(tes_table$nearest_DMR_distance + 1),
     col = ifelse(tes_table$has_overlapping_DMR, "#E41A1C", "#CCCCCC"),
     pch = 16, cex = 0.5,
     main = "TES Peaks: Distance to Nearest DMR",
     xlab = "Peak index", ylab = "log10(Distance + 1) [bp]")
abline(h = log10(5001), col = "blue", lty = 2)
legend("topright", c("Overlapping DMR", "No overlap", "5kb threshold"),
       col = c("#E41A1C", "#CCCCCC", "blue"), pch = c(16, 16, NA), lty = c(NA, NA, 2))

# TEAD1 peaks
plot(1:nrow(tead1_table), log10(tead1_table$nearest_DMR_distance + 1),
     col = ifelse(tead1_table$has_overlapping_DMR, "#377EB8", "#CCCCCC"),
     pch = 16, cex = 0.5,
     main = "TEAD1 Peaks: Distance to Nearest DMR",
     xlab = "Peak index", ylab = "log10(Distance + 1) [bp]")
abline(h = log10(5001), col = "blue", lty = 2)
legend("topright", c("Overlapping DMR", "No overlap", "5kb threshold"),
       col = c("#377EB8", "#CCCCCC", "blue"), pch = c(16, 16, NA), lty = c(NA, NA, 2))

dev.off()

cat("  Created: peak_DMR_scatter.pdf\n")

# Plot 7: Cumulative distribution of distances
cat("  Creating cumulative distance plot...\n")

pdf(file.path(OUTPUT_BASE, "plots", "cumulative_distance.pdf"), width = 8, height = 6)

tes_dist <- sort(tes_peaks$nearest_dmr_distance[!is.na(tes_peaks$nearest_dmr_distance)])
tead1_dist <- sort(tead1_peaks$nearest_dmr_distance[!is.na(tead1_peaks$nearest_dmr_distance)])

plot(tes_dist, (1:length(tes_dist))/length(tes_dist), type = "l", col = "#E41A1C", lwd = 2,
     xlim = c(0, 100000), ylim = c(0, 1),
     main = "Cumulative Distribution: Peak to Nearest DMR Distance",
     xlab = "Distance to nearest DMR (bp)", ylab = "Cumulative fraction of peaks")
lines(tead1_dist, (1:length(tead1_dist))/length(tead1_dist), col = "#377EB8", lwd = 2)
abline(v = 5000, col = "gray", lty = 2)
abline(v = 10000, col = "gray", lty = 2)
legend("bottomright", c("TES peaks", "TEAD1 peaks", "5kb", "10kb"),
       col = c("#E41A1C", "#377EB8", "gray", "gray"), lty = c(1, 1, 2, 2), lwd = 2)

dev.off()

cat("  Created: cumulative_distance.pdf\n")

cat("\n")

# =============================================================================
# PHASE 7: SUMMARY STATISTICS
# =============================================================================

cat("=== PHASE 7: Summary Statistics ===\n\n")

summary_stats <- list(
    "TES_peaks_total" = length(tes_peaks),
    "TEAD1_peaks_total" = length(tead1_peaks),
    "DMRs_total" = length(dmrs),
    "DMRs_hyper" = sum(dmr_df$direction == "Hyper"),
    "DMRs_hypo" = sum(dmr_df$direction == "Hypo"),
    "TES_peaks_overlapping_DMR" = sum(tes_table$has_overlapping_DMR),
    "TES_peaks_overlapping_DMR_pct" = round(100 * mean(tes_table$has_overlapping_DMR), 2),
    "TEAD1_peaks_overlapping_DMR" = sum(tead1_table$has_overlapping_DMR),
    "TEAD1_peaks_overlapping_DMR_pct" = round(100 * mean(tead1_table$has_overlapping_DMR), 2),
    "TES_peaks_nearby_DMR_5kb" = sum(tes_table$has_nearby_DMR_5kb),
    "TES_peaks_nearby_DMR_5kb_pct" = round(100 * mean(tes_table$has_nearby_DMR_5kb), 2),
    "TEAD1_peaks_nearby_DMR_5kb" = sum(tead1_table$has_nearby_DMR_5kb),
    "TEAD1_peaks_nearby_DMR_5kb_pct" = round(100 * mean(tead1_table$has_nearby_DMR_5kb), 2),
    "DMRs_at_TES_peaks" = sum(dmr_table$has_TES_peak),
    "DMRs_at_TEAD1_peaks" = sum(dmr_table$has_TEAD1_peak),
    "median_TES_to_DMR_distance" = median(tes_peaks$nearest_dmr_distance, na.rm = TRUE),
    "median_TEAD1_to_DMR_distance" = median(tead1_peaks$nearest_dmr_distance, na.rm = TRUE),
    "genes_with_TES_peak" = length(tes_genes),
    "genes_with_TEAD1_peak" = length(tead1_genes),
    "genes_with_DMR" = length(dmr_genes),
    "genes_with_TES_and_DMR" = length(intersect(tes_genes, dmr_genes)),
    "genes_with_TEAD1_and_DMR" = length(intersect(tead1_genes, dmr_genes)),
    "genes_with_all_three" = length(Reduce(intersect, list(tes_genes, tead1_genes, dmr_genes)))
)

# Print summary
cat("FEATURE COUNTS:\n")
cat(sprintf("  TES peaks: %d\n", summary_stats$TES_peaks_total))
cat(sprintf("  TEAD1 peaks: %d\n", summary_stats$TEAD1_peaks_total))
cat(sprintf("  DMRs: %d (Hyper: %d, Hypo: %d)\n",
            summary_stats$DMRs_total, summary_stats$DMRs_hyper, summary_stats$DMRs_hypo))

cat("\nOVERLAP STATISTICS:\n")
cat(sprintf("  TES peaks overlapping DMR: %d (%.1f%%)\n",
            summary_stats$TES_peaks_overlapping_DMR, summary_stats$TES_peaks_overlapping_DMR_pct))
cat(sprintf("  TEAD1 peaks overlapping DMR: %d (%.1f%%)\n",
            summary_stats$TEAD1_peaks_overlapping_DMR, summary_stats$TEAD1_peaks_overlapping_DMR_pct))
cat(sprintf("  TES peaks with DMR within 5kb: %d (%.1f%%)\n",
            summary_stats$TES_peaks_nearby_DMR_5kb, summary_stats$TES_peaks_nearby_DMR_5kb_pct))
cat(sprintf("  TEAD1 peaks with DMR within 5kb: %d (%.1f%%)\n",
            summary_stats$TEAD1_peaks_nearby_DMR_5kb, summary_stats$TEAD1_peaks_nearby_DMR_5kb_pct))

cat("\nDISTANCE STATISTICS:\n")
cat(sprintf("  Median distance TES peak to nearest DMR: %.0f bp\n", summary_stats$median_TES_to_DMR_distance))
cat(sprintf("  Median distance TEAD1 peak to nearest DMR: %.0f bp\n", summary_stats$median_TEAD1_to_DMR_distance))

cat("\nGENE-LEVEL OVERLAP:\n")
cat(sprintf("  Genes with TES peak: %d\n", summary_stats$genes_with_TES_peak))
cat(sprintf("  Genes with TEAD1 peak: %d\n", summary_stats$genes_with_TEAD1_peak))
cat(sprintf("  Genes with DMR: %d\n", summary_stats$genes_with_DMR))
cat(sprintf("  Genes with TES peak + DMR: %d\n", summary_stats$genes_with_TES_and_DMR))
cat(sprintf("  Genes with TEAD1 peak + DMR: %d\n", summary_stats$genes_with_TEAD1_and_DMR))
cat(sprintf("  Genes with TES + TEAD1 + DMR: %d\n", summary_stats$genes_with_all_three))

# Save summary
write.csv(as.data.frame(summary_stats),
          file.path(OUTPUT_BASE, "tables", "summary_statistics.csv"), row.names = FALSE)

cat("\n")
cat("==========================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("==========================================================\n")
cat("Completed:", as.character(Sys.time()), "\n")
cat(sprintf("Output directory: %s\n", OUTPUT_BASE))
cat("\nGenerated files:\n")
cat("  Tables:\n")
cat("    - TES_peaks_with_DMR_info.csv\n")
cat("    - TEAD1_peaks_with_DMR_info.csv\n")
cat("    - DMRs_with_peak_info.csv\n")
cat("    - gene_level_summary.csv\n")
cat("    - summary_statistics.csv\n")
cat("  Plots:\n")
cat("    - venn_genes_peaks_DMRs.pdf\n")
cat("    - distance_distributions.pdf\n")
cat("    - overlap_summary_barplot.pdf\n")
cat("    - DMR_direction_at_peaks.pdf\n")
cat("    - genomic_annotations.pdf\n")
cat("    - peak_DMR_scatter.pdf\n")
cat("    - cumulative_distance.pdf\n")
