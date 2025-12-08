#!/usr/bin/env Rscript

library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(ggplot2)
library(ggpubr)

setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top")

message("Loading data...")

# Load binding classification
load("SRF_Eva_integrated_analysis/scripts/analysis_3/results/01_binding_classification/binding_classification_data.RData")

# Load meDIP DMRs
medip_dmrs_file <- "meDIP/results/07_differential_MEDIPS/TES_vs_GFP_DMRs_FDR05_FC2.csv"

if (!file.exists(medip_dmrs_file)) {
  message("WARNING: meDIP DMRs not found, skipping DMR overlap analysis")
  dmrs_available <- FALSE
} else {
  dmrs <- read.csv(medip_dmrs_file)
  message("  Loaded ", nrow(dmrs), " DMRs")

  # Convert to GRanges (add "chr" prefix to match binding sites)
  dmrs_gr <- GRanges(
    seqnames = paste0("chr", dmrs$chr),
    ranges = IRanges(start = dmrs$start, end = dmrs$stop),
    logFC = dmrs$logFC,
    FDR = dmrs$FDR
  )
  dmrs_available <- TRUE
}

# Separate hypermethylated and hypomethylated DMRs
if (dmrs_available) {
  hyper_dmrs <- dmrs_gr[dmrs_gr$logFC > 0]
  hypo_dmrs <- dmrs_gr[dmrs_gr$logFC < 0]

  message("  Hypermethylated DMRs: ", length(hyper_dmrs))
  message("  Hypomethylated DMRs: ", length(hypo_dmrs))
}

# Analyze overlap with binding sites
message("\nAnalyzing DMR-peak overlaps...")

OUTPUT_DIR <- "SRF_Eva_integrated_analysis/scripts/analysis_3/results/07_methylation_binding"

if (dmrs_available) {
  # Find overlaps
  hyper_overlaps <- findOverlaps(all_peaks_unique, hyper_dmrs)
  hypo_overlaps <- findOverlaps(all_peaks_unique, hypo_dmrs)

  # Add methylation info to peaks
  all_peaks_unique$has_hyper_dmr <- FALSE
  all_peaks_unique$has_hypo_dmr <- FALSE
  all_peaks_unique$has_hyper_dmr[queryHits(hyper_overlaps)] <- TRUE
  all_peaks_unique$has_hypo_dmr[queryHits(hypo_overlaps)] <- TRUE

  # Count DMRs per category
  dmr_counts <- as.data.frame(all_peaks_unique) %>%
    group_by(category) %>%
    summarise(
      n_peaks = n(),
      n_with_hyper = sum(has_hyper_dmr),
      pct_hyper = round(sum(has_hyper_dmr) / n() * 100, 2),
      n_with_hypo = sum(has_hypo_dmr),
      pct_hypo = round(sum(has_hypo_dmr) / n() * 100, 2),
      .groups = "drop"
    )

  write.csv(dmr_counts,
            file.path(OUTPUT_DIR, "dmr_overlap_by_category.csv"),
            row.names = FALSE)

  # Statistical test: Are TES sites enriched for hypermethylation?
  tes_related <- all_peaks_unique$category %in% c("TES_unique", "Shared_TES_dominant")
  tead1_related <- all_peaks_unique$category %in% c("TEAD1_unique", "Shared_TEAD1_dominant")

  # Fisher's exact test (with error handling for insufficient data)
  tes_hyper_test <- tryCatch({
    tbl <- table(tes_related, as.data.frame(all_peaks_unique)$has_hyper_dmr)
    if (nrow(tbl) >= 2 && ncol(tbl) >= 2) {
      fisher.test(tbl)
    } else {
      list(p.value = NA, estimate = NA, conf.int = c(NA, NA),
           method = "Fisher's Exact Test (insufficient data)")
    }
  }, error = function(e) {
    list(p.value = NA, estimate = NA, conf.int = c(NA, NA),
         method = paste("Fisher's Exact Test failed:", e$message))
  })

  tead1_hyper_test <- tryCatch({
    tbl <- table(tead1_related, as.data.frame(all_peaks_unique)$has_hyper_dmr)
    if (nrow(tbl) >= 2 && ncol(tbl) >= 2) {
      fisher.test(tbl)
    } else {
      list(p.value = NA, estimate = NA, conf.int = c(NA, NA),
           method = "Fisher's Exact Test (insufficient data)")
    }
  }, error = function(e) {
    list(p.value = NA, estimate = NA, conf.int = c(NA, NA),
         method = paste("Fisher's Exact Test failed:", e$message))
  })

  # Visualizations
  message("Creating visualizations...")

  # DMR overlap barplot
  pdf(file.path(OUTPUT_DIR, "dmr_overlap_by_category.pdf"), width = 10, height = 6)
  dmr_counts_long <- dmr_counts %>%
    tidyr::pivot_longer(cols = c(pct_hyper, pct_hypo),
                        names_to = "dmr_type",
                        values_to = "percentage") %>%
    mutate(dmr_type = ifelse(dmr_type == "pct_hyper", "Hypermethylated", "Hypomethylated"))

  p <- ggplot(dmr_counts_long, aes(x = category, y = percentage, fill = dmr_type)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_text(aes(label = paste0(round(percentage, 1), "%")),
              position = position_dodge(width = 0.9), vjust = -0.5, size = 3) +
    scale_fill_manual(values = c("Hypermethylated" = "red", "Hypomethylated" = "blue")) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = "DMR Overlap with Binding Sites",
      subtitle = "Percentage of peaks overlapping differentially methylated regions",
      x = "Binding Category",
      y = "% of Peaks with DMR",
      fill = "DMR Type"
    )
  print(p)
  dev.off()

  # Summary report
  sink(file.path(OUTPUT_DIR, "PHASE3_1_SUMMARY.txt"))
  cat("=== Phase 3.1: Methylation at Binding Sites ===\n")
  cat("Date:", as.character(Sys.time()), "\n\n")

  cat("DMR Overlap Counts:\n")
  cat("===================\n")
  print(dmr_counts)

  cat("\n\nStatistical Tests:\n")
  cat("==================\n")
  cat("TES-related peaks vs hypermethylation:\n")
  print(tes_hyper_test)
  cat("\n")
  cat("TEAD1-related peaks vs hypermethylation:\n")
  print(tead1_hyper_test)

  sink()
}

message("\n=== Phase 3.1 Complete ===")
