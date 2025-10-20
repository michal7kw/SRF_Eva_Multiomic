#!/usr/bin/env Rscript

cat("=================================================================\n")
cat("APPROACH 5: Differential Binding + Differential Expression\n")
cat("=================================================================\n")

# Load workspace and helper functions
source("SRF_Eva_integrated_analysis/focused_enrichment_analysis/scripts/00_helper_functions.R")
load_libraries()
setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top")
load("SRF_Eva_integrated_analysis/focused_enrichment_analysis/data/loaded_data.RData")

# Output directory
base_dir <- "SRF_Eva_integrated_analysis/focused_enrichment_analysis"
approach5_dir <- file.path(base_dir, "approach5_diffbind")

# Create directories
dir.create(file.path(approach5_dir, "results"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(approach5_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(approach5_dir, "gene_lists"), recursive = TRUE, showWarnings = FALSE)

################################################################################
# ANALYSIS
################################################################################

# Filter for TES-specific differential peaks (TES > TESmut)
if ("FDR" %in% colnames(diffbind_tes_vs_tesmut)) {
  tes_specific_diffbind <- diffbind_tes_vs_tesmut %>%
    filter(FDR < 0.05, Fold > 1.5)

  cat("Total differential sites:", nrow(diffbind_tes_vs_tesmut), "\n")
  cat("TES-specific sites (vs TESmut):", nrow(tes_specific_diffbind), "\n")

  # DiffBind results don't have gene annotations
  # Overlap with annotated peaks to get genes
  cat("Overlapping differential peaks with annotated peaks...\n")

  # Create GRanges for differential peaks
  # Add "chr" prefix if missing for chromosome name compatibility
  diffbind_seqnames <- tes_specific_diffbind$seqnames
  if (!grepl("^chr", diffbind_seqnames[1])) {
    diffbind_seqnames <- paste0("chr", diffbind_seqnames)
  }

  diffbind_gr <- GRanges(
    seqnames = diffbind_seqnames,
    ranges = IRanges(start = tes_specific_diffbind$start,
                     end = tes_specific_diffbind$end)
  )

  # Create GRanges for annotated peaks
  annotated_gr <- GRanges(
    seqnames = tes_peaks_annotated$seqnames,
    ranges = IRanges(start = tes_peaks_annotated$start,
                     end = tes_peaks_annotated$end)
  )

  # Find overlaps
  overlaps <- findOverlaps(diffbind_gr, annotated_gr)
  overlapping_annotated_idx <- subjectHits(overlaps)

  cat("  DiffBind peaks:", length(diffbind_gr), "\n")
  cat("  Annotated peaks:", length(annotated_gr), "\n")
  cat("  Overlaps found:", length(overlaps), "\n")

  # Get genes from overlapping annotated peaks
  overlapping_peaks <- tes_peaks_annotated[overlapping_annotated_idx, ]
  tes_diffbind_genes <- extract_genes_from_peaks(overlapping_peaks)

  if (length(tes_diffbind_genes) > 0) {
    cat("Unique genes at differential sites:", length(tes_diffbind_genes), "\n")

    # Intersect with downregulated DEGs (since TES is repressor)
    # Load tes_direct_down_genes from approach2
    tes_direct_down <- tes_direct %>%
      filter(log2FoldChange < 0, padj < 0.05)
    tes_direct_down_genes <- tes_direct_down$gene_symbol

    tes_diffbind_down <- intersect(tes_diffbind_genes, tes_direct_down_genes)
    cat("Differential binding + downregulated:", length(tes_diffbind_down), "\n")

    # Save gene lists
    if (length(tes_diffbind_genes) > 0) {
      write.table(tes_diffbind_genes,
                  file.path(approach5_dir, "gene_lists", "TES_diffbind_genes.txt"),
                  row.names = FALSE, col.names = FALSE, quote = FALSE)
    }

    if (length(tes_diffbind_down) > 0) {
      write.table(tes_diffbind_down,
                  file.path(approach5_dir, "gene_lists", "TES_diffbind_down.txt"),
                  row.names = FALSE, col.names = FALSE, quote = FALSE)
    }

    # GO enrichment - try with all diffbind genes if down is too few
    genes_for_enrichment <- if (length(tes_diffbind_down) >= 10) {
      cat("Using differential binding + downregulated genes for enrichment\n")
      tes_diffbind_down
    } else if (length(tes_diffbind_genes) >= 10) {
      cat("Too few downregulated genes, using all differential binding genes instead\n")
      tes_diffbind_genes
    } else {
      cat("Too few genes for enrichment analysis\n")
      NULL
    }

    if (!is.null(genes_for_enrichment) && length(genes_for_enrichment) >= 10) {
      cat("Running GO enrichment on", length(genes_for_enrichment), "genes...\n")
      ego_approach5 <- perform_GO_enrichment(
        gene_list = genes_for_enrichment,
        background_genes = all_expressed_genes,
        ont = "BP"
      )

      results_approach5 <- save_enrichment_results(ego_approach5, "approach5_diffbind", approach5_dir)

      # Extract migration terms
      if (!is.null(results_approach5)) {
        migration_approach5 <- extract_migration_terms(results_approach5)
        if (!is.null(migration_approach5) && nrow(migration_approach5) > 0) {
          write_csv(migration_approach5,
                    file.path(approach5_dir, "results", "approach5_migration_terms.csv"))
          cat("  Found", nrow(migration_approach5), "migration-related terms\n")
          cat("  Top migration term rank:", which(results_approach5$ID == migration_approach5$ID[1]), "\n")
        }
      }
    } else {
      cat("Skipping enrichment due to insufficient genes\n")
    }
  }
} else {
  cat("DiffBind results not in expected format, skipping Approach 5\n")
}

cat("\n=================================================================\n")
cat("APPROACH 5 COMPLETE\n")
cat("=================================================================\n")
