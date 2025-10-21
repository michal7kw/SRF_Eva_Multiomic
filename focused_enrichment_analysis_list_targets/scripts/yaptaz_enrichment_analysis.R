#!/usr/bin/env Rscript

cat("=================================================================\n")
cat("YAP/TAZ Target Enrichment Analysis\n")
cat("Using known YAP/TAZ targets instead of RNA-seq DEGs\n")
cat("=================================================================\n")

# Load workspace and helper functions
source("SRF_Eva_integrated_analysis/focused_enrichment_analysis/scripts/00_helper_functions.R")
load_libraries()
setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top")
load("SRF_Eva_integrated_analysis/focused_enrichment_analysis/data/loaded_data.RData")

################################################################################
# CUSTOM PLOTTING FUNCTION FOR THIS ANALYSIS
################################################################################

save_enrichment_results_custom <- function(enrich_obj, prefix, output_dir) {

  if (is.null(enrich_obj) || nrow(enrich_obj) == 0) {
    cat("  No significant enrichment found for", prefix, "\n")
    return(NULL)
  }

  cat("  Processing", nrow(enrich_obj), "enriched terms for", prefix, "\n")

  # Save full results
  results_df <- as.data.frame(enrich_obj)
  write_csv(results_df, file.path(output_dir, "results", paste0(prefix, "_GO_enrichment.csv")))
  cat("  Saved CSV results\n")

  # Save top 20 results
  if (nrow(results_df) > 0) {
    top20 <- head(results_df, 20)
    write_csv(top20, file.path(output_dir, "results", paste0(prefix, "_top20_terms.csv")))
  }

  # Generate plots
  if (nrow(enrich_obj) >= 10) {
    cat("  Generating plots for", prefix, "...\n")

    # Determine number of categories to show
    n_categories <- min(20, nrow(enrich_obj))

    # Calculate dynamic height based on number of terms
    # Base height + 0.6 inches per term (ensures adequate spacing)
    plot_height <- max(16, 6 + (n_categories * 0.6))

    cat("  Plot dimensions: width=16, height=", plot_height, "\n")

    # Dot plot
    p1 <- dotplot(enrich_obj, showCategory = n_categories) +
      ggtitle(paste(prefix, "- Top", n_categories, "GO Terms")) +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, lineheight = 1.3, margin = margin(r = 8)),
        axis.text.x = element_text(size = 11),
        axis.title = element_text(size = 12, face = "bold"),
        plot.margin = margin(10, 15, 10, 10)
      ) +
      scale_y_discrete(labels = function(x) {
        # Wrap long labels at 50 characters
        sapply(x, function(label) {
          if (nchar(label) > 50) {
            paste(strwrap(label, width = 50), collapse = "\n")
          } else {
            label
          }
        })
      })

    # Save dotplot as PDF
    pdf_file <- file.path(output_dir, "plots", paste0(prefix, "_dotplot.pdf"))
    cat("  Saving dotplot PDF:", basename(pdf_file), "\n")
    ggsave(pdf_file, p1, width = 16, height = plot_height, device = "pdf")

    # Save dotplot as PNG
    png_file <- file.path(output_dir, "plots", paste0(prefix, "_dotplot.png"))
    cat("  Saving dotplot PNG:", basename(png_file), "\n")
    ggsave(png_file, p1, width = 16, height = plot_height, device = "png", dpi = 300)

    # Bar plot
    p2 <- barplot(enrich_obj, showCategory = n_categories) +
      ggtitle(paste(prefix, "- Top", n_categories, "GO Terms")) +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, lineheight = 1.3, margin = margin(r = 8)),
        axis.text.x = element_text(size = 11),
        axis.title = element_text(size = 12, face = "bold"),
        plot.margin = margin(10, 15, 10, 10)
      ) +
      scale_y_discrete(labels = function(x) {
        # Wrap long labels at 50 characters
        sapply(x, function(label) {
          if (nchar(label) > 50) {
            paste(strwrap(label, width = 50), collapse = "\n")
          } else {
            label
          }
        })
      })

    # Save barplot as PDF
    pdf_file <- file.path(output_dir, "plots", paste0(prefix, "_barplot.pdf"))
    cat("  Saving barplot PDF:", basename(pdf_file), "\n")
    ggsave(pdf_file, p2, width = 16, height = plot_height, device = "pdf")

    # Save barplot as PNG
    png_file <- file.path(output_dir, "plots", paste0(prefix, "_barplot.png"))
    cat("  Saving barplot PNG:", basename(png_file), "\n")
    ggsave(png_file, p2, width = 16, height = plot_height, device = "png", dpi = 300)

    cat("  All plots saved successfully!\n")
  } else {
    cat("  Too few enriched terms for plotting (", nrow(enrich_obj), ")\n")
  }

  return(results_df)
}

# Output directory
base_dir <- "SRF_Eva_integrated_analysis/focused_enrichment_analysis_list_targets"

# Create directories (should already exist but ensure)
dir.create(file.path(base_dir, "results"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(base_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(base_dir, "gene_lists"), recursive = TRUE, showWarnings = FALSE)

################################################################################
# LOAD YAP/TAZ TARGET GENES
################################################################################

cat("\n=== Loading YAP/TAZ target genes ===\n")

# Load YAP/TAZ target gene list from TES_degs.txt
yaptaz_targets_raw <- read.table(
  "SRF_Eva_integrated_analysis/data/TES_degs.txt",
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE
)

# Extract gene symbols (remove row numbers and arrow)
yaptaz_targets <- gsub("^.*→", "", yaptaz_targets_raw$V1)
yaptaz_targets <- unique(yaptaz_targets[yaptaz_targets != ""])

cat("Total YAP/TAZ target genes loaded:", length(yaptaz_targets), "\n")

# Save the cleaned gene list
write.table(yaptaz_targets,
            file.path(base_dir, "gene_lists", "yaptaz_targets_all.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)

################################################################################
# INTERSECT YAP/TAZ TARGETS WITH TES/TEAD1 BINDING
################################################################################

cat("\n=== Intersection with TES/TEAD1 binding data ===\n")

# Get genes bound by TES from Cut&Tag data
tes_bound_genes <- extract_genes_from_peaks(tes_peaks_annotated)
cat("Genes with TES binding peaks:", length(tes_bound_genes), "\n")

# Get genes bound by TEAD1
tead1_bound_genes <- extract_genes_from_peaks(tead1_peaks_annotated)
cat("Genes with TEAD1 binding peaks:", length(tead1_bound_genes), "\n")

# Intersections
yaptaz_tes_bound <- intersect(yaptaz_targets, tes_bound_genes)
yaptaz_tead1_bound <- intersect(yaptaz_targets, tead1_bound_genes)
yaptaz_both_bound <- intersect(yaptaz_tes_bound, yaptaz_tead1_bound)
yaptaz_tes_only <- setdiff(yaptaz_tes_bound, tead1_bound_genes)
yaptaz_tead1_only <- setdiff(yaptaz_tead1_bound, tes_bound_genes)
yaptaz_no_binding <- setdiff(yaptaz_targets, union(tes_bound_genes, tead1_bound_genes))

cat("\nYAP/TAZ targets bound by TES:", length(yaptaz_tes_bound), "\n")
cat("YAP/TAZ targets bound by TEAD1:", length(yaptaz_tead1_bound), "\n")
cat("YAP/TAZ targets bound by both:", length(yaptaz_both_bound), "\n")
cat("YAP/TAZ targets bound by TES only:", length(yaptaz_tes_only), "\n")
cat("YAP/TAZ targets bound by TEAD1 only:", length(yaptaz_tead1_only), "\n")
cat("YAP/TAZ targets with no TES/TEAD1 binding:", length(yaptaz_no_binding), "\n")

# Save gene lists
write.table(yaptaz_tes_bound,
            file.path(base_dir, "gene_lists", "yaptaz_TES_bound.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(yaptaz_tead1_bound,
            file.path(base_dir, "gene_lists", "yaptaz_TEAD1_bound.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(yaptaz_both_bound,
            file.path(base_dir, "gene_lists", "yaptaz_both_bound.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(yaptaz_tes_only,
            file.path(base_dir, "gene_lists", "yaptaz_TES_only_bound.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(yaptaz_tead1_only,
            file.path(base_dir, "gene_lists", "yaptaz_TEAD1_only_bound.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(yaptaz_no_binding,
            file.path(base_dir, "gene_lists", "yaptaz_no_binding.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)

################################################################################
# HYPERGEOMETRIC ENRICHMENT TESTS
################################################################################

cat("\n=== Hypergeometric enrichment tests ===\n")

# Universe: all genes with TES or TEAD1 peaks
universe_genes <- union(tes_bound_genes, tead1_bound_genes)
universe_size <- length(universe_genes)

cat("Universe size (all genes with TES or TEAD1 peaks):", universe_size, "\n")

# Test 1: Are YAP/TAZ targets enriched in TES-bound genes?
if (length(yaptaz_tes_bound) > 0) {
  cat("\n--- TES binding enrichment ---\n")

  # Number of YAP/TAZ targets in universe
  yaptaz_in_universe <- length(intersect(yaptaz_targets, universe_genes))
  tes_bound_size <- length(tes_bound_genes)
  overlap <- length(yaptaz_tes_bound)

  phyper_tes <- phyper(q = overlap - 1,
                       m = yaptaz_in_universe,
                       n = universe_size - yaptaz_in_universe,
                       k = tes_bound_size,
                       lower.tail = FALSE)

  cat("YAP/TAZ targets in TES-bound genes:", overlap, "\n")
  cat("Expected by chance:", round((tes_bound_size * yaptaz_in_universe) / universe_size, 2), "\n")
  cat("Fold enrichment:", round(overlap / ((tes_bound_size * yaptaz_in_universe) / universe_size), 2), "\n")
  cat("Hypergeometric p-value:", format(phyper_tes, scientific = TRUE), "\n")
}

# Test 2: Are YAP/TAZ targets enriched in TEAD1-bound genes?
if (length(yaptaz_tead1_bound) > 0) {
  cat("\n--- TEAD1 binding enrichment ---\n")

  yaptaz_in_universe <- length(intersect(yaptaz_targets, universe_genes))
  tead1_bound_size <- length(tead1_bound_genes)
  overlap <- length(yaptaz_tead1_bound)

  phyper_tead1 <- phyper(q = overlap - 1,
                         m = yaptaz_in_universe,
                         n = universe_size - yaptaz_in_universe,
                         k = tead1_bound_size,
                         lower.tail = FALSE)

  cat("YAP/TAZ targets in TEAD1-bound genes:", overlap, "\n")
  cat("Expected by chance:", round((tead1_bound_size * yaptaz_in_universe) / universe_size, 2), "\n")
  cat("Fold enrichment:", round(overlap / ((tead1_bound_size * yaptaz_in_universe) / universe_size), 2), "\n")
  cat("Hypergeometric p-value:", format(phyper_tead1, scientific = TRUE), "\n")
}

################################################################################
# DETAILED ANNOTATION TABLES
################################################################################

cat("\n=== Creating detailed annotation tables ===\n")

# For YAP/TAZ targets bound by TES, get expression data if available
if (length(yaptaz_tes_bound) > 0) {
  yaptaz_tes_detailed <- data.frame(
    gene_symbol = yaptaz_tes_bound,
    stringsAsFactors = FALSE
  )

  # Add expression data if gene is in RNA-seq results
  yaptaz_tes_detailed <- yaptaz_tes_detailed %>%
    left_join(deseq_results %>% select(gene_symbol, baseMean, log2FoldChange, padj),
              by = "gene_symbol") %>%
    mutate(
      in_rnaseq_degs = !is.na(padj) & padj < 0.05,
      bound_by_tead1 = gene_symbol %in% tead1_bound_genes
    ) %>%
    arrange(desc(in_rnaseq_degs), padj)

  write_csv(yaptaz_tes_detailed,
            file.path(base_dir, "results", "yaptaz_TES_bound_detailed.csv"))

  cat("Saved detailed table for YAP/TAZ targets bound by TES\n")
}

# For YAP/TAZ targets bound by TEAD1
if (length(yaptaz_tead1_bound) > 0) {
  yaptaz_tead1_detailed <- data.frame(
    gene_symbol = yaptaz_tead1_bound,
    stringsAsFactors = FALSE
  )

  yaptaz_tead1_detailed <- yaptaz_tead1_detailed %>%
    left_join(deseq_results %>% select(gene_symbol, baseMean, log2FoldChange, padj),
              by = "gene_symbol") %>%
    mutate(
      in_rnaseq_degs = !is.na(padj) & padj < 0.05,
      bound_by_tes = gene_symbol %in% tes_bound_genes
    ) %>%
    arrange(desc(in_rnaseq_degs), padj)

  write_csv(yaptaz_tead1_detailed,
            file.path(base_dir, "results", "yaptaz_TEAD1_bound_detailed.csv"))

  cat("Saved detailed table for YAP/TAZ targets bound by TEAD1\n")
}

################################################################################
# GO ENRICHMENT ANALYSIS
################################################################################

cat("\n=== GO enrichment analysis ===\n")

# Enrichment for YAP/TAZ targets bound by TES
if (length(yaptaz_tes_bound) >= 10) {
  cat("\nRunning GO enrichment for YAP/TAZ targets bound by TES...\n")

  go_yaptaz_tes <- perform_GO_enrichment(
    gene_list = yaptaz_tes_bound,
    background_genes = all_expressed_genes,
    ont = "BP",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05
  )

  if (!is.null(go_yaptaz_tes) && nrow(go_yaptaz_tes) > 0) {
    cat("Found", nrow(go_yaptaz_tes), "enriched GO terms\n")
    save_enrichment_results_custom(go_yaptaz_tes, "yaptaz_TES_bound", base_dir)
  } else {
    cat("No significant GO enrichment found\n")
  }
} else {
  cat("Too few genes for GO enrichment (", length(yaptaz_tes_bound), ")\n")
}

# Enrichment for YAP/TAZ targets bound by TEAD1
if (length(yaptaz_tead1_bound) >= 10) {
  cat("\nRunning GO enrichment for YAP/TAZ targets bound by TEAD1...\n")

  go_yaptaz_tead1 <- perform_GO_enrichment(
    gene_list = yaptaz_tead1_bound,
    background_genes = all_expressed_genes,
    ont = "BP",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05
  )

  if (!is.null(go_yaptaz_tead1) && nrow(go_yaptaz_tead1) > 0) {
    cat("Found", nrow(go_yaptaz_tead1), "enriched GO terms\n")
    save_enrichment_results_custom(go_yaptaz_tead1, "yaptaz_TEAD1_bound", base_dir)
  } else {
    cat("No significant GO enrichment found\n")
  }
} else {
  cat("Too few genes for GO enrichment (", length(yaptaz_tead1_bound), ")\n")
}

# Enrichment for YAP/TAZ targets bound by both TES and TEAD1
if (length(yaptaz_both_bound) >= 10) {
  cat("\nRunning GO enrichment for YAP/TAZ targets bound by both TES and TEAD1...\n")

  go_yaptaz_both <- perform_GO_enrichment(
    gene_list = yaptaz_both_bound,
    background_genes = all_expressed_genes,
    ont = "BP",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05
  )

  if (!is.null(go_yaptaz_both) && nrow(go_yaptaz_both) > 0) {
    cat("Found", nrow(go_yaptaz_both), "enriched GO terms\n")
    save_enrichment_results_custom(go_yaptaz_both, "yaptaz_both_bound", base_dir)
  } else {
    cat("No significant GO enrichment found\n")
  }
} else {
  cat("Too few genes for GO enrichment (", length(yaptaz_both_bound), ")\n")
}

################################################################################
# VISUALIZATIONS
################################################################################

cat("\n=== Creating visualizations ===\n")

# Plot 1: Venn diagram showing overlap
pdf(file.path(base_dir, "plots", "yaptaz_tes_tead1_venn.pdf"), width = 8, height = 8)

venn_list <- list(
  "YAP/TAZ targets\n(n=224)" = yaptaz_targets,
  "TES-bound" = tes_bound_genes,
  "TEAD1-bound" = tead1_bound_genes
)

venn.diagram(
  x = venn_list,
  filename = NULL,
  category.names = names(venn_list),
  fill = c("#E69F00", "#56B4E9", "#009E73"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.2,
  cat.fontface = "bold",
  margin = 0.1
) %>% grid.draw()

dev.off()

# Plot 2: Bar plot showing binding categories
binding_summary <- data.frame(
  Category = c("TES only", "TEAD1 only", "Both TES & TEAD1", "No binding"),
  Count = c(
    length(yaptaz_tes_only),
    length(yaptaz_tead1_only),
    length(yaptaz_both_bound),
    length(yaptaz_no_binding)
  ),
  Percentage = c(
    length(yaptaz_tes_only) / length(yaptaz_targets) * 100,
    length(yaptaz_tead1_only) / length(yaptaz_targets) * 100,
    length(yaptaz_both_bound) / length(yaptaz_targets) * 100,
    length(yaptaz_no_binding) / length(yaptaz_targets) * 100
  )
)

pdf(file.path(base_dir, "plots", "yaptaz_binding_categories.pdf"), width = 10, height = 6)

p1 <- ggplot(binding_summary, aes(x = reorder(Category, -Count), y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = sprintf("%d\n(%.1f%%)", Count, Percentage)),
            vjust = -0.5, size = 4) +
  labs(title = "YAP/TAZ Target Genes: TES/TEAD1 Binding Categories",
       subtitle = sprintf("Total YAP/TAZ targets: %d", length(yaptaz_targets)),
       x = "Binding Category",
       y = "Number of Genes") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 11),
        plot.title = element_text(size = 14, face = "bold")) +
  scale_fill_brewer(palette = "Set2") +
  ylim(0, max(binding_summary$Count) * 1.15)

print(p1)
dev.off()

# Plot 3: Enrichment comparison (if both tests were significant)
if (exists("phyper_tes") && exists("phyper_tead1")) {
  enrichment_data <- data.frame(
    Factor = c("TES", "TEAD1"),
    Overlap = c(length(yaptaz_tes_bound), length(yaptaz_tead1_bound)),
    PValue = c(phyper_tes, phyper_tead1),
    NegLog10P = c(-log10(phyper_tes), -log10(phyper_tead1))
  )

  pdf(file.path(base_dir, "plots", "yaptaz_enrichment_comparison.pdf"), width = 8, height = 6)

  p2 <- ggplot(enrichment_data, aes(x = Factor, y = NegLog10P, fill = Factor)) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "darkred") +
    geom_text(aes(label = sprintf("p = %.2e\nn = %d", PValue, Overlap)),
              vjust = -0.5, size = 4) +
    labs(title = "YAP/TAZ Target Enrichment in TES vs TEAD1 Binding",
       subtitle = "Hypergeometric test",
       x = "Transcription Factor",
       y = "-log10(p-value)") +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 12),
          plot.title = element_text(size = 14, face = "bold")) +
    scale_fill_manual(values = c("TES" = "#56B4E9", "TEAD1" = "#009E73"))

  print(p2)
  dev.off()
}

# Plot 4: If we have expression data, show log2FC distribution
yaptaz_with_expression <- yaptaz_targets[yaptaz_targets %in% deseq_results$gene_symbol]

if (length(yaptaz_with_expression) > 0) {
  expression_data <- deseq_results %>%
    filter(gene_symbol %in% yaptaz_targets) %>%
    mutate(
      binding_status = case_when(
        gene_symbol %in% yaptaz_both_bound ~ "Both TES & TEAD1",
        gene_symbol %in% yaptaz_tes_only ~ "TES only",
        gene_symbol %in% yaptaz_tead1_only ~ "TEAD1 only",
        TRUE ~ "No binding"
      ),
      is_deg = !is.na(padj) & padj < 0.05
    )

  pdf(file.path(base_dir, "plots", "yaptaz_expression_by_binding.pdf"), width = 10, height = 6)

  p3 <- ggplot(expression_data, aes(x = binding_status, y = log2FoldChange, fill = binding_status)) +
    geom_boxplot(alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    labs(title = "Expression Changes of YAP/TAZ Targets by TES/TEAD1 Binding Status",
         subtitle = "TES vs GFP (RNA-seq)",
         x = "Binding Status",
         y = "log2 Fold Change (TES vs GFP)") +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 11),
          plot.title = element_text(size = 14, face = "bold")) +
    scale_fill_brewer(palette = "Set2")

  print(p3)
  dev.off()
}

################################################################################
# SUMMARY REPORT
################################################################################

cat("\n=== Generating summary report ===\n")

summary_text <- paste0(
  "=================================================================\n",
  "YAP/TAZ TARGET ENRICHMENT ANALYSIS SUMMARY\n",
  "=================================================================\n\n",
  "Total YAP/TAZ target genes: ", length(yaptaz_targets), "\n",
  "YAP/TAZ targets with expression data: ", length(yaptaz_with_expression), "\n\n",
  "--- Binding Analysis ---\n",
  "Total genes with TES peaks: ", length(tes_bound_genes), "\n",
  "Total genes with TEAD1 peaks: ", length(tead1_bound_genes), "\n\n",
  "YAP/TAZ targets bound by TES: ", length(yaptaz_tes_bound), " (",
  round(length(yaptaz_tes_bound)/length(yaptaz_targets)*100, 1), "%)\n",
  "YAP/TAZ targets bound by TEAD1: ", length(yaptaz_tead1_bound), " (",
  round(length(yaptaz_tead1_bound)/length(yaptaz_targets)*100, 1), "%)\n",
  "YAP/TAZ targets bound by both: ", length(yaptaz_both_bound), " (",
  round(length(yaptaz_both_bound)/length(yaptaz_targets)*100, 1), "%)\n",
  "YAP/TAZ targets bound by TES only: ", length(yaptaz_tes_only), "\n",
  "YAP/TAZ targets bound by TEAD1 only: ", length(yaptaz_tead1_only), "\n",
  "YAP/TAZ targets with no binding: ", length(yaptaz_no_binding), "\n\n"
)

if (exists("phyper_tes")) {
  summary_text <- paste0(
    summary_text,
    "--- Enrichment Statistics ---\n",
    "TES binding enrichment p-value: ", format(phyper_tes, scientific = TRUE), "\n"
  )
}

if (exists("phyper_tead1")) {
  summary_text <- paste0(
    summary_text,
    "TEAD1 binding enrichment p-value: ", format(phyper_tead1, scientific = TRUE), "\n\n"
  )
}

# Check overlap with RNA-seq DEGs
yaptaz_in_degs <- intersect(yaptaz_targets, deseq_results$gene_symbol[!is.na(deseq_results$padj) & deseq_results$padj < 0.05])

summary_text <- paste0(
  summary_text,
  "--- Overlap with RNA-seq DEGs ---\n",
  "YAP/TAZ targets that are also DEGs (padj < 0.05): ", length(yaptaz_in_degs), " (",
  round(length(yaptaz_in_degs)/length(yaptaz_targets)*100, 1), "%)\n",
  "YAP/TAZ DEGs bound by TES: ", length(intersect(yaptaz_in_degs, tes_bound_genes)), "\n",
  "YAP/TAZ DEGs bound by TEAD1: ", length(intersect(yaptaz_in_degs, tead1_bound_genes)), "\n",
  "YAP/TAZ DEGs bound by both: ", length(intersect(yaptaz_in_degs, intersect(tes_bound_genes, tead1_bound_genes))), "\n\n",
  "=================================================================\n",
  "Analysis complete: ", date(), "\n",
  "=================================================================\n"
)

cat(summary_text)

writeLines(summary_text, file.path(base_dir, "results", "analysis_summary.txt"))

cat("\n=================================================================\n")
cat("YAP/TAZ TARGET ANALYSIS COMPLETE\n")
cat("Results saved to:", base_dir, "\n")
cat("=================================================================\n")
