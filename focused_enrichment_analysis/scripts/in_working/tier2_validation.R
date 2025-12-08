#!/usr/bin/env Rscript

cat("=================================================================\n")
cat("TIER 2: Validation with Orthogonal Approach\n")
cat("=================================================================\n")

# Load workspace and helper functions
source("SRF_Eva_integrated_analysis/focused_enrichment_analysis/scripts/00_helper_functions.R")
load_libraries()
setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top")
load("SRF_Eva_integrated_analysis/focused_enrichment_analysis/data/loaded_data.RData")

# Output directory
base_dir <- "SRF_Eva_integrated_analysis/focused_enrichment_analysis"
tier2_dir <- file.path(base_dir, "tier2_validation")

# Create directories
dir.create(file.path(tier2_dir, "results"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(tier2_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(tier2_dir, "gene_lists"), recursive = TRUE, showWarnings = FALSE)

################################################################################
# LOAD GENE LISTS
################################################################################

cat("\n=== Loading gene lists ===\n")

# Load direct downregulated genes
tes_direct_down_genes <- tes_direct %>%
  filter(log2FoldChange < 0, padj < 0.05) %>%
  pull(gene_symbol)

# Load diffbind down genes
diffbind_file <- file.path(base_dir, "approach5_diffbind", "gene_lists", "TES_diffbind_down.txt")
if (file.exists(diffbind_file)) {
  tes_diffbind_down <- read.table(diffbind_file, stringsAsFactors = FALSE)$V1
  cat("Loaded", length(tes_diffbind_down), "DiffBind + Down genes\n")
} else {
  cat("WARNING: DiffBind + Down gene list not found. Run approach5 first.\n")
  tes_diffbind_down <- character(0)
}

################################################################################
# VALIDATION ANALYSIS
################################################################################

if (length(tes_diffbind_down) > 0) {

  cat("\n=== Validation Analysis ===\n")
  cat("Direct downregulated genes:", length(tes_direct_down_genes), "\n")
  cat("DiffBind + Down genes:", length(tes_diffbind_down), "\n")

  # Venn diagram: Direct Down vs DiffBind Down
  venn_data <- list(
    "Direct\nDownregulated" = tes_direct_down_genes,
    "DiffBind\n+ Down" = tes_diffbind_down
  )

  pdf(file.path(tier2_dir, "plots", "validation_venn.pdf"), width = 8, height = 8)
  venn.plot <- venn.diagram(
    x = venn_data,
    filename = NULL,
    fill = c("lightblue", "pink"),
    alpha = 0.5,
    category.names = names(venn_data),
    main = "Validation: Direct vs DiffBind Approaches"
  )
  grid.draw(venn.plot)
  dev.off()

  # High-confidence genes (in both)
  high_conf_validated <- intersect(tes_direct_down_genes, tes_diffbind_down)
  cat("High-confidence validated genes:", length(high_conf_validated), "\n")

  write.table(high_conf_validated,
              file.path(tier2_dir, "gene_lists", "validated_highconf_genes.txt"),
              row.names = FALSE, col.names = FALSE, quote = FALSE)

  # Detailed table
  validated_detailed <- tes_direct %>%
    filter(gene_symbol %in% high_conf_validated) %>%
    arrange(padj)

  write_csv(validated_detailed,
            file.path(tier2_dir, "gene_lists", "validated_highconf_detailed.csv"))

  # Summary statistics
  summary_stats <- data.frame(
    Category = c("Direct Downregulated Only",
                 "DiffBind + Down Only",
                 "High-Confidence (Both)"),
    N_genes = c(
      length(setdiff(tes_direct_down_genes, tes_diffbind_down)),
      length(setdiff(tes_diffbind_down, tes_direct_down_genes)),
      length(high_conf_validated)
    )
  )

  write_csv(summary_stats,
            file.path(tier2_dir, "results", "validation_summary.csv"))

  cat("\n=== Validation Summary ===\n")
  print(summary_stats)

  # Enrichment on validated genes
  if (length(high_conf_validated) >= 10) {
    cat("\nRunning GO enrichment on validated high-confidence genes...\n")
    ego_tier2 <- perform_GO_enrichment(
      gene_list = high_conf_validated,
      background_genes = all_expressed_genes,
      ont = "BP"
    )

    results_tier2 <- save_enrichment_results(ego_tier2, "tier2_validated", tier2_dir)

    # Extract migration terms
    if (!is.null(results_tier2)) {
      migration_tier2 <- extract_migration_terms(results_tier2)
      if (!is.null(migration_tier2) && nrow(migration_tier2) > 0) {
        write_csv(migration_tier2,
                  file.path(tier2_dir, "results", "tier2_migration_terms.csv"))
        cat("  Found", nrow(migration_tier2), "migration-related terms\n")
        cat("  Top migration term rank:", which(results_tier2$ID == migration_tier2$ID[1]), "\n")
      }
    }
  } else {
    cat("\nToo few validated genes for GO enrichment (need >=10)\n")
  }

  ################################################################################
  # COMPARISON PLOTS
  ################################################################################

  cat("\n=== Creating comparison plots ===\n")

  # Plot 1: Overlap proportions
  overlap_data <- data.frame(
    Approach = c("Direct\nDownregulated", "DiffBind\n+ Down"),
    Total = c(length(tes_direct_down_genes), length(tes_diffbind_down)),
    Validated = c(length(high_conf_validated), length(high_conf_validated)),
    Unique = c(
      length(setdiff(tes_direct_down_genes, tes_diffbind_down)),
      length(setdiff(tes_diffbind_down, tes_direct_down_genes))
    )
  ) %>%
    mutate(
      PercentValidated = (Validated / Total) * 100,
      PercentUnique = (Unique / Total) * 100
    )

  # Stacked bar plot
  overlap_long <- overlap_data %>%
    select(Approach, Validated, Unique) %>%
    tidyr::pivot_longer(cols = c(Validated, Unique),
                        names_to = "Category",
                        values_to = "Count")

  pdf(file.path(tier2_dir, "plots", "validation_overlap_proportions.pdf"), width = 8, height = 6)
  p1 <- ggplot(overlap_long, aes(x = Approach, y = Count, fill = Category)) +
    geom_bar(stat = "identity", position = "stack") +
    geom_text(data = overlap_data,
              aes(x = Approach, y = Total, label = Total, fill = NULL),
              vjust = -0.5, size = 4) +
    scale_fill_manual(values = c("Validated" = "darkgreen", "Unique" = "lightgray"),
                      labels = c("Validated" = "High-Confidence (Both)",
                                "Unique" = "Approach-Specific")) +
    labs(title = "Validation Analysis: Overlap Between Approaches",
         subtitle = "Numbers show total genes per approach",
         x = "Approach",
         y = "Number of Genes",
         fill = "") +
    theme_bw() +
    theme(plot.title = element_text(size = 14, face = "bold"),
          legend.position = "bottom")
  print(p1)
  dev.off()

  # Plot 2: Expression distribution comparison
  if (length(high_conf_validated) > 0) {
    expression_data <- tes_direct %>%
      filter(gene_symbol %in% c(tes_direct_down_genes, tes_diffbind_down)) %>%
      mutate(Category = case_when(
        gene_symbol %in% high_conf_validated ~ "High-Confidence (Both)",
        gene_symbol %in% tes_direct_down_genes ~ "Direct Down Only",
        gene_symbol %in% tes_diffbind_down ~ "DiffBind Down Only",
        TRUE ~ "Other"
      ))

    pdf(file.path(tier2_dir, "plots", "validation_expression_distribution.pdf"),
        width = 10, height = 6)
    p2 <- ggplot(expression_data, aes(x = Category, y = log2FoldChange, fill = Category)) +
      geom_boxplot() +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
      scale_fill_brewer(palette = "Set2") +
      labs(title = "Expression Changes by Validation Category",
           subtitle = "High-confidence genes identified by both approaches",
           x = "",
           y = "log2 Fold Change (TES vs GFP)") +
      theme_bw() +
      theme(plot.title = element_text(size = 14, face = "bold"),
            legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1))
    print(p2)
    dev.off()

    # Plot 3: Top validated genes heatmap
    if (nrow(validated_detailed) >= 10) {
      top_validated <- validated_detailed %>%
        head(20)

      pdf(file.path(tier2_dir, "plots", "top_validated_genes_heatmap.pdf"),
          width = 6, height = 8)
      p3 <- ggplot(top_validated, aes(x = 1, y = reorder(gene_symbol, -log2FoldChange),
                                      fill = log2FoldChange)) +
        geom_tile(color = "white") +
        geom_text(aes(label = sprintf("%.2f", log2FoldChange)), size = 3) +
        scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                             name = "log2FC") +
        labs(title = "Top 20 Validated High-Confidence Genes",
             subtitle = "Confirmed by both Direct and DiffBind approaches",
             x = "", y = "") +
        theme_minimal() +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              plot.title = element_text(size = 14, face = "bold"))
      print(p3)
      dev.off()
    }
  }

  cat("Plots created successfully!\n")

} else {
  cat("\nWARNING: Cannot perform validation analysis without DiffBind results.\n")
  cat("Please run approach5_diffbind.R first.\n")
}

cat("\n=================================================================\n")
cat("TIER 2 COMPLETE\n")
cat("=================================================================\n")
