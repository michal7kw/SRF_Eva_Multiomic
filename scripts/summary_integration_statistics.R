#!/usr/bin/env Rscript
################################################################################
# Multi-Omics Integration Summary Statistics
# Combines Binding + Methylation + Expression data summaries
#
# Author: Generated for SRF_Eva project
# Date: 2024-12
#
# Usage: Rscript summary_integration_statistics.R
# Output: results/summary/integration_summary_statistics.csv
#         results/summary/integration_top_genes_by_category.csv
################################################################################

suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
})

cat("========================================\n")
cat("Multi-Omics Integration Summary Generator\n")
cat("========================================\n\n")

# Define paths
BASE_DIR <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top"
INTEGRATION_DIR <- file.path(BASE_DIR, "SRF_Eva_integrated_analysis")
MEDIP_DIR <- file.path(BASE_DIR, "meDIP")
RNA_DIR <- file.path(BASE_DIR, "SRF_Eva_RNA")
CUTNTAG_DIR <- file.path(BASE_DIR, "SRF_Eva_CUTandTAG")

OUTPUT_DIR <- file.path(INTEGRATION_DIR, "results/summary")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

################################################################################
# 1. Load All Data Sources
################################################################################

cat("1. Loading data sources...\n")

# --- Gene Classification (Binding + Expression) ---
classification_file <- file.path(INTEGRATION_DIR,
    "results/20_methylation_binding_expression/gene_classification_summary.csv")

if (file.exists(classification_file)) {
    classification <- read.csv(classification_file, stringsAsFactors = FALSE)
    cat(sprintf("   Loaded gene classification: %d genes\n", nrow(classification)))
} else {
    classification <- NULL
    cat("   WARNING: Gene classification file not found\n")
}

# --- meDIP Integration Data ---
medip_integration_file <- file.path(MEDIP_DIR,
    "results/16_advanced_visualization/integrated_medip_rnaseq_data.csv")

if (file.exists(medip_integration_file)) {
    medip_data <- read.csv(medip_integration_file, stringsAsFactors = FALSE)
    cat(sprintf("   Loaded meDIP integration: %d genes\n", nrow(medip_data)))
} else {
    medip_data <- NULL
    cat("   WARNING: meDIP integration file not found\n")
}

# --- DMR Data ---
dmr_file <- file.path(MEDIP_DIR, "results/07_differential_MEDIPS/TES_vs_GFP_DMRs_FDR05.csv")
if (file.exists(dmr_file)) {
    dmrs <- read.csv(dmr_file, stringsAsFactors = FALSE)
    cat(sprintf("   Loaded DMRs: %d regions\n", nrow(dmrs)))
} else {
    dmrs <- NULL
}

# --- Annotated DMRs ---
annotated_dmr_file <- file.path(MEDIP_DIR, "results/08_annotation/TES_vs_GFP_annotated.csv")
if (file.exists(annotated_dmr_file)) {
    annotated_dmrs <- read.csv(annotated_dmr_file, stringsAsFactors = FALSE)
    cat(sprintf("   Loaded annotated DMRs: %d regions\n", nrow(annotated_dmrs)))
} else {
    annotated_dmrs <- NULL
}

# --- Gene Set Summary ---
geneset_file <- file.path(MEDIP_DIR, "results/12_gene_sets/gene_set_summary.csv")
if (file.exists(geneset_file)) {
    geneset_summary <- read.csv(geneset_file, stringsAsFactors = FALSE)
    cat(sprintf("   Loaded gene set summary\n"))
} else {
    geneset_summary <- NULL
}

################################################################################
# 2. Initialize Summary Statistics
################################################################################

cat("\n2. Calculating summary statistics...\n")

summary_stats <- data.frame(
    Category = character(),
    Metric = character(),
    Value = character(),
    stringsAsFactors = FALSE
)

add_stat <- function(category, metric, value) {
    summary_stats <<- rbind(summary_stats, data.frame(
        Category = category,
        Metric = metric,
        Value = as.character(value),
        stringsAsFactors = FALSE
    ))
}

################################################################################
# 3. Binding + Expression Classification Statistics
################################################################################

if (!is.null(classification)) {
    cat("   Calculating binding + expression statistics...\n")

    # Total counts
    total_genes <- nrow(classification)
    add_stat("Overview", "Total_genes_analyzed", total_genes)

    # Expression status
    expr_counts <- table(classification$expression_status)
    for (status in names(expr_counts)) {
        add_stat("Expression_Status", paste0("N_", status), expr_counts[status])
        add_stat("Expression_Status", paste0("Pct_", status),
                 sprintf("%.1f%%", 100 * expr_counts[status] / total_genes))
    }

    # Binding status
    binding_counts <- table(classification$binding_status)
    for (status in names(binding_counts)) {
        add_stat("Binding_Status", paste0("N_", status), binding_counts[status])
        add_stat("Binding_Status", paste0("Pct_", status),
                 sprintf("%.1f%%", 100 * binding_counts[status] / total_genes))
    }

    # Combined classification
    combined_counts <- table(classification$combined_class)
    for (combo in names(combined_counts)) {
        add_stat("Combined_Classification", combo, combined_counts[combo])
    }

    # Direct target summary (bound + DE)
    degs <- classification %>% filter(expression_status %in% c("upregulated", "downregulated"))
    n_degs <- nrow(degs)
    add_stat("DEG_Summary", "Total_DEGs", n_degs)

    if (n_degs > 0) {
        # DEGs with binding
        tes_bound_degs <- sum(degs$tes_bound, na.rm = TRUE)
        tead1_bound_degs <- sum(degs$tead1_bound, na.rm = TRUE)
        both_bound_degs <- sum(degs$tes_bound & degs$tead1_bound, na.rm = TRUE)
        any_bound_degs <- sum(degs$tes_bound | degs$tead1_bound, na.rm = TRUE)
        unbound_degs <- n_degs - any_bound_degs

        add_stat("Direct_Targets", "DEGs_with_TES_binding", tes_bound_degs)
        add_stat("Direct_Targets", "DEGs_with_TEAD1_binding", tead1_bound_degs)
        add_stat("Direct_Targets", "DEGs_with_TES_and_TEAD1_binding", both_bound_degs)
        add_stat("Direct_Targets", "DEGs_with_any_binding", any_bound_degs)
        add_stat("Direct_Targets", "DEGs_without_binding_indirect", unbound_degs)

        add_stat("Direct_Targets", "Pct_DEGs_TES_bound",
                 sprintf("%.1f%%", 100 * tes_bound_degs / n_degs))
        add_stat("Direct_Targets", "Pct_DEGs_any_bound",
                 sprintf("%.1f%%", 100 * any_bound_degs / n_degs))
        add_stat("Direct_Targets", "Pct_DEGs_indirect",
                 sprintf("%.1f%%", 100 * unbound_degs / n_degs))

        # Direction breakdown
        up_degs <- degs %>% filter(expression_status == "upregulated")
        down_degs <- degs %>% filter(expression_status == "downregulated")

        add_stat("Direction_Summary", "N_upregulated", nrow(up_degs))
        add_stat("Direction_Summary", "N_downregulated", nrow(down_degs))
        add_stat("Direction_Summary", "Up_to_Down_ratio",
                 sprintf("%.2f", nrow(up_degs) / max(nrow(down_degs), 1)))

        # Bound + direction
        add_stat("Direction_by_Binding", "TES_bound_upregulated",
                 sum(up_degs$tes_bound, na.rm = TRUE))
        add_stat("Direction_by_Binding", "TES_bound_downregulated",
                 sum(down_degs$tes_bound, na.rm = TRUE))
        add_stat("Direction_by_Binding", "TEAD1_bound_upregulated",
                 sum(up_degs$tead1_bound, na.rm = TRUE))
        add_stat("Direction_by_Binding", "TEAD1_bound_downregulated",
                 sum(down_degs$tead1_bound, na.rm = TRUE))
    }
}

################################################################################
# 4. Methylation Statistics
################################################################################

cat("   Calculating methylation statistics...\n")

if (!is.null(dmrs)) {
    add_stat("DMR_Overview", "Total_DMRs", nrow(dmrs))

    hyper <- sum(dmrs$logFC > 0, na.rm = TRUE)
    hypo <- sum(dmrs$logFC < 0, na.rm = TRUE)

    add_stat("DMR_Overview", "Hypermethylated_DMRs", hyper)
    add_stat("DMR_Overview", "Hypomethylated_DMRs", hypo)
    add_stat("DMR_Overview", "Pct_hypermethylated", sprintf("%.1f%%", 100 * hyper / nrow(dmrs)))
    add_stat("DMR_Overview", "Hyper_to_Hypo_ratio", sprintf("%.1f:1", hyper / max(hypo, 1)))

    # Effect sizes
    add_stat("DMR_Effect_Size", "Mean_logFC", sprintf("%.3f", mean(dmrs$logFC, na.rm = TRUE)))
    add_stat("DMR_Effect_Size", "Median_logFC", sprintf("%.3f", median(dmrs$logFC, na.rm = TRUE)))
    add_stat("DMR_Effect_Size", "Max_fold_change", sprintf("%.1f", max(dmrs$fold_change, na.rm = TRUE)))
}

# Genomic distribution
if (!is.null(annotated_dmrs)) {
    promoter_dmrs <- sum(grepl("Promoter", annotated_dmrs$annotation), na.rm = TRUE)
    genebody_dmrs <- sum(grepl("Intron|Exon|UTR", annotated_dmrs$annotation), na.rm = TRUE)
    intergenic_dmrs <- sum(grepl("Intergenic", annotated_dmrs$annotation), na.rm = TRUE)

    add_stat("DMR_Genomic_Distribution", "Promoter_DMRs", promoter_dmrs)
    add_stat("DMR_Genomic_Distribution", "Gene_body_DMRs", genebody_dmrs)
    add_stat("DMR_Genomic_Distribution", "Intergenic_DMRs", intergenic_dmrs)
    add_stat("DMR_Genomic_Distribution", "Pct_Promoter",
             sprintf("%.1f%%", 100 * promoter_dmrs / nrow(annotated_dmrs)))

    # Unique genes with DMRs
    unique_dmr_genes <- length(unique(na.omit(annotated_dmrs$SYMBOL)))
    add_stat("DMR_Gene_Coverage", "Unique_genes_with_DMRs", unique_dmr_genes)
}

################################################################################
# 5. Methylation-Expression Correlation
################################################################################

if (!is.null(medip_data)) {
    cat("   Calculating methylation-expression correlation...\n")

    # Filter for valid data
    valid_data <- medip_data %>%
        filter(!is.na(meDIP_log2FC) & !is.na(log2FoldChange))

    if (nrow(valid_data) > 10) {
        # Overall correlation
        pearson <- cor.test(valid_data$meDIP_log2FC, valid_data$log2FoldChange, method = "pearson")
        spearman <- cor.test(valid_data$meDIP_log2FC, valid_data$log2FoldChange, method = "spearman")

        add_stat("Meth_Expr_Correlation", "Pearson_r", sprintf("%.4f", pearson$estimate))
        add_stat("Meth_Expr_Correlation", "Pearson_pvalue", sprintf("%.2e", pearson$p.value))
        add_stat("Meth_Expr_Correlation", "Spearman_rho", sprintf("%.4f", spearman$estimate))
        add_stat("Meth_Expr_Correlation", "Spearman_pvalue", sprintf("%.2e", spearman$p.value))

        # Methylation by expression category
        if ("padj" %in% colnames(valid_data)) {
            up_genes <- valid_data %>% filter(padj < 0.05 & log2FoldChange > 1)
            down_genes <- valid_data %>% filter(padj < 0.05 & log2FoldChange < -1)
            unchanged <- valid_data %>% filter(is.na(padj) | padj >= 0.05 | abs(log2FoldChange) <= 1)

            if (nrow(up_genes) > 5 && nrow(down_genes) > 5) {
                add_stat("Meth_by_Expression", "Mean_meDIP_delta_upregulated",
                         sprintf("%.3f", mean(up_genes$meDIP_delta, na.rm = TRUE)))
                add_stat("Meth_by_Expression", "Mean_meDIP_delta_downregulated",
                         sprintf("%.3f", mean(down_genes$meDIP_delta, na.rm = TRUE)))
                add_stat("Meth_by_Expression", "Mean_meDIP_delta_unchanged",
                         sprintf("%.3f", mean(unchanged$meDIP_delta, na.rm = TRUE)))

                # Wilcoxon test between up and down
                wtest <- wilcox.test(up_genes$meDIP_delta, down_genes$meDIP_delta)
                add_stat("Meth_by_Expression", "Wilcoxon_Up_vs_Down_pvalue",
                         sprintf("%.2e", wtest$p.value))
            }
        }
    }
}

################################################################################
# 6. DMR-DEG Enrichment Analysis
################################################################################

cat("   Calculating DMR-DEG enrichment...\n")

if (!is.null(classification) && !is.null(annotated_dmrs)) {
    # Get genes with DMRs
    dmr_genes <- unique(na.omit(annotated_dmrs$SYMBOL))

    # Get DEGs and non-DEGs
    degs_symbols <- classification %>%
        filter(expression_status %in% c("upregulated", "downregulated")) %>%
        pull(gene_symbol) %>%
        unique() %>%
        na.omit()

    non_degs_symbols <- classification %>%
        filter(!expression_status %in% c("upregulated", "downregulated")) %>%
        pull(gene_symbol) %>%
        unique() %>%
        na.omit()

    # Calculate overlap
    degs_with_dmrs <- length(intersect(degs_symbols, dmr_genes))
    degs_without_dmrs <- length(setdiff(degs_symbols, dmr_genes))
    non_degs_with_dmrs <- length(intersect(non_degs_symbols, dmr_genes))
    non_degs_without_dmrs <- length(setdiff(non_degs_symbols, dmr_genes))

    add_stat("DMR_DEG_Enrichment", "DEGs_with_DMRs", degs_with_dmrs)
    add_stat("DMR_DEG_Enrichment", "DEGs_without_DMRs", degs_without_dmrs)
    add_stat("DMR_DEG_Enrichment", "NonDEGs_with_DMRs", non_degs_with_dmrs)
    add_stat("DMR_DEG_Enrichment", "NonDEGs_without_DMRs", non_degs_without_dmrs)

    pct_degs_with_dmrs <- 100 * degs_with_dmrs / (degs_with_dmrs + degs_without_dmrs)
    pct_nondegs_with_dmrs <- 100 * non_degs_with_dmrs / (non_degs_with_dmrs + non_degs_without_dmrs)

    add_stat("DMR_DEG_Enrichment", "Pct_DEGs_with_DMRs", sprintf("%.1f%%", pct_degs_with_dmrs))
    add_stat("DMR_DEG_Enrichment", "Pct_NonDEGs_with_DMRs", sprintf("%.1f%%", pct_nondegs_with_dmrs))

    # Fisher's exact test
    contingency <- matrix(c(degs_with_dmrs, degs_without_dmrs,
                           non_degs_with_dmrs, non_degs_without_dmrs),
                         nrow = 2, byrow = TRUE)

    fisher_result <- fisher.test(contingency)
    add_stat("DMR_DEG_Enrichment", "Fisher_odds_ratio", sprintf("%.2f", fisher_result$estimate))
    add_stat("DMR_DEG_Enrichment", "Fisher_pvalue", sprintf("%.2e", fisher_result$p.value))
    add_stat("DMR_DEG_Enrichment", "Fisher_95CI_low", sprintf("%.2f", fisher_result$conf.int[1]))
    add_stat("DMR_DEG_Enrichment", "Fisher_95CI_high", sprintf("%.2f", fisher_result$conf.int[2]))
}

################################################################################
# 7. Generate Top Genes by Category
################################################################################

cat("\n3. Generating top genes by category...\n")

top_genes_list <- list()

if (!is.null(classification)) {
    # Most significant DEGs by binding category
    categories <- c("TES_only_bound", "TEAD1_only_bound", "TES_TEAD1_bound", "Neither_bound")

    for (cat in categories) {
        # Upregulated
        up_genes <- classification %>%
            filter(binding_status == cat & expression_status == "upregulated") %>%
            arrange(padj) %>%
            head(50) %>%
            mutate(Category = paste0(cat, "_upregulated"),
                   Direction = "Up")

        # Downregulated
        down_genes <- classification %>%
            filter(binding_status == cat & expression_status == "downregulated") %>%
            arrange(padj) %>%
            head(50) %>%
            mutate(Category = paste0(cat, "_downregulated"),
                   Direction = "Down")

        if (nrow(up_genes) > 0) top_genes_list[[paste0(cat, "_up")]] <- up_genes
        if (nrow(down_genes) > 0) top_genes_list[[paste0(cat, "_down")]] <- down_genes
    }

    # Combine all top genes
    if (length(top_genes_list) > 0) {
        top_genes <- bind_rows(top_genes_list) %>%
            select(gene_symbol, ensembl_id, baseMean, log2FoldChange, padj,
                   binding_status, expression_status, Category, Direction)
        cat(sprintf("   Generated top gene lists: %d entries\n", nrow(top_genes)))
    }
}

################################################################################
# 8. Save Results
################################################################################

cat("\n4. Saving results...\n")

# Save summary statistics
summary_file <- file.path(OUTPUT_DIR, "integration_summary_statistics.csv")
write.csv(summary_stats, summary_file, row.names = FALSE)
cat(sprintf("   Saved: %s\n", summary_file))

# Save top genes
if (exists("top_genes") && nrow(top_genes) > 0) {
    top_genes_file <- file.path(OUTPUT_DIR, "integration_top_genes_by_category.csv")
    write.csv(top_genes, top_genes_file, row.names = FALSE)
    cat(sprintf("   Saved: %s\n", top_genes_file))
}

# Create wide format summary
summary_wide <- summary_stats %>%
    unite(Full_Metric, Category, Metric, sep = "__") %>%
    pivot_wider(names_from = Full_Metric, values_from = Value)

summary_wide_file <- file.path(OUTPUT_DIR, "integration_summary_wide.csv")
write.csv(summary_wide, summary_wide_file, row.names = FALSE)
cat(sprintf("   Saved: %s\n", summary_wide_file))

################################################################################
# 9. Print Summary to Console
################################################################################

cat("\n========================================\n")
cat("KEY INTEGRATION STATISTICS\n")
cat("========================================\n\n")

key_categories <- c("Overview", "Direct_Targets", "DMR_Overview",
                   "Meth_Expr_Correlation", "DMR_DEG_Enrichment")

for (cat in key_categories) {
    cat_stats <- summary_stats %>% filter(Category == cat)
    if (nrow(cat_stats) > 0) {
        cat(sprintf("\n%s:\n", cat))
        for (i in 1:min(nrow(cat_stats), 10)) {
            cat(sprintf("  %s: %s\n", cat_stats$Metric[i], cat_stats$Value[i]))
        }
    }
}

cat("\n========================================\n")
cat("Integration summary complete!\n")
cat(sprintf("Output directory: %s\n", OUTPUT_DIR))
cat("========================================\n")
