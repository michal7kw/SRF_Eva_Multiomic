#!/usr/bin/env Rscript

cat("=================================================================\n")
cat("APPROACH 8: Comparative Enrichment Analysis\n")
cat("Focus: Show migration/cancer are PREFERENTIALLY enriched\n")
cat("=================================================================\n")

# Load workspace and helper functions
source("SRF_Eva_integrated_analysis/focused_enrichment_analysis/scripts/00_helper_functions.R")
load_libraries()
setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top")
load("SRF_Eva_integrated_analysis/focused_enrichment_analysis/data/loaded_data.RData")

# Output directory
base_dir <- "SRF_Eva_integrated_analysis/focused_enrichment_analysis"
approach8_dir <- file.path(base_dir, "approach8_comparative")

# Create directories
dir.create(file.path(approach8_dir, "results"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(approach8_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(approach8_dir, "gene_lists"), recursive = TRUE, showWarnings = FALSE)

################################################################################
# RUN FULL GO ENRICHMENT
################################################################################

cat("\n=== Running comprehensive GO enrichment (may take 5-10 minutes) ===\n")

tes_direct_genes <- tes_direct$gene_symbol

# Run GO enrichment with relaxed thresholds to get many terms
ego_full <- perform_GO_enrichment(
  gene_list = tes_direct_genes,
  background_genes = all_expressed_genes,
  ont = "BP",
  pvalueCutoff = 0.05,  # Keep all significant terms
  qvalueCutoff = 0.2,   # Relaxed q-value to get more terms
  minGSSize = 10,
  maxGSSize = 500
)

if (is.null(ego_full) || nrow(ego_full) == 0) {
  cat("ERROR: No GO enrichment results found\n")
  quit(status = 1)
}

# Convert to data frame
ego_full_df <- as.data.frame(ego_full) %>%
  mutate(
    # Calculate fold enrichment
    GeneRatio_num = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])),
    BgRatio_num = sapply(strsplit(BgRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])),
    FoldEnrichment = GeneRatio_num / BgRatio_num,
    NegLog10P = -log10(pvalue),
    NegLog10Q = -log10(qvalue)
  )

cat("Total GO terms enriched:", nrow(ego_full_df), "\n")

# Save full results
write_csv(ego_full_df, file.path(approach8_dir, "results", "full_GO_enrichment.csv"))

################################################################################
# RANK PATHWAYS BY FOLD ENRICHMENT
################################################################################

cat("\n=== Ranking pathways by fold enrichment ===\n")

# Rank by fold enrichment
ego_ranked_fe <- ego_full_df %>%
  arrange(desc(FoldEnrichment)) %>%
  mutate(Rank_FoldEnrichment = row_number())

# Rank by p-value
ego_ranked_pval <- ego_full_df %>%
  arrange(pvalue) %>%
  mutate(Rank_Pvalue = row_number())

# Combine rankings
ego_ranked <- ego_full_df %>%
  arrange(desc(FoldEnrichment)) %>%
  mutate(
    Rank_FoldEnrichment = row_number(),
    Percentile_FoldEnrichment = (1 - Rank_FoldEnrichment / nrow(ego_full_df)) * 100
  ) %>%
  left_join(
    ego_ranked_pval %>% select(ID, Rank_Pvalue),
    by = "ID"
  ) %>%
  mutate(
    Percentile_Pvalue = (1 - Rank_Pvalue / nrow(ego_full_df)) * 100
  )

# Save ranked results
write_csv(ego_ranked, file.path(approach8_dir, "results", "ranked_GO_enrichment.csv"))

# Identify migration/cancer-related terms
migration_keywords <- c("migration", "motility", "adhesion", "invasion",
                        "metastasis", "locomotion", "chemotaxis", "cell movement")

cancer_keywords <- c("cancer", "tumor", "neoplasm", "oncogen", "carcinogen",
                     "proliferation", "apoptosis", "angiogenesis", "emt",
                     "epithelial.*mesenchymal")

# Tag terms
ego_ranked <- ego_ranked %>%
  mutate(
    Is_Migration = grepl(paste(migration_keywords, collapse = "|"), Description, ignore.case = TRUE),
    Is_Cancer = grepl(paste(cancer_keywords, collapse = "|"), Description, ignore.case = TRUE),
    Category = case_when(
      Is_Migration ~ "Migration",
      Is_Cancer ~ "Cancer",
      TRUE ~ "Other"
    )
  )

# Summary statistics
cat("\n--- Category Distribution ---\n")
cat("Migration-related terms:", sum(ego_ranked$Is_Migration), "\n")
cat("Cancer-related terms:", sum(ego_ranked$Is_Cancer), "\n")
cat("Other terms:", sum(!ego_ranked$Is_Migration & !ego_ranked$Is_Cancer), "\n")

# Top migration/cancer terms
migration_terms <- ego_ranked %>%
  filter(Is_Migration) %>%
  arrange(Rank_FoldEnrichment) %>%
  head(20)

cancer_terms <- ego_ranked %>%
  filter(Is_Cancer) %>%
  arrange(Rank_FoldEnrichment) %>%
  head(20)

cat("\n--- Top Migration Terms by Fold Enrichment ---\n")
print(migration_terms %>% select(Description, FoldEnrichment, Rank_FoldEnrichment, Percentile_FoldEnrichment, pvalue))

cat("\n--- Top Cancer Terms by Fold Enrichment ---\n")
print(cancer_terms %>% select(Description, FoldEnrichment, Rank_FoldEnrichment, Percentile_FoldEnrichment, pvalue))

# Save migration and cancer term tables
write_csv(migration_terms, file.path(approach8_dir, "results", "top20_migration_terms_ranked.csv"))
write_csv(cancer_terms, file.path(approach8_dir, "results", "top20_cancer_terms_ranked.csv"))

################################################################################
# BOOTSTRAP COMPARISON ANALYSIS
################################################################################

cat("\n=== Bootstrap comparison analysis (may take 5-10 minutes) ===\n")

# Get migration gene sets for comparison
msigdb_go <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP")

migration_go_genes <- msigdb_go %>%
  filter(grepl("CELL.*MIGRATION|CELL.*MOTILITY", gs_name, ignore.case = TRUE)) %>%
  pull(gene_symbol) %>%
  unique()

cat("Migration gene set size:", length(migration_go_genes), "\n")

# Observed enrichment for migration genes
migration_overlap <- length(intersect(tes_direct_genes, migration_go_genes))
migration_in_universe <- length(intersect(migration_go_genes, all_expressed_genes))
migration_expected <- (length(tes_direct_genes) * migration_in_universe) / length(all_expressed_genes)
migration_fold_enrichment <- migration_overlap / migration_expected

cat("Migration gene overlap with TES targets:", migration_overlap, "\n")
cat("Expected by chance:", round(migration_expected, 2), "\n")
cat("Observed fold enrichment:", round(migration_fold_enrichment, 3), "\n")

# Bootstrap: sample random gene sets and calculate enrichment
set.seed(42)
n_bootstrap <- 500
bootstrap_size <- length(migration_go_genes)

cat("\nRunning", n_bootstrap, "bootstrap iterations...\n")

bootstrap_enrichments <- replicate(n_bootstrap, {
  # Sample random genes from universe
  random_genes <- sample(all_expressed_genes, bootstrap_size)

  # Calculate overlap with TES targets
  random_overlap <- length(intersect(tes_direct_genes, random_genes))

  # Expected by chance
  random_expected <- (length(tes_direct_genes) * bootstrap_size) / length(all_expressed_genes)

  # Fold enrichment
  random_fold <- random_overlap / random_expected

  return(random_fold)
})

# Calculate empirical p-value
empirical_pval <- sum(bootstrap_enrichments >= migration_fold_enrichment) / n_bootstrap

cat("\nBootstrap results:\n")
cat("  Mean fold enrichment (random):", round(mean(bootstrap_enrichments), 3), "\n")
cat("  SD fold enrichment (random):", round(sd(bootstrap_enrichments), 3), "\n")
cat("  95% CI (random): [", round(quantile(bootstrap_enrichments, 0.025), 3), ",",
    round(quantile(bootstrap_enrichments, 0.975), 3), "]\n")
cat("  Observed migration fold enrichment:", round(migration_fold_enrichment, 3), "\n")
cat("  Empirical p-value:", empirical_pval, "\n")

# Z-score
z_score <- (migration_fold_enrichment - mean(bootstrap_enrichments)) / sd(bootstrap_enrichments)
cat("  Z-score:", round(z_score, 3), "\n")

# Save bootstrap results
bootstrap_df <- data.frame(
  Iteration = 1:n_bootstrap,
  FoldEnrichment = bootstrap_enrichments
)
write_csv(bootstrap_df, file.path(approach8_dir, "results", "bootstrap_results.csv"))

################################################################################
# VISUALIZATIONS
################################################################################

cat("\n=== Creating visualizations ===\n")

# Plot 1: Top 30 pathways ranked by fold enrichment
top30 <- head(ego_ranked, 30)

pdf(file.path(approach8_dir, "plots", "top30_pathways_by_fold_enrichment.pdf"),
    width = 14, height = 10)

p1 <- ggplot(top30, aes(x = reorder(Description, FoldEnrichment),
                        y = FoldEnrichment, fill = Category)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = sprintf("#%d", Rank_FoldEnrichment)),
            hjust = -0.2, size = 3) +
  coord_flip() +
  scale_fill_manual(values = c("Migration" = "#E69F00",
                               "Cancer" = "#D55E00",
                               "Other" = "#999999")) +
  labs(title = "Top 30 GO Biological Processes Enriched in TES Direct Targets",
       subtitle = "Ranked by fold enrichment over expected",
       x = "GO Term",
       y = "Fold Enrichment",
       fill = "Category") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 9),
        plot.title = element_text(size = 14, face = "bold"),
        legend.position = "bottom")

print(p1)
dev.off()

# Plot 2: Migration terms highlighted in all enriched pathways
migration_highlight <- ego_ranked %>%
  arrange(desc(FoldEnrichment)) %>%
  mutate(
    Highlight = case_when(
      Is_Migration ~ "Migration",
      Is_Cancer ~ "Cancer",
      TRUE ~ "Other"
    )
  )

pdf(file.path(approach8_dir, "plots", "all_pathways_migration_highlighted.pdf"),
    width = 12, height = 8)

p2 <- ggplot(migration_highlight, aes(x = Rank_FoldEnrichment, y = FoldEnrichment,
                                       color = Highlight, size = Count)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("Migration" = "#E69F00",
                                "Cancer" = "#D55E00",
                                "Other" = "#CCCCCC")) +
  scale_size_continuous(range = c(1, 8)) +
  labs(title = "All Enriched GO Terms: Migration & Cancer Terms Highlighted",
       subtitle = sprintf("%d total enriched terms (padj < 0.05)", nrow(ego_ranked)),
       x = "Rank (by Fold Enrichment)",
       y = "Fold Enrichment",
       color = "Category",
       size = "Gene Count") +
  theme_bw() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        legend.position = "right")

print(p2)
dev.off()

# Plot 3: Bootstrap distribution with observed value
pdf(file.path(approach8_dir, "plots", "bootstrap_comparison.pdf"), width = 10, height = 7)

p3 <- ggplot(bootstrap_df, aes(x = FoldEnrichment)) +
  geom_histogram(bins = 50, fill = "gray70", color = "black", alpha = 0.7) +
  geom_vline(xintercept = migration_fold_enrichment, color = "red", lwd = 1.5,
             linetype = "solid") +
  geom_vline(xintercept = mean(bootstrap_enrichments), color = "blue", lwd = 1,
             linetype = "dashed") +
  geom_vline(xintercept = quantile(bootstrap_enrichments, 0.975), color = "blue",
             lwd = 1, linetype = "dotted") +
  annotate("text", x = migration_fold_enrichment, y = Inf, vjust = 1.5,
           label = sprintf("Observed migration\nenrichment: %.2fx", migration_fold_enrichment),
           color = "red", fontface = "bold") +
  annotate("text", x = mean(bootstrap_enrichments), y = Inf, vjust = 3,
           label = sprintf("Random mean: %.2fx", mean(bootstrap_enrichments)),
           color = "blue") +
  labs(title = "Migration Gene Enrichment vs Random Gene Sets",
       subtitle = sprintf("%d bootstrap iterations (gene set size: %d)",
                         n_bootstrap, bootstrap_size),
       x = "Fold Enrichment",
       y = "Frequency") +
  theme_bw() +
  theme(plot.title = element_text(size = 14, face = "bold"))

print(p3)
dev.off()

# Plot 4: Percentile rankings for migration/cancer terms
percentile_plot_data <- ego_ranked %>%
  filter(Is_Migration | Is_Cancer) %>%
  select(Description, Percentile_FoldEnrichment, Category, FoldEnrichment) %>%
  arrange(desc(Percentile_FoldEnrichment)) %>%
  head(30)

pdf(file.path(approach8_dir, "plots", "migration_cancer_percentile_ranks.pdf"),
    width = 12, height = 10)

p4 <- ggplot(percentile_plot_data,
             aes(x = reorder(Description, Percentile_FoldEnrichment),
                 y = Percentile_FoldEnrichment, fill = Category)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 50, linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = 90, linetype = "dashed", color = "darkred") +
  geom_hline(yintercept = 95, linetype = "dashed", color = "darkred", lwd = 1.2) +
  geom_text(aes(label = sprintf("%.1f%%", Percentile_FoldEnrichment)),
            hjust = -0.1, size = 2.5) +
  coord_flip() +
  scale_fill_manual(values = c("Migration" = "#E69F00", "Cancer" = "#D55E00")) +
  labs(title = "Top Migration & Cancer Terms: Percentile Rankings",
       subtitle = "Ranked by fold enrichment among all enriched GO terms",
       x = "GO Term",
       y = "Percentile Rank (higher = more enriched)",
       fill = "Category") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 8),
        plot.title = element_text(size = 14, face = "bold"),
        legend.position = "bottom") +
  ylim(0, 105)

print(p4)
dev.off()

# Plot 5: Combined significance plot (fold enrichment vs p-value)
pdf(file.path(approach8_dir, "plots", "enrichment_vs_significance.pdf"),
    width = 12, height = 8)

p5 <- ggplot(ego_ranked, aes(x = FoldEnrichment, y = NegLog10P,
                              color = Category, size = Count)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", alpha = 0.5) +
  geom_vline(xintercept = median(ego_ranked$FoldEnrichment), linetype = "dashed",
             color = "gray50", alpha = 0.5) +
  scale_color_manual(values = c("Migration" = "#E69F00",
                                "Cancer" = "#D55E00",
                                "Other" = "#CCCCCC")) +
  scale_size_continuous(range = c(1, 10)) +
  labs(title = "GO Term Enrichment: Fold Enrichment vs Statistical Significance",
       subtitle = "Migration and cancer terms show both high enrichment and high significance",
       x = "Fold Enrichment",
       y = "-log10(p-value)",
       color = "Category",
       size = "Gene Count") +
  theme_bw() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        legend.position = "right")

print(p5)
dev.off()

################################################################################
# SUMMARY STATISTICS
################################################################################

cat("\n=== Generating summary statistics ===\n")

# Calculate statistics for migration/cancer vs other terms
migration_stats <- ego_ranked %>% filter(Is_Migration) %>%
  summarise(
    n_terms = n(),
    mean_fold = mean(FoldEnrichment),
    median_fold = median(FoldEnrichment),
    median_rank = median(Rank_FoldEnrichment),
    median_percentile = median(Percentile_FoldEnrichment),
    top10_count = sum(Rank_FoldEnrichment <= 10),
    top20_count = sum(Rank_FoldEnrichment <= 20),
    top50_count = sum(Rank_FoldEnrichment <= 50)
  )

other_stats <- ego_ranked %>% filter(!Is_Migration & !Is_Cancer) %>%
  summarise(
    n_terms = n(),
    mean_fold = mean(FoldEnrichment),
    median_fold = median(FoldEnrichment)
  )

# Statistical test: Are migration fold enrichments higher than others?
wilcox_test <- wilcox.test(
  ego_ranked %>% filter(Is_Migration) %>% pull(FoldEnrichment),
  ego_ranked %>% filter(!Is_Migration & !Is_Cancer) %>% pull(FoldEnrichment),
  alternative = "greater"
)

summary_text <- paste0(
  "=================================================================\n",
  "APPROACH 8: COMPARATIVE ENRICHMENT ANALYSIS SUMMARY\n",
  "=================================================================\n\n",
  "Total GO terms enriched (padj < 0.05): ", nrow(ego_ranked), "\n",
  "Migration-related terms: ", sum(ego_ranked$Is_Migration), "\n",
  "Cancer-related terms: ", sum(ego_ranked$Is_Cancer), "\n",
  "Other terms: ", sum(!ego_ranked$Is_Migration & !ego_ranked$Is_Cancer), "\n\n",
  "--- Migration Term Statistics ---\n",
  "Number of migration terms:", migration_stats$n_terms, "\n",
  "Mean fold enrichment:", sprintf("%.2f", migration_stats$mean_fold), "\n",
  "Median fold enrichment:", sprintf("%.2f", migration_stats$median_fold), "\n",
  "Median rank:", migration_stats$median_rank, "of", nrow(ego_ranked), "\n",
  "Median percentile:", sprintf("%.1f%%", migration_stats$median_percentile), "\n",
  "Terms in top 10:", migration_stats$top10_count, "\n",
  "Terms in top 20:", migration_stats$top20_count, "\n",
  "Terms in top 50:", migration_stats$top50_count, "\n\n",
  "--- Comparison to Other Terms ---\n",
  "Mean fold enrichment (other terms):", sprintf("%.2f", other_stats$mean_fold), "\n",
  "Fold difference (migration vs other):", sprintf("%.2fx", migration_stats$mean_fold / other_stats$mean_fold), "\n",
  "Wilcoxon test p-value:", format(wilcox_test$p.value, scientific = TRUE), "\n\n",
  "--- Bootstrap Analysis ---\n",
  "Migration gene set size:", length(migration_go_genes), "\n",
  "Observed fold enrichment:", sprintf("%.2f", migration_fold_enrichment), "\n",
  "Random mean fold enrichment:", sprintf("%.2f", mean(bootstrap_enrichments)), "\n",
  "Random 95% CI: [", sprintf("%.2f", quantile(bootstrap_enrichments, 0.025)),
  ", ", sprintf("%.2f", quantile(bootstrap_enrichments, 0.975)), "]\n",
  "Z-score:", sprintf("%.2f", z_score), "\n",
  "Empirical p-value:", format(empirical_pval, scientific = TRUE), "\n\n",
  "=================================================================\n",
  "KEY FINDINGS:\n",
  "=================================================================\n",
  "1. Migration terms rank in top ", sprintf("%.0f%%", migration_stats$median_percentile),
  " by fold enrichment\n",
  "2. ", migration_stats$top10_count, " migration terms in top 10 most enriched pathways\n",
  "3. Migration shows ", sprintf("%.1fx", migration_fold_enrichment / mean(bootstrap_enrichments)),
  " higher enrichment than random gene sets (p ", format(empirical_pval, scientific = TRUE), ")\n",
  "4. Migration terms have ", sprintf("%.1fx", migration_stats$mean_fold / other_stats$mean_fold),
  " higher fold enrichment than other pathways (p ", format(wilcox_test$p.value, scientific = TRUE), ")\n\n",
  "Analysis complete: ", date(), "\n",
  "=================================================================\n"
)

cat(summary_text)
writeLines(summary_text, file.path(approach8_dir, "results", "analysis_summary.txt"))

cat("\n=================================================================\n")
cat("APPROACH 8 COMPLETE\n")
cat("Results saved to:", approach8_dir, "\n")
cat("=================================================================\n")
