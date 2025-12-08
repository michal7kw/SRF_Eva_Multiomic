#!/usr/bin/env Rscript

cat("=================================================================\n")
cat("TIER 1: Show Progression of Specificity\n")
cat("=================================================================\n")

# Load workspace and helper functions
source("SRF_Eva_integrated_analysis/focused_enrichment_analysis/scripts/00_helper_functions.R")
load_libraries()
setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top")
load("SRF_Eva_integrated_analysis/focused_enrichment_analysis/data/loaded_data.RData")

# Output directory
base_dir <- "SRF_Eva_integrated_analysis/focused_enrichment_analysis"
tier1_dir <- file.path(base_dir, "tier1_progression")

# Create directories
dir.create(file.path(tier1_dir, "results"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(tier1_dir, "plots"), recursive = TRUE, showWarnings = FALSE)

################################################################################
# LOAD APPROACH RESULTS
################################################################################

cat("\n=== Loading approach results ===\n")

# Helper function to safely load migration results
load_migration_results <- function(approach_dir, approach_name) {
  file_path <- file.path(base_dir, approach_dir, "results",
                         paste0(approach_name, "_migration_terms.csv"))
  if (file.exists(file_path)) {
    return(read_csv(file_path, show_col_types = FALSE))
  }
  return(NULL)
}

# Helper function to safely load GO enrichment results
load_go_results <- function(approach_dir, prefix) {
  file_path <- file.path(base_dir, approach_dir, "results",
                         paste0(prefix, "_GO_enrichment.csv"))
  if (file.exists(file_path)) {
    return(read_csv(file_path, show_col_types = FALSE))
  }
  return(NULL)
}

# Load results from each approach
results_approach1 <- load_go_results("approach1_direct_targets", "approach1_TES_direct")
results_approach2 <- load_go_results("approach2_downregulated", "approach2_TES_direct_down")
results_approach3 <- load_go_results("approach3_promoter_peaks", "approach3_TES_promoter_DEGs")
results_approach4 <- load_go_results("approach4_high_confidence", "approach4_TES_highconf_DEGs")
results_approach5 <- load_go_results("approach5_diffbind", "approach5_diffbind")

# Load migration terms
migration_approach1 <- load_migration_results("approach1_direct_targets", "approach1")
migration_approach2 <- load_migration_results("approach2_downregulated", "approach2")
migration_approach3 <- load_migration_results("approach3_promoter_peaks", "approach3")
migration_approach4 <- load_migration_results("approach4_high_confidence", "approach4")
migration_approach5 <- load_migration_results("approach5_diffbind", "approach5")

# Load gene lists
tes_direct_genes <- tes_direct$gene_symbol
tes_direct_down_genes <- tes_direct %>%
  filter(log2FoldChange < 0, padj < 0.05) %>%
  pull(gene_symbol)

# Load promoter DEGs
promoter_file <- file.path(base_dir, "approach3_promoter_peaks", "gene_lists", "TES_promoter_DEGs.txt")
if (file.exists(promoter_file)) {
  tes_promoter_degs <- read.table(promoter_file, stringsAsFactors = FALSE)$V1
} else {
  tes_promoter_degs <- character(0)
}

# Load highconf DEGs
highconf_file <- file.path(base_dir, "approach4_high_confidence", "gene_lists", "TES_highconf_DEGs.txt")
if (file.exists(highconf_file)) {
  tes_highconf_degs <- read.table(highconf_file, stringsAsFactors = FALSE)$V1
} else {
  tes_highconf_degs <- character(0)
}

# Load diffbind down
diffbind_file <- file.path(base_dir, "approach5_diffbind", "gene_lists", "TES_diffbind_down.txt")
if (file.exists(diffbind_file)) {
  tes_diffbind_down <- read.table(diffbind_file, stringsAsFactors = FALSE)$V1
} else {
  tes_diffbind_down <- character(0)
}

################################################################################
# COMPILE PROGRESSION DATA
################################################################################

cat("\n=== Compiling progression data ===\n")

progression_data <- data.frame(
  Approach = c("All TES peaks",
               "Direct targets (any)",
               "Direct downregulated",
               "Promoter DEGs",
               "High-conf DEGs",
               "DiffBind + Down"),
  N_genes = c(
    length(unique(tes_peaks_annotated$SYMBOL[!is.na(tes_peaks_annotated$SYMBOL)])),
    length(tes_direct_genes),
    length(tes_direct_down_genes),
    length(tes_promoter_degs),
    length(tes_highconf_degs),
    length(tes_diffbind_down)
  ),
  N_sig_terms = c(
    NA,  # Would need to load original analysis
    ifelse(!is.null(results_approach1), sum(results_approach1$qvalue < 0.05), 0),
    ifelse(!is.null(results_approach2), sum(results_approach2$qvalue < 0.05), 0),
    ifelse(!is.null(results_approach3), sum(results_approach3$qvalue < 0.05), 0),
    ifelse(!is.null(results_approach4), sum(results_approach4$qvalue < 0.05), 0),
    ifelse(!is.null(results_approach5), sum(results_approach5$qvalue < 0.05), NA)
  ),
  Top_migration_rank = c(
    NA,
    ifelse(!is.null(migration_approach1) && nrow(migration_approach1) > 0,
           which(results_approach1$ID == migration_approach1$ID[1])[1], NA),
    ifelse(!is.null(migration_approach2) && nrow(migration_approach2) > 0,
           which(results_approach2$ID == migration_approach2$ID[1])[1], NA),
    ifelse(!is.null(migration_approach3) && nrow(migration_approach3) > 0,
           which(results_approach3$ID == migration_approach3$ID[1])[1], NA),
    ifelse(!is.null(migration_approach4) && nrow(migration_approach4) > 0,
           which(results_approach4$ID == migration_approach4$ID[1])[1], NA),
    ifelse(!is.null(migration_approach5) && nrow(migration_approach5) > 0,
           which(results_approach5$ID == migration_approach5$ID[1])[1], NA)
  )
)

write_csv(progression_data,
          file.path(tier1_dir, "results", "progression_summary.csv"))

cat("\n=== PROGRESSION SUMMARY ===\n")
print(progression_data)

################################################################################
# CREATE VISUALIZATIONS
################################################################################

cat("\n=== Creating visualizations ===\n")

# Plot 1: Migration term ranking progression
if (nrow(progression_data %>% filter(!is.na(Top_migration_rank))) > 0) {
  p_progression <- ggplot(progression_data %>% filter(!is.na(Top_migration_rank)),
                          aes(x = reorder(Approach, -Top_migration_rank),
                              y = Top_migration_rank,
                              fill = N_genes)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = N_genes), hjust = -0.2) +
    coord_flip() +
    scale_fill_gradient(low = "lightblue", high = "darkblue") +
    labs(title = "Progression of Specificity: Migration Term Ranking",
         subtitle = "Lower rank = stronger enrichment",
         x = "Approach",
         y = "Rank of Top Migration Term",
         fill = "N genes") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5))

  ggsave(file.path(tier1_dir, "plots", "progression_comparison.pdf"),
         p_progression, width = 10, height = 6)
}

# Plot 2: Number of significant terms
if (nrow(progression_data %>% filter(!is.na(N_sig_terms))) > 0) {
  p_sig_terms <- ggplot(progression_data %>% filter(!is.na(N_sig_terms)),
                        aes(x = reorder(Approach, -N_sig_terms),
                            y = N_sig_terms,
                            fill = Approach)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = N_sig_terms), vjust = -0.5) +
    labs(title = "Number of Significant GO Terms by Approach",
         subtitle = "q-value < 0.05",
         x = "Approach",
         y = "Number of Significant GO Terms") +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(size = 14, face = "bold"))

  ggsave(file.path(tier1_dir, "plots", "sig_terms_comparison.pdf"),
         p_sig_terms, width = 10, height = 6)
}

# Plot 3: Gene count progression
p_gene_counts <- ggplot(progression_data,
                        aes(x = reorder(Approach, -N_genes),
                            y = N_genes,
                            fill = Approach)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = N_genes), vjust = -0.5) +
  labs(title = "Gene Count by Approach",
       subtitle = "Progression from broad to specific",
       x = "Approach",
       y = "Number of Genes") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 14, face = "bold"))

ggsave(file.path(tier1_dir, "plots", "gene_counts_progression.pdf"),
       p_gene_counts, width = 10, height = 6)

cat("Visualizations created successfully!\n")

cat("\n=================================================================\n")
cat("TIER 1 COMPLETE\n")
cat("=================================================================\n")
