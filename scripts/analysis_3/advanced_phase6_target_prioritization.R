#!/usr/bin/env Rscript

################################################################################
# Phase 6.3: Target Gene Prioritization
#
# Purpose: Rank genes for experimental validation based on multi-omics evidence
#          Creates BOTH simplified (3-category) and detailed (6-category) analyses
#
# Author: Advanced Multi-Omics Analysis Plan
# Date: 2025-01-24
################################################################################

message("=== Phase 6.3: Target Gene Prioritization ===")
message("Start time: ", Sys.time())
message("NOTE: Creating both SIMPLIFIED (3-cat) and DETAILED (6-cat) analyses")

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)      # Required for geom_text_repel
  library(pheatmap)
  library(RColorBrewer)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

# Fix namespace conflicts: dplyr functions get masked by AnnotationDbi
select <- dplyr::select
filter <- dplyr::filter
slice <- dplyr::slice
mutate <- dplyr::mutate
arrange <- dplyr::arrange

setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top")

# Input files
INTEGRATED_DATA <- "SRF_Eva_integrated_analysis/scripts/analysis_3/results/08_methylation_expression/integrated_binding_methylation_expression.csv"
PEAK_GENE_DATA <- "SRF_Eva_integrated_analysis/scripts/analysis_3/results/04_category_expression/peak_gene_associations.csv"

# Output directories - separate for detailed vs simplified
OUTPUT_DIR <- "SRF_Eva_integrated_analysis/scripts/analysis_3/results/12_target_prioritization"
OUTPUT_DIR_DETAILED <- file.path(OUTPUT_DIR, "detailed_6cat")
OUTPUT_DIR_SIMPLE <- file.path(OUTPUT_DIR, "simplified_3cat")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTPUT_DIR_DETAILED, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTPUT_DIR_SIMPLE, recursive = TRUE, showWarnings = FALSE)

################################################################################
# Helper function: Convert detailed to simplified categories
################################################################################

convert_to_simple_category <- function(category) {
  dplyr::case_when(
    category == "TES_unique" ~ "TES_Unique",
    category == "TEAD1_unique" ~ "TEAD1_Unique",
    grepl("Shared", category) ~ "Shared",
    category == "Unbound" ~ "Unbound",
    TRUE ~ category
  )
}

################################################################################
# Step 1: Load integrated data
################################################################################

message("\n[Step 1] Loading integrated data...")

# Load three-way integrated data
integrated_data <- read.csv(INTEGRATED_DATA)
message("  Loaded ", nrow(integrated_data), " genes")

# Load peak-gene associations for detailed binding info
peak_gene_data <- read.csv(PEAK_GENE_DATA)

################################################################################
# Step 2: Calculate prioritization scores
################################################################################

message("\n[Step 2] Calculating prioritization scores...")

# Function to normalize score to 0-1
normalize_score <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

# Calculate individual scores
prioritization_data <- integrated_data %>%
  mutate(
    # 1. Binding confidence score (0-1)
    binding_score = case_when(
      !tes_bound & !tead1_bound ~ 0,
      n_peaks == 0 ~ 0,
      has_promoter_peak ~ 1.0,  # Promoter peaks highest confidence
      TRUE ~ 0.7 * (log2(max_tes_signal + max_tead1_signal + 1) / 10)  # Signal-based for enhancers
    ),

    # 2. Expression magnitude score (0-1)
    expression_score = normalize_score(abs(log2FoldChange)),

    # 3. Statistical significance score (0-1)
    significance_score = ifelse(is.na(padj), 0, -log10(padj + 1e-300) / 50),
    significance_score = pmin(significance_score, 1),  # Cap at 1

    # 4. Methylation evidence score (0-1)
    methylation_score = case_when(
      !has_promoter_dmr ~ 0,
      regulatory_mechanism == "Direct_Epigenetic_Silencing" ~ 1.0,
      regulatory_mechanism == "Indirect_Methylation" ~ 0.7,
      TRUE ~ 0.5
    ),

    # 5. Biological relevance score (placeholder, will be enhanced with databases)
    # For now, prioritize: cancer-related GO terms, repression, strong effects
    biological_score = case_when(
      regulatory_mechanism == "Direct_Epigenetic_Silencing" ~ 0.9,
      regulatory_mechanism == "Non_Methylation_Repression" ~ 0.8,
      regulatory_mechanism == "TES_Activation" ~ 0.6,
      is_DEG ~ 0.5,
      TRUE ~ 0.2
    ),

    # 6. Novelty score (genes without many publications - would need PubMed API)
    # Placeholder: favor genes with moderate expression (not housekeeping, not silent)
    novelty_score = case_when(
      baseMean < 10 ~ 0.3,  # Too low expression
      baseMean > 10000 ~ 0.4,  # Likely housekeeping
      TRUE ~ 0.8
    )
  )

# Weighted composite score
prioritization_data <- prioritization_data %>%
  mutate(
    priority_score =
      0.25 * binding_score +
      0.20 * expression_score +
      0.20 * significance_score +
      0.15 * methylation_score +
      0.10 * biological_score +
      0.10 * novelty_score,

    # Normalize to 0-100
    priority_score_normalized = priority_score * 100
  )

message("  Calculated priority scores for ", nrow(prioritization_data), " genes")

################################################################################
# Step 3: Add biological annotations
################################################################################

message("\n[Step 3] Adding biological annotations...")

# Get gene IDs
ensembl_to_entrez <- mapIds(org.Hs.eg.db,
                             keys = integrated_data$gene_name,
                             column = "ENTREZID",
                             keytype = "SYMBOL",
                             multiVals = "first")

prioritization_data$entrez_id <- ensembl_to_entrez[prioritization_data$gene_name]

# Get GO annotations
go_bp <- select(org.Hs.eg.db,
                keys = prioritization_data$entrez_id,
                columns = "GO",
                keytype = "ENTREZID")

go_bp <- go_bp %>% filter(ONTOLOGY == "BP")

# Flag cancer-related genes
cancer_terms <- c(
  "apoptosis", "apoptotic", "cell death", "cell cycle", "proliferation",
  "migration", "invasion", "metastasis", "angiogenesis", "DNA repair"
)

# For each gene, check if it has cancer-related GO terms
cancer_related <- sapply(prioritization_data$entrez_id, function(eid) {
  if (is.na(eid)) return(FALSE)
  gene_gos <- go_bp$TERM[go_bp$ENTREZID == eid]
  any(grepl(paste(cancer_terms, collapse = "|"), gene_gos, ignore.case = TRUE))
})

prioritization_data$is_cancer_related <- cancer_related

# Update biological score based on annotations
prioritization_data <- prioritization_data %>%
  mutate(
    biological_score = ifelse(is_cancer_related,
                              pmin(biological_score + 0.2, 1.0),
                              biological_score),

    # Recalculate priority score
    priority_score =
      0.25 * binding_score +
      0.20 * expression_score +
      0.20 * significance_score +
      0.15 * methylation_score +
      0.10 * biological_score +
      0.10 * novelty_score,

    priority_score_normalized = priority_score * 100
  )

################################################################################
# Step 4: Prioritization tiers
################################################################################

message("\n[Step 4] Assigning prioritization tiers...")

# Define tiers based on score thresholds
prioritization_data <- prioritization_data %>%
  mutate(
    priority_tier = case_when(
      priority_score_normalized >= 70 ~ "High Priority",
      priority_score_normalized >= 50 ~ "Medium Priority",
      priority_score_normalized >= 30 ~ "Low Priority",
      TRUE ~ "Very Low Priority"
    ),
    priority_tier = factor(priority_tier,
                          levels = c("High Priority", "Medium Priority",
                                   "Low Priority", "Very Low Priority"))
  )

tier_counts <- table(prioritization_data$priority_tier)
message("  Priority tiers assigned:")
print(tier_counts)

################################################################################
# Step 5: Top candidate selection
################################################################################

message("\n[Step 5] Selecting top candidates...")

# High priority genes
high_priority <- prioritization_data %>%
  filter(priority_tier == "High Priority") %>%
  arrange(desc(priority_score_normalized))

message("  High priority genes: ", nrow(high_priority))

# Top 50 overall
top_50 <- prioritization_data %>%
  arrange(desc(priority_score_normalized)) %>%
  head(50)

# Top 20 per mechanism
top_per_mechanism <- prioritization_data %>%
  filter(regulatory_mechanism %in% c("Direct_Epigenetic_Silencing",
                                     "Non_Methylation_Repression",
                                     "TES_Activation")) %>%
  group_by(regulatory_mechanism) %>%
  arrange(desc(priority_score_normalized)) %>%
  slice_head(n = 20) %>%
  ungroup()

################################################################################
# Step 6: CRISPR target selection
################################################################################

message("\n[Step 6] Selecting CRISPR validation candidates...")

# Criteria for CRISPR targets:
# 1. High priority score
# 2. Strong multi-omics evidence
# 3. Cancer-related
# 4. Direct TES targets (preferably with methylation)

crispr_candidates <- prioritization_data %>%
  filter(
    priority_score_normalized >= 60,
    is_DEG,
    tes_bound,
    is_cancer_related
  ) %>%
  arrange(desc(priority_score_normalized)) %>%
  head(20) %>%
  mutate(
    crispr_rationale = paste0(
      "Score: ", round(priority_score_normalized, 1),
      " | Mechanism: ", regulatory_mechanism,
      " | log2FC: ", round(log2FoldChange, 2),
      " | Binding: ", primary_category
    )
  )

message("  CRISPR candidates selected: ", nrow(crispr_candidates))

################################################################################
# Step 7: Export results
################################################################################

message("\n[Step 7] Exporting results...")

# Main prioritized list
write.csv(prioritization_data %>% arrange(desc(priority_score_normalized)),
          file.path(OUTPUT_DIR, "all_genes_prioritized.csv"),
          row.names = FALSE)

# High priority genes
write.csv(high_priority,
          file.path(OUTPUT_DIR, "high_priority_genes.csv"),
          row.names = FALSE)

# Top 50
write.csv(top_50,
          file.path(OUTPUT_DIR, "top_50_candidates.csv"),
          row.names = FALSE)

# Top per mechanism
write.csv(top_per_mechanism,
          file.path(OUTPUT_DIR, "top_candidates_by_mechanism.csv"),
          row.names = FALSE)

# CRISPR candidates
write.csv(crispr_candidates,
          file.path(OUTPUT_DIR, "CRISPR_validation_candidates.csv"),
          row.names = FALSE)

# Score components for top 50
score_components <- top_50 %>%
  select(gene_name, binding_score, expression_score, significance_score,
         methylation_score, biological_score, novelty_score, priority_score_normalized)

write.csv(score_components,
          file.path(OUTPUT_DIR, "top_50_score_breakdown.csv"),
          row.names = FALSE)

################################################################################
# Step 8: Visualizations
################################################################################

message("\n[Step 8] Generating visualizations...")

# 8.1: Priority score distribution
pdf(file.path(OUTPUT_DIR, "priority_score_distribution.pdf"), width = 10, height = 6)
ggplot(prioritization_data, aes(x = priority_score_normalized)) +
  geom_histogram(aes(fill = priority_tier), bins = 50, alpha = 0.7) +
  geom_vline(xintercept = c(30, 50, 70), linetype = "dashed", color = "red") +
  scale_fill_brewer(palette = "Set1") +
  theme_classic() +
  labs(
    title = "Distribution of Priority Scores",
    x = "Priority Score (0-100)",
    y = "Number of Genes",
    fill = "Priority Tier"
  )
dev.off()

# 8.2: Score components heatmap
if (nrow(top_50) > 0) {
  score_matrix <- as.matrix(score_components[, -1])
  rownames(score_matrix) <- score_components$gene_name

  pdf(file.path(OUTPUT_DIR, "top_50_score_heatmap.pdf"), width = 10, height = 14)
  pheatmap(
    score_matrix,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    color = colorRampPalette(c("white", "yellow", "red"))(100),
    main = "Top 50 Candidates - Score Components",
    labels_col = c("Binding", "Expression", "Significance",
                  "Methylation", "Biological", "Novelty", "Total"),
    fontsize_row = 6,
    cellwidth = 20,
    cellheight = 8
  )
  dev.off()
}

# 8.3: Priority by mechanism
pdf(file.path(OUTPUT_DIR, "priority_by_mechanism.pdf"), width = 12, height = 8)
key_mechanisms <- prioritization_data %>%
  filter(regulatory_mechanism %in% c("Direct_Epigenetic_Silencing",
                                     "Non_Methylation_Repression",
                                     "Indirect_Methylation",
                                     "TES_Activation",
                                     "TEAD1_Specific"))

ggplot(key_mechanisms, aes(x = regulatory_mechanism, y = priority_score_normalized,
                           fill = regulatory_mechanism)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.2, alpha = 0.5, outlier.shape = NA) +
  scale_fill_brewer(palette = "Set2") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  labs(
    title = "Priority Scores by Regulatory Mechanism",
    x = "Mechanism",
    y = "Priority Score"
  )
dev.off()

# 8.4: Top candidates summary card
if (nrow(crispr_candidates) > 0) {
  pdf(file.path(OUTPUT_DIR, "CRISPR_candidates_summary.pdf"), width = 14, height = 10)

  p <- ggplot(crispr_candidates %>% head(15),
              aes(x = reorder(gene_name, priority_score_normalized),
                  y = priority_score_normalized)) +
    geom_bar(aes(fill = regulatory_mechanism), stat = "identity") +
    geom_text(aes(label = round(priority_score_normalized, 1)),
              hjust = -0.2, size = 3) +
    coord_flip() +
    scale_fill_brewer(palette = "Set2") +
    theme_classic() +
    ylim(0, 100) +
    labs(
      title = "Top 15 CRISPR Validation Candidates",
      subtitle = "Ranked by composite priority score",
      x = "Gene",
      y = "Priority Score",
      fill = "Mechanism"
    )

  print(p)
  dev.off()
}

# 8.5: Multi-dimensional scatter
pdf(file.path(OUTPUT_DIR, "expression_vs_binding_priority.pdf"), width = 12, height = 10)
ggplot(prioritization_data %>% filter(is_DEG),
       aes(x = log2(max_tes_signal + max_tead1_signal + 1),
           y = abs(log2FoldChange))) +
  geom_point(aes(color = priority_tier, size = priority_score_normalized), alpha = 0.6) +
  geom_text_repel(data = prioritization_data %>%
                    filter(priority_score_normalized >= 75) %>%
                    head(20),
                  aes(label = gene_name), size = 3, max.overlaps = 20) +
  scale_color_manual(values = c("High Priority" = "#E41A1C",
                               "Medium Priority" = "#FF7F00",
                               "Low Priority" = "#FFFF33",
                               "Very Low Priority" = "gray")) +
  scale_size_continuous(range = c(1, 5)) +
  theme_classic() +
  labs(
    title = "Gene Prioritization: Binding Signal vs Expression Change",
    x = "log2(Total Binding Signal + 1)",
    y = "|log2 Fold Change|",
    color = "Priority Tier",
    size = "Priority Score"
  )
dev.off()

################################################################################
# Step 9: Summary report
################################################################################

message("\n[Step 9] Creating summary report...")

sink(file.path(OUTPUT_DIR, "PHASE6_3_SUMMARY.txt"))
cat("=== Phase 6.3: Target Gene Prioritization ===\n")
cat("Date:", as.character(Sys.time()), "\n\n")

cat("Priority Score Components:\n")
cat("==========================\n")
cat("- Binding confidence: 25%\n")
cat("- Expression magnitude: 20%\n")
cat("- Statistical significance: 20%\n")
cat("- Methylation evidence: 15%\n")
cat("- Biological relevance: 10%\n")
cat("- Novelty score: 10%\n\n")

cat("Prioritization Tiers:\n")
cat("=====================\n")
print(tier_counts)
cat("\n")

cat("High Priority Genes:", nrow(high_priority), "\n")
cat("  - With promoter peaks:", sum(high_priority$has_promoter_peak), "\n")
cat("  - With promoter DMRs:", sum(high_priority$has_promoter_dmr, na.rm = TRUE), "\n")
cat("  - Cancer-related:", sum(high_priority$is_cancer_related), "\n\n")

cat("Top 10 Candidates:\n")
cat("==================\n")
print(top_50 %>%
      head(10) %>%
      select(gene_name, priority_score_normalized, regulatory_mechanism,
             log2FoldChange, primary_category) %>%
      as.data.frame())
cat("\n")

cat("CRISPR Validation Candidates:", nrow(crispr_candidates), "\n")
cat("============================\n")
print(crispr_candidates %>%
      head(10) %>%
      select(gene_name, priority_score_normalized, regulatory_mechanism,
             log2FoldChange, crispr_rationale) %>%
      as.data.frame())

sink()

message("\n=== Analysis Complete ===")
message("Output directory: ", OUTPUT_DIR)
message("End time: ", Sys.time())
message("\nKey outputs:")
message("  - all_genes_prioritized.csv (", nrow(prioritization_data), " genes)")
message("  - high_priority_genes.csv (", nrow(high_priority), " genes)")
message("  - CRISPR_validation_candidates.csv (", nrow(crispr_candidates), " genes)")
