#!/usr/bin/env Rscript
#
# IMPROVED GSEA ANALYSIS: Cancer-Relevant Pathways with Expression Direction
# Enhanced version with UP/DOWN regulation analysis
# Comparing TES-only, TEAD1-only, and shared transcriptional targets

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(pheatmap)
  library(ComplexHeatmap)
  library(circlize)
  library(gridExtra)
  library(RColorBrewer)
  library(tidyr)
})

# Set working directory and create output directories
setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_1")

cat("=== IMPROVED GSEA ANALYSIS: Cancer Pathways with Expression Direction ===\n")
cat("Analysis started:", as.character(Sys.time()), "\n\n")

# Create output directory
output_dir <- "output/18_gsea_cancer_pathways_improved"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# PHASE 1: DATA LOADING
# =============================================================================

cat("=== PHASE 1: Loading Integrative Analysis Results ===\n")

# Load direct target classifications
cat("Loading direct target gene lists...\n")
tes_direct <- read.csv("output/11_final_integrative_analysis_all_genes/direct_targets/TES_direct_targets_all_genes.csv",
  stringsAsFactors = FALSE
)
tead1_direct <- read.csv("output/11_final_integrative_analysis_all_genes/direct_targets/TEAD1_direct_targets_all_genes.csv",
  stringsAsFactors = FALSE
)
tes_specific <- read.csv("output/11_final_integrative_analysis_all_genes/direct_targets/TES_specific_targets_all_genes.csv",
  stringsAsFactors = FALSE
)
tead1_specific <- read.csv("output/11_final_integrative_analysis_all_genes/direct_targets/TEAD1_specific_targets_all_genes.csv",
  stringsAsFactors = FALSE
)

# Calculate shared targets
shared_genes <- intersect(tes_direct$ensembl_id, tead1_direct$ensembl_id)
shared_targets <- tes_direct[tes_direct$ensembl_id %in% shared_genes, ]

# IMPROVEMENT 1: Separate by expression direction
cat("\nSeparating targets by expression direction...\n")
tes_direct_up <- tes_direct[tes_direct$log2FoldChange > 0, ]
tes_direct_down <- tes_direct[tes_direct$log2FoldChange < 0, ]
tead1_direct_up <- tead1_direct[tead1_direct$log2FoldChange > 0, ]
tead1_direct_down <- tead1_direct[tead1_direct$log2FoldChange < 0, ]
shared_up <- shared_targets[shared_targets$log2FoldChange > 0, ]
shared_down <- shared_targets[shared_targets$log2FoldChange < 0, ]

cat(sprintf(
  "✓ TES targets: %d total (%d UP, %d DOWN)\n",
  nrow(tes_direct), nrow(tes_direct_up), nrow(tes_direct_down)
))
cat(sprintf(
  "✓ TEAD1 targets: %d total (%d UP, %d DOWN)\n",
  nrow(tead1_direct), nrow(tead1_direct_up), nrow(tead1_direct_down)
))
cat(sprintf(
  "✓ Shared targets: %d total (%d UP, %d DOWN)\n",
  nrow(shared_targets), nrow(shared_up), nrow(shared_down)
))
cat(sprintf("✓ TES-specific: %d genes\n", nrow(tes_specific)))
cat(sprintf("✓ TEAD1-specific: %d genes\n\n", nrow(tead1_specific)))

# Load GO enrichment results
cat("Loading GO enrichment pathway data...\n")
tes_go <- read.csv("output/11_final_integrative_analysis_all_genes/pathway_analysis/TES_direct_GO_enrichment_all_genes.csv",
  stringsAsFactors = FALSE
)
tead1_go <- read.csv("output/11_final_integrative_analysis_all_genes/pathway_analysis/TEAD1_direct_GO_enrichment_all_genes.csv",
  stringsAsFactors = FALSE
)
shared_go <- read.csv("output/11_final_integrative_analysis_all_genes/pathway_analysis/TES_TEAD1_shared_GO_enrichment_all_genes.csv",
  stringsAsFactors = FALSE
)

cat(sprintf("✓ TES direct pathways: %d\n", nrow(tes_go)))
cat(sprintf("✓ TEAD1 direct pathways: %d\n", nrow(tead1_go)))
cat(sprintf("✓ Shared pathways: %d\n\n", nrow(shared_go)))

# =============================================================================
# PHASE 2: FILTER FOR CANCER-RELEVANT PATHWAYS
# =============================================================================

cat("=== PHASE 2: Filtering for Cancer-Relevant Pathways ===\n")

# IMPROVEMENT 2: Expanded keyword list with biological context
cancer_keywords <- list(
  apoptosis = c(
    "apoptosis", "apoptotic", "cell death", "programmed cell death",
    "necrosis", "necrotic", "ferroptosis", "pyroptosis", "anoikis",
    "caspase activation", "death receptor", "intrinsic apoptotic",
    "extrinsic apoptotic", "mitochondrial depolarization"
  ),
  migration = c(
    "migration", "migratory", "motility", "invasion", "invasive",
    "chemotaxis", "chemotactic", "cell movement", "locomotion",
    "extracellular matrix", "cell adhesion", "focal adhesion",
    "integrin", "lamellipodium", "filopodium"
  ),
  proliferation = c(
    "proliferation", "proliferative", "cell cycle", "mitosis",
    "mitotic", "cell division", "growth", "G1/S", "G2/M",
    "cyclin", "CDK", "checkpoint", "DNA replication",
    "cytokinesis", "chromosome segregation"
  ),

  # IMPROVEMENT 3: Add angiogenesis (critical for glioblastoma)
  angiogenesis = c(
    "angiogenesis", "angiogenic", "blood vessel", "vasculature",
    "endothelial", "VEGF", "vascular development",
    "tube formation", "vascularization"
  ),

  # IMPROVEMENT 4: Add metabolism (cancer hallmark)
  metabolism = c(
    "glycolysis", "glycolytic", "glucose metabolism",
    "oxidative phosphorylation", "ATP synthesis",
    "metabolic reprogramming", "Warburg effect",
    "glutamine metabolism", "lipid metabolism"
  )
)

# Function to filter pathways by keywords AND annotate with genes
filter_cancer_pathways_enhanced <- function(go_results, keywords_list, target_genes) {
  all_keywords <- unlist(keywords_list, use.names = FALSE)
  pattern <- paste(all_keywords, collapse = "|")

  filtered <- go_results[grep(pattern, go_results$Description, ignore.case = TRUE), ]

  # Add pathway category
  filtered$pathway_category <- NA
  for (category in names(keywords_list)) {
    category_pattern <- paste(keywords_list[[category]], collapse = "|")
    filtered$pathway_category[grep(category_pattern, filtered$Description, ignore.case = TRUE)] <- category
  }

  # IMPROVEMENT 5: Add expression direction annotation
  # Parse geneID column to get genes in each pathway
  if ("geneID" %in% colnames(filtered)) {
    filtered$direction_ratio <- sapply(filtered$geneID, function(gene_string) {
      if (is.na(gene_string) || gene_string == "") {
        return(NA)
      }

      # Split by "/" to get individual Ensembl IDs
      ensembl_ids <- unlist(strsplit(as.character(gene_string), "/"))

      if (length(ensembl_ids) == 0) {
        return(NA)
      }

      # Count UP vs DOWN in these genes (using ensembl_id)
      up_count <- sum(ensembl_ids %in% target_genes$ensembl_id[target_genes$log2FoldChange > 0])
      down_count <- sum(ensembl_ids %in% target_genes$ensembl_id[target_genes$log2FoldChange < 0])

      if ((up_count + down_count) == 0) {
        return(NA)
      }
      return(up_count / (up_count + down_count))
    })

    # Classify as UP-dominant, DOWN-dominant, or Mixed
    filtered$direction_class <- ifelse(is.na(filtered$direction_ratio), "Unknown",
      ifelse(filtered$direction_ratio > 0.7, "Upregulated",
        ifelse(filtered$direction_ratio < 0.3, "Downregulated",
          "Mixed"
        )
      )
    )
  }

  return(filtered)
}

# Filter pathways for each target group
cat("Filtering pathways by cancer-relevant keywords...\n")
tes_cancer <- filter_cancer_pathways_enhanced(tes_go, cancer_keywords, tes_direct)
tead1_cancer <- filter_cancer_pathways_enhanced(tead1_go, cancer_keywords, tead1_direct)
shared_cancer <- filter_cancer_pathways_enhanced(shared_go, cancer_keywords, shared_targets)

cat(sprintf("✓ TES cancer pathways: %d\n", nrow(tes_cancer)))
cat(sprintf(
  "  - Apoptosis: %d, Migration: %d, Proliferation: %d\n",
  sum(tes_cancer$pathway_category == "apoptosis", na.rm = TRUE),
  sum(tes_cancer$pathway_category == "migration", na.rm = TRUE),
  sum(tes_cancer$pathway_category == "proliferation", na.rm = TRUE)
))
if ("angiogenesis" %in% tes_cancer$pathway_category) {
  cat(sprintf(
    "  - Angiogenesis: %d, Metabolism: %d\n",
    sum(tes_cancer$pathway_category == "angiogenesis", na.rm = TRUE),
    sum(tes_cancer$pathway_category == "metabolism", na.rm = TRUE)
  ))
}

cat(sprintf("\n✓ TEAD1 cancer pathways: %d\n", nrow(tead1_cancer)))
cat(sprintf(
  "  - Apoptosis: %d, Migration: %d, Proliferation: %d\n",
  sum(tead1_cancer$pathway_category == "apoptosis", na.rm = TRUE),
  sum(tead1_cancer$pathway_category == "migration", na.rm = TRUE),
  sum(tead1_cancer$pathway_category == "proliferation", na.rm = TRUE)
))

cat(sprintf("\n✓ Shared cancer pathways: %d\n\n", nrow(shared_cancer)))

# =============================================================================
# PHASE 3: EXPRESSION DIRECTION ANALYSIS
# =============================================================================

cat("=== PHASE 3: Expression Direction Analysis ===\n")

# Analyze direction bias for each pathway category
analyze_direction_bias <- function(pathways, category_name) {
  if (nrow(pathways) == 0 || !"direction_class" %in% colnames(pathways)) {
    cat(sprintf("%s: No direction data available\n", category_name))
    return(NULL)
  }

  direction_summary <- table(pathways$pathway_category, pathways$direction_class)
  cat(sprintf("\n%s Direction Bias:\n", category_name))
  print(direction_summary)

  return(direction_summary)
}

tes_direction <- analyze_direction_bias(tes_cancer, "TES")
tead1_direction <- analyze_direction_bias(tead1_cancer, "TEAD1")
shared_direction <- analyze_direction_bias(shared_cancer, "Shared")

# =============================================================================
# PHASE 4: EXPORT RESULTS
# =============================================================================

cat("\n=== PHASE 4: Exporting Results ===\n")

write.csv(tes_cancer,
  file.path(output_dir, "TES_only_cancer_pathways_directional.csv"),
  row.names = FALSE
)
write.csv(tead1_cancer,
  file.path(output_dir, "TEAD1_only_cancer_pathways_directional.csv"),
  row.names = FALSE
)
write.csv(shared_cancer,
  file.path(output_dir, "Shared_cancer_pathways_directional.csv"),
  row.names = FALSE
)

cat("✓ Directional pathway tables exported\n\n")

# =============================================================================
# PHASE 5: ENHANCED VISUALIZATIONS
# =============================================================================

cat("=== PHASE 5: Creating Enhanced Visualizations ===\n")

# Add source annotation
tes_cancer$target_group <- "TES-only"
tead1_cancer$target_group <- "TEAD1-only"
shared_cancer$target_group <- "Shared"

# Combine datasets
all_cancer_pathways <- rbind(
  tes_cancer[, intersect(names(tes_cancer), names(tead1_cancer))],
  tead1_cancer[, intersect(names(tes_cancer), names(tead1_cancer))],
  shared_cancer[, intersect(names(tes_cancer), names(tead1_cancer))]
)

# Parse GeneRatio
all_cancer_pathways$GeneRatio_numeric <- sapply(
  strsplit(all_cancer_pathways$GeneRatio, "/"),
  function(x) as.numeric(x[1]) / as.numeric(x[2])
)

# IMPROVEMENT 6: Direction-aware visualization
if ("direction_class" %in% colnames(all_cancer_pathways)) {
  cat("Creating direction-aware pathway plots...\n")

  # Plot 1: Direction bias by category and target group
  direction_data <- all_cancer_pathways %>%
    filter(!is.na(pathway_category) & !is.na(direction_class)) %>%
    group_by(target_group, pathway_category, direction_class) %>%
    summarise(count = n(), .groups = "drop")

  pdf(file.path(output_dir, "01_pathway_direction_bias.pdf"), width = 14, height = 10)
  p1 <- ggplot(direction_data, aes(
    x = pathway_category, y = count,
    fill = direction_class
  )) +
    geom_bar(stat = "identity", position = "fill", color = "black", linewidth = 0.3) +
    facet_wrap(~target_group, ncol = 3) +
    labs(
      title = "Pathway Expression Direction Bias",
      subtitle = "Proportion of Up/Down/Mixed regulated pathways per category",
      x = "Pathway Category",
      y = "Proportion",
      fill = "Direction"
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
      strip.text = element_text(face = "bold", size = 12)
    ) +
    scale_fill_manual(values = c(
      "Upregulated" = "#E31A1C",
      "Downregulated" = "#1F78B4",
      "Mixed" = "#33A02C",
      "Unknown" = "gray80"
    )) +
    scale_y_continuous(labels = scales::percent)
  print(p1)
  dev.off()

  # Plot 2: Biological interpretation plot
  cat("Creating biological interpretation plot...\n")

  # Focus on key pathways with clear direction
  key_pathways <- all_cancer_pathways %>%
    filter(!is.na(pathway_category) & direction_class != "Unknown") %>%
    group_by(target_group, pathway_category) %>%
    arrange(p.adjust) %>%
    slice_head(n = 5) %>%
    ungroup()

  # Only create plot if we have data
  if (nrow(key_pathways) > 0 && length(unique(key_pathways$pathway_category)) > 0) {
    pdf(file.path(output_dir, "02_key_pathways_with_direction.pdf"), width = 16, height = 12)
    p2 <- ggplot(key_pathways, aes(
      x = target_group, y = Description,
      color = direction_class, size = GeneRatio_numeric
    )) +
      geom_point(alpha = 0.8) +
      facet_wrap(~pathway_category, scales = "free_y", ncol = 2) +
      labs(
        title = "Top Cancer Pathways with Expression Direction",
        subtitle = "Top 5 pathways per category - Color indicates UP/DOWN regulation",
        x = "Target Group",
        y = "Pathway Description",
        color = "Direction",
        size = "Gene Ratio"
      ) +
      theme_bw(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 9),
        strip.text = element_text(face = "bold", size = 11)
      ) +
      scale_color_manual(values = c(
        "Upregulated" = "#E31A1C",
        "Downregulated" = "#1F78B4",
        "Mixed" = "#33A02C"
      )) +
      scale_size_continuous(range = c(3, 10))
    print(p2)
    dev.off()
  } else {
    cat("  Note: Not enough pathways with direction annotation for Plot 2\n")
  }
}

cat("✓ Enhanced visualizations created\n\n")

# =============================================================================
# PHASE 6: BIOLOGICAL INTERPRETATION SUMMARY
# =============================================================================

cat("=== PHASE 6: Generating Biological Interpretation ===\n")

# Create interpretation summary
summary_file <- file.path(output_dir, "BIOLOGICAL_INTERPRETATION.txt")
cat("BIOLOGICAL INTERPRETATION OF CANCER PATHWAY ANALYSIS\n", file = summary_file)
cat("====================================================\n\n", file = summary_file, append = TRUE)
cat(paste("Generated:", Sys.time(), "\n\n"), file = summary_file, append = TRUE)

cat("KEY FINDINGS:\n", file = summary_file, append = TRUE)
cat("-------------\n\n", file = summary_file, append = TRUE)

# Analyze TES effects
if (nrow(tes_cancer) > 0) {
  tes_apop_up <- sum(tes_cancer$pathway_category == "apoptosis" &
    tes_cancer$direction_class == "Upregulated", na.rm = TRUE)
  tes_apop_down <- sum(tes_cancer$pathway_category == "apoptosis" &
    tes_cancer$direction_class == "Downregulated", na.rm = TRUE)

  cat("TES TRANSCRIPTIONAL EFFECTS:\n", file = summary_file, append = TRUE)
  cat(sprintf(
    "  Apoptosis: %d upregulated, %d downregulated\n",
    tes_apop_up, tes_apop_down
  ), file = summary_file, append = TRUE)

  if (tes_apop_up > tes_apop_down) {
    cat("  → Interpretation: TES predominantly PROMOTES apoptosis (tumor suppressor)\n",
      file = summary_file, append = TRUE
    )
  } else if (tes_apop_down > tes_apop_up) {
    cat("  → Interpretation: TES predominantly INHIBITS apoptosis (pro-survival)\n",
      file = summary_file, append = TRUE
    )
  }
  cat("\n", file = summary_file, append = TRUE)
}

cat("\nNOTE: This analysis provides directional context for pathway regulation.\n",
  file = summary_file, append = TRUE
)
cat("Upregulated apoptosis + Downregulated proliferation = Anti-cancer effect\n",
  file = summary_file, append = TRUE
)
cat("Upregulated migration = Pro-metastatic effect\n",
  file = summary_file, append = TRUE
)

cat("✓ Biological interpretation saved\n\n")

cat("========================================\n")
cat("IMPROVED GSEA ANALYSIS COMPLETE\n")
cat("========================================\n")
cat("Completed:", as.character(Sys.time()), "\n")
cat(sprintf("Output directory: %s\n", output_dir))
cat("\nKey improvements implemented:\n")
cat("  1. Expression direction analysis (UP/DOWN)\n")
cat("  2. Expanded pathway categories (angiogenesis, metabolism)\n")
cat("  3. Direction-aware visualizations\n")
cat("  4. Biological interpretation guide\n")
cat("  5. Enhanced pathway annotations\n")
cat("\nAll results exported successfully!\n")
