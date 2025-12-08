#!/usr/bin/env Rscript
#
# TRUE GSEA ANALYSIS: Rank-based Gene Set Enrichment Analysis
# Uses fgsea package for proper ranking-based enrichment
# Analyzes ALL genes ranked by fold change (not just significant ones)

# %%
suppressPackageStartupMessages({
  library(fgsea)
  library(dplyr)
  library(ggplot2)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(readr)
  library(stringr)
  library(tidyr)
})

setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_1")

cat("=== TRUE GSEA ANALYSIS: Rank-Based Enrichment ===\n")
cat("Analysis started:", as.character(Sys.time()), "\n\n")

# Create output directory
output_dir <- "output/01_true_gsea_analysis"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# PHASE 1: LOAD RNA-SEQ DATA AND CREATE RANKED GENE LIST
# =============================================================================

cat("=== PHASE 1: Creating Ranked Gene List ===\n")

# %%
# Load complete RNA-seq results (ALL genes, not just significant)
rna_all <- read.delim("../../../SRF_Eva_RNA/results/05_deseq2/deseq2_results_TES_vs_GFP.txt",
  stringsAsFactors = FALSE
)

cat(sprintf("✓ Loaded %d genes from RNA-seq\n", nrow(rna_all)))

# Clean and prepare for ranking
rna_all$ensembl_id <- gsub("\\..*", "", rna_all$gene_id)

head(rna_all)

# Get gene symbols
rna_all$gene_symbol <- mapIds(org.Hs.eg.db,
  keys = rna_all$ensembl_id,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

head(rna_all)

# %%
# Remove genes without symbols or log2FC
rna_ranked <- rna_all %>%
  filter(!is.na(gene_symbol) & !is.na(log2FoldChange)) %>%
  arrange(desc(log2FoldChange))

cat(sprintf("✓ Ranked %d genes with valid symbols and fold changes\n", nrow(rna_ranked)))

# Create ranked list (gene symbol → log2FC)
gene_ranks <- setNames(rna_ranked$log2FoldChange, rna_ranked$gene_symbol)

# Remove duplicates (keep first = highest ranking)
gene_ranks <- gene_ranks[!duplicated(names(gene_ranks))]

cat(sprintf("✓ Final ranked list: %d unique genes\n", length(gene_ranks)))
cat(sprintf("  Range: %.2f to %.2f log2FC\n", min(gene_ranks), max(gene_ranks)))
cat(sprintf("  Most upregulated: %s (%.2f)\n", names(gene_ranks)[1], gene_ranks[1]))
cat(sprintf(
  "  Most downregulated: %s (%.2f)\n\n", names(gene_ranks)[length(gene_ranks)],
  gene_ranks[length(gene_ranks)]
))

# %%

# =============================================================================
# PHASE 2: PREPARE GENE SETS FOR GSEA
# =============================================================================


cat("=== PHASE 2: Preparing Gene Sets ===\n")

# Load GO gene sets from org.Hs.eg.db
cat("Loading GO Biological Process gene sets...\n")

go_bp <- AnnotationDbi::select(org.Hs.eg.db,
  keys = keys(org.Hs.eg.db, keytype = "GOALL"),
  columns = c("SYMBOL", "ONTOLOGYALL"),
  keytype = "GOALL"
) %>%
  filter(ONTOLOGYALL == "BP") %>%
  filter(!is.na(SYMBOL))

# Convert to list format for fgsea
go_gene_sets <- split(go_bp$SYMBOL, go_bp$GOALL)

# Remove very small and very large gene sets
go_gene_sets <- go_gene_sets[sapply(go_gene_sets, length) >= 10 &
  sapply(go_gene_sets, length) <= 500]

cat(sprintf("✓ Loaded %d GO BP gene sets (10-500 genes each)\n", length(go_gene_sets)))

# Create cancer-focused gene sets manually
cancer_gene_sets <- list(
  APOPTOSIS = c(
    "BAX", "BAK1", "BID", "BIM", "PUMA", "NOXA", "CASP3", "CASP8", "CASP9",
    "FAS", "FASL", "TNFRSF10A", "TNFRSF10B", "TP53", "BCL2", "BCL2L1",
    "MCL1", "APAF1", "CYCS", "DIABLO"
  ),
  ANTI_APOPTOSIS = c(
    "BCL2", "BCL2L1", "BCL2L2", "MCL1", "BCL2A1", "BCLW",
    "BIRC2", "BIRC3", "BIRC5", "XIAP", "NAIP", "CFLIP"
  ),
  CELL_CYCLE = c(
    "CDK1", "CDK2", "CDK4", "CDK6", "CCNA1", "CCNA2", "CCNB1", "CCNB2",
    "CCND1", "CCND2", "CCND3", "CCNE1", "CCNE2", "E2F1", "E2F2", "E2F3",
    "RB1", "CDKN1A", "CDKN1B", "CDKN2A", "CDKN2B"
  ),
  MIGRATION = c(
    "ITGB1", "ITGB3", "ITGB5", "ITGA5", "ITGAV", "CDH1", "CDH2",
    "VIM", "SNAI1", "SNAI2", "TWIST1", "ZEB1", "ZEB2", "MMP2", "MMP9",
    "TIMP1", "TIMP2", "CXCR4", "CXCL12"
  ),
  HIPPO_YAP = c(
    "YAP1", "WWTR1", "TEAD1", "TEAD2", "TEAD3", "TEAD4", "LATS1", "LATS2",
    "STK3", "STK4", "MOB1A", "MOB1B", "SAV1", "NF2", "AMOT", "AMOTL1",
    "AMOTL2", "PTPN14", "FAT1", "FAT2", "FAT3", "FAT4"
  ),
  ANGIOGENESIS = c(
    "VEGFA", "VEGFB", "VEGFC", "VEGFD", "FLT1", "KDR", "FLT4",
    "ANGPT1", "ANGPT2", "TEK", "PDGFA", "PDGFB", "PDGFRA", "PDGFRB",
    "FGF2", "FGFR1", "FGFR2", "HIF1A", "EPAS1"
  ),
  EMT = c(
    "CDH1", "CDH2", "VIM", "FN1", "SNAI1", "SNAI2", "SLUG", "TWIST1", "TWIST2",
    "ZEB1", "ZEB2", "GSC", "FOXC2", "TCF3", "TCF4"
  ),
  GLIOBLASTOMA_CORE = c(
    "EGFR", "PTEN", "TP53", "CDKN2A", "CDKN2B", "NF1", "RB1",
    "PIK3CA", "PIK3R1", "PDGFRA", "MET", "BRAF", "IDH1", "IDH2",
    "ATRX", "H3F3A", "TERT", "MDM2", "MDM4"
  )
)

# Combine GO and custom gene sets
all_gene_sets <- c(go_gene_sets, cancer_gene_sets)

cat(sprintf("✓ Total gene sets for GSEA: %d\n", length(all_gene_sets)))
cat(sprintf("  - GO BP: %d\n", length(go_gene_sets)))
cat(sprintf("  - Custom cancer sets: %d\n\n", length(cancer_gene_sets)))

# %%

# =============================================================================
# PHASE 3: RUN FGSEA
# =============================================================================

cat("=== PHASE 3: Running Gene Set Enrichment Analysis ===\n")
cat("This may take several minutes...\n\n")

set.seed(42) # For reproducibility

# Run fgsea
fgsea_results <- fgsea(
  pathways = all_gene_sets,
  stats = gene_ranks,
  minSize = 10,
  maxSize = 500,
  nproc = 8, # Parallel processing
  nPermSimple = 10000
) # Number of permutations

# Filter significant results
fgsea_sig <- fgsea_results %>%
  filter(padj < 0.05) %>%
  arrange(padj)

cat(sprintf("✓ GSEA complete!\n"))
cat(sprintf("  Total gene sets tested: %d\n", nrow(fgsea_results)))
cat(sprintf("  Significant gene sets (FDR < 0.05): %d\n", nrow(fgsea_sig)))
cat(sprintf("  Upregulated (positive NES): %d\n", sum(fgsea_sig$NES > 0)))
cat(sprintf("  Downregulated (negative NES): %d\n\n", sum(fgsea_sig$NES < 0)))

# %%

# =============================================================================
# PHASE 4: FILTER FOR CANCER-RELEVANT PATHWAYS
# =============================================================================

cat("=== PHASE 4: Filtering Cancer-Relevant Pathways ===\n")

# Define keywords - EXPANDED for comprehensive cancer pathway capture
cancer_keywords <- c(
  # Cell death pathways
  "apoptosis", "apoptotic", "cell death", "programmed cell death",
  "necrosis", "necrotic", "ferroptosis", "pyroptosis", "anoikis",
  "autophagy", "autophagic", "survival", "viability", "senescence",

  # Cell cycle and proliferation
  "proliferation", "proliferative", "cell cycle", "mitosis", "mitotic",
  "cell division", "growth", "G1/S", "G2/M", "S phase", "M phase",
  "DNA replication", "chromosome", "spindle", "cytokinesis",
  "cyclin", "checkpoint", "DNA repair", "DNA damage",

  # Migration and invasion
  "migration", "migratory", "motility", "invasion", "invasive",
  "chemotaxis", "chemotactic", "cell movement", "locomotion",
  "adhesion", "cell adhesion", "focal adhesion", "cell junction",
  "cytoskeleton", "actin", "tubulin", "microtubule",

  # Angiogenesis and vasculature
  "angiogenesis", "angiogenic", "blood vessel", "vasculature",
  "endothelial", "VEGF", "vascular",

  # Signaling pathways
  "signaling", "signal transduction", "kinase", "phosphorylation",
  "growth factor", "receptor", "activation", "cascade",

  # Metabolism
  "glycolysis", "metabolism", "metabolic", "glucose", "ATP",
  "oxidative", "respiration", "biosynthesis", "catabolic",

  # Transcription and chromatin
  "transcription", "gene expression", "chromatin", "histone",
  "RNA processing", "splicing", "translation",

  # Cancer-specific terms
  "tumor", "cancer", "oncogenic", "transformation",
  "EMT", "epithelial", "mesenchymal", "stemness",
  "Hippo", "YAP", "TEAD", "Wnt", "Notch", "Hedgehog",

  # Stress response
  "stress", "oxidative stress", "hypoxia", "ER stress",
  "unfolded protein", "heat shock", "inflammatory"
)

pattern <- paste(cancer_keywords, collapse = "|")

fgsea_cancer <- fgsea_sig %>%
  filter(grepl(pattern, pathway, ignore.case = TRUE))

cat(sprintf("✓ Cancer-relevant pathways: %d\n", nrow(fgsea_cancer)))
cat(sprintf("  Upregulated: %d\n", sum(fgsea_cancer$NES > 0)))
cat(sprintf("  Downregulated: %d\n\n", sum(fgsea_cancer$NES < 0)))

# Categorize pathways into functional groups (priority order matters - first match wins)
fgsea_cancer$category <- NA

# Cell death (highest priority for apoptosis-related terms)
fgsea_cancer$category[grepl("apoptosis|apoptotic|death|necrosis|ferroptosis|pyroptosis|anoikis|autophagy|survival|senescence",
  fgsea_cancer$pathway,
  ignore.case = TRUE
)] <- "Cell Death"

# Cell cycle and proliferation
fgsea_cancer$category[grepl("proliferation|cell cycle|mitosis|mitotic|division|growth|G1|G2|S phase|M phase|replication|chromosome|spindle|cyclin|checkpoint",
  fgsea_cancer$pathway,
  ignore.case = TRUE
)] <- "Proliferation"

# Migration and invasion
fgsea_cancer$category[grepl("migration|invasion|motility|chemotaxis|locomotion|adhesion|cytoskeleton|actin|tubulin",
  fgsea_cancer$pathway,
  ignore.case = TRUE
)] <- "Migration"

# Angiogenesis
fgsea_cancer$category[grepl("angiogenesis|angiogenic|blood vessel|vascular|vasculature|endothelial|VEGF",
  fgsea_cancer$pathway,
  ignore.case = TRUE
)] <- "Angiogenesis"

# Metabolism
fgsea_cancer$category[grepl("glycolysis|metabolism|metabolic|glucose|ATP|oxidative|respiration|biosynthesis",
  fgsea_cancer$pathway,
  ignore.case = TRUE
)] <- "Metabolism"

# Signaling pathways
fgsea_cancer$category[grepl("signaling|signal transduction|kinase|phosphorylation|growth factor|receptor|Hippo|YAP|TEAD|Wnt|Notch",
  fgsea_cancer$pathway,
  ignore.case = TRUE
)] <- "Signaling"

# Transcription and chromatin regulation
fgsea_cancer$category[grepl("transcription|gene expression|chromatin|histone|RNA processing|splicing",
  fgsea_cancer$pathway,
  ignore.case = TRUE
)] <- "Transcription"

# Stress response
fgsea_cancer$category[grepl("stress|hypoxia|ER stress|unfolded protein|heat shock|inflammatory",
  fgsea_cancer$pathway,
  ignore.case = TRUE
)] <- "Stress Response"

# EMT and transformation
fgsea_cancer$category[grepl("EMT|epithelial|mesenchymal|transformation|tumor|cancer|oncogenic",
  fgsea_cancer$pathway,
  ignore.case = TRUE
)] <- "EMT/Transformation"

# %%

# =============================================================================
# PHASE 5: EXPORT RESULTS
# =============================================================================

cat("=== PHASE 5: Exporting Results ===\n")

# Convert list columns to character strings for CSV export
fgsea_results_export <- fgsea_results %>%
  mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";"))

fgsea_sig_export <- fgsea_sig %>%
  mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";"))

fgsea_cancer_export <- fgsea_cancer %>%
  mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";"))

# Export all results
write.csv(fgsea_results_export, file.path(output_dir, "fgsea_all_pathways.csv"), row.names = FALSE)
write.csv(fgsea_sig_export, file.path(output_dir, "fgsea_significant_pathways.csv"), row.names = FALSE)
write.csv(fgsea_cancer_export, file.path(output_dir, "fgsea_cancer_pathways.csv"), row.names = FALSE)

cat("✓ Results exported\n\n")

# %%

# =============================================================================
# PHASE 6: VISUALIZATIONS
# =============================================================================

cat("=== PHASE 6: Creating Visualizations ===\n")

# Plot 1: Enrichment score plot for top pathways
cat("Creating enrichment score plots...\n")

top_pathways_up <- fgsea_cancer %>%
  filter(NES > 0) %>%
  arrange(desc(NES)) %>%
  head(10)

top_pathways_down <- fgsea_cancer %>%
  filter(NES < 0) %>%
  arrange(NES) %>%
  head(10)

top_pathways <- rbind(top_pathways_up, top_pathways_down)

if (nrow(top_pathways) > 0) {
  pdf(file.path(output_dir, "01_top_pathways_barplot.pdf"), width = 14, height = 10)
  p1 <- ggplot(top_pathways, aes(x = reorder(pathway, NES), y = NES, fill = NES > 0)) +
    geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
    coord_flip() +
    labs(
      title = "Top 20 Cancer Pathways from GSEA",
      subtitle = "Normalized Enrichment Score (NES)",
      x = "Pathway",
      y = "NES",
      fill = "Direction"
    ) +
    scale_fill_manual(
      values = c("FALSE" = "#1F78B4", "TRUE" = "#E31A1C"),
      labels = c("Downregulated", "Upregulated")
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      axis.text.y = element_text(size = 9)
    )
  print(p1)
  dev.off()
}

# Plot 2: Enrichment plots for specific pathways
cat("Creating detailed enrichment plots for key pathways...\n")

# Select interesting pathways
key_pathway_names <- c()
if ("APOPTOSIS" %in% names(all_gene_sets)) key_pathway_names <- c(key_pathway_names, "APOPTOSIS")
if ("CELL_CYCLE" %in% names(all_gene_sets)) key_pathway_names <- c(key_pathway_names, "CELL_CYCLE")
if ("MIGRATION" %in% names(all_gene_sets)) key_pathway_names <- c(key_pathway_names, "MIGRATION")
if ("HIPPO_YAP" %in% names(all_gene_sets)) key_pathway_names <- c(key_pathway_names, "HIPPO_YAP")

if (length(key_pathway_names) > 0) {
  pdf(file.path(output_dir, "02_detailed_enrichment_plots.pdf"), width = 12, height = 8)
  for (pathway_name in key_pathway_names) {
    p <- plotEnrichment(all_gene_sets[[pathway_name]], gene_ranks) +
      labs(title = paste("GSEA Enrichment Plot:", pathway_name)) +
      theme_bw(base_size = 14)
    print(p)
  }
  dev.off()
}

# Plot 3: Bubble plot of cancer pathways
cat("Creating bubble plot...\n")

if (nrow(fgsea_cancer) > 0 && !all(is.na(fgsea_cancer$category))) {
  pdf(file.path(output_dir, "03_cancer_pathways_bubble.pdf"), width = 14, height = 10)
  p3 <- ggplot(
    fgsea_cancer %>% filter(!is.na(category)),
    aes(x = NES, y = -log10(padj), color = category, size = size)
  ) +
    geom_point(alpha = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
    labs(
      title = "Cancer Pathway Enrichment (GSEA)",
      subtitle = "Bubble size = number of genes in pathway",
      x = "Normalized Enrichment Score (NES)",
      y = "-log10(Adjusted p-value)",
      color = "Category",
      size = "Gene Set Size"
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5)
    ) +
    scale_color_brewer(palette = "Set1")
  print(p3)
  dev.off()
}

# Plot 4: Category summary
cat("Creating category summary plot...\n")

if (!all(is.na(fgsea_cancer$category))) {
  category_summary <- fgsea_cancer %>%
    filter(!is.na(category)) %>%
    mutate(direction = ifelse(NES > 0, "Up", "Down")) %>%
    group_by(category, direction) %>%
    summarise(count = n(), .groups = "drop")

  pdf(file.path(output_dir, "04_category_summary.pdf"), width = 10, height = 7)
  p4 <- ggplot(category_summary, aes(x = category, y = count, fill = direction)) +
    geom_bar(stat = "identity", position = "dodge", color = "black", linewidth = 0.3) +
    labs(
      title = "Cancer Pathway Categories (GSEA)",
      subtitle = "Number of significantly enriched pathways per category",
      x = "Category",
      y = "Number of Pathways",
      fill = "Direction"
    ) +
    scale_fill_manual(values = c("Up" = "#E31A1C", "Down" = "#1F78B4")) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5)
    )
  print(p4)
  dev.off()
}

cat("✓ All visualizations created\n\n")

# %%

# =============================================================================
# PHASE 7: SUMMARY REPORT
# =============================================================================

cat("=== PHASE 7: Generating Summary Report ===\n")

summary_file <- file.path(output_dir, "GSEA_SUMMARY.txt")
cat("TRUE GSEA ANALYSIS SUMMARY\n", file = summary_file)
cat("==========================\n\n", file = summary_file, append = TRUE)
cat(paste("Generated:", Sys.time(), "\n\n"), file = summary_file, append = TRUE)

cat("METHOD:\n", file = summary_file, append = TRUE)
cat("  - Algorithm: fgsea (Fast Gene Set Enrichment Analysis)\n", file = summary_file, append = TRUE)
cat("  - Ranking metric: log2FoldChange (TES vs GFP)\n", file = summary_file, append = TRUE)
cat(sprintf("  - Total genes ranked: %d\n", length(gene_ranks)), file = summary_file, append = TRUE)
cat(sprintf("  - Gene sets tested: %d\n", nrow(fgsea_results)), file = summary_file, append = TRUE)
cat(sprintf("  - Permutations: 10,000\n\n"), file = summary_file, append = TRUE)

cat("RESULTS:\n", file = summary_file, append = TRUE)
cat(sprintf("  Significant pathways (FDR < 0.05): %d\n", nrow(fgsea_sig)), file = summary_file, append = TRUE)
cat(sprintf("    - Positively enriched (NES > 0): %d\n", sum(fgsea_sig$NES > 0)), file = summary_file, append = TRUE)
cat(sprintf("    - Negatively enriched (NES < 0): %d\n\n", sum(fgsea_sig$NES < 0)), file = summary_file, append = TRUE)

cat("CANCER-RELEVANT PATHWAYS:\n", file = summary_file, append = TRUE)
cat(sprintf("  Total cancer pathways: %d\n", nrow(fgsea_cancer)), file = summary_file, append = TRUE)

if (!all(is.na(fgsea_cancer$category))) {
  for (cat_name in unique(fgsea_cancer$category[!is.na(fgsea_cancer$category)])) {
    cat_pathways <- fgsea_cancer %>% filter(category == cat_name)
    cat(
      sprintf(
        "    %s: %d (%d up, %d down)\n",
        cat_name,
        nrow(cat_pathways),
        sum(cat_pathways$NES > 0),
        sum(cat_pathways$NES < 0)
      ),
      file = summary_file, append = TRUE
    )
  }
}

cat("\n\nTOP 10 UPREGULATED PATHWAYS:\n", file = summary_file, append = TRUE)
if (nrow(fgsea_sig) > 0) {
  top_up <- fgsea_sig %>%
    filter(NES > 0) %>%
    arrange(desc(NES)) %>%
    head(10)
  for (i in 1:min(nrow(top_up), 10)) {
    cat(sprintf("  %d. %s (NES=%.2f, FDR=%.2e)\n", i, top_up$pathway[i], top_up$NES[i], top_up$padj[i]),
      file = summary_file, append = TRUE
    )
  }
}

cat("\n\nTOP 10 DOWNREGULATED PATHWAYS:\n", file = summary_file, append = TRUE)
if (nrow(fgsea_sig) > 0) {
  top_down <- fgsea_sig %>%
    filter(NES < 0) %>%
    arrange(NES) %>%
    head(10)
  for (i in 1:min(nrow(top_down), 10)) {
    cat(sprintf("  %d. %s (NES=%.2f, FDR=%.2e)\n", i, top_down$pathway[i], top_down$NES[i], top_down$padj[i]),
      file = summary_file, append = TRUE
    )
  }
}

cat("\n✓ Summary report saved\n\n")

cat("========================================\n")
cat("TRUE GSEA ANALYSIS COMPLETE\n")
cat("========================================\n")
cat("Completed:", as.character(Sys.time()), "\n")
cat(sprintf("Output directory: %s\n\n", output_dir))
cat("Key files:\n")
cat("  - fgsea_all_pathways.csv (all tested pathways)\n")
cat("  - fgsea_significant_pathways.csv (FDR < 0.05)\n")
cat("  - fgsea_cancer_pathways.csv (cancer-relevant subset)\n")
cat("  - 01_top_pathways_barplot.pdf\n")
cat("  - 02_detailed_enrichment_plots.pdf\n")
cat("  - 03_cancer_pathways_bubble.pdf\n")
cat("  - 04_category_summary.pdf\n")
cat("  - GSEA_SUMMARY.txt\n\n")
cat("This is TRUE GSEA using ranked gene lists!\n")
cat("All genes contribute to enrichment score calculation.\n")
