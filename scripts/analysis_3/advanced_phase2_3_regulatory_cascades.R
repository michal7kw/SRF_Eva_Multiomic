#!/usr/bin/env Rscript

################################################################################
# Phase 2.3: Regulatory Cascades and Secondary Targets
#
# Purpose: Infer regulatory cascades and identify secondary targets
#
# Author: Advanced Multi-Omics Analysis Plan
# Date: 2025-01-24
################################################################################

message("=== Phase 2.3: Regulatory Cascades and Secondary Targets ===")
message("Start time: ", Sys.time())

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(igraph)
  library(visNetwork)
  library(msigdbr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggrepel)  # For geom_text_repel
})

# Fix namespace conflicts: dplyr functions get masked by AnnotationDbi
select <- dplyr::select
filter <- dplyr::filter
slice <- dplyr::slice
mutate <- dplyr::mutate
arrange <- dplyr::arrange

setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top")

# Input files
RNA_SEQ_FILE <- "SRF_Eva_RNA/results/05_deseq2/deseq2_results_TES_vs_GFP.txt"
GENE_BINDING_EXPR <- "SRF_Eva_integrated_analysis/scripts/analysis_3/results/04_category_expression/genes_with_binding_and_expression.csv"
# NOTE: Master regulators file is OPTIONAL - Phase 5.2 output
# If not available, this script will still work but without master regulator integration
MASTER_REGULATORS <- "SRF_Eva_integrated_analysis/scripts/analysis_3/results/09_tf_networks/master_regulator_analysis.csv"

# Output directories - separate for detailed vs simplified
OUTPUT_DIR <- "SRF_Eva_integrated_analysis/scripts/analysis_3/results/06_regulatory_cascades"
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
# Step 1: Load data
################################################################################

message("\n[Step 1] Loading data...")

# Load RNA-seq results
rna_results <- read.delim(RNA_SEQ_FILE, stringsAsFactors = FALSE)
message("  Loaded RNA-seq: ", nrow(rna_results), " genes")

# Load gene-binding data
gene_data <- read.csv(GENE_BINDING_EXPR)
message("  Loaded gene-binding: ", nrow(gene_data), " genes")

# Load master regulators (OPTIONAL - may not exist if Phase 5.2 hasn't run)
master_regs <- NULL
if (file.exists(MASTER_REGULATORS)) {
  master_regs <- read.csv(MASTER_REGULATORS)
  message("  Loaded master regulators: ", nrow(master_regs), " TFs")
} else {
  message("  NOTE: Master regulators file not found (Phase 5.2 not yet run)")
  message("        This script will proceed without master regulator integration")
  message("        Re-run after Phase 5.2 for full analysis")
}

################################################################################
# Step 2: Identify primary transcription factor targets
################################################################################

message("\n[Step 2] Identifying primary TF targets...")

# Get TF annotation from GO
tf_go_terms <- c("GO:0003700", "GO:0000981", "GO:0001228")

tf_annotations <- lapply(tf_go_terms, function(go_term) {
  genes <- get(go_term, org.Hs.egGO2ALLEGS)
  symbols <- mapIds(org.Hs.eg.db, keys = genes, column = "SYMBOL",
                   keytype = "ENTREZID", multiVals = "first")
  symbols[!is.na(symbols)]
})

all_tfs <- unique(unlist(tf_annotations))
message("  Total TFs annotated: ", length(all_tfs))

# Identify primary TF targets (TES/TEAD1 bound + DEG + is a TF)
primary_tf_targets <- gene_data %>%
  filter(
    gene_name %in% all_tfs,
    primary_category != "Unbound",
    is_DEG
  )

message("  Primary TF targets (bound + DE): ", nrow(primary_tf_targets))

# Classify by binding category
primary_tf_summary <- primary_tf_targets %>%
  group_by(primary_category) %>%
  summarise(
    n_tfs = n(),
    n_up = sum(regulation == "Up"),
    n_down = sum(regulation == "Down"),
    .groups = "drop"
  )

message("\n  Primary TF targets by category:")
print(primary_tf_summary)

################################################################################
# Step 3: Get TF target gene sets
################################################################################

message("\n[Step 3] Loading TF target gene sets...")

# Get TF targets from MSigDB
msigdb_tfs <- msigdbr(species = "Homo sapiens", category = "C3", subcategory = "TFT:GTRD")

# Filter for our primary TF targets
primary_tf_symbols <- primary_tf_targets$gene_name

relevant_tf_sets <- msigdb_tfs %>%
  mutate(tf_name = toupper(gsub("_.*", "", gs_name))) %>%
  filter(tf_name %in% toupper(primary_tf_symbols))

message("  Found gene sets for ", length(unique(relevant_tf_sets$gs_name)), " primary TFs")

################################################################################
# Step 4: Identify secondary targets
################################################################################

message("\n[Step 4] Identifying secondary targets...")

# For each primary TF, identify its targets among DEGs
cascade_results <- list()

for (tf in primary_tf_symbols) {
  # Get TF's expression
  tf_data <- primary_tf_targets %>% filter(gene_name == tf)

  if (nrow(tf_data) == 0) next

  # Get this TF's target genes
  tf_targets <- relevant_tf_sets %>%
    filter(toupper(gsub("_.*", "", gs_name)) == toupper(tf)) %>%
    pull(gene_symbol) %>%
    unique()

  if (length(tf_targets) < 10) next  # Need sufficient targets

  # Find which targets are DEGs but NOT directly bound by TES/TEAD1
  secondary_targets <- rna_results %>%
    filter(
      gene_symbol %in% tf_targets,
      !is.na(padj),
      padj < 0.05
    ) %>%
    left_join(gene_data %>% select(gene_name, primary_category),
             by = c("gene_symbol" = "gene_name")) %>%
    mutate(
      is_indirect = is.na(primary_category) | primary_category == "Unbound"
    )

  n_indirect <- sum(secondary_targets$is_indirect)

  if (n_indirect < 5) next  # Need sufficient indirect targets

  cascade_results[[tf]] <- list(
    tf_name = tf,
    tf_log2FC = tf_data$log2FoldChange,
    tf_padj = tf_data$padj,
    tf_category = tf_data$primary_category,
    n_total_targets = length(tf_targets),
    n_targets_deg = nrow(secondary_targets),
    n_indirect_targets = n_indirect,
    pct_indirect = round(n_indirect / nrow(secondary_targets) * 100, 2),
    indirect_targets = secondary_targets %>%
      filter(is_indirect) %>%
      pull(gene_symbol)
  )
}

message("  Found ", length(cascade_results), " TFs with secondary targets")

# Convert to data frame
cascade_df <- do.call(rbind, lapply(cascade_results, function(x) {
  data.frame(
    TF = x$tf_name,
    TF_log2FC = x$tf_log2FC,
    TF_padj = x$tf_padj,
    TF_binding_category = x$tf_category,
    n_total_targets = x$n_total_targets,
    n_targets_DEG = x$n_targets_deg,
    n_indirect_targets = x$n_indirect_targets,
    pct_indirect = x$pct_indirect,
    stringsAsFactors = FALSE
  )
}))

# Sort by number of indirect targets
cascade_df <- cascade_df %>% arrange(desc(n_indirect_targets))

################################################################################
# Step 5: Build regulatory network
################################################################################

message("\n[Step 5] Building regulatory network...")

# Create network edges
network_edges <- list()

for (tf in names(cascade_results)) {
  cascade <- cascade_results[[tf]]

  # TES/TEAD1 -> TF edge
  network_edges[[length(network_edges) + 1]] <- data.frame(
    from = "TES/TEAD1",
    to = tf,
    edge_type = "direct",
    regulation = ifelse(cascade$tf_log2FC > 0, "activation", "repression"),
    stringsAsFactors = FALSE
  )

  # TF -> secondary target edges (top 10 per TF to keep manageable)
  secondary <- head(cascade$indirect_targets, 10)

  for (target in secondary) {
    # Get target's expression
    target_expr <- rna_results %>%
      filter(gene_symbol == target) %>%
      slice(1)

    if (nrow(target_expr) > 0) {
      network_edges[[length(network_edges) + 1]] <- data.frame(
        from = tf,
        to = target,
        edge_type = "indirect",
        regulation = ifelse(target_expr$log2FoldChange > 0, "activation", "repression"),
        stringsAsFactors = FALSE
      )
    }
  }
}

network_edges_df <- do.call(rbind, network_edges)

# Create network nodes
all_nodes <- unique(c(network_edges_df$from, network_edges_df$to))

network_nodes <- data.frame(
  id = all_nodes,
  label = all_nodes,
  node_type = case_when(
    all_nodes == "TES/TEAD1" ~ "source",
    all_nodes %in% cascade_df$TF ~ "primary_TF",
    TRUE ~ "secondary_target"
  ),
  stringsAsFactors = FALSE
)

# Add expression data to nodes
network_nodes <- network_nodes %>%
  left_join(
    rna_results %>% select(gene_symbol, log2FoldChange, padj),
    by = c("id" = "gene_symbol")
  )

################################################################################
# Step 6: Export results
################################################################################

message("\n[Step 6] Exporting results...")

write.csv(cascade_df,
          file.path(OUTPUT_DIR, "regulatory_cascades_summary.csv"),
          row.names = FALSE)

write.csv(primary_tf_targets,
          file.path(OUTPUT_DIR, "primary_TF_targets.csv"),
          row.names = FALSE)

write.csv(network_edges_df,
          file.path(OUTPUT_DIR, "regulatory_network_edges.csv"),
          row.names = FALSE)

write.csv(network_nodes,
          file.path(OUTPUT_DIR, "regulatory_network_nodes.csv"),
          row.names = FALSE)

# Export detailed secondary targets for each TF
for (tf in names(cascade_results)) {
  cascade <- cascade_results[[tf]]

  secondary_genes <- rna_results %>%
    filter(gene_symbol %in% cascade$indirect_targets) %>%
    select(gene_symbol, baseMean, log2FoldChange, padj)

  write.csv(secondary_genes,
            file.path(OUTPUT_DIR, paste0(tf, "_secondary_targets.csv")),
            row.names = FALSE)
}

################################################################################
# Step 7: Visualizations
################################################################################

message("\n[Step 7] Generating visualizations...")

# 7.1: Cascade overview barplot
pdf(file.path(OUTPUT_DIR, "cascade_TF_overview.pdf"), width = 12, height = 8)
top_cascades <- head(cascade_df, 15)

ggplot(top_cascades, aes(x = reorder(TF, n_indirect_targets), y = n_indirect_targets)) +
  geom_bar(aes(fill = TF_binding_category), stat = "identity") +
  geom_text(aes(label = n_indirect_targets), hjust = -0.2, size = 3) +
  coord_flip() +
  scale_fill_brewer(palette = "Set2") +
  theme_classic() +
  labs(
    title = "Regulatory Cascades: Primary TFs and Their Secondary Targets",
    subtitle = "Number of indirect DEGs regulated by each TF",
    x = "Transcription Factor",
    y = "Number of Secondary Targets (Indirect DEGs)",
    fill = "TF Binding\nCategory"
  )
dev.off()

# 7.2: TF expression vs cascade size
pdf(file.path(OUTPUT_DIR, "TF_expression_vs_cascade_size.pdf"), width = 10, height = 8)
ggplot(cascade_df, aes(x = TF_log2FC, y = n_indirect_targets)) +
  geom_point(aes(color = TF_binding_category, size = n_targets_DEG), alpha = 0.7) +
  geom_text_repel(aes(label = TF), size = 3, max.overlaps = 15) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_brewer(palette = "Set2") +
  scale_size_continuous(name = "Total DEG\nTargets") +
  theme_classic() +
  labs(
    title = "TF Expression Change vs Cascade Size",
    subtitle = "Do more strongly regulated TFs have larger cascades?",
    x = "TF log2 Fold Change (TES vs GFP)",
    y = "Number of Secondary Targets",
    color = "TF Binding\nCategory"
  )
dev.off()

# 7.3: Network visualization with igraph
if (nrow(network_edges_df) > 0) {
  # Create igraph object
  g <- graph_from_data_frame(network_edges_df, directed = TRUE, vertices = network_nodes)

  # Set node colors
  V(g)$color <- case_when(
    V(g)$node_type == "source" ~ "gold",
    V(g)$node_type == "primary_TF" ~ "skyblue",
    TRUE ~ "lightgray"
  )

  # Set node sizes
  V(g)$size <- case_when(
    V(g)$node_type == "source" ~ 20,
    V(g)$node_type == "primary_TF" ~ 10,
    TRUE ~ 5
  )

  # Set edge colors
  E(g)$color <- ifelse(E(g)$edge_type == "direct", "red", "gray")
  E(g)$width <- ifelse(E(g)$edge_type == "direct", 2, 0.5)

  # Plot network (simplified - show only top cascades)
  top_tfs <- head(cascade_df$TF, 5)
  nodes_to_keep <- c("TES/TEAD1",
                     top_tfs,
                     network_edges_df$to[network_edges_df$from %in% top_tfs])

  g_subset <- induced_subgraph(g, V(g)$name %in% nodes_to_keep)

  pdf(file.path(OUTPUT_DIR, "regulatory_cascade_network.pdf"), width = 14, height = 14)
  plot(g_subset,
       layout = layout_with_fr(g_subset),
       vertex.label.cex = 0.7,
       vertex.label.color = "black",
       edge.arrow.size = 0.3,
       main = "Regulatory Cascade Network (Top 5 TFs)")
  dev.off()
}

# 7.4: Interactive network with visNetwork (HTML output)
if (nrow(network_edges_df) > 0) {
  # Prepare for visNetwork
  vis_nodes <- network_nodes %>%
    mutate(
      title = paste0(id, "\nlog2FC: ", round(log2FoldChange, 2)),
      group = node_type,
      size = case_when(
        node_type == "source" ~ 30,
        node_type == "primary_TF" ~ 20,
        TRUE ~ 10
      )
    )

  vis_edges <- network_edges_df %>%
    mutate(
      color = ifelse(edge_type == "direct", "red", "gray"),
      width = ifelse(edge_type == "direct", 3, 1)
    )

  # Create interactive network (top cascades only)
  vis_nodes_subset <- vis_nodes %>%
    filter(id %in% c("TES/TEAD1", top_tfs) |
           id %in% vis_edges$to[vis_edges$from %in% top_tfs])

  vis_edges_subset <- vis_edges %>%
    filter(from %in% vis_nodes_subset$id & to %in% vis_nodes_subset$id)

  network_plot <- visNetwork(vis_nodes_subset, vis_edges_subset, width = "100%") %>%
    visGroups(groupname = "source", color = "gold") %>%
    visGroups(groupname = "primary_TF", color = "skyblue") %>%
    visGroups(groupname = "secondary_target", color = "lightgray") %>%
    visLegend() %>%
    visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
    visLayout(randomSeed = 42)

  visSave(network_plot, file.path(OUTPUT_DIR, "interactive_cascade_network.html"))
}

################################################################################
# Step 8: Functional enrichment of secondary targets
################################################################################

message("\n[Step 8] Enriching secondary target pathways...")

# For top TFs, perform GO enrichment on their secondary targets
enrichment_results <- list()

for (tf in head(cascade_df$TF, 10)) {
  cascade <- cascade_results[[tf]]

  if (length(cascade$indirect_targets) < 10) next

  # Convert to Entrez IDs
  entrez_ids <- mapIds(org.Hs.eg.db,
                       keys = cascade$indirect_targets,
                       column = "ENTREZID",
                       keytype = "SYMBOL",
                       multiVals = "first")

  entrez_ids <- entrez_ids[!is.na(entrez_ids)]

  if (length(entrez_ids) < 10) next

  # GO enrichment
  ego <- enrichGO(
    gene = entrez_ids,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05,
    readable = TRUE
  )

  if (!is.null(ego) && nrow(ego@result) > 0) {
    enrichment_results[[tf]] <- ego@result %>%
      head(10) %>%
      mutate(TF = tf)
  }
}

if (length(enrichment_results) > 0) {
  all_enrichment <- do.call(rbind, enrichment_results)

  write.csv(all_enrichment,
            file.path(OUTPUT_DIR, "secondary_target_pathway_enrichment.csv"),
            row.names = FALSE)
}

################################################################################
# Step 9: Summary report
################################################################################

message("\n[Step 9] Creating summary report...")

sink(file.path(OUTPUT_DIR, "PHASE2_3_SUMMARY.txt"))
cat("=== Phase 2.3: Regulatory Cascades and Secondary Targets ===\n")
cat("Date:", as.character(Sys.time()), "\n\n")

cat("Primary TF Targets (Direct TES/TEAD1 targets that are TFs):\n")
cat("============================================================\n")
cat("Total primary TF targets:", nrow(primary_tf_targets), "\n\n")
print(primary_tf_summary)

cat("\n\nRegulatory Cascades Identified:\n")
cat("================================\n")
cat("TFs with secondary targets:", nrow(cascade_df), "\n\n")

cat("Top 10 Cascades:\n")
print(head(cascade_df %>%
           select(TF, TF_log2FC, TF_binding_category,
                  n_indirect_targets, pct_indirect), 10))

cat("\n\nRegulatory Network Statistics:\n")
cat("===============================\n")
cat("Total nodes:", nrow(network_nodes), "\n")
cat("  - Source (TES/TEAD1):", sum(network_nodes$node_type == "source"), "\n")
cat("  - Primary TFs:", sum(network_nodes$node_type == "primary_TF"), "\n")
cat("  - Secondary targets:", sum(network_nodes$node_type == "secondary_target"), "\n")
cat("Total edges:", nrow(network_edges_df), "\n")
cat("  - Direct (TES/TEAD1 -> TF):", sum(network_edges_df$edge_type == "direct"), "\n")
cat("  - Indirect (TF -> secondary):", sum(network_edges_df$edge_type == "indirect"), "\n")

sink()

message("\n=== Analysis Complete ===")
message("Output directory: ", OUTPUT_DIR)
message("End time: ", Sys.time())
message("\nKey findings:")
message("  - ", nrow(primary_tf_targets), " TFs are direct TES/TEAD1 targets")
message("  - ", nrow(cascade_df), " TFs have identifiable regulatory cascades")
message("  - Total secondary targets: ", sum(cascade_df$n_indirect_targets))
