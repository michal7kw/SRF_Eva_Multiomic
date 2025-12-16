#!/usr/bin/env Rscript

################################################################################
# Phase 5.3: Pathway Crosstalk Analysis
# TES vs TEAD1 Comparative Study (Excluding TESmut)
#
# Purpose: Identify hub pathways and pathway interactions
#
# Analysis:
# 1. Load enriched pathways from all previous analyses
# 2. Build pathway-gene network
# 3. Calculate pathway similarity (gene overlap)
# 4. Identify hub pathways (high connectivity)
# 5. Community detection (pathway modules)
# 6. Pathway-pathway correlation analysis
# 7. Visualize pathway network
#
# Author: Advanced Multi-Omics Analysis Pipeline
# Date: 2025-01-24
################################################################################

library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(ggplot2)
library(igraph)
library(ggraph)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(visNetwork)

# Fix namespace conflicts: dplyr functions get masked by AnnotationDbi
select <- dplyr::select
filter <- dplyr::filter
slice <- dplyr::slice
mutate <- dplyr::mutate
arrange <- dplyr::arrange

# Set working directory
setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top")

# Create output directories - separate for detailed vs simplified
output_dir <- "SRF_Eva_integrated_analysis/scripts/analysis_3/results/11_pathway_crosstalk"
output_dir_detailed <- file.path(output_dir, "detailed_6cat")
output_dir_simple <- file.path(output_dir, "simplified_3cat")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir_detailed, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir_simple, recursive = TRUE, showWarnings = FALSE)

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

# Logging function
log_message <- function(msg) {
  cat(paste0("[", Sys.time(), "] ", msg, "\n"))
}

log_message("Starting Phase 5.3: Pathway Crosstalk Analysis")

################################################################################
# 1. Load Pathway Enrichment Results
################################################################################

log_message("Loading pathway enrichment results...")

# Try multiple possible locations for GO enrichment results
possible_go_dirs <- c(
  "SRF_Eva_integrated_analysis/scripts/analysis_3/results/04_category_expression",
  "SRF_Eva_integrated_analysis/scripts/analysis_1/output/10_final_integrative_analysis/pathway_analysis",
  "SRF_Eva_integrated_analysis/scripts/analysis_1/output/11_final_integrative_analysis_all_genes/pathway_analysis"
)

all_go_results <- list()

for (go_dir in possible_go_dirs) {
  if (!dir.exists(go_dir)) next

  # Try multiple patterns
  go_files <- list.files(go_dir, pattern = "GO_enrichment.*\\.csv$", full.names = TRUE)

  if (length(go_files) == 0) next

  log_message(sprintf("  Found %d GO files in %s", length(go_files), go_dir))

  for (file in go_files) {
    # Extract category from filename
    fname <- basename(file)
    category <- gsub("_GO_enrichment.*\\.csv$", "", fname)

    tryCatch({
      go_data <- read_csv(file, show_col_types = FALSE)

      # Check for required columns
      if (!"ID" %in% colnames(go_data) || !"geneID" %in% colnames(go_data)) {
        log_message(sprintf("    Skipping %s: missing required columns", fname))
        next
      }

      # Add source and category columns if not present
      go_data <- go_data %>%
        mutate(source = "GO", category = category)

      all_go_results[[paste0(go_dir, "_", category)]] <- go_data
    }, error = function(e) {
      log_message(sprintf("    Error reading %s: %s", fname, e$message))
    })
  }
}

go_combined <- bind_rows(all_go_results)

log_message(sprintf("  Loaded %d GO enrichment results from %d files",
                   nrow(go_combined), length(all_go_results)))

# Check if we have any data to work with
if (nrow(go_combined) == 0) {
  log_message("WARNING: No GO enrichment data found. Creating minimal outputs...")

  # Create empty/minimal output files
  write_csv(data.frame(pathway = character(), degree = integer()),
            file.path(output_dir, "hub_pathways_ranked.csv"))
  write_csv(data.frame(community = integer(), n_pathways = integer()),
            file.path(output_dir, "pathway_communities.csv"))
  write_csv(data.frame(ID = character(), n_categories = integer()),
            file.path(output_dir, "pathway_category_overlap.csv"))

  # Create summary report
  sink(file.path(output_dir, "PHASE5_3_SUMMARY.txt"))
  cat("================================================================================\n")
  cat("Phase 5.3: Pathway Crosstalk Analysis - Analysis Summary\n")
  cat("================================================================================\n\n")
  cat("WARNING: No GO enrichment data found.\n")
  cat("Please run Phase 2.1 (expression by category) or Phase 10 (integrative analysis) first.\n\n")
  cat("Expected locations checked:\n")
  for (d in possible_go_dirs) {
    cat(sprintf("  - %s\n", d))
  }
  sink()

  log_message("Phase 5.3 completed (no data to analyze)")
  quit(save = "no", status = 0)
}

# Check if qvalue column exists (might be named differently)
if (!"qvalue" %in% colnames(go_combined)) {
  if ("p.adjust" %in% colnames(go_combined)) {
    go_combined$qvalue <- go_combined$p.adjust
  } else if ("pvalue" %in% colnames(go_combined)) {
    # Use raw p-value if no adjusted available
    go_combined$qvalue <- go_combined$pvalue
    log_message("  WARNING: Using raw p-values as qvalue (no adjusted p-values found)")
  } else {
    # Create a placeholder qvalue column with 0 (all pass filter)
    go_combined$qvalue <- 0.01
    log_message("  WARNING: No p-value column found, using all pathways")
  }
}

# Filter for significant pathways
sig_pathways <- go_combined %>%
  filter(qvalue < 0.05)

log_message(sprintf("  %d significant pathways (qvalue < 0.05)", nrow(sig_pathways)))

################################################################################
# 2. Build Pathway-Gene Bipartite Network
################################################################################

log_message("Building pathway-gene bipartite network...")

# Parse geneID column (format: "GENE1/GENE2/GENE3")
pathway_gene_edges <- sig_pathways %>%
  select(ID, Description, category, qvalue, geneID) %>%
  mutate(genes = str_split(geneID, "/")) %>%
  unnest(genes) %>%
  select(pathway = ID, pathway_name = Description, gene = genes, category, qvalue)

log_message(sprintf("  Created %d pathway-gene edges", nrow(pathway_gene_edges)))

# Get unique pathways and genes
unique_pathways <- sig_pathways %>%
  select(ID, Description, category, qvalue, Count) %>%
  distinct()

unique_genes <- pathway_gene_edges %>%
  distinct(gene)

log_message(sprintf("  %d unique pathways, %d unique genes",
                   nrow(unique_pathways), nrow(unique_genes)))

################################################################################
# 3. Calculate Pathway Similarity (Jaccard Index)
################################################################################

log_message("Calculating pathway similarity...")

# Create pathway-gene membership matrix
pathway_list <- split(pathway_gene_edges$gene, pathway_gene_edges$pathway)

# Function to calculate Jaccard similarity
jaccard_similarity <- function(set1, set2) {
  intersection <- length(intersect(set1, set2))
  union <- length(union(set1, set2))
  return(intersection / union)
}

# Calculate pairwise similarity for all pathway pairs
pathway_ids <- names(pathway_list)
n_pathways <- length(pathway_ids)

# Create similarity matrix (only for top pathways to save memory)
# Select top 100 pathways by significance
top_pathways <- unique_pathways %>%
  arrange(qvalue) %>%
  head(100)

top_pathway_ids <- top_pathways$ID
top_pathway_list <- pathway_list[top_pathway_ids]

log_message(sprintf("  Computing similarity for top %d pathways...", length(top_pathway_ids)))

similarity_matrix <- matrix(0, nrow = length(top_pathway_ids), ncol = length(top_pathway_ids))
rownames(similarity_matrix) <- top_pathway_ids
colnames(similarity_matrix) <- top_pathway_ids

for (i in 1:length(top_pathway_ids)) {
  for (j in i:length(top_pathway_ids)) {
    sim <- jaccard_similarity(top_pathway_list[[i]], top_pathway_list[[j]])
    similarity_matrix[i, j] <- sim
    similarity_matrix[j, i] <- sim
  }

  if (i %% 10 == 0) {
    log_message(sprintf("    Progress: %d/%d pathways", i, length(top_pathway_ids)))
  }
}

log_message("  Similarity matrix complete")

################################################################################
# 4. Build Pathway-Pathway Network
################################################################################

log_message("Building pathway-pathway similarity network...")

# Convert similarity matrix to edge list (threshold: Jaccard > 0.2)
SIMILARITY_THRESHOLD <- 0.2

edges <- data.frame()
for (i in 1:(nrow(similarity_matrix)-1)) {
  for (j in (i+1):ncol(similarity_matrix)) {
    sim <- similarity_matrix[i, j]
    if (sim > SIMILARITY_THRESHOLD && sim < 1.0) {  # Exclude self-loops
      edges <- rbind(edges, data.frame(
        from = rownames(similarity_matrix)[i],
        to = colnames(similarity_matrix)[j],
        weight = sim
      ))
    }
  }
}

log_message(sprintf("  Created %d pathway-pathway edges (similarity > %.2f)",
                   nrow(edges), SIMILARITY_THRESHOLD))

# Create igraph network
if (nrow(edges) > 0) {
  # Ensure no duplicate pathway IDs (can happen when combining multiple GO files)
  top_pathways_unique <- top_pathways %>%
    dplyr::distinct(ID, .keep_all = TRUE) %>%
    dplyr::select(ID, Description, category, qvalue)

  pathway_network <- graph_from_data_frame(edges, directed = FALSE,
                                          vertices = top_pathways_unique)

  # Calculate network statistics
  # Use igraph:: prefix to avoid circlize masking
  degree_centrality <- as.vector(igraph::degree(pathway_network))
  betweenness_centrality <- as.vector(igraph::betweenness(pathway_network))
  eigenvector_centrality <- as.vector(igraph::eigen_centrality(pathway_network)$vector)

  # Add to vertex attributes
  V(pathway_network)$degree <- degree_centrality
  V(pathway_network)$betweenness <- betweenness_centrality
  V(pathway_network)$eigenvector <- eigenvector_centrality

  log_message("  Network statistics calculated")
} else {
  log_message("  WARNING: No pathway edges created. Check similarity threshold.")
  pathway_network <- NULL
}

################################################################################
# 5. Identify Hub Pathways
################################################################################

log_message("Identifying hub pathways...")

if (!is.null(pathway_network)) {
  # Extract vertex attributes safely (some may be lists after merging)
  # Build vectors first, then create data.frame
  v_pathway <- V(pathway_network)$name
  v_desc <- sapply(V(pathway_network)$Description, function(x) if(is.null(x)) NA_character_ else as.character(x[1]))
  v_cat <- sapply(V(pathway_network)$category, function(x) if(is.null(x)) NA_character_ else as.character(x[1]))
  v_qval <- sapply(V(pathway_network)$qvalue, function(x) if(is.null(x)) NA_real_ else as.numeric(x[1]))

  # Use the already-converted numeric vectors
  v_degree <- degree_centrality
  v_betweenness <- betweenness_centrality
  v_eigenvector <- eigenvector_centrality

  hub_pathways <- data.frame(
    pathway = v_pathway,
    description = v_desc,
    category = v_cat,
    qvalue = v_qval,
    degree = v_degree,
    betweenness = v_betweenness,
    eigenvector = v_eigenvector,
    stringsAsFactors = FALSE
  )
  hub_pathways <- hub_pathways[order(-hub_pathways$degree), ]

  write_csv(hub_pathways, file.path(output_dir, "hub_pathways_ranked.csv"))

  log_message(sprintf("  Top 10 hub pathways (by degree centrality):"))
  for (i in 1:min(10, nrow(hub_pathways))) {
    log_message(sprintf("    %d. %s (degree = %d)",
                       i, hub_pathways$description[i], hub_pathways$degree[i]))
  }
} else {
  hub_pathways <- data.frame()
}

################################################################################
# 6. Community Detection
################################################################################

log_message("Performing community detection...")

if (!is.null(pathway_network) && vcount(pathway_network) > 0) {
  # Louvain community detection
  communities <- cluster_louvain(pathway_network, weights = E(pathway_network)$weight)

  V(pathway_network)$community <- membership(communities)

  log_message(sprintf("  Identified %d pathway communities", max(membership(communities))))

  # Summarize communities
  community_summary <- hub_pathways %>%
    mutate(community = V(pathway_network)$community[match(pathway, V(pathway_network)$name)]) %>%
    group_by(community) %>%
    summarise(
      n_pathways = n(),
      avg_degree = mean(degree),
      top_pathway = description[which.max(degree)],
      categories = paste(unique(category), collapse = ", "),
      .groups = "drop"
    ) %>%
    arrange(desc(n_pathways))

  write_csv(community_summary, file.path(output_dir, "pathway_communities.csv"))

  # Export detailed community membership
  community_membership <- hub_pathways %>%
    mutate(community = V(pathway_network)$community[match(pathway, V(pathway_network)$name)]) %>%
    arrange(community, desc(degree))

  write_csv(community_membership, file.path(output_dir, "pathway_community_membership.csv"))
} else {
  communities <- NULL
  community_summary <- data.frame()
}

################################################################################
# 7. Pathway-Pathway Correlation (Expression-Based)
################################################################################

log_message("Calculating expression-based pathway correlations...")

# Load expression data
expr_file <- "SRF_Eva_integrated_analysis/scripts/analysis_3/results/04_category_expression/genes_with_binding_and_expression.csv"

if (file.exists(expr_file)) {
  expr_data <- read_csv(expr_file, show_col_types = FALSE)

  # For each pathway, calculate average expression change of member genes
  # Note: Column name is gene_name (not gene_symbol)
  pathway_scores <- pathway_gene_edges %>%
    left_join(expr_data %>% dplyr::select(gene_name, log2FoldChange),
              by = c("gene" = "gene_name")) %>%
    filter(!is.na(log2FoldChange)) %>%
    group_by(pathway, pathway_name, category) %>%
    summarise(
      n_genes = n(),
      mean_log2fc = mean(log2FoldChange, na.rm = TRUE),
      median_log2fc = median(log2FoldChange, na.rm = TRUE),
      .groups = "drop"
    )

  write_csv(pathway_scores, file.path(output_dir, "pathway_expression_scores.csv"))

  # Calculate correlation between pathway scores
  top_pathway_scores <- pathway_scores %>%
    filter(pathway %in% top_pathway_ids, n_genes >= 5)

  if (nrow(top_pathway_scores) > 1) {
    score_matrix <- top_pathway_scores %>%
      select(pathway, mean_log2fc) %>%
      pivot_wider(names_from = pathway, values_from = mean_log2fc, values_fill = 0) %>%
      as.data.frame()

    # Can't correlate pathways directly, but can cluster by expression pattern
    log_message("  Pathway expression scores calculated")
  }
} else {
  pathway_scores <- data.frame()
  log_message("  Expression data not found, skipping expression-based analysis")
}

################################################################################
# 8. Visualizations
################################################################################

log_message("Creating pathway network visualizations...")

# Plot 1: Pathway similarity heatmap
if (nrow(similarity_matrix) > 0) {
  pdf(file.path(output_dir, "pathway_similarity_heatmap.pdf"), width = 14, height = 14)

  # Cluster pathways
  pathway_clusters <- hclust(as.dist(1 - similarity_matrix), method = "complete")

  # Annotate by category (deduplicate IDs first)
  pathway_annotations <- top_pathways %>%
    dplyr::distinct(ID, .keep_all = TRUE) %>%
    dplyr::select(ID, category) %>%
    as.data.frame()
  rownames(pathway_annotations) <- pathway_annotations$ID
  pathway_annotations$ID <- NULL

  # Handle case where fewer than 3 categories exist
  n_categories <- length(unique(pathway_annotations$category))
  category_colors <- brewer.pal(max(n_categories, 3), "Set2")[1:n_categories]
  names(category_colors) <- unique(pathway_annotations$category)

  row_ha <- rowAnnotation(
    Category = pathway_annotations[rownames(similarity_matrix), "category"],
    col = list(Category = category_colors)
  )

  col_fun <- colorRamp2(c(0, 0.5, 1), c("white", "yellow", "red"))

  ht <- Heatmap(
    similarity_matrix,
    name = "Jaccard\nSimilarity",
    col = col_fun,
    cluster_rows = pathway_clusters,
    cluster_columns = pathway_clusters,
    show_row_names = FALSE,
    show_column_names = FALSE,
    row_split = 5,
    column_split = 5,
    right_annotation = row_ha,
    heatmap_legend_param = list(
      title = "Similarity\n(Jaccard Index)"
    )
  )

  draw(ht, heatmap_legend_side = "right")

  dev.off()
}

# Plot 2: Pathway network visualization
if (!is.null(pathway_network) && vcount(pathway_network) > 0) {
  pdf(file.path(output_dir, "pathway_network_graph.pdf"), width = 16, height = 16)

  # Set layout
  set.seed(42)
  layout <- layout_with_fr(pathway_network)

  # Color by community
  if (!is.null(communities)) {
    community_colors <- setNames(
      brewer.pal(max(3, max(membership(communities))), "Set3"),
      1:max(membership(communities))
    )
    vertex_colors <- community_colors[V(pathway_network)$community]
  } else {
    vertex_colors <- "lightblue"
  }

  # Size by degree
  vertex_sizes <- scales::rescale(V(pathway_network)$degree, to = c(3, 15))

  plot(pathway_network,
       layout = layout,
       vertex.size = vertex_sizes,
       vertex.color = vertex_colors,
       vertex.label = ifelse(V(pathway_network)$degree > quantile(V(pathway_network)$degree, 0.9),
                            V(pathway_network)$Description, NA),
       vertex.label.cex = 0.5,
       vertex.label.dist = 0.5,
       edge.width = E(pathway_network)$weight * 2,
       edge.color = rgb(0.5, 0.5, 0.5, 0.3),
       main = "Pathway Similarity Network\n(nodes = pathways, edges = shared genes)",
       sub = sprintf("%d pathways, %d communities", vcount(pathway_network),
                    ifelse(is.null(communities), 1, max(membership(communities)))))

  # Add legend
  if (!is.null(communities)) {
    legend("bottomleft",
           legend = paste("Community", 1:length(community_colors)),
           fill = community_colors,
           cex = 0.7,
           title = "Pathway Communities")
  }

  dev.off()
}

# Plot 3: Hub pathway bar chart
if (nrow(hub_pathways) > 0) {
  pdf(file.path(output_dir, "hub_pathways_barchart.pdf"), width = 12, height = 8)

  top_hubs <- hub_pathways %>%
    head(20) %>%
    mutate(description = str_trunc(description, 50))

  p1 <- ggplot(top_hubs, aes(x = reorder(description, degree), y = degree)) +
    geom_bar(stat = "identity", aes(fill = category)) +
    coord_flip() +
    scale_fill_brewer(palette = "Set2") +
    labs(
      title = "Top 20 Hub Pathways",
      subtitle = "Ranked by degree centrality (number of connections)",
      x = "Pathway",
      y = "Degree Centrality",
      fill = "Category"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10)
    )

  print(p1)

  dev.off()
}

# Plot 4: Pathway expression scores
if (nrow(pathway_scores) > 0) {
  pdf(file.path(output_dir, "pathway_expression_scores.pdf"), width = 12, height = 8)

  top_expr_pathways <- pathway_scores %>%
    arrange(desc(abs(mean_log2fc))) %>%
    head(30) %>%
    mutate(pathway_name = str_trunc(pathway_name, 50))

  p2 <- ggplot(top_expr_pathways, aes(x = reorder(pathway_name, mean_log2fc), y = mean_log2fc)) +
    geom_bar(stat = "identity", aes(fill = category)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    coord_flip() +
    scale_fill_brewer(palette = "Set2") +
    labs(
      title = "Top 30 Pathways by Expression Change",
      subtitle = "Mean log2 fold change of pathway member genes",
      x = "Pathway",
      y = "Mean log2 Fold Change",
      fill = "Category"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10)
    )

  print(p2)

  dev.off()
}

# Plot 5: Interactive network (if network exists)
if (!is.null(pathway_network) && vcount(pathway_network) > 0) {
  log_message("Creating interactive network visualization...")

  # Prepare nodes
  nodes <- data.frame(
    id = V(pathway_network)$name,
    label = str_trunc(V(pathway_network)$Description, 40),
    title = V(pathway_network)$Description,
    group = ifelse(!is.null(communities), V(pathway_network)$community, 1),
    value = V(pathway_network)$degree,
    category = V(pathway_network)$category
  )

  # Prepare edges
  edges_vis <- data.frame(
    from = as_edgelist(pathway_network)[,1],
    to = as_edgelist(pathway_network)[,2],
    value = E(pathway_network)$weight,
    title = sprintf("Similarity: %.2f", E(pathway_network)$weight)
  )

  # Create interactive network
  vis_net <- visNetwork(nodes, edges_vis, width = "100%", height = "800px") %>%
    visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE),
              selectedBy = "category",
              nodesIdSelection = TRUE) %>%
    visPhysics(stabilization = TRUE) %>%
    visLayout(randomSeed = 42) %>%
    visInteraction(navigationButtons = TRUE, zoomView = TRUE)

  visSave(vis_net, file.path(output_dir, "pathway_network_interactive.html"))

  log_message("  Interactive network saved")
}

################################################################################
# 9. Cross-Category Pathway Overlap Analysis
################################################################################

log_message("Analyzing cross-category pathway overlap...")

# Which pathways are enriched in multiple categories?
pathway_category_counts <- sig_pathways %>%
  group_by(ID, Description) %>%
  summarise(
    n_categories = n_distinct(category),
    categories = paste(unique(category), collapse = "; "),
    min_qvalue = min(qvalue),
    .groups = "drop"
  ) %>%
  arrange(desc(n_categories), min_qvalue)

write_csv(pathway_category_counts, file.path(output_dir, "pathway_category_overlap.csv"))

shared_pathways <- pathway_category_counts %>%
  filter(n_categories >= 2)

log_message(sprintf("  %d pathways enriched in multiple categories", nrow(shared_pathways)))

if (nrow(shared_pathways) > 0) {
  pdf(file.path(output_dir, "shared_pathways_upset.pdf"), width = 12, height = 8)

  # Simple bar chart of multi-category pathways
  top_shared <- shared_pathways %>%
    head(20) %>%
    mutate(Description = str_trunc(Description, 50))

  p3 <- ggplot(top_shared, aes(x = reorder(Description, n_categories), y = n_categories)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    labs(
      title = "Pathways Enriched Across Multiple Categories",
      subtitle = "Top 20 pathways by number of categories",
      x = "Pathway",
      y = "Number of Categories"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10)
    )

  print(p3)

  dev.off()
}

################################################################################
# 10. Generate Summary Report
################################################################################

log_message("Generating summary report...")

summary_file <- file.path(output_dir, "PHASE5_3_SUMMARY.txt")

sink(summary_file)

cat("================================================================================\n")
cat("Phase 5.3: Pathway Crosstalk Analysis - Analysis Summary\n")
cat("TES vs TEAD1 Comparative Study (Excluding TESmut)\n")
cat("================================================================================\n\n")

cat("Analysis Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("INPUT DATA:\n")
cat(sprintf("  - Total GO enrichment results: %d\n", nrow(go_combined)))
cat(sprintf("  - Significant pathways (qvalue < 0.05): %d\n", nrow(sig_pathways)))
cat(sprintf("  - Unique pathways: %d\n", nrow(unique_pathways)))
cat(sprintf("  - Unique genes in pathways: %d\n", nrow(unique_genes)))
cat("\n")

cat("PATHWAY NETWORK:\n")
if (!is.null(pathway_network)) {
  cat(sprintf("  - Network nodes (pathways): %d\n", vcount(pathway_network)))
  cat(sprintf("  - Network edges (similarities): %d\n", ecount(pathway_network)))
  cat(sprintf("  - Network density: %.4f\n", edge_density(pathway_network)))
  cat(sprintf("  - Average degree: %.2f\n", mean(degree(pathway_network))))
  if (!is.null(communities)) {
    cat(sprintf("  - Number of communities: %d\n", max(membership(communities))))
  }
} else {
  cat("  - Network could not be constructed (insufficient similarity)\n")
}
cat("\n")

cat("HUB PATHWAYS:\n")
if (nrow(hub_pathways) > 0) {
  cat("  Top 10 hub pathways (by degree centrality):\n")
  for (i in 1:min(10, nrow(hub_pathways))) {
    cat(sprintf("    %d. %s\n", i, hub_pathways$description[i]))
    cat(sprintf("       Category: %s, Degree: %d, Betweenness: %.1f\n",
                hub_pathways$category[i], hub_pathways$degree[i],
                hub_pathways$betweenness[i]))
  }
} else {
  cat("  No hub pathways identified\n")
}
cat("\n")

cat("PATHWAY COMMUNITIES:\n")
if (nrow(community_summary) > 0) {
  for (i in 1:nrow(community_summary)) {
    cat(sprintf("  Community %d:\n", community_summary$community[i]))
    cat(sprintf("    - Pathways: %d\n", community_summary$n_pathways[i]))
    cat(sprintf("    - Avg degree: %.1f\n", community_summary$avg_degree[i]))
    cat(sprintf("    - Top pathway: %s\n", community_summary$top_pathway[i]))
    cat(sprintf("    - Categories: %s\n", community_summary$categories[i]))
  }
} else {
  cat("  No communities detected\n")
}
cat("\n")

cat("CROSS-CATEGORY OVERLAP:\n")
cat(sprintf("  - Pathways in multiple categories: %d\n", nrow(shared_pathways)))
if (nrow(shared_pathways) > 0) {
  cat("  Top shared pathways:\n")
  for (i in 1:min(5, nrow(shared_pathways))) {
    cat(sprintf("    %d. %s (%d categories: %s)\n",
                i, shared_pathways$Description[i],
                shared_pathways$n_categories[i],
                shared_pathways$categories[i]))
  }
}
cat("\n")

cat("KEY FINDINGS:\n")
cat("  1. Hub pathways represent core biological processes affected by TES/TEAD1\n")
cat("  2. Pathway communities reveal functional modules\n")
cat("  3. Cross-category pathways indicate shared regulatory mechanisms\n")
cat("  4. High similarity pathways may represent redundant or hierarchical processes\n")
cat("\n")

cat("OUTPUT FILES:\n")
cat("  1. hub_pathways_ranked.csv - All pathways ranked by centrality\n")
cat("  2. pathway_communities.csv - Community summary\n")
cat("  3. pathway_community_membership.csv - Detailed community assignments\n")
cat("  4. pathway_expression_scores.csv - Expression-based pathway scores\n")
cat("  5. pathway_category_overlap.csv - Cross-category pathway analysis\n")
cat("  6. pathway_similarity_heatmap.pdf - Clustered similarity matrix\n")
cat("  7. pathway_network_graph.pdf - Network visualization\n")
cat("  8. hub_pathways_barchart.pdf - Top hub pathways\n")
cat("  9. pathway_expression_scores.pdf - Expression-based ranking\n")
cat(" 10. pathway_network_interactive.html - Interactive network (if generated)\n")
cat(" 11. shared_pathways_upset.pdf - Multi-category pathway overlap\n")
cat("\n")

cat("BIOLOGICAL INTERPRETATION:\n")
cat("  - Hub pathways: Central processes coordinating multiple functions\n")
cat("  - Communities: Functional modules operating in coordination\n")
cat("  - High similarity: Overlapping gene sets, possible redundancy\n")
cat("  - Cross-category pathways: Core processes affected by both TES and TEAD1\n")
cat("\n")

cat("NEXT STEPS:\n")
cat("  - Validate hub pathways with orthogonal datasets\n")
cat("  - Test functional dependency of hub pathway genes\n")
cat("  - Examine temporal dynamics of pathway activation\n")
cat("  - Integrate with protein-protein interaction networks\n")
cat("\n")

cat("================================================================================\n")
cat("Analysis Complete\n")
cat("================================================================================\n")

sink()

log_message("Phase 5.3 Analysis Complete!")
log_message(sprintf("Results saved to: %s", output_dir))
log_message(sprintf("Summary report: %s", summary_file))

cat("\n")
cat("================================================================================\n")
cat("Phase 5.3: Pathway Crosstalk Analysis - COMPLETE\n")
cat("================================================================================\n")
