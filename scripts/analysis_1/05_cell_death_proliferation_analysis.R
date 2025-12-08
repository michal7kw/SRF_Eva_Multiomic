#!/usr/bin/env Rscript

# ===============================================================================
# SCRIPT: cell_death_proliferation_analysis.R
# PURPOSE: Focused analysis of cell death and proliferation pathways
#
# DESCRIPTION:
# This script performs targeted enrichment analysis focusing on cell death,
# apoptosis, and proliferation pathways. This is particularly relevant for
# TES/TEAD1 research in glioblastoma cells, as these transcription factors
# are known to regulate cell survival and proliferation pathways.
#
# KEY ANALYSES:
# 1. GO enrichment focused on cell death/proliferation terms
# 2. KEGG pathway analysis for cancer-related pathways
# 3. Reactome pathway analysis
# 4. Custom gene set analysis for oncology-relevant pathways
# 5. Comparative analysis between TES and TEAD1
# ===============================================================================

library(clusterProfiler)
library(org.Hs.eg.db)
library(KEGGREST)
library(ReactomePA)
library(ggplot2)
library(dplyr)
library(VennDiagram)
library(ComplexHeatmap)
library(RColorBrewer)

# Fix namespace conflicts: clusterProfiler masks dplyr::filter and stats::filter
# Explicitly use dplyr functions
filter <- dplyr::filter
select <- dplyr::select
mutate <- dplyr::mutate
arrange <- dplyr::arrange

# Set working directory
setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_1")

# Create output directory
dir.create("output/05_cell_death_proliferation", recursive = TRUE, showWarnings = FALSE)

cat("=== CELL DEATH AND PROLIFERATION PATHWAY ANALYSIS ===\n")
cat("Focus: TES/TEAD1 regulation of cell survival in glioblastoma\n\n")

# Load annotated peak data (from CUTandTAG pipeline)
annotation_dir <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_CUTandTAG/results/07_analysis_narrow"
tes_file <- file.path(annotation_dir, "TES_peaks_annotated.csv")
tead1_file <- file.path(annotation_dir, "TEAD1_peaks_annotated.csv")

# Function to extract gene IDs from annotation files
extract_gene_ids <- function(file_path, sample_name) {
    if (!file.exists(file_path)) {
        cat("Warning: File not found:", file_path, "\n")
        return(character(0))
    }

    data <- read.csv(file_path, stringsAsFactors = FALSE)

    if ("geneId" %in% colnames(data)) {
        gene_ids <- unique(data$geneId[!is.na(data$geneId)])
        cat("Loaded", length(gene_ids), "genes for", sample_name, "\n")
        return(gene_ids)
    } else {
        cat("Warning: No geneId column found in", file_path, "\n")
        return(character(0))
    }
}

# Extract gene IDs for each condition
tes_genes <- extract_gene_ids(tes_file, "TES")
tead1_genes <- extract_gene_ids(tead1_file, "TEAD1")

# Create gene lists for analysis
gene_lists <- list(
    TES = tes_genes,
    TEAD1 = tead1_genes
)

# Remove empty lists
gene_lists <- gene_lists[sapply(gene_lists, length) > 0]

if (length(gene_lists) == 0) {
    stop("No gene lists available for analysis")
}

cat("Gene list summary:\n")
for (name in names(gene_lists)) {
    cat(sprintf("  %s: %d genes\n", name, length(gene_lists[[name]])))
}
cat("\n")

# ===============================================================================
# 1. FOCUSED GO ANALYSIS - CELL DEATH AND PROLIFERATION
# ===============================================================================

cat("1. GO ENRICHMENT ANALYSIS - CELL DEATH & PROLIFERATION\n")
cat("======================================================\n")

# Define cell death and proliferation related GO terms
cell_death_terms <- c(
    "GO:0008219", # cell death
    "GO:0012501", # programmed cell death
    "GO:0006915", # apoptotic process
    "GO:0043065", # positive regulation of apoptotic process
    "GO:0043066", # negative regulation of apoptotic process
    "GO:0097190", # apoptotic signaling pathway
    "GO:0006917", # induction of apoptosis
    "GO:0043523", # regulation of neuron apoptotic process
    "GO:0070059", # intrinsic apoptotic signaling pathway
    "GO:0097191", # extrinsic apoptotic signaling pathway
    "GO:0001844", # autophagy of mitochondrion
    "GO:0016264" # gap junction-mediated intercellular transport
)

proliferation_terms <- c(
    "GO:0008283", # cell proliferation
    "GO:0042127", # regulation of cell proliferation
    "GO:0008284", # positive regulation of cell proliferation
    "GO:0008285", # negative regulation of cell proliferation
    "GO:0051301", # cell division
    "GO:0000278", # mitotic cell cycle
    "GO:0007049", # cell cycle
    "GO:0045787", # positive regulation of cell cycle
    "GO:0045786", # negative regulation of cell cycle
    "GO:0051726", # regulation of cell cycle
    "GO:0006270", # DNA replication initiation
    "GO:0000082", # G1/S transition of mitotic cell cycle
    "GO:0000086", # G2/M transition of mitotic cell cycle
    "GO:0044772", # mitotic cell cycle phase transition
    "GO:0007067" # mitosis
)

# Function to perform focused GO analysis
perform_focused_go <- function(gene_ids, sample_name, focus_terms, focus_name) {
    if (length(gene_ids) < 10) {
        cat("Skipping", sample_name, focus_name, "- insufficient genes\n")
        return(NULL)
    }

    cat("Analyzing", focus_name, "for", sample_name, "...\n")

    # Standard GO enrichment
    ego <- enrichGO(
        gene = gene_ids,
        OrgDb = org.Hs.eg.db,
        keyType = "ENTREZID",
        ont = "BP",
        pAdjustMethod = "BH",
        qvalueCutoff = 0.05,
        readable = TRUE
    )

    if (is.null(ego) || nrow(ego@result) == 0) {
        cat("No significant GO terms found for", sample_name, focus_name, "\n")
        return(NULL)
    }

    # Filter for terms of interest
    focused_results <- ego@result[ego@result$ID %in% focus_terms, ]

    # Also include terms containing key words
    keywords <- if (focus_name == "cell_death") {
        c("death", "apoptosis", "apoptotic", "necrosis", "autophagy")
    } else {
        c("proliferation", "division", "cycle", "mitosis", "growth")
    }

    keyword_results <- ego@result[grepl(paste(keywords, collapse = "|"),
        ego@result$Description,
        ignore.case = TRUE
    ), ]

    # Combine results
    combined_results <- rbind(focused_results, keyword_results)
    combined_results <- combined_results[!duplicated(combined_results$ID), ]

    if (nrow(combined_results) > 0) {
        # Create plots
        plot_title <- paste(sample_name, focus_name, "GO Enrichment")

        # Dot plot
        p1 <- ggplot(head(combined_results, 15), aes(x = Count, y = reorder(Description, Count))) +
            geom_point(aes(color = p.adjust, size = Count)) +
            scale_color_gradient(low = "red", high = "blue", name = "p.adjust") +
            labs(title = plot_title, x = "Gene Count", y = "GO Terms") +
            theme_minimal() +
            theme(axis.text.y = element_text(size = 8))

        # Save plots
        filename_base <- paste0("output/05_cell_death_proliferation/", sample_name, "_", focus_name, "_GO")

        ggsave(paste0(filename_base, ".pdf"), p1, width = 12, height = 8)
        ggsave(paste0(filename_base, ".png"), p1, width = 12, height = 8, dpi = 300)

        # Save results
        write.csv(combined_results, paste0(filename_base, "_results.csv"), row.names = FALSE)

        return(combined_results)
    }

    return(NULL)
}

# Perform focused GO analysis for each sample
go_results <- list()
for (sample_name in names(gene_lists)) {
    genes <- gene_lists[[sample_name]]

    # Cell death analysis
    go_results[[paste0(sample_name, "_death")]] <-
        perform_focused_go(genes, sample_name, cell_death_terms, "cell_death")

    # Proliferation analysis
    go_results[[paste0(sample_name, "_prolif")]] <-
        perform_focused_go(genes, sample_name, proliferation_terms, "proliferation")
}

# ===============================================================================
# 2. KEGG PATHWAY ANALYSIS - CANCER PATHWAYS
# ===============================================================================

cat("\n2. KEGG PATHWAY ANALYSIS - CANCER RELATED\n")
cat("=========================================\n")

# Function to perform KEGG analysis
perform_kegg_analysis <- function(gene_ids, sample_name) {
    if (length(gene_ids) < 10) {
        cat("Skipping KEGG for", sample_name, "- insufficient genes\n")
        return(NULL)
    }

    cat("Running KEGG analysis for", sample_name, "...\n")

    kegg_result <- enrichKEGG(
        gene = gene_ids,
        organism = "hsa",
        keyType = "kegg",
        pAdjustMethod = "BH",
        qvalueCutoff = 0.05
    )

    if (is.null(kegg_result) || nrow(kegg_result@result) == 0) {
        cat("No significant KEGG pathways found for", sample_name, "\n")
        return(NULL)
    }

    # Focus on cancer-related pathways
    cancer_keywords <- c(
        "cancer", "carcinoma", "apoptosis", "cell cycle",
        "proliferation", "tumor", "oncogene", "suppressor",
        "DNA repair", "p53", "PI3K", "MAPK", "Wnt", "Hippo"
    )

    cancer_pathways <- kegg_result@result[
        grepl(paste(cancer_keywords, collapse = "|"),
            kegg_result@result$Description,
            ignore.case = TRUE
        ),
    ]

    if (nrow(cancer_pathways) > 0) {
        # Create visualization
        p1 <- ggplot(
            head(cancer_pathways, 10),
            aes(x = Count, y = reorder(Description, Count))
        ) +
            geom_point(aes(color = p.adjust, size = Count)) +
            scale_color_gradient(low = "red", high = "blue", name = "p.adjust") +
            labs(
                title = paste(sample_name, "Cancer-Related KEGG Pathways"),
                x = "Gene Count", y = "KEGG Pathways"
            ) +
            theme_minimal() +
            theme(axis.text.y = element_text(size = 8))

        # Save results
        filename_base <- paste0("output/05_cell_death_proliferation/", sample_name, "_KEGG_cancer")

        ggsave(paste0(filename_base, ".pdf"), p1, width = 14, height = 8)
        ggsave(paste0(filename_base, ".png"), p1, width = 14, height = 8, dpi = 300)
        write.csv(cancer_pathways, paste0(filename_base, "_results.csv"), row.names = FALSE)

        return(cancer_pathways)
    }

    return(NULL)
}

# Perform KEGG analysis for each sample
kegg_results <- list()
for (sample_name in names(gene_lists)) {
    kegg_results[[sample_name]] <- perform_kegg_analysis(gene_lists[[sample_name]], sample_name)
}

# ===============================================================================
# 3. REACTOME PATHWAY ANALYSIS
# ===============================================================================

cat("\n3. REACTOME PATHWAY ANALYSIS\n")
cat("============================\n")

# Function to perform Reactome analysis
perform_reactome_analysis <- function(gene_ids, sample_name) {
    if (length(gene_ids) < 10) {
        cat("Skipping Reactome for", sample_name, "- insufficient genes\n")
        return(NULL)
    }

    cat("Running Reactome analysis for", sample_name, "...\n")

    # Try to run Reactome analysis with error handling
    tryCatch(
        {
            reactome_result <- enrichPathway(
                gene = gene_ids,
                organism = "human",
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05,
                readable = TRUE
            )

            if (is.null(reactome_result) || nrow(reactome_result@result) == 0) {
                cat("No significant Reactome pathways found for", sample_name, "\n")
                return(NULL)
            }

            # Focus on cell death and proliferation pathways
            focus_keywords <- c(
                "apoptosis", "death", "proliferation", "cycle",
                "division", "DNA repair", "checkpoint", "tumor",
                "oncogene", "suppressor", "survival"
            )

            focused_pathways <- reactome_result@result[
                grepl(paste(focus_keywords, collapse = "|"),
                    reactome_result@result$Description,
                    ignore.case = TRUE
                ),
            ]

            if (nrow(focused_pathways) > 0) {
                # Create visualization
                p1 <- ggplot(
                    head(focused_pathways, 10),
                    aes(x = Count, y = reorder(Description, Count))
                ) +
                    geom_point(aes(color = p.adjust, size = Count)) +
                    scale_color_gradient(low = "red", high = "blue", name = "p.adjust") +
                    labs(
                        title = paste(sample_name, "Cell Death/Proliferation Reactome Pathways"),
                        x = "Gene Count", y = "Reactome Pathways"
                    ) +
                    theme_minimal() +
                    theme(axis.text.y = element_text(size = 7))

                # Save results
                filename_base <- paste0("output/05_cell_death_proliferation/", sample_name, "_Reactome")

                ggsave(paste0(filename_base, ".pdf"), p1, width = 16, height = 8)
                ggsave(paste0(filename_base, ".png"), p1, width = 16, height = 8, dpi = 300)
                write.csv(focused_pathways, paste0(filename_base, "_results.csv"), row.names = FALSE)

                return(focused_pathways)
            }

            return(NULL)
        },
        error = function(e) {
            cat("Error in Reactome analysis for", sample_name, ":", e$message, "\n")
            return(NULL)
        }
    )
}

# Perform Reactome analysis for each sample
reactome_results <- list()
for (sample_name in names(gene_lists)) {
    reactome_results[[sample_name]] <- perform_reactome_analysis(gene_lists[[sample_name]], sample_name)
}

# ===============================================================================
# 4. COMPARATIVE ANALYSIS
# ===============================================================================

cat("\n4. COMPARATIVE ANALYSIS\n")
cat("=======================\n")

# Compare gene sets between conditions
if (length(gene_lists) >= 2) {
    # Create Venn diagram for gene overlap
    if (length(gene_lists) == 2) {
        venn_data <- gene_lists
        venn.diagram(venn_data,
            filename = "output/05_cell_death_proliferation/gene_overlap_venn.png",
            category.names = names(venn_data),
            output = TRUE,
            imagetype = "png",
            height = 2000, width = 2000, resolution = 300,
            fill = c("lightblue", "lightcoral")
        )
    } else if (length(gene_lists) == 3) {
        venn_data <- gene_lists
        venn.diagram(venn_data,
            filename = "output/05_cell_death_proliferation/gene_overlap_venn.png",
            category.names = names(venn_data),
            output = TRUE,
            imagetype = "png",
            height = 2000, width = 2000, resolution = 300,
            fill = c("lightblue", "lightcoral", "lightgreen")
        )
    }

    # Calculate overlap statistics
    overlap_stats <- data.frame(
        comparison = character(),
        overlap_count = numeric(),
        jaccard_index = numeric(),
        stringsAsFactors = FALSE
    )

    sample_names <- names(gene_lists)
    for (i in 1:(length(sample_names) - 1)) {
        for (j in (i + 1):length(sample_names)) {
            name1 <- sample_names[i]
            name2 <- sample_names[j]
            genes1 <- gene_lists[[name1]]
            genes2 <- gene_lists[[name2]]

            overlap <- length(intersect(genes1, genes2))
            union_size <- length(union(genes1, genes2))
            jaccard <- overlap / union_size

            overlap_stats <- rbind(overlap_stats, data.frame(
                comparison = paste(name1, "vs", name2),
                overlap_count = overlap,
                jaccard_index = jaccard,
                stringsAsFactors = FALSE
            ))
        }
    }

    write.csv(overlap_stats, "output/05_cell_death_proliferation/gene_overlap_stats.csv", row.names = FALSE)
    cat("Gene overlap analysis completed\n")
}

# ===============================================================================
# 5. SUMMARY REPORT
# ===============================================================================

cat("\n5. GENERATING SUMMARY REPORT\n")
cat("============================\n")

# Create comprehensive summary
summary_report <- list(
    analysis_date = Sys.Date(),
    samples_analyzed = names(gene_lists),
    gene_counts = sapply(gene_lists, length)
)

# Count significant results
go_counts <- sapply(go_results, function(x) if (!is.null(x)) nrow(x) else 0)
kegg_counts <- sapply(kegg_results, function(x) if (!is.null(x)) nrow(x) else 0)
reactome_counts <- sapply(reactome_results, function(x) if (!is.null(x)) nrow(x) else 0)

# Create summary table
summary_table <- data.frame(
    Sample = rep(names(gene_lists), each = 3),
    Analysis = rep(c("GO_cell_death", "GO_proliferation", "KEGG"), length(gene_lists)),
    Significant_Terms = 0,
    stringsAsFactors = FALSE
)

# Fill in counts
for (sample in names(gene_lists)) {
    death_key <- paste0(sample, "_death")
    prolif_key <- paste0(sample, "_prolif")

    summary_table[summary_table$Sample == sample & summary_table$Analysis == "GO_cell_death", "Significant_Terms"] <-
        if (death_key %in% names(go_counts)) go_counts[[death_key]] else 0

    summary_table[summary_table$Sample == sample & summary_table$Analysis == "GO_proliferation", "Significant_Terms"] <-
        if (prolif_key %in% names(go_counts)) go_counts[[prolif_key]] else 0

    summary_table[summary_table$Sample == sample & summary_table$Analysis == "KEGG", "Significant_Terms"] <-
        if (sample %in% names(kegg_counts)) kegg_counts[[sample]] else 0
}

write.csv(summary_table, "output/05_cell_death_proliferation/analysis_summary.csv", row.names = FALSE)

# Print summary
cat("\n=== ANALYSIS SUMMARY ===\n")
print(summary_table)

cat("\n=== OUTPUT FILES ===\n")
cat("Results saved in: output/05_cell_death_proliferation/\n")
cat("- GO enrichment plots and results (PDF/PNG/CSV)\n")
cat("- KEGG pathway analysis (PDF/PNG/CSV)\n")
cat("- Reactome pathway analysis (PDF/PNG/CSV)\n")
cat("- Gene overlap analysis and Venn diagrams\n")
cat("- Summary statistics and reports\n")

cat("\nCell death and proliferation analysis complete!\n")
