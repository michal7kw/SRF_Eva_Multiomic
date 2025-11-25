# Integrative Analysis Phase 1: GSEA and Candidate Genes

This directory contains scripts for the first phase of integrative analysis, focusing on Gene Set Enrichment Analysis (GSEA), candidate gene profiling, and promoter binding visualization.

## Overview

The analysis is divided into several key modules:
1.  **TRUE GSEA**: Rank-based enrichment analysis using MSigDB collections.
2.  **Candidate Gene Analysis**: Heatmaps and profiles for specific splicing-related genes.
3.  **Promoter Heatmaps**: Visualizing TES and TEAD1 binding at promoter regions.
4.  **Integrative Plots**: Combining RNA-seq and Cut&Tag data.

## Pipeline Steps

### 1. TRUE GSEA Analysis (`true_gsea_analysis.R`)
-   **Description**: Performs GSEA using `fgsea` on RNA-seq ranked lists against MSigDB collections (Hallmark, C2, C5).
-   **Output Directory**: `output/true_gsea_analysis/`
-   **Key Outputs**:
    -   `fgsea_all_pathways.csv`: Complete GSEA results.
    -   `fgsea_significant_pathways.csv`: Significant pathways (padj < 0.05).
    -   `01_top_pathways_barplot.pdf`: Barplot of top enriched pathways.
    -   `02_detailed_enrichment_plots.pdf`: Enrichment plots for top pathways.
    -   `03_cancer_pathways_bubble.pdf`: Bubble plot of cancer-related pathways.

### 2. Candidate Gene Heatmaps (`candidate_gene_heatmap.R`)
-   **Description**: Generates expression and binding heatmaps for a defined list of candidate genes (e.g., splicing factors).
-   **Output Directory**: `output/candidate_gene_heatmap/`
-   **Key Outputs**:
    -   `candidate_genes_with_peaks.csv`: Binding status of candidate genes.
    -   `01_binding_comparison.pdf`: Comparison of binding at candidate loci.
    -   `02_peak_count_distribution.pdf`: Distribution of peak counts.
    -   `03_expression_vs_binding.pdf`: Correlation plot.

### 3. Promoter Heatmaps (`binding_heatmap.R`)
-   **Description**: Creates global signal heatmaps centered on TSS for TES and TEAD1.
-   **Output Directory**: `output/binding_heatmap/`
-   **Key Outputs**:
    -   `promoter_binding_heatmap.pdf`: Signal intensity at promoters.
    -   `promoter_binding_heatmap_with_fc.pdf`: Sorted by fold-change.

### 4. Integrative Analysis (`final_integrative_analysis.R`)
-   **Description**: Integrates differential expression with binding data to identify direct targets.
-   **Output Directory**: `output/only_degs/results/` (for DEGs) or `output/results/` (for all genes)
-   **Key Outputs**:
    -   `direct_targets/TES_direct_targets.csv`: Genes bound and regulated by TES.
    -   `plots/01_gene_classification_pie.pdf`: Proportion of direct vs indirect targets.
    -   `plots/04_target_overlap_venn.pdf`: Overlap of TES and TEAD1 targets.

## Generated Outputs

### GSEA Results
![Top Pathways](output/true_gsea_analysis/01_top_pathways_barplot.png)
*Top enriched pathways from GSEA analysis.*

![Enrichment Plots](output/true_gsea_analysis/02_detailed_enrichment_plots.png)
*Detailed enrichment plots for key pathways.*

### Candidate Gene Analysis
![Expression vs Binding](output/candidate_gene_heatmap/03_expression_vs_binding.png)
*Relationship between TES binding and gene expression for candidate genes.*

### Integrative Plots
![Target Overlap](output/only_degs/plots/04_target_overlap_venn.png)
*Overlap of TES and TEAD1 direct targets.*

![Gene Classification](output/only_degs/plots/01_gene_classification_pie.png)
*Classification of DEGs into direct and indirect targets.*

## Usage

To run the complete analysis suite:

```bash
sbatch 0_run_all_analyses.sh
```

Or run individual modules:

```bash
# GSEA
sbatch true_gsea_analysis.sh

# Candidate Genes
sbatch candidate_gene_heatmap.sh

# Integrative Analysis
sbatch final_integrative_analysis.sh
```
