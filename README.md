### Direct Target Identification

- **`integrative_analysis/results/direct_targets/TES_direct_targets.csv`**: This file lists the direct target genes of TES, which are genes that are both bound by TES and are differentially expressed.
- **`integrative_analysis/results/direct_targets/TEAD1_direct_targets.csv`**: Similar to the above, this file lists the direct target genes of TEAD1.
- **`integrative_analysis/results/direct_targets/TES_specific_targets.csv`**: This file contains genes that are direct targets of TES but not TEAD1.
- **`integrative_analysis/results/direct_targets/TEAD1_specific_targets.csv`**: This file contains genes that are direct targets of TEAD1 but not TES.

### Gene Classification

- **`integrative_analysis/results/regulatory_networks/gene_classification.csv`**: This file provides a classification for each gene based on its binding and expression status. The classes include:
    - `TES_direct`: Bound by TES and differentially expressed.
    - `TEAD1_direct`: Bound by TEAD1 and differentially expressed.
    - `TES_TEAD1_shared`: Bound by both and differentially expressed.
    - `Indirect`: Differentially expressed but not bound by either.
    - `Unbound`: Not differentially expressed and not bound.

### Pathway Enrichment Analysis

- **`integrative_analysis/results/pathway_analysis/*_GO_enrichment.csv`**: These files contain the results of Gene Ontology (GO) enrichment analysis for different gene sets (e.g., `TES_direct`, `TEAD1_direct`). They identify biological processes that are over-represented in the list of target genes, providing insights into the functional roles of TES and TEAD1.

### Visualization

- **`integrative_analysis/plots/*.pdf`**: This directory contains several plots that visualize the results of the analysis:
    - `01_gene_classification_pie.pdf`: A pie chart showing the proportion of genes in each regulatory class.
    - `02_expression_histograms.pdf`: Histograms showing the distribution of expression changes for each regulatory class.
    - `03_expression_boxplots.pdf`: Boxplots comparing the expression changes between different regulatory classes.
    - `04_target_overlap_venn.pdf`: A Venn diagram showing the overlap between TES and TEAD1 direct targets and all differentially expressed genes.
    - `05_summary_barchart.pdf`: A bar chart summarizing the number of genes in each category.
    - `*_pathway_enrichment.pdf`: Dot plots visualizing the results of the GO enrichment analysis.

### Summary Report

- **`integrative_analysis/results/FINAL_ANALYSIS_SUMMARY.txt`**: A text file summarizing the key statistics from the analysis, such as the number of peaks, bound genes, direct targets, and enriched pathways.