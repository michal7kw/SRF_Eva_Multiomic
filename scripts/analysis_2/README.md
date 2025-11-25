# Multiomic Analysis of TES and TEAD1 Regulation

This pipeline integrates Cut&Tag (TES, TEAD1), RNA-seq (TES vs GFP), and meDIP-seq (TES vs GFP) data to dissect the regulatory logic of TES and TEAD1.

## Overview

The analysis classifies binding sites into Shared (TES+TEAD1) and Unique categories, then characterizes them by:
1.  **Binding Logic**: Overlap and differential binding.
2.  **Motif Enrichment**: Identifying co-factors.
3.  **Epigenetic State**: DNA methylation levels (meDIP).
4.  **Transcriptional Impact**: Correlation with gene expression changes.
5.  **Functional Pathways**: GO/KEGG enrichment.

## Pipeline Steps

### Module 1: Refined Peak Classification (`01_peak_classification.R`)
-   **Input**: Consensus peaks (`TES_consensus_peaks.bed`, `TEAD1_consensus_peaks.bed`)
-   **Output**:
    -   `Shared_TES_TEAD1.bed`: Intersected peaks
    -   `TES_Unique.bed`: TES peaks without TEAD1 overlap
    -   `TEAD1_Unique.bed`: TEAD1 peaks without TES overlap
    -   `Master_Peak_Annotation.csv`: Detailed annotation with gene symbols
    -   **Plots**:
        -   `Peak_Overlap_Venn.pdf`: Counts of shared vs unique peaks
        -   `Peak_Annotation_Distribution.pdf`: Genomic distribution (Promoter, Intron, etc.)
        -   `Distance_to_TSS.pdf`: Binding distance relative to TSS

### Module 2: Motif Enrichment (`02_motif_analysis.sh`)
-   **Input**: BED files from Module 1
-   **Output**: HOMER motif enrichment results in `results/02_motif_enrichment/{Category}/`
-   **Method**: `findMotifsGenome.pl` looking for known and de novo motifs.

### Module 3: meDIP Integration (`03_medip_integration.R`)
-   **Input**: Peak categories, meDIP BigWigs, DMRs
-   **Output**:
    -   `Methylation_at_Binding_Sites.csv`: Mean methylation signal per peak
    -   `DMR_Overlap_Stats.csv`: Overlap counts with DMRs
    -   **Plots**:
        -   `Methylation_Boxplot.pdf`: Methylation levels by binding category
        -   `DMR_Overlap_Barplot.pdf`: % of peaks overlapping DMRs

### Module 4: Gene Regulatory Logic (`04_gene_regulatory_logic.R`)
-   **Input**: Annotated peaks, RNA-seq DEGs (`deseq2_results_TES_vs_GFP.txt`)
-   **Output**:
    -   `Peaks_Linked_to_Expression.csv`: Merged binding and expression data
    -   **Plots**:
        -   `Expression_by_Binding_Category.pdf`: Log2FC of target genes by category
        -   `Expression_Promoter_vs_Enhancer.pdf`: Log2FC split by regulatory context

### Module 5: Functional Enrichment (`05_functional_enrichment.R`)
-   **Input**: Annotated peaks
-   **Output**:
    -   `GO_BP_{Category}.csv`: GO Biological Process enrichment tables
    -   **Plots**:
        -   `Dotplot_{Category}.pdf`: Enrichment dotplots per category
        -   `Comparison_Dotplot.pdf`: Comparative enrichment across categories

## Generated Outputs

### Peak Classification
![Peak Overlap](results/01_peak_classification/Peak_Overlap_Venn.png)
*Note: PDF output converted to PNG for display*

### Methylation Analysis
![Methylation Levels](results/03_medip_integration/Methylation_Boxplot.png)
*DNA methylation levels at TES binding sites, stratified by TEAD1 co-binding.*

### Gene Expression Impact
![Expression Changes](results/04_gene_regulatory_logic/Expression_by_Binding_Category.png)
*Transcriptional changes (Log2FC) of genes bound by TES, TEAD1, or both.*

### Functional Pathways
![Pathway Enrichment](results/05_functional_enrichment/Comparison_Dotplot.png)
*Comparative functional enrichment of target genes.*

## Usage

Run the master script to execute all modules (except Motif Analysis which requires manual HOMER setup):

```bash
sbatch 00_run_all_new_modules.sh
```
