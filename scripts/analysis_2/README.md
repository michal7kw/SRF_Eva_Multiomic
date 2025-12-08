# Analysis 2: Multiomic Regulation by TES and TEAD1
**SRF_Eva_integrated_analysis**

## Overview
This directory contains the integrated analysis pipeline for dissecting the regulatory logic of TES and TEAD1. It combines Cut&Tag (binding), RNA-seq (expression), and meDIP-seq (methylation) data to classify binding sites, explore co-factor motifs, and test mechanistic hypotheses about DNA methylation changes.

## Master Script
**`00_run_all_new_modules.sh`**
This is the main execution script. It manages the entire pipeline using SLURM job dependencies.
*   **Usage:** `sbatch 00_run_all_new_modules.sh`
*   **Function:** Submits jobs for all 6 analysis modules (Module 2 requires manual HOMER setup and is skipped by default).

## Pipeline Modules

### Module 1: Peak Classification
*   **Script:** `01_peak_classification.R`
*   **Description:** Generates consensus peaks from replicates (requiring 2/3 overlap) and classifies them into Shared (TES+TEAD1), TES-Unique, and TEAD1-Unique categories.
*   **Output:** `results/01_peak_classification/` (BED files, Venn diagrams, annotation plots).

### Module 2: Motif Analysis
*   **Script:** `02_motif_analysis.sh`
*   **Description:** Performs de novo and known motif enrichment analysis using HOMER on the classified peak sets.
*   **Output:** `results/02_motif_enrichment/` (HOMER HTML reports).

### Module 3: meDIP Integration
*   **Script:** `03_medip_integration.R`
*   **Description:** Quantifies DNA methylation levels at binding sites and correlates them with Differentially Methylated Regions (DMRs).
*   **Output:** `results/03_medip_integration/` (Methylation boxplots, DMR overlap stats).

### Module 4: Gene Regulatory Logic
*   **Script:** `04_gene_regulatory_logic.R`
*   **Description:** Links binding sites to target genes and integrates RNA-seq data to determine the transcriptional effect (Activation vs Repression) of each binding category. Distinguishes between Promoter and Enhancer regulation.
*   **Output:** `results/04_gene_regulatory_logic/` (Expression boxplots, regulatory category analysis).

### Module 5: Functional Enrichment
*   **Script:** `05_functional_enrichment.R`
*   **Description:** Performs GO Biological Process enrichment analysis for genes associated with each binding category.
*   **Output:** `results/05_functional_enrichment/` (Dotplots, enrichment tables).

### Module 6: Methylation at Regulated Genes
*   **Script:** `06_methylation_at_regulated_genes.R`
*   **Description:** Tests the "Indirect Methylation Hypothesis". Instead of looking at binding sites, it examines whether TES-mediated repression is associated with broad methylation changes across the *gene bodies* of downregulated genes (potentially via DNMT3A/3L recruitment).
*   **Output:** `results/06_methylation_at_regulated_genes/` (Gene body methylation profiles, correlation plots).

## Directory Structure
```
analysis_2/
├── 00_run_all_new_modules.sh         # MASTER SCRIPT
├── *.R                               # R analysis scripts (Modules 1, 3, 4, 5, 6)
├── *.sh                              # SLURM wrapper scripts
├── results/                          # Main output directory
│   ├── 01_peak_classification/
│   ├── ...
│   └── 06_methylation_at_regulated_genes/
├── README.md
└── RESULTS.md                        # Summary of key findings
```

## Bioinformatics Notes (Review Findings)
A comprehensive review of the scripts in this directory was conducted to ensure logic, consistency, and adherence to best practices.

1.  **Robust Peak Calling:** Module 1 uses a replicate-based consensus approach (requiring peaks to be present in at least 2 of 3 replicates) rather than simple pooling. This significantly reduces false positives and ensures that the "Shared" vs "Unique" classification is biologically meaningful.
2.  **Integrative Logic:** The pipeline correctly links data across three modalities (Cut&Tag, RNA-seq, meDIP). ID conversion is handled consistently using `org.Hs.eg.db`.
3.  **Mechanistic Hypotheses:**
    *   **Direct Regulation:** Module 4 tests for direct transcriptional effects at binding sites.
    *   **Indirect/Spreading Methylation:** Module 6 adds a critical layer by testing for methylation changes at *gene bodies*. This is essential for understanding TES function, as it is known to recruit DNMTs which may methylate gene bodies to enforce silencing, even if the binding site itself is at a distal enhancer.
4.  **Statistical Rigor:** The scripts employ appropriate statistical tests (Fisher's exact test for overlaps, Wilcoxon tests for signal differences) to validate findings.

## Usage
To run the full pipeline:
1.  Ensure you are in the `analysis_2` directory.
2.  Submit the master script:
    ```bash
    sbatch 00_run_all_new_modules.sh
    ```
3.  Monitor progress using `squeue -u <username>`.