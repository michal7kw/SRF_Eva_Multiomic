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
![Target Overlap](output/10_final_integrative_analysis/plots/04_target_overlap_venn.png)
*Overlap of TES and TEAD1 direct targets.*

![Gene Classification](output/10_final_integrative_analysis/plots/01_gene_classification_pie.png)
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

---

# Overview

Phase 1 integrative analysis scripts combining Cut&Tag chromatin profiling with RNA-seq differential expression data for TES/TEAD1 regulatory network characterization. Identifies direct transcriptional targets (TF-bound + differentially expressed) vs indirect targets.

## Running Analyses

All scripts are submitted via SLURM. Working directory is this folder (`scripts/analysis_1`).

```bash
# Activate R environment first
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate r_chipseq_env

# Run complete pipeline
bash 0_run_all_analyses.sh

# Or run individual analyses
sbatch 10_final_integrative_analysis.sh      # RECOMMENDED: DEGs-only integrative analysis
sbatch 11_final_integrative_analysis_all_genes.sh  # All genes version
sbatch 12_msigdb_gsea_by_collection.sh       # MSigDB GSEA (~2-4 hours)
sbatch true_gsea_analysis.sh              # Rank-based GSEA
sbatch candidate_gene_heatmap.sh          # Candidate gene vs random controls
sbatch binding_heatmap.sh                 # Promoter binding signal heatmaps
sbatch enhanced_integrative_plots.sh      # Publication-quality figures
```

Monitor jobs: `squeue -u $USER`
View logs: `tail -f logs/<job_name>.out`

## Directory Structure

```
analysis_1/
├── *.R, *.sh          # Analysis scripts (paired R + SLURM wrapper)
├── logs/              # SLURM job logs
├── output/            # All generated results
│   ├── only_degs/     # DEG-focused analysis (RECOMMENDED)
│   │   ├── results/10_direct_targets/    # TES/TEAD1 target gene lists
│   │   ├── results/10_pathway_analysis/  # GO enrichment
│   │   ├── results/10_regulatory_networks/  # gene_classification.csv
│   │   └── plots/
│   ├── 12_msigdb_by_collection/  # GSEA by MSigDB collection
│   ├── candidate_gene_heatmap/
│   ├── binding_heatmap/
│   └── enhanced_integrative_plots/
```

Shared data accessed via relative paths:
- `../../data/` - Cached RData files, TES_degs.txt candidate gene list
- `../../GENE_SETS/` - 10,000+ MSigDB .grp gene set files
- `../../../SRF_Eva_RNA/results/05_deseq2/` - RNA-seq DESeq2 results
- `../../../SRF_Eva_CUTandTAG/results/` - Cut&Tag peaks and BigWig files

## Key Scripts

| Script | Purpose |
|--------|---------|
| `final_integrative_analysis.R` | Core integration: maps peaks to DEGs, classifies direct vs indirect targets, GO enrichment |
| `msigdb_gsea_by_collection.R` | Comprehensive GSEA using fgsea against all MSigDB collections |
| `true_gsea_analysis.R` | Rank-based GSEA using log2FC-ranked gene list |
| `candidate_gene_heatmap.R` | Compare TES/TEAD1 binding at 223 candidate genes vs matched random controls |
| `binding_heatmap.R` | TSS-centered binding signal heatmap from BigWig files |
| `cutandtag_density_plot.R` / `cutntag_density_heatmap.R` | Density profiles at promoters |
| `enhanced_integrative_plots.R` | Publication-quality multi-panel figures |

## Gene Classification

Scripts classify genes into regulatory categories:
- **TES_direct**: TES-bound + differentially expressed
- **TEAD1_direct**: TEAD1-bound + differentially expressed
- **TES_TEAD1_shared**: Both TF-bound + DE
- **Indirect**: DE but not TF-bound (downstream effects)
- **Unbound**: Not DE and not bound (all-genes analysis only)

## Gene ID Conversion

RNA-seq uses Ensembl IDs; Cut&Tag uses Entrez IDs. Scripts handle conversion via `org.Hs.eg.db`:
```r
# Ensembl → Symbol
mapIds(org.Hs.eg.db, keys = ensembl_id, column = "SYMBOL", keytype = "ENSEMBL")

# Entrez → Ensembl
mapIds(org.Hs.eg.db, keys = entrez_id, column = "ENSEMBL", keytype = "ENTREZID")
```

## Required R Packages

Genomics: `GenomicRanges`, `rtracklayer`, `ChIPseeker`, `TxDb.Hsapiens.UCSC.hg38.knownGene`
Enrichment: `clusterProfiler`, `fgsea`, `org.Hs.eg.db`
Visualization: `ComplexHeatmap`, `ggplot2`, `VennDiagram`, `pheatmap`, `EnrichedHeatmap`
Data: `dplyr`, `readr`, `tidyr`, `stringr`

