# YAP/TAZ Target Enrichment Analysis

## Overview

This analysis investigates the overlap between known YAP/TAZ target genes and TES/TEAD1 binding sites from Cut&Tag data. Unlike other analyses in this directory that focus on differentially expressed genes (DEGs) from RNA-seq, this analysis uses a curated list of 224 known YAP/TAZ target genes from the literature.

## Input Data

### YAP/TAZ Target Genes
- **Source**: `/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/data/TES_degs.txt`
- **Description**: List of 224 known YAP/TAZ target genes
- **Format**: Text file with gene symbols

### Cut&Tag Binding Data
- **TES peaks**: From `SRF_Eva_CUTandTAG/results/07_analysis_narrow/TES_peaks_annotated.csv`
- **TEAD1 peaks**: From `SRF_Eva_CUTandTAG/results/07_analysis_narrow/TEAD1_peaks_annotated.csv`

### RNA-seq Expression Data (Optional)
- **DESeq2 results**: From `SRF_Eva_RNA/results/05_deseq2/deseq2_results_TES_vs_GFP.txt`
- Used to check if YAP/TAZ targets are also differentially expressed in our system

## Analysis Strategy

### 1. Gene Set Intersection
- Identify YAP/TAZ targets bound by TES
- Identify YAP/TAZ targets bound by TEAD1
- Classify into categories:
  - Bound by TES only
  - Bound by TEAD1 only
  - Bound by both TES and TEAD1
  - Not bound by either

### 2. Enrichment Testing
- Hypergeometric test: Are YAP/TAZ targets enriched in TES-bound genes?
- Hypergeometric test: Are YAP/TAZ targets enriched in TEAD1-bound genes?
- Universe: All genes with TES or TEAD1 peaks

### 3. GO Enrichment Analysis
- Run Gene Ontology enrichment on:
  - YAP/TAZ targets bound by TES
  - YAP/TAZ targets bound by TEAD1
  - YAP/TAZ targets bound by both

### 4. Expression Analysis
- Cross-reference with RNA-seq data to identify:
  - YAP/TAZ targets that are also DEGs
  - Expression direction (up/down) by binding status

## Running the Analysis

### Prerequisites
- Loaded data workspace: `SRF_Eva_integrated_analysis/focused_enrichment_analysis/data/loaded_data.RData`
- Ensure this file exists by running: `SRF_Eva_integrated_analysis/focused_enrichment_analysis/scripts/00_load_data.R`

### Execution

```bash
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top

# Submit to SLURM
sbatch SRF_Eva_integrated_analysis/focused_enrichment_analysis/yaptaz_targets_analysis/scripts/submit_yaptaz_analysis.sh

# Check job status
squeue -u $USER

# Monitor logs
tail -f SRF_Eva_integrated_analysis/focused_enrichment_analysis/yaptaz_targets_analysis/logs/yaptaz_analysis.out
```

### Alternative: Run locally
```bash
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate r-bio

Rscript SRF_Eva_integrated_analysis/focused_enrichment_analysis/yaptaz_targets_analysis/scripts/yaptaz_enrichment_analysis.R
```

## Output Files

### Gene Lists (`gene_lists/`)
- `yaptaz_targets_all.txt` - All 224 YAP/TAZ target genes (cleaned)
- `yaptaz_TES_bound.txt` - YAP/TAZ targets with TES binding
- `yaptaz_TEAD1_bound.txt` - YAP/TAZ targets with TEAD1 binding
- `yaptaz_both_bound.txt` - YAP/TAZ targets bound by both
- `yaptaz_TES_only_bound.txt` - YAP/TAZ targets bound only by TES
- `yaptaz_TEAD1_only_bound.txt` - YAP/TAZ targets bound only by TEAD1
- `yaptaz_no_binding.txt` - YAP/TAZ targets with no TES/TEAD1 binding

### Results Tables (`results/`)
- `yaptaz_TES_bound_detailed.csv` - Detailed annotation for TES-bound YAP/TAZ targets
- `yaptaz_TEAD1_bound_detailed.csv` - Detailed annotation for TEAD1-bound YAP/TAZ targets
- `yaptaz_*_GO_enrichment.csv` - GO enrichment results for each category
- `analysis_summary.txt` - Text summary of key findings

### Plots (`plots/`)
- `yaptaz_tes_tead1_venn.pdf` - Venn diagram showing three-way overlap
- `yaptaz_binding_categories.pdf` - Bar plot of binding categories
- `yaptaz_enrichment_comparison.pdf` - Enrichment significance comparison
- `yaptaz_expression_by_binding.pdf` - Expression changes by binding status
- `yaptaz_*_dotplot.pdf` - GO enrichment dot plots
- `yaptaz_*_barplot.pdf` - GO enrichment bar plots

## Key Differences from Approach 6

This analysis differs from the migration-focused Approach 6 in several ways:

1. **Gene Set Source**: Uses curated YAP/TAZ targets instead of GO/Hallmark gene sets
2. **Focus**: Binding pattern analysis rather than functional category enrichment
3. **Universe**: Uses all genes with peaks (not all expressed genes) for enrichment tests
4. **Direction**: Tests if known YAP/TAZ targets are enriched in our binding data (vs testing if bound genes are enriched for pathways)

## Biological Questions Addressed

1. **Conservation**: Do endogenous TEAD1 and synthetic TES bind to the same YAP/TAZ target genes?
2. **Specificity**: Are YAP/TAZ targets preferentially bound by TES vs TEAD1?
3. **Functional regulation**: Are YAP/TAZ targets also differentially expressed in our TES overexpression system?
4. **Pathway insight**: What biological processes are enriched among bound YAP/TAZ targets?

## Expected Results

- **High overlap expected**: Since TES contains a TEAD1 binding domain, we expect significant overlap with known YAP/TAZ targets
- **Enrichment p-values**: Should be highly significant if our Cut&Tag worked well
- **Expression correlation**: May not be perfect - binding doesn't always mean regulation in this context
- **GO enrichment**: Should reflect canonical YAP/TAZ biology (proliferation, migration, etc.)

## Computational Requirements

- **Runtime**: ~30-60 minutes
- **Memory**: 32GB (conservative, likely needs less)
- **CPUs**: 8 cores
- **Storage**: ~50MB for results

## Notes

- The YAP/TAZ target list contains some duplicate gene symbols (e.g., GADD45B, STMN1, CDC20) which are automatically handled
- Some gene symbols may not map to our RNA-seq data due to naming differences
- GO enrichment requires at least 10 genes, so smaller categories may not have enrichment results
