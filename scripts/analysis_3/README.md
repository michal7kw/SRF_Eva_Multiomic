# Advanced Integrative Analysis Pipeline (Analysis 3)

## Overview

Advanced multi-omics integrative analysis pipeline for **TES vs TEAD1** comparative study in SNB19 glioblastoma cells. This pipeline integrates:

- **Cut&Tag**: TES/TEAD1 binding site classification and signal quantification
- **RNA-seq**: Differential expression (TES vs GFP)
- **meDIP-seq**: DNA methylation at regulatory regions

The pipeline classifies binding sites using a **dual classification system**, discovers motifs, maps regulatory cascades, integrates methylation data, and prioritizes targets for experimental validation.

---

## Dual Peak Classification System

### Background: Why Two Classification Schemes?

The pipeline provides two complementary classification approaches:

| Scheme | Categories | Purpose |
|--------|-----------|---------|
| **Detailed (6-category)** | Signal-based subcategorization | Fine-grained mechanistic analysis |
| **Simplified (3-category)** | Overlap-based only | Direct comparison with analysis_2 |

### Detailed Classification (6 Categories)

Based on peak overlap AND relative signal intensity:

| Category | Definition | Biological Interpretation |
|----------|------------|---------------------------|
| **TES_unique** | TES peak with no TEAD1 overlap | Novel TES binding sites |
| **TEAD1_unique** | TEAD1 peak with no TES overlap | Endogenous TEAD1-only sites |
| **Shared_high** | Overlapping, both signals ≥75th percentile | High-confidence co-bound sites |
| **Shared_TES_dominant** | Overlapping, TES signal >2× TEAD1 | TES preferentially binds |
| **Shared_TEAD1_dominant** | Overlapping, TEAD1 signal >2× TES | TEAD1 preferentially binds |
| **Shared_equivalent** | Overlapping, similar signal (within 2×) | Equal co-occupancy |

**Classification Parameters:**
- Overlap threshold: 500bp
- Signal ratio threshold: 0.2 (log2, ~1.15-fold)
- High signal percentile: 75th percentile

### Simplified Classification (3 Categories)

Direct overlap-based classification (comparable to analysis_2):

| Category | Definition | Maps From |
|----------|------------|-----------|
| **TES_Unique** | TES peak with no TEAD1 overlap | TES_unique |
| **TEAD1_Unique** | TEAD1 peak with no TES overlap | TEAD1_unique |
| **Shared** | Any overlapping peaks | All Shared_* categories |

### Comparison with Analysis_2 Pipeline

| Aspect | Analysis_2 | Analysis_3 |
|--------|-----------|------------|
| **Input peaks** | Replicate consensus (≥2/3 support) | Pooled MACS2 peak calls |
| **Classification** | 3 categories (overlap-only) | 6 categories (signal-based) + 3-category simplified |
| **Signal quantification** | None | BigWig mean signal at peaks |
| **Output** | Single classification | Dual classification system |

**Expected count differences:**
- Analysis_3 has more peaks (pooled calling captures more sites)
- Simplified 3-category counts enable direct pipeline comparison

---

## Pipeline Architecture

### Phase Overview

```
Phase 1: Binding & Signal Analysis
├── 1.1: Binding Classification (dual 6-cat + 3-cat)
├── 1.2: Motif Discovery (HOMER)
└── 1.3: Signal Dynamics (deepTools)

Phase 2: Regulatory Logic
├── 2.1: Expression by Binding Category
├── 2.2: Promoter vs Enhancer Effects
└── 2.3: Regulatory Cascades (TF networks)

Phase 3: Methylation Integration
├── 3.1: Methylation at Binding Sites
└── 3.2: Three-way Integration (binding + methylation + expression)

Phase 5: Network Analysis
├── 5.1: Coregulator Discovery
├── 5.2: Master Regulator Identification
└── 5.3: Pathway Crosstalk

Phase 6: Target Prioritization
└── Multi-metric scoring for experimental validation

Phase 8: Publication Figures
└── Final publication-ready visualizations
```

*Note: Phase 4 is intentionally skipped to align with top-level project numbering.*

### Dependency Graph

```
Phase 1.1 (binding classification) ─┬─> Phase 1.2 (motifs)
                                    ├─> Phase 1.3 (signal dynamics)
                                    └─> Phase 2.1 (expression by category)
                                            │
                                            └─> Phase 2.2 (promoter/enhancer)
                                                    │
                                                    └─> Phase 3 (methylation integration)
                                                            │
                                                            └─> Phase 5.2 (master regulators) ──┐
                                                                    │                           │
                                                                    ├─> Phase 5.1, 5.3          │
                                                                    │                           │
                                                                    └─> Phase 2.3 (cascades)*───┘
                                                                            │
                                                                            └─> Phase 6 (prioritization)
                                                                                    │
                                                                                    └─> Phase 8 (figures)

* Phase 2.3 requires Phase 5.2 output - run Phase 5 before Phase 2.3
```

---

## Phase 1: Binding & Signal Analysis

### Phase 1.1: Advanced Binding Classification

**Script:** `advanced_01_binding_classification.R`

**Purpose:** Classify TES and TEAD1 peaks into binding categories using both overlap-based and signal-based criteria.

**Methodology:**
1. Load TES and TEAD1 narrowPeak files (pooled MACS2 calls)
2. Extract mean CPM signal from BigWig files at each peak
3. Identify overlapping peaks (500bp threshold)
4. Classify based on signal ratios into 6 detailed categories
5. Generate simplified 3-category classification for comparison
6. Annotate genomic context using ChIPseeker

**Input Files:**
- `SRF_Eva_CUTandTAG/results/05_peaks_narrow/TES_peaks.narrowPeak`
- `SRF_Eva_CUTandTAG/results/05_peaks_narrow/TEAD1_peaks.narrowPeak`
- `SRF_Eva_CUTandTAG/results/06_bigwig/TES-{1,2,3}_CPM.bw`
- `SRF_Eva_CUTandTAG/results/06_bigwig/TEAD1-{1,2,3}_CPM.bw`

**Output Files:**

| File | Description |
|------|-------------|
| `binding_site_classification.csv` | Main classification table with both `category` (6-cat) and `category_simple` (3-cat) columns |
| `binding_category_summary.txt` | Summary statistics for both classification schemes |
| `binding_classification_data.RData` | R objects for downstream analysis |

**Detailed Category BED Files (6 categories):**
| File | Description |
|------|-------------|
| `TES_unique.bed` | TES-only binding sites |
| `TEAD1_unique.bed` | TEAD1-only binding sites |
| `Shared_high.bed` | High-signal co-bound sites |
| `Shared_TES_dominant.bed` | TES-dominant shared sites |
| `Shared_TEAD1_dominant.bed` | TEAD1-dominant shared sites |
| `Shared_equivalent.bed` | Equivalent co-bound sites |

**Simplified Category BED Files (3 categories):**
| File | Description |
|------|-------------|
| `TES_Unique_simple.bed` | All TES-only peaks |
| `Shared_simple.bed` | All shared peaks (combined) |
| `TEAD1_Unique_simple.bed` | All TEAD1-only peaks |

**Visualization Outputs:**

| File | Description |
|------|-------------|
| `category_counts_barplot.pdf` | Detailed 6-category distribution |
| `category_counts_barplot_simple.pdf` | Simplified 3-category distribution |
| `category_comparison_both_schemes.pdf` | Side-by-side comparison of both schemes |
| `shared_category_breakdown.pdf` | Pie chart: how Shared splits into subcategories |
| `signal_comparison_scatterplot.pdf` | TES vs TEAD1 signal (6 categories) |
| `signal_comparison_simple.pdf` | TES vs TEAD1 signal (3 categories) |
| `binding_category_venn.pdf` | Venn diagram of peak overlap |
| `genomic_context_enrichment.pdf` | Genomic distribution by category |
| `peak_width_distribution.pdf` | Peak width violin plots |

### Phase 1.2: Motif Analysis

**Script:** `advanced_02_motif_analysis.sh`

**Purpose:** De novo motif discovery using HOMER for each binding category.

**Prerequisites:** HOMER hg38 genome must be installed:
```bash
perl /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/genomics_env/share/homer/configureHomer.pl -install hg38
```

**Output:** `results/02_motif_analysis/<category>_motifs/`

### Phase 1.3: Signal Dynamics

**Script:** `advanced_03_signal_dynamics.sh`

**Purpose:** Quantify binding signal intensity and dynamics using deepTools.

**Output:** `results/03_signal_dynamics/`

---

## Phase 2: Regulatory Logic

### Phase 2.1: Expression by Binding Category

**Script:** `advanced_phase2_expression_by_category.R`

**Purpose:** Correlate binding category with expression magnitude.

**Key Outputs:**
- `results/04_category_expression/expression_by_category_boxplots.pdf`
- `results/04_category_expression/binding_signal_vs_expression.pdf`

### Phase 2.2: Promoter vs Enhancer Effects

**Script:** `advanced_phase2_promoter_enhancer.R`

**Purpose:** Compare regulatory effects of promoter-proximal vs distal binding.

**Key Outputs:**
- `results/05_promoter_enhancer/promoter_enhancer_effects.pdf`
- `results/05_promoter_enhancer/distance_decay_plot.pdf`

### Phase 2.3: Regulatory Cascades

**Script:** `advanced_phase2_3_regulatory_cascades.R`

**Purpose:** Construct TF-to-TF regulatory networks.

**Dependency:** Requires Phase 5.2 (master regulators) output.

**Key Outputs:**
- `results/06_regulatory_cascades/regulatory_cascade_network.pdf`
- `results/06_regulatory_cascades/cascade_TF_overview.pdf`

---

## Phase 3: Methylation Integration

### Phase 3.1: Methylation at Binding Sites

**Script:** `advanced_phase3_1_methylation_binding.R`

**Purpose:** Analyze overlap between differentially methylated regions (DMRs) and binding site categories.

**Key Findings (from latest run):**
- TEAD1-unique peaks are significantly enriched for hypermethylation (OR=1.89, p=0.003)
- TES-related peaks show reduced hypermethylation (OR=0.53, p=0.006)
- TES_unique peaks uniquely show hypomethylation (0.29% vs 0% for other categories)

**Output Files:**
- `results/07_methylation_binding/PHASE3_1_SUMMARY.txt`
- `results/07_methylation_binding/dmr_overlap_by_category.csv`
- `results/07_methylation_binding/methylation_binding_barplot.pdf`

### Phase 3.2: Three-Way Integration

**Script:** `advanced_phase3_promoter_methylation_expression.R`

**Purpose:** Integrate binding, methylation, and expression for mechanistic classification.

**Key Outputs:**
- `results/08_methylation_expression/three_way_integration_heatmap.pdf`
- `results/08_methylation_expression/mechanistic_classification_pie.pdf`

---

## Phase 5: Network Analysis

### Phase 5.1: Coregulator Discovery

**Script:** `advanced_phase5_1_coregulators.R`

**Purpose:** Predict cooperative transcription factors based on motif enrichment and expression correlation.

**Output:** `results/10_coregulators/`

### Phase 5.2: Master Regulator Identification

**Script:** `advanced_phase5_master_regulators.R`

**Purpose:** Identify key transcription factor drivers of the regulatory program.

**Output:** `results/09_tf_networks/master_regulator_analysis.csv`

### Phase 5.3: Pathway Crosstalk

**Script:** `advanced_phase5_3_pathway_crosstalk.R`

**Purpose:** Map interactions between enriched pathways.

**Output:** `results/11_pathway_crosstalk/`

---

## Phase 6: Target Prioritization

**Script:** `advanced_phase6_target_prioritization.R`

**Purpose:** Rank genes for CRISPR/experimental validation using multi-metric scoring.

**Scoring Metrics:**
- Binding signal strength
- Expression fold change
- Statistical significance
- Methylation status
- Pathway centrality
- Druggability

**Output:** `results/12_target_prioritization/`

---

## Phase 8: Publication Figures

**Script:** `advanced_phase8_publication_figures.R`

**Purpose:** Generate publication-ready figures summarizing the entire study.

**Output:** `results/13_publication_figures/`

---

## Prerequisites

### Required Data

Ensure upstream pipelines have completed:

```
SRF_Eva_CUTandTAG/results/05_peaks_narrow/
  ├── TES_peaks.narrowPeak       # Pooled MACS2 peak calls
  └── TEAD1_peaks.narrowPeak

SRF_Eva_CUTandTAG/results/06_bigwig/
  ├── TES-{1,2,3}_CPM.bw         # CPM-normalized coverage
  └── TEAD1-{1,2,3}_CPM.bw

SRF_Eva_RNA/results/05_deseq2/
  └── deseq2_results_TES_vs_GFP.txt

meDIP/results/07_differential_MEDIPS/
  └── TES_vs_GFP_DMRs_FDR05_FC2.csv  # Optional, for Phase 3
```

### Conda Environments

```bash
# Activate base conda
source /opt/common/tools/ric.cosr/miniconda3/bin/activate

# Primary environment for R analysis
conda activate r_chipseq_env

# For HOMER motif analysis (Phase 1.2)
conda activate genomics_env

# For deepTools signal dynamics (Phase 1.3)
conda activate signal_dynamics
```

### HOMER hg38 Setup (for Phase 1.2)

```bash
perl /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/genomics_env/share/homer/configureHomer.pl -install hg38
```

---

## Quick Start Guide

### Run Complete Pipeline

```bash
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis

# ============================================================
# Phase 1: Binding Classification and Motif Analysis
# ============================================================
sbatch scripts/analysis_3/advanced_01_binding_classification.sh  # START HERE
# Wait for completion
sbatch scripts/analysis_3/advanced_02_motif_analysis.sh
sbatch scripts/analysis_3/advanced_03_signal_dynamics.sh

# ============================================================
# Phase 2.1-2.2: Regulatory Logic (first part)
# ============================================================
sbatch scripts/analysis_3/advanced_phase2_expression_by_category.sh
sbatch scripts/analysis_3/advanced_phase2_promoter_enhancer.sh

# ============================================================
# Phase 3: Epigenetic Mechanisms
# ============================================================
sbatch scripts/analysis_3/advanced_phase3_1_methylation_binding.sh
sbatch scripts/analysis_3/advanced_phase3_promoter_methylation_expression.sh

# ============================================================
# Phase 5: Network Analysis (BEFORE Phase 2.3!)
# ============================================================
sbatch scripts/analysis_3/advanced_phase5_master_regulators.sh
sbatch scripts/analysis_3/advanced_phase5_1_coregulators.sh
sbatch scripts/analysis_3/advanced_phase5_3_pathway_crosstalk.sh

# ============================================================
# Phase 2.3: Regulatory Cascades (requires Phase 5 output)
# ============================================================
sbatch scripts/analysis_3/advanced_phase2_3_regulatory_cascades.sh

# ============================================================
# Phase 6: Target Prioritization
# ============================================================
sbatch scripts/analysis_3/advanced_phase6_target_prioritization.sh

# ============================================================
# Phase 8: Publication Figures (run last)
# ============================================================
sbatch scripts/analysis_3/advanced_phase8_publication_figures.sh
```

### Monitor Jobs

```bash
squeue -u $USER
tail -f scripts/analysis_3/logs/advanced_*.out
```

---

## Key Analysis Parameters

### Peak Classification (Phase 1.1)

| Parameter | Value | Description |
|-----------|-------|-------------|
| Overlap threshold | 500bp | Maximum distance for peak overlap |
| Signal ratio threshold | 0.2 (log2) | ~1.15-fold for equivalent classification |
| High signal percentile | 75th | Threshold for Shared_high category |

### Peak-to-Gene Mapping

| Parameter | Value | Description |
|-----------|-------|-------------|
| Promoter annotation | ChIPseeker | Based on TxDb annotations |
| Enhancer threshold | 50kb | Maximum distance from TSS |

### Statistical Thresholds

| Parameter | Value | Description |
|-----------|-------|-------------|
| DEG significance | padj < 0.05 | FDR-corrected |
| Pathway enrichment | qvalue < 0.05 | FDR-corrected |
| log2FC threshold | None | Captures all significant changes |

### Gene ID Harmonization

```
RNA-seq (Ensembl) → Cut&Tag (Entrez) → meDIP (Symbol)
```
Conversion via `org.Hs.eg.db` R package.

---

## SLURM Configuration

Standard job parameters:
```bash
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it
```

Typical resources:
| Phase | Memory | CPUs | Time |
|-------|--------|------|------|
| 1.1 Binding | 32GB | 8 | 2h |
| 1.2 Motifs | 16GB | 8 | 2h |
| 2.x Regulatory | 32GB | 8 | 2h |
| 3.x Methylation | 32GB | 8 | 2h |
| 5.x Networks | 32GB | 8 | 2h |
| 6 Prioritization | 32GB | 8 | 2h |
| 8 Figures | 16GB | 4 | 1h |

---

## Key R Libraries

| Category | Packages |
|----------|----------|
| **Genomics** | GenomicRanges, rtracklayer, ChIPseeker, TxDb.Hsapiens.UCSC.hg38.knownGene |
| **Enrichment** | clusterProfiler, org.Hs.eg.db, msigdbr, fgsea |
| **Visualization** | ggplot2, ComplexHeatmap, pheatmap, VennDiagram, ggrepel |
| **Networks** | igraph, visNetwork |
| **Data manipulation** | dplyr, tidyr, readr, stringr |

---

## Troubleshooting

### Phase 1.1: Signal values are NaN

**Cause:** Chromosome naming mismatch between peaks and BigWig files.

**Solution:** The script automatically handles both UCSC (chr1) and Ensembl (1) naming conventions. If issues persist:
1. Check BigWig files exist: `ls SRF_Eva_CUTandTAG/results/06_bigwig/*.bw`
2. Verify chromosome names: `bigWigInfo TES-1_CPM.bw | head`

### Phase 1.2: HOMER fails

**Cause:** hg38 genome not installed for HOMER.

**Solution:**
```bash
perl /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/genomics_env/share/homer/configureHomer.pl -install hg38
```

### Phase 2.3: Master regulators file not found

**Cause:** Phase 5.2 not completed before Phase 2.3.

**Solution:** Run Phase 5.2 first:
```bash
sbatch scripts/analysis_3/advanced_phase5_master_regulators.sh
# Wait for completion, then:
sbatch scripts/analysis_3/advanced_phase2_3_regulatory_cascades.sh
```

### Phase 3: No DMR overlaps

**Cause:** meDIP analysis not completed or no significant DMRs.

**Solution:** Verify meDIP output exists:
```bash
ls meDIP/results/07_differential_MEDIPS/TES_vs_GFP_DMRs_FDR05_FC2.csv
```

---

## Output Directory Structure

```
results/
├── 01_binding_classification/
│   ├── binding_site_classification.csv   # Main classification table
│   ├── binding_category_summary.txt      # Summary statistics
│   ├── *.bed                             # Category-specific BED files
│   └── *.pdf                             # Visualizations
├── 02_motif_analysis/
│   └── <category>_motifs/                # HOMER results per category
├── 03_signal_dynamics/
├── 04_category_expression/
├── 05_promoter_enhancer/
├── 06_regulatory_cascades/
├── 07_methylation_binding/
├── 08_methylation_expression/
├── 09_tf_networks/
├── 10_coregulators/
├── 11_pathway_crosstalk/
├── 12_target_prioritization/
└── 13_publication_figures/
```

---

## Notes

- **TESmut is excluded** from this comparative analysis (TES vs TEAD1 only)
- All scripts set working directory to `/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top`
- Output paths are relative to this working directory
- Phase 4 is intentionally skipped to maintain consistency with top-level project numbering

---

## Related Documentation

- **CLAUDE.md**: AI assistant guidance for this pipeline
- **INTERPRETATION_GUIDE.md**: Biological interpretation of results
- **RESULTS.md**: Summary of analysis findings
- **Parent CLAUDE.md**: `SRF_Eva_integrated_analysis/CLAUDE.md`
- **Project CLAUDE.md**: `/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/CLAUDE.md`

---

*Last updated: 2025-11-27*
