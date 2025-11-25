# Advanced Integrative Analysis: Regulatory Networks and Mechanisms

This directory contains scripts for the advanced phase of the multiomic analysis, focusing on regulatory cascades, coregulator discovery, pathway crosstalk, and target prioritization.

## Prerequisites

### Required Data
- Completed Cut&Tag peak calling (`SRF_Eva_CUTandTAG/results/05_peaks_narrow/`)
- BigWig files (`SRF_Eva_CUTandTAG/results/06_bigwig/`)
- RNA-seq differential expression results (`SRF_Eva_RNA/results/05_deseq2/`)
- meDIP methylation data (optional, for Phase 3)

### Conda Environments
- **r_chipseq_env**: For R-based analysis (GenomicRanges, ChIPseeker, clusterProfiler)
- **genomics_env**: For HOMER motif analysis

### HOMER Setup (for Phase 1.2 Motif Analysis)
HOMER requires the hg38 genome to be installed. Install it with:

```bash
perl /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/genomics_env/share/homer/configureHomer.pl -install hg38
```

This downloads the hg38 genome sequences and prepares them for motif analysis.

## Overview

The advanced analysis pipeline builds upon the initial integration to:
1.  **Map Regulatory Cascades**: Identify downstream transcription factors regulated by TES/TEAD1.
2.  **Dissect Binding Context**: Analyze promoter vs. enhancer effects and signal dynamics.
3.  **Integrate Methylation Mechanisms**: Classify targets by methylation-dependent vs. independent regulation.
4.  **Identify Coregulators**: Discover cooperative TFs using motif and expression data.
5.  **Prioritize Targets**: Rank genes for validation using a multi-metric scoring system.

## Pipeline Steps

### Phase 1: Advanced Binding & Signal Analysis
-   **`advanced_01_binding_classification.sh/R`**: Classify peaks into TES-unique, TEAD1-unique, and shared categories with signal quantification.
    -   *Outputs*: `results/01_binding_classification/binding_site_classification.csv`, `<category>.bed`, visualization PDFs.
-   **`advanced_02_motif_analysis.sh`**: HOMER de novo motif discovery for each binding category.
    -   *Requires*: HOMER hg38 genome installed (see Prerequisites).
    -   *Outputs*: `results/02_motif_analysis/<category>_motifs/`, `motif_enrichment_heatmap.pdf`, `top_motifs_by_category.pdf`.
-   **`advanced_03_signal_dynamics.sh`**: Analysis of binding signal intensity and width.
    -   *Outputs*: `signal_comparison_scatterplot.pdf`, `peak_width_distribution.pdf`.

### 2. Regulatory Logic (Phase 2)
-   **`advanced_phase2_expression_by_category.R`**: Correlates binding category with expression magnitude.
    -   *Outputs*: `expression_by_category_boxplots.pdf`, `binding_signal_vs_expression.pdf`.
-   **`advanced_phase2_promoter_enhancer.R`**: Compares efficacy of promoter vs. enhancer binding.
    -   *Outputs*: `promoter_enhancer_effects.pdf`, `distance_decay_plot.pdf`.
-   **`advanced_phase2_3_regulatory_cascades.R`**: Constructs TF-to-TF regulatory networks.
    -   *Outputs*: `regulatory_cascade_network.pdf`, `cascade_TF_overview.pdf`.

### 3. Epigenetic Mechanisms (Phase 3)
-   **`advanced_phase3_promoter_methylation_expression.R`**: Integrates meDIP, RNA-seq, and binding to classify mechanisms.
    -   *Outputs*: `three_way_integration_heatmap.pdf`, `mechanistic_classification_pie.pdf`.

### 4. Network & Pathway Analysis (Phase 5)
-   **`advanced_phase5_master_regulators.R`**: Identifies key drivers of the observed phenotype.
    -   *Outputs*: `master_regulator_volcano.pdf`, `top_master_regulators_barplot.pdf`.
-   **`advanced_phase5_1_coregulators.R`**: Predicts co-factors based on motif enrichment and expression correlation.
    -   *Outputs*: `coregulator_enrichment_heatmap.pdf`, `motif_expression_bubble.pdf`.
-   **`advanced_phase5_3_pathway_crosstalk.R`**: Maps interactions between enriched pathways.
    -   *Outputs*: `pathway_network_graph.pdf`, `pathway_similarity_heatmap.pdf`.

### 5. Target Prioritization (Phase 6)
-   **`advanced_phase6_target_prioritization.R`**: Ranks targets for CRISPR/experimental validation.
    -   *Outputs*: `top_50_score_heatmap.pdf`, `CRISPR_candidates_summary.pdf`, `priority_score_distribution.pdf`.

### 6. Publication Figures (Phase 8)
-   **`advanced_phase8_publication_figures.R`**: Generates polished, publication-ready figures summarizing the entire study.
    -   *Outputs*: `Figure_1_Overview.pdf`, `Figure_2_Mechanisms.pdf`, etc.

## Generated Outputs

### Regulatory Networks
![Cascade Network](results/advanced_analysis/phase2_regulatory_cascades/regulatory_cascade_network.png)
*Network of transcription factors directly regulated by TES/TEAD1.*

### Epigenetic Mechanisms
![Three-Way Heatmap](results/advanced_analysis/phase3_methylation/three_way_integration_heatmap.png)
*Integration of Binding, Methylation, and Expression to define regulatory modes.*

### Coregulator Discovery
![Coregulator Heatmap](results/advanced_analysis/phase5_coregulators/coregulator_enrichment_heatmap.png)
*Enrichment of co-factor motifs in TES/TEAD1 binding sites.*

### Target Prioritization
![Top Targets](results/advanced_analysis/phase6_prioritization/top_50_score_heatmap.png)
*Top 50 high-priority targets ranked by multi-omic evidence.*

## Usage

### IMPORTANT: Execution Order and Dependencies

The scripts have specific dependencies and must be run in the correct order. **Phase 4 was intentionally skipped** in the numbering to maintain consistency with the overall project phase structure.

### Dependency Graph
```
Phase 1.1 (Binding Classification)
    │
    ├── Phase 1.2 (Motif Analysis)
    │
    └── Phase 2.1 (Expression by Category)
            │
            └── Phase 2.2 (Promoter/Enhancer)
                    │
                    └── Phase 3 (Methylation Integration)
                            │
                            └── Phase 5.2 (Master Regulators) ─────┐
                                    │                              │
                                    ├── Phase 5.1 (Coregulators)   │
                                    │                              │
                                    └── Phase 2.3 (Regulatory Cascades)*
                                            │
                                            └── Phase 6 (Target Prioritization)
                                                    │
                                                    └── Phase 8 (Publication Figures)

* Phase 2.3 requires Phase 5.2 output (master regulators) - this is intentional
```

Run the scripts sequentially from the `SRF_Eva_integrated_analysis` directory:

```bash
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis

# ============================================================
# Phase 1: Binding Classification and Motif Analysis
# ============================================================
sbatch scripts/analysis_3/advanced_01_binding_classification.sh  # START HERE
# Wait for completion before continuing
sbatch scripts/analysis_3/advanced_02_motif_analysis.sh          # Requires Phase 1.1 BED files
sbatch scripts/analysis_3/advanced_03_signal_dynamics.sh         # Can run parallel with 1.2

# ============================================================
# Phase 2.1-2.2: Regulatory Logic (first part)
# ============================================================
sbatch scripts/analysis_3/advanced_phase2_expression_by_category.sh
sbatch scripts/analysis_3/advanced_phase2_promoter_enhancer.sh

# ============================================================
# Phase 3: Epigenetic Mechanisms (requires meDIP data)
# ============================================================
sbatch scripts/analysis_3/advanced_phase3_promoter_methylation_expression.sh

# ============================================================
# Phase 5: Network Analysis (BEFORE Phase 2.3!)
# ============================================================
sbatch scripts/analysis_3/advanced_phase5_master_regulators.sh   # Run BEFORE Phase 2.3
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

### Key Notes:
- **Phase 4 is not missing** - it was intentionally skipped to align with the top-level project numbering
- **Phase 2.3 (Regulatory Cascades)** requires Phase 5.2 (Master Regulators) output - run Phase 5 before 2.3
- **Phase 8 (Publication Figures)** should be run last as it aggregates all outputs

Monitor job progress with:
```bash
squeue -u $USER
tail -f scripts/analysis_3/logs/advanced_*.out
```
