# Analysis 1: Integrated Transcriptional Regulation Analysis
**SRF_Eva_integrated_analysis**

**Maintainer:** kubacki.michal@hsr.it

## Overview
This directory contains the core analysis pipeline for integrating Cut&Tag (TES, TEAD1) and RNA-seq (TES vs GFP) data. The analysis aims to define the direct transcriptional targets of TES and TEAD1, explore their regulatory logic, and identify downstream biological pathways relevant to cancer.

## 🚀 Quick Start
To run the full pipeline:
```bash
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_1
sbatch 00_run_all_analysis_1_modules.sh
```
This master script manages the entire pipeline using SLURM job dependencies.

## 📂 Pipeline Modules

### Phase 1: Foundation (Basic Analysis)
| Module | Script | Output |
|:---|:---|:---|
| **01** | `01_true_gsea_analysis.R` | Rank-based GSEA (fgsea) | `output/01_true_gsea_analysis/` |
| **02** | `02_directional_go_enrichment.R` | GO Enrichment (Up vs Down) | `output/02_directional_go_enrichment/` |
| **03** | `03_candidate_gene_heatmap.R` | Candidate gene binding profiles | `output/03_candidate_gene_heatmap/` |
| **05** | `05_cell_death_proliferation...` | Targeted pathway analysis | `output/05_cell_death_proliferation/` |

### Phase 2: Integrative Analysis (Core Logic)
| Module | Script | Output |
|:---|:---|:---|
| **10** | `10_final_integrative_analysis.R` | **Core Integration**: Binding + Expression | `output/10_final_integrative_analysis/` |
| **11** | `11_final_integrative_..._all_genes.R` | Genome-wide integration (incl. non-DEGs) | `output/11_final_integrative_analysis_all_genes/` |

### Phase 3: MSigDB & Advanced GSEA
| Module | Script | Output |
|:---|:---|:---|
| **12** | `12_msigdb_gsea_by_collection.R` | GSEA across all MSigDB collections | `output/12_msigdb_by_collection/` |
| **13** | `13_msigdb_gsea_selected.R` | Focused GSEA on select gene sets | `output/13_msigdb_gsea_selected/` |
| **17** | `17_gsea_cancer_pathways.R` | Cancer-specific pathway extraction | `output/17_gsea_cancer_pathways/` |
| **18** | `18_gsea_cancer_..._IMPROVED.R` | Directional Cancer GSEA (Up/Down) | `output/18_gsea_cancer_pathways_improved/` |

### Phase 4: Visualization & Deep Dives
| Module | Script | Output |
|:---|:---|:---|
| **04** | `04_generate_promoter_heatmaps.R` | Promoter region definitions | `output/04_generate_promoter_heatmaps/` |
| **14** | `14_binding_heatmap.R` | Signal intensity heatmaps | `output/14_binding_heatmap/` |
| **15** | `15_cutandtag_density_plot.R` | TSS Density Plots | `output/15_cutandtag_density_plot/` |
| **16** | `16_enhanced_integrative_plots.R` | Publication-ready summary figures | `output/16_enhanced_integrative_plots/` |
| **16** | `16_enhanced_directional_go.R` | Enhanced GO visualization | `output/16_enhanced_directional_go/` |
| **16** | `16_enhanced_gsea_visualizations.R` | Enhanced GSEA plots | `output/16_enhanced_gsea_visualizations/` |

---

## 🔬 Bioinformatics Notes
A comprehensive review of the scripts in this directory was conducted to ensure logic, consistency, and adherence to best practices.

1.  **GSEA Methodology:** The pipeline correctly uses `fgsea` for rank-based enrichment analysis (Scripts 01, 12, 13, 17, 18), which is superior to simple overlap-based methods (ORA) for detecting subtle but coordinated pathway changes.
2.  **Integrative Logic:** The definition of "Direct Targets" (Script 10) is robust, requiring both significant differential expression (RNA-seq) and proximal binding (Cut&Tag peaks within promoter/enhancer regions).
3.  **Directionality:** The "Improved" cancer pathway analysis (Script 18) adds critical biological context by distinguishing between pathways driven by up-regulated vs. down-regulated genes.
4.  **Control Selection:** In `03_candidate_gene_heatmap.R`, the selection of random control genes is based on expression matching.
5.  **Gene ID Conversion:** All scripts consistently handle Gene ID conversions (Ensembl <-> Entrez <-> Symbol) using `org.Hs.eg.db`.

