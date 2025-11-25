#!/bin/bash

echo "=========================================="
echo "MASTER ANALYSIS EXECUTION SCRIPT"
echo "=========================================="
echo "Started: $(date)"
echo ""
echo "This script will submit all remaining analyses:"
echo "  1. TRUE GSEA (rank-based enrichment)"
echo "  2. RNA-seq directional GO enrichment"
echo "  3. Candidate gene analysis (TES_degs.txt)"
echo "  4. Promoter heatmaps (TES/TEAD1 at UP/DOWN genes)"
echo "  5. Cell death & proliferation pathway analysis"
echo ""
echo "Jobs will be submitted to SLURM queue."
echo "Monitor with: squeue -u $USER"
echo ""

# Base directory
BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top"

# =============================================================================
# PHASE 1: TRUE GSEA (already completed, but can re-run if needed)
# =============================================================================

echo "=== Phase 1: TRUE GSEA ==="
echo "Status: Already completed"
echo "Output: ${BASE_DIR}/SRF_Eva_integrated_analysis/output/true_gsea_analysis/"
echo "  Skipping (already completed - remove this block to re-run)"
echo ""

# To re-run TRUE GSEA, uncomment these lines:
cd ${BASE_DIR}/SRF_Eva_integrated_analysis
JOB1=$(sbatch --parsable scripts/analysis_1/true_gsea_analysis.sh)
echo "  Submitted TRUE GSEA: Job ID $JOB1"
GSEA_DEP="--dependency=afterok:$JOB1"

GSEA_DEP=""
echo ""

# === SRF_Eva_RNA ==========================================================================
# PHASE 2: RNA-seq Directional GO Enrichment
# =============================================================================

echo "=== Phase 2: RNA-seq Directional GO Enrichment ==="
cd ${BASE_DIR}/SRF_Eva_RNA
JOB2=$(sbatch --parsable scripts/directional_go_enrichment.sh)
echo "  Submitted Directional GO: Job ID $JOB2"
echo ""

# ==== SRF_Eva_integrated_analysis =========================================================================
# PHASE 3: Candidate Gene Heatmap
# =============================================================================

echo "=== Phase 3: Candidate Gene Analysis ==="
cd ${BASE_DIR}/SRF_Eva_integrated_analysis
JOB3=$(sbatch --parsable scripts/analysis_1/candidate_gene_heatmap.sh)
echo "  Submitted Candidate Analysis: Job ID $JOB3"
echo ""

# ==== SRF_Eva_CUTandTAG =========================================================================
# PHASE 4: Promoter Heatmaps
# NOTE: Requires integrative analysis results (direct_targets/*_all_genes.csv)
# =============================================================================

echo "=== Phase 4: Promoter Heatmaps ==="
cd ${BASE_DIR}/SRF_Eva_CUTandTAG

# Check if integrative analysis results exist
if [ -f "${BASE_DIR}/SRF_Eva_integrated_analysis/output/results/direct_targets/TES_direct_targets_all_genes.csv" ]; then
    echo "  Integrative analysis results found - proceeding..."
    JOB4=$(sbatch --parsable scripts/generate_promoter_heatmaps.sh)
    echo "  Submitted Promoter Heatmaps: Job ID $JOB4"
else
    echo "  WARNING: Integrative analysis results not found!"
    echo "  Please run integrative analysis first: sbatch SRF_Eva_integrated_analysis/scripts/analysis_1/final_integrative_analysis_all_genes.sh"
    JOB4=""
fi
echo ""

# ==== SRF_Eva_CUTandTAG =========================================================================
# PHASE 5: Cell Death & Proliferation Pathway Analysis
# NOTE: Requires peak annotation from Cut&Tag pipeline
# =============================================================================

echo "=== Phase 5: Cell Death & Proliferation Analysis ==="
cd ${BASE_DIR}/SRF_Eva_CUTandTAG

# Check if peak annotations exist
if [ -f "${BASE_DIR}/SRF_Eva_CUTandTAG/results/07_analysis_narrow/TES_peaks_annotated.csv" ]; then
    echo "  Peak annotations found - proceeding..."
    JOB5=$(sbatch --parsable cell_death_proliferation_analysis.sh)
    echo "  Submitted Cell Death Analysis: Job ID $JOB5"
else
    echo "  WARNING: Peak annotations not found!"
    echo "  Please run peak annotation first: sbatch 8_annotate_narrow.sh"
    JOB5=""
fi
echo ""

# =============================================================================
# SUMMARY
# =============================================================================

echo "=========================================="
echo "ALL JOBS SUBMITTED!"
echo "=========================================="
echo ""
echo "Job IDs:"
if [ ! -z "$JOB1" ]; then
    echo "  TRUE GSEA: $JOB1"
fi
echo "  Directional GO: $JOB2"
echo "  Candidate Analysis: $JOB3"
if [ ! -z "$JOB4" ]; then
    echo "  Promoter Heatmaps: $JOB4"
fi
if [ ! -z "$JOB5" ]; then
    echo "  Cell Death Analysis: $JOB5"
fi
echo ""
echo "Monitor progress:"
echo "  squeue -u $USER"
echo ""
echo "Check logs:"
echo "  tail -f ${BASE_DIR}/SRF_Eva_RNA/logs/directional_go.out"
echo "  tail -f ${BASE_DIR}/SRF_Eva_integrated_analysis/logs/candidate_heatmap.out"
if [ ! -z "$JOB4" ]; then
    echo "  tail -f ${BASE_DIR}/SRF_Eva_CUTandTAG/logs/promoter_heatmaps.out"
fi
if [ ! -z "$JOB5" ]; then
    echo "  tail -f ${BASE_DIR}/SRF_Eva_CUTandTAG/logs/cell_death_analysis.out"
fi
echo ""
echo "Output directories:"
echo "  ${BASE_DIR}/SRF_Eva_integrated_analysis/output/true_gsea_analysis/"
echo "  ${BASE_DIR}/SRF_Eva_RNA/results/06_directional_go/"
echo "  ${BASE_DIR}/SRF_Eva_integrated_analysis/output/candidate_gene_heatmap/"
if [ ! -z "$JOB4" ]; then
    echo "  ${BASE_DIR}/SRF_Eva_CUTandTAG/results/08_promoter_heatmaps/"
fi
if [ ! -z "$JOB5" ]; then
    echo "  ${BASE_DIR}/SRF_Eva_CUTandTAG/results/cell_death_proliferation/"
fi
echo ""
echo "Estimated completion: 2-4 hours"
echo ""
