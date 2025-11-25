#!/bin/bash

echo "=========================================="
echo "MASTER ANALYSIS EXECUTION SCRIPT"
echo "=========================================="
echo "Started: $(date)"
echo ""
echo "This script will submit all analyses:"
echo "  1. TRUE GSEA (rank-based enrichment)"
echo "  2. Directional GO enrichment (UP vs DOWN)"
echo "  3. Candidate gene analysis (TES_degs.txt)"
echo "  4. Promoter heatmaps (TES/TEAD1 at UP/DOWN genes)"
echo "  5. Cell death & proliferation pathway analysis"
echo ""
echo "All output will be saved to: output/<script_name>/"
echo ""
echo "Jobs will be submitted to SLURM queue."
echo "Monitor with: squeue -u $USER"
echo ""

# Base directory
BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_1"
LOGS_DIR="${BASE_DIR}/logs"

# Change to working directory
cd ${BASE_DIR}

# Create logs directory
mkdir -p ${LOGS_DIR}

# =============================================================================
# PHASE 1: TRUE GSEA
# =============================================================================

echo "=== Phase 1: TRUE GSEA ==="
JOB1=$(sbatch --parsable --output=${LOGS_DIR}/true_gsea.out --error=${LOGS_DIR}/true_gsea.err 01_true_gsea_analysis.sh)
echo "  Submitted TRUE GSEA: Job ID $JOB1"
echo ""

# =============================================================================
# PHASE 2: Directional GO Enrichment
# =============================================================================

echo "=== Phase 2: Directional GO Enrichment ==="
JOB2=$(sbatch --parsable --output=${LOGS_DIR}/directional_go.out --error=${LOGS_DIR}/directional_go.err 02_directional_go_enrichment.sh)
echo "  Submitted Directional GO: Job ID $JOB2"
echo ""

# =============================================================================
# PHASE 3: Candidate Gene Analysis
# =============================================================================

echo "=== Phase 3: Candidate Gene Analysis ==="
JOB3=$(sbatch --parsable --output=${LOGS_DIR}/candidate_gene.out --error=${LOGS_DIR}/candidate_gene.err 03_candidate_gene_heatmap.sh)
echo "  Submitted Candidate Analysis: Job ID $JOB3"
echo ""

# =============================================================================
# PHASE 4: Promoter Heatmaps
# NOTE: Requires integrative analysis results (direct_targets/*_all_genes.csv)
# =============================================================================

echo "=== Phase 4: Promoter Heatmaps ==="

# Check if integrative analysis results exist
INTEG_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis"
if [ -f "${INTEG_DIR}/output/results/10_direct_targets/TES_direct_targets_all_genes.csv" ]; then
    echo "  Integrative analysis results found - proceeding..."
    JOB4=$(sbatch --parsable --output=${LOGS_DIR}/promoter_heatmaps.out --error=${LOGS_DIR}/promoter_heatmaps.err 04_generate_promoter_heatmaps.sh)
    echo "  Submitted Promoter Heatmaps: Job ID $JOB4"
else
    echo "  WARNING: Integrative analysis results not found!"
    echo "  Please run: sbatch final_integrative_analysis_all_genes.sh"
    JOB4=""
fi
echo ""

# =============================================================================
# PHASE 5: Cell Death & Proliferation Pathway Analysis
# NOTE: Requires peak annotation from Cut&Tag pipeline
# =============================================================================

echo "=== Phase 5: Cell Death & Proliferation Analysis ==="

# Check if peak annotations exist
CUTANDTAG_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_CUTandTAG"
if [ -f "${CUTANDTAG_DIR}/results/07_analysis_narrow/TES_peaks_annotated.csv" ]; then
    echo "  Peak annotations found - proceeding..."
    JOB5=$(sbatch --parsable --output=${LOGS_DIR}/cell_death.out --error=${LOGS_DIR}/cell_death.err 05_cell_death_proliferation_analysis.sh)
    echo "  Submitted Cell Death Analysis: Job ID $JOB5"
else
    echo "  WARNING: Peak annotations not found!"
    echo "  Please run Cut&Tag annotation pipeline first."
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
echo "  TRUE GSEA: $JOB1"
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
echo "  All logs are saved in: ${LOGS_DIR}/"
echo "  tail -f ${LOGS_DIR}/true_gsea_*.out"
echo "  tail -f ${LOGS_DIR}/directional_go_*.out"
echo "  tail -f ${LOGS_DIR}/candidate_gene_*.out"
if [ ! -z "$JOB4" ]; then
    echo "  tail -f ${LOGS_DIR}/promoter_heatmaps_*.out"
fi
if [ ! -z "$JOB5" ]; then
    echo "  tail -f ${LOGS_DIR}/cell_death_*.out"
fi
echo ""
echo "Output directories:"
echo "  ${BASE_DIR}/output/true_gsea_analysis/"
echo "  ${BASE_DIR}/output/directional_go_enrichment/"
echo "  ${BASE_DIR}/output/candidate_gene_heatmap/"
if [ ! -z "$JOB4" ]; then
    echo "  ${BASE_DIR}/output/promoter_heatmaps/"
fi
if [ ! -z "$JOB5" ]; then
    echo "  ${BASE_DIR}/output/cell_death_proliferation/"
fi
echo ""
echo "Estimated completion: 2-4 hours"
echo ""
