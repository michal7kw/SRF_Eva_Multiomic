#!/bin/bash

# Master script to run all analysis_1 modules with proper dependency management
# Executes modules in logical order with parallel execution where possible

set -e

# Navigate to the correct directory
SCRIPT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_1"
cd "$SCRIPT_DIR"

# Create logs directory if it doesn't exist
mkdir -p logs

echo "=============================================="
echo "Starting Analysis_1 Pipeline"
echo "Working directory: $(pwd)"
echo "Date: $(date)"
echo "=============================================="

# Define module arrays for parallel execution
# Phase 1: Foundation analyses (can run in parallel)
PHASE_1_MODULES=(
    "01_true_gsea_analysis.sh"
    "02_directional_go_enrichment.sh"
    "03_candidate_gene_heatmap.sh"
    "05_cell_death_proliferation_analysis.sh"
)

# Phase 2: Integrative analyses (depend on Phase 1)
PHASE_2_MODULES=(
    "10_final_integrative_analysis.sh"
    "11_final_integrative_analysis_all_genes.sh"
)

# Phase 3: MSigDB analyses (depend on Phase 1)
PHASE_3_MODULES=(
    "12_msigdb_gsea_by_collection.sh"
    "13_msigdb_gsea_selected.sh"
)

# Phase 4: Visualization modules (depend on their respective analyses)
PHASE_4A_MODULES=(
    "04_generate_promoter_heatmaps.sh"  # depends on integrative analysis
)

PHASE_4B_MODULES=(
    "14_binding_heatmap.sh"
    "15_cutandtag_density_plot.sh"
)

PHASE_4C_MODULES=(
    "16_enhanced_integrative_plots.sh"
)

# NOTE: PHASE_4D_MODULES removed - these R scripts are already run by 16_enhanced_integrative_plots.sh

PHASE_4E_MODULES=(
    "17_gsea_cancer_pathways.sh"
    "18_gsea_cancer_pathways_IMPROVED.sh"
)

# Phase 1: Submit foundation modules in parallel
echo ""
echo "[PHASE 1] Submitting Foundation Analyses (parallel execution)..."

# Submit Phase 1 modules
PHASE_1_JOB_IDS=()
for module in "${PHASE_1_MODULES[@]}"; do
    if [ -f "$module" ]; then
        echo "  Submitting: $module"
        JOB_ID=$(sbatch --parsable "$module")
        PHASE_1_JOB_IDS+=("$JOB_ID")
        echo "    Job ID: $JOB_ID"
    else
        echo "  Warning: $module not found, skipping"
    fi
done

# Create dependency string for Phase 1
PHASE_1_DEPS=$(IFS=:; echo "${PHASE_1_JOB_IDS[*]}")
echo "  Phase 1 dependency string: $PHASE_1_DEPS"

# Phase 2: Submit integrative analyses (depend on all Phase 1)
echo ""
echo "[PHASE 2] Submitting Integrative Analyses (depends on Phase 1)..."

PHASE_2_JOB_IDS=()
for module in "${PHASE_2_MODULES[@]}"; do
    if [ -f "$module" ]; then
        echo "  Submitting: $module"
        JOB_ID=$(sbatch --parsable --dependency=afterok:$PHASE_1_DEPS "$module")
        PHASE_2_JOB_IDS+=("$JOB_ID")
        echo "    Job ID: $JOB_ID"
    else
        echo "  Warning: $module not found, skipping"
    fi
done

# Create dependency string for Phase 2
PHASE_2_DEPS=$(IFS=:; echo "${PHASE_2_JOB_IDS[*]}")
echo "  Phase 2 dependency string: $PHASE_2_DEPS"

# Phase 3: Submit MSigDB analyses (depend on Phase 1)
echo ""
echo "[PHASE 3] Submitting MSigDB Analyses (depends on Phase 1)..."

PHASE_3_JOB_IDS=()
for module in "${PHASE_3_MODULES[@]}"; do
    if [ -f "$module" ]; then
        echo "  Submitting: $module"
        JOB_ID=$(sbatch --parsable --dependency=afterok:$PHASE_1_DEPS "$module")
        PHASE_3_JOB_IDS+=("$JOB_ID")
        echo "    Job ID: $JOB_ID"
    else
        echo "  Warning: $module not found, skipping"
    fi
done

# Create dependency string for Phase 3
PHASE_3_DEPS=$(IFS=:; echo "${PHASE_3_JOB_IDS[*]}")
echo "  Phase 3 dependency string: $PHASE_3_DEPS"

# Phase 4A: Submit promoter heatmap generation (depends on Phase 2)
echo ""
echo "[PHASE 4A] Submitting Promoter Heatmap Generation (depends on Phase 2)..."

PHASE_4A_JOB_IDS=()
for module in "${PHASE_4A_MODULES[@]}"; do
    if [ -f "$module" ]; then
        echo "  Submitting: $module"
        JOB_ID=$(sbatch --parsable --dependency=afterok:$PHASE_2_DEPS "$module")
        PHASE_4A_JOB_IDS+=("$JOB_ID")
        echo "    Job ID: $JOB_ID"
    else
        echo "  Warning: $module not found, skipping"
    fi
done

# Phase 4B: Submit binding and density analyses (depends on Phase 1)
echo ""
echo "[PHASE 4B] Submitting Binding and Density Analyses (depends on Phase 1)..."

PHASE_4B_JOB_IDS=()
for module in "${PHASE_4B_MODULES[@]}"; do
    if [ -f "$module" ]; then
        echo "  Submitting: $module"
        JOB_ID=$(sbatch --parsable --dependency=afterok:$PHASE_1_DEPS "$module")
        PHASE_4B_JOB_IDS+=("$JOB_ID")
        echo "    Job ID: $JOB_ID"
    else
        echo "  Warning: $module not found, skipping"
    fi
done

# Phase 4C: Submit enhanced integrative plots (depends on Phase 2)
echo ""
echo "[PHASE 4C] Submitting Enhanced Integrative Plots (depends on Phase 2)..."

PHASE_4C_JOB_IDS=()
for module in "${PHASE_4C_MODULES[@]}"; do
    if [ -f "$module" ]; then
        echo "  Submitting: $module"
        JOB_ID=$(sbatch --parsable --dependency=afterok:$PHASE_2_DEPS "$module")
        PHASE_4C_JOB_IDS+=("$JOB_ID")
        echo "    Job ID: $JOB_ID"
    else
        echo "  Warning: $module not found, skipping"
    fi
done

# NOTE: Phase 4D removed - 16_enhanced_integrative_plots.sh (Phase 4C) already runs
# 16_enhanced_directional_go.R and 16_enhanced_gsea_visualizations.R internally.
# The previous implementation had bugs:
# 1. Re-submitted Phase 1 jobs that were already running
# 2. Used invalid sbatch syntax (sbatch ... Rscript X.R)

# Phase 4E: Submit cancer pathway analyses (depends on Phase 2)
echo ""
echo "[PHASE 4E] Submitting Cancer Pathway Analyses (depends on Phase 2)..."

PHASE_4E_JOB_IDS=()
for module in "${PHASE_4E_MODULES[@]}"; do
    if [ -f "$module" ]; then
        echo "  Submitting: $module"
        JOB_ID=$(sbatch --parsable --dependency=afterok:$PHASE_2_DEPS "$module")
        PHASE_4E_JOB_IDS+=("$JOB_ID")
        echo "    Job ID: $JOB_ID"
    else
        echo "  Warning: $module not found, skipping"
    fi
done

# Combine all job IDs for final summary
ALL_JOB_IDS=("${PHASE_1_JOB_IDS[@]}" "${PHASE_2_JOB_IDS[@]}" "${PHASE_3_JOB_IDS[@]}" "${PHASE_4A_JOB_IDS[@]}" "${PHASE_4B_JOB_IDS[@]}" "${PHASE_4C_JOB_IDS[@]}" "${PHASE_4E_JOB_IDS[@]}")

echo ""
echo "=============================================="
echo "All jobs submitted!"
echo ""
echo "Job dependency structure:"
echo "  Phase 1 (Foundation): ${#PHASE_1_JOB_IDS[@]} jobs"
for i in "${!PHASE_1_JOB_IDS[@]}"; do
    echo "    ${PHASE_1_MODULES[$i]}: ${PHASE_1_JOB_IDS[$i]}"
done
echo ""
echo "  Phase 2 (Integrative): ${#PHASE_2_JOB_IDS[@]} jobs (depends on all Phase 1)"
for i in "${!PHASE_2_JOB_IDS[@]}"; do
    echo "    ${PHASE_2_MODULES[$i]}: ${PHASE_2_JOB_IDS[$i]}"
done
echo ""
echo "  Phase 3 (MSigDB): ${#PHASE_3_JOB_IDS[@]} jobs (depends on Phase 1)"
for i in "${!PHASE_3_JOB_IDS[@]}"; do
    echo "    ${PHASE_3_MODULES[$i]}: ${PHASE_3_JOB_IDS[$i]}"
done
echo ""
echo "  Phase 4A (Promoter Heatmaps): ${#PHASE_4A_JOB_IDS[@]} jobs (depends on Phase 2)"
for i in "${!PHASE_4A_JOB_IDS[@]}"; do
    echo "    ${PHASE_4A_MODULES[$i]}: ${PHASE_4A_JOB_IDS[$i]}"
done
echo ""
echo "  Phase 4B (Binding/Density): ${#PHASE_4B_JOB_IDS[@]} jobs (depends on Phase 1)"
for i in "${!PHASE_4B_JOB_IDS[@]}"; do
    echo "    ${PHASE_4B_MODULES[$i]}: ${PHASE_4B_JOB_IDS[$i]}"
done
echo ""
echo "  Phase 4C (Enhanced Integrative): ${#PHASE_4C_JOB_IDS[@]} jobs (depends on Phase 2)"
for i in "${!PHASE_4C_JOB_IDS[@]}"; do
    echo "    ${PHASE_4C_MODULES[$i]}: ${PHASE_4C_JOB_IDS[$i]}"
done
echo ""
echo "  Phase 4E (Cancer Pathways): ${#PHASE_4E_JOB_IDS[@]} jobs (depends on Phase 2)"
for i in "${!PHASE_4E_JOB_IDS[@]}"; do
    echo "    ${PHASE_4E_MODULES[$i]}: ${PHASE_4E_JOB_IDS[$i]}"
done
echo ""
echo "Total jobs submitted: ${#ALL_JOB_IDS[@]}"
echo ""
echo "Monitor progress with:"
echo "  squeue -u \$USER"
echo "  tail -f logs/*.out"
echo ""
echo "Output directories:"
echo "  - output/01_true_gsea_analysis/"
echo "  - output/02_directional_go_enrichment/"
echo "  - output/03_candidate_gene_heatmap/"
echo "  - output/05_cell_death_proliferation/"
echo "  - output/10_final_integrative_analysis/"
echo "  - output/11_final_integrative_analysis_all_genes/"
echo "  - output/12_msigdb_by_collection/"
echo "  - output/13_msigdb_gsea_selected/"
echo "  - output/04_generate_promoter_heatmaps/"
echo "  - output/14_binding_heatmap/"
echo "  - output/15_cutandtag_density_plot/"
echo "  - output/16_enhanced_integrative_plots/"
echo "  - output/16_enhanced_directional_go/"
echo "  - output/16_enhanced_gsea_visualizations/"
echo "  - output/17_gsea_cancer_pathways/"
echo "  - output/18_gsea_cancer_pathways_improved/"
echo "=============================================="