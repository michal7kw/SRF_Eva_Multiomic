#!/bin/bash

# Master script to run the advanced multiomic analysis modules with parallel execution
# Optimized dependency management for maximum parallelization

set -e

# Navigate to the correct directory
SCRIPT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_3"
cd "$SCRIPT_DIR"

# Create logs directory if it doesn't exist
mkdir -p logs

echo "=============================================="
echo "Starting Advanced Multiomic Analysis Pipeline"
echo "Working directory: $(pwd)"
echo "Date: $(date)"
echo "=============================================="

# Phase 1.1: Binding Classification (must run first - no dependencies)
echo ""
echo "[Phase 1.1] Submitting Binding Classification..."
PHASE1_1_JOB=$(sbatch --parsable advanced_01_binding_classification.sh)
echo "  Job ID: $PHASE1_1_JOB"

# Phase 1.2: Motif Analysis (depends on Phase 1.1)
echo ""
echo "[Phase 1.2] Submitting Motif Analysis (depends on Job $PHASE1_1_JOB)..."
PHASE1_2_JOB=$(sbatch --parsable --dependency=afterok:$PHASE1_1_JOB advanced_02_motif_analysis.sh)
echo "  Job ID: $PHASE1_2_JOB"

# Phase 1.3: Signal Dynamics (depends on Phase 1.1)
echo ""
echo "[Phase 1.3] Submitting Signal Dynamics (depends on Job $PHASE1_1_JOB)..."
PHASE1_3_JOB=$(sbatch --parsable --dependency=afterok:$PHASE1_1_JOB advanced_03_signal_dynamics.sh)
echo "  Job ID: $PHASE1_3_JOB"

# Phase 2.1: Expression Analysis (depends on Phase 1.1)
echo ""
echo "[Phase 2.1] Submitting Expression Analysis (depends on Job $PHASE1_1_JOB)..."
PHASE2_1_JOB=$(sbatch --parsable --dependency=afterok:$PHASE1_1_JOB advanced_phase2_expression_by_category.sh)
echo "  Job ID: $PHASE2_1_JOB"

# Phase 2.2: Promoter/Enhancer Effects (depends on Phase 1.1)
echo ""
echo "[Phase 2.2] Submitting Promoter/Enhancer Effects (depends on Job $PHASE1_1_JOB)..."
PHASE2_2_JOB=$(sbatch --parsable --dependency=afterok:$PHASE1_1_JOB advanced_phase2_promoter_enhancer.sh)
echo "  Job ID: $PHASE2_2_JOB"

# Phase 5.1: Master Regulators (depends on Phase 2.1)
echo ""
echo "[Phase 5.1] Submitting Master Regulators (depends on Job $PHASE2_1_JOB)..."
PHASE5_1_JOB=$(sbatch --parsable --dependency=afterok:$PHASE2_1_JOB advanced_phase5_master_regulators.sh)
echo "  Job ID: $PHASE5_1_JOB"

# Phase 5.2: Co-regulators (depends on Phase 2.1)
echo ""
echo "[Phase 5.2] Submitting Co-regulators (depends on Job $PHASE2_1_JOB)..."
PHASE5_2_JOB=$(sbatch --parsable --dependency=afterok:$PHASE2_1_JOB advanced_phase5_1_coregulators.sh)
echo "  Job ID: $PHASE5_2_JOB"

# Phase 5.3: Pathway Crosstalk (depends on Phase 2.1)
echo ""
echo "[Phase 5.3] Submitting Pathway Crosstalk (depends on Job $PHASE2_1_JOB)..."
PHASE5_3_JOB=$(sbatch --parsable --dependency=afterok:$PHASE2_1_JOB advanced_phase5_3_pathway_crosstalk.sh)
echo "  Job ID: $PHASE5_3_JOB"

# Phase 3.1: Methylation at Binding (depends on Phase 1.1)
echo ""
echo "[Phase 3.1] Submitting Methylation at Binding (depends on Job $PHASE1_1_JOB)..."
PHASE3_1_JOB=$(sbatch --parsable --dependency=afterok:$PHASE1_1_JOB advanced_phase3_methylation_at_binding.sh)
echo "  Job ID: $PHASE3_1_JOB"

# Phase 3.2: Promoter Methylation & Expression (depends on Phase 3.1 AND Phase 2.1)
echo ""
echo "[Phase 3.2] Submitting Promoter Methylation & Expression (depends on Jobs $PHASE3_1_JOB AND $PHASE2_1_JOB)..."
PHASE3_2_JOB=$(sbatch --parsable --dependency=afterok:$PHASE3_1_JOB:$PHASE2_1_JOB advanced_phase3_promoter_methylation_expression.sh)
echo "  Job ID: $PHASE3_2_JOB"

# Phase 2.3: Regulatory Cascades (depends on Phase 5.2 - Co-regulators)
echo ""
echo "[Phase 2.3] Submitting Regulatory Cascades (depends on Job $PHASE5_2_JOB)..."
PHASE2_3_JOB=$(sbatch --parsable --dependency=afterok:$PHASE5_2_JOB advanced_phase2_3_regulatory_cascades.sh)
echo "  Job ID: $PHASE2_3_JOB"

# Phase 6: Target Prioritization (depends on Phase 2.3 AND Phase 3.2)
echo ""
echo "[Phase 6] Submitting Target Prioritization (depends on Jobs $PHASE2_3_JOB AND $PHASE3_2_JOB)..."
PHASE6_JOB=$(sbatch --parsable --dependency=afterok:$PHASE2_3_JOB:$PHASE3_2_JOB advanced_phase6_target_prioritization.sh)
echo "  Job ID: $PHASE6_JOB"

# Phase 8: Publication Figures (depends on ALL previous phases completing)
# This is the final step that depends on all other phases
echo ""
echo "[Phase 8] Submitting Publication Figures (depends on ALL previous phases)..."
# Create dependency list for all critical jobs
ALL_JOBS_DEPS="afterok:$PHASE1_1_JOB:$PHASE1_2_JOB:$PHASE1_3_JOB:$PHASE2_1_JOB:$PHASE2_2_JOB:$PHASE5_1_JOB:$PHASE5_2_JOB:$PHASE5_3_JOB:$PHASE3_1_JOB:$PHASE3_2_JOB:$PHASE2_3_JOB:$PHASE6_JOB"
PHASE8_JOB=$(sbatch --parsable --dependency=$ALL_JOBS_DEPS advanced_phase8_publication_figures.sh)
echo "  Job ID: $PHASE8_JOB"

echo ""
echo "=============================================="
echo "All advanced analysis jobs submitted!"
echo ""
echo "Job dependency structure:"
echo "  Phase 1.1 (Binding Classification):    $PHASE1_1_JOB"
echo "    ├─ Phase 1.2 (Motif Analysis):         $PHASE1_2_JOB"
echo "    ├─ Phase 1.3 (Signal Dynamics):        $PHASE1_3_JOB"
echo "    ├─ Phase 2.1 (Expression Analysis):    $PHASE2_1_JOB"
echo "    │   ├─ Phase 5.1 (Master Regulators):   $PHASE5_1_JOB"
echo "    │   ├─ Phase 5.2 (Co-regulators):       $PHASE5_2_JOB"
echo "    │   │   └─ Phase 2.3 (Regulatory Cascades): $PHASE2_3_JOB"
echo "    │   │       └─ Phase 6 (Target Prioritization): $PHASE6_JOB"
echo "    │   └─ Phase 5.3 (Pathway Crosstalk):  $PHASE5_3_JOB"
echo "    ├─ Phase 2.2 (Promoter/Enhancer Effects): $PHASE2_2_JOB"
echo "    └─ Phase 3.1 (Methylation at Binding): $PHASE3_1_JOB"
echo "        └─ Phase 3.2 (Promoter Methylation): $PHASE3_2_JOB"
echo "            └─ Phase 6 (Target Prioritization): $PHASE6_JOB"
echo "  Phase 8 (Publication Figures):           $PHASE8_JOB (depends on ALL)"
echo ""
echo "Parallel execution groups:"
echo "  Group 1 (after Phase 1.1): Phase 1.2, 1.3, 2.1, 2.2, 3.1 run in parallel"
echo "  Group 2 (after Phase 2.1): Phase 5.1, 5.2, 5.3 run in parallel"
echo "  Group 3 (after Phase 3.1): Phase 3.2 runs"
echo ""
echo "Monitor progress with:"
echo "  squeue -u \$USER"
echo "  tail -f logs/*.out"
echo ""
echo "Check job dependencies:"
echo "  scontrol show job <jobid>"
echo "=============================================="