#!/bin/bash

# Master script to run the new multiomic analysis modules
# FIXED: Added job dependencies to ensure proper execution order

set -e

# Navigate to the correct directory
SCRIPT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_2"
cd "$SCRIPT_DIR"

# Create logs directory if it doesn't exist
mkdir -p logs

echo "=============================================="
echo "Starting Multiomic Analysis Pipeline"
echo "Working directory: $(pwd)"
echo "Date: $(date)"
echo "=============================================="

# Module 1: Peak Classification (no dependencies - runs first)
echo ""
echo "[1/5] Submitting Module 1: Peak Classification..."
JOB1=$(sbatch --parsable 01_peak_classification.sh)
echo "  Job ID: $JOB1"

# Module 2: Motif Analysis (depends on Module 1 for BED files) - SKIPPED
echo ""
echo "[2/5] Skipping Module 2: Motif Analysis..."
# JOB2=$(sbatch --parsable --dependency=afterok:$JOB1 02_motif_analysis.sh)
# echo "  Job ID: $JOB2"

# Module 3: meDIP Integration (depends on Module 1 for peak categories)
echo ""
echo "[3/5] Submitting Module 3: meDIP Integration (depends on Job $JOB1)..."
JOB3=$(sbatch --parsable --dependency=afterok:$JOB1 03_medip_integration.sh)
echo "  Job ID: $JOB3"

# Module 4: Gene Regulatory Logic (depends on Module 1 for peak categories)
echo ""
echo "[4/5] Submitting Module 4: Gene Regulatory Logic (depends on Job $JOB1)..."
JOB4=$(sbatch --parsable --dependency=afterok:$JOB1 04_gene_regulatory_logic.sh)
echo "  Job ID: $JOB4"

# Module 5: Functional Enrichment (depends on Module 1 for annotations)
echo ""
echo "[5/5] Submitting Module 5: Functional Enrichment (depends on Job $JOB1)..."
JOB5=$(sbatch --parsable --dependency=afterok:$JOB1 05_functional_enrichment.sh)
echo "  Job ID: $JOB5"

# Module 6: Methylation at Regulated Genes (depends on Module 1 for annotations)
echo ""
echo "[6/6] Submitting Module 6: Methylation at Regulated Genes (depends on Job $JOB1)..."
JOB6=$(sbatch --parsable --dependency=afterok:$JOB1 06_methylation_at_regulated_genes.sh)
echo "  Job ID: $JOB6"

echo ""
echo "=============================================="
echo "All jobs submitted!"
echo ""
echo "Job dependency chain:"
echo "  Module 1 (Peak Classification):    $JOB1"
echo "  Module 2 (Motif Analysis):         SKIPPED"
echo "  Module 3 (meDIP Integration):      $JOB3 -> depends on $JOB1"
echo "  Module 4 (Gene Regulatory Logic):  $JOB4 -> depends on $JOB1"
echo "  Module 5 (Functional Enrichment):  $JOB5 -> depends on $JOB1"
echo "  Module 6 (Methylation Analysis):   $JOB6 -> depends on $JOB1"
echo ""
echo "Monitor progress with:"
echo "  squeue -u \$USER"
echo "  tail -f logs/*.out"
echo "=============================================="
