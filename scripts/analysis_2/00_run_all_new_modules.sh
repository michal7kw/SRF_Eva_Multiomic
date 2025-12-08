#!/bin/bash

# Master script to run the new multiomic analysis modules
# Job dependencies ensure proper execution order:
#   - Modules 3, 4, 5 depend on Module 1 (need peak classification)
#   - Module 6 is INDEPENDENT (only needs RNA-seq + meDIP data)

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
echo "[1/6] Submitting Module 1: Peak Classification..."
JOB1=$(sbatch --parsable 01_peak_classification.sh)
echo "  Job ID: $JOB1"

# Module 2: Motif Analysis (depends on Module 1 for BED files) - SKIPPED
echo ""
echo "[2/6] Skipping Module 2: Motif Analysis..."
# JOB2=$(sbatch --parsable --dependency=afterok:$JOB1 02_motif_analysis.sh)
# echo "  Job ID: $JOB2"

# Module 3: meDIP Integration (depends on Module 1 for peak categories)
echo ""
echo "[3/6] Submitting Module 3: meDIP Integration (depends on Job $JOB1)..."
JOB3=$(sbatch --parsable --dependency=afterok:$JOB1 03_medip_integration.sh)
echo "  Job ID: $JOB3"

# Module 4: Gene Regulatory Logic (depends on Module 1 for peak categories)
echo ""
echo "[4/6] Submitting Module 4: Gene Regulatory Logic (depends on Job $JOB1)..."
JOB4=$(sbatch --parsable --dependency=afterok:$JOB1 04_gene_regulatory_logic.sh)
echo "  Job ID: $JOB4"

# Module 5: Functional Enrichment (depends on Module 1 for annotations)
echo ""
echo "[5/6] Submitting Module 5: Functional Enrichment (depends on Job $JOB1)..."
JOB5=$(sbatch --parsable --dependency=afterok:$JOB1 05_functional_enrichment.sh)
echo "  Job ID: $JOB5"

# Module 6: Methylation at Regulated Genes (INDEPENDENT - only needs RNA-seq + meDIP)
# Can run in parallel with Module 1
echo ""
echo "[6/6] Submitting Module 6: Methylation at Regulated Genes (INDEPENDENT)..."
JOB6=$(sbatch --parsable 06_methylation_at_regulated_genes.sh)
echo "  Job ID: $JOB6"

echo ""
echo "=============================================="
echo "All jobs submitted!"
echo ""
echo "Job dependency structure:"
echo "  Module 1 (Peak Classification):    $JOB1 (runs immediately)"
echo "  Module 2 (Motif Analysis):         SKIPPED"
echo "  Module 3 (meDIP Integration):      $JOB3 -> depends on $JOB1"
echo "  Module 4 (Gene Regulatory Logic):  $JOB4 -> depends on $JOB1"
echo "  Module 5 (Functional Enrichment):  $JOB5 -> depends on $JOB1"
echo "  Module 6 (Methylation Analysis):   $JOB6 (runs immediately, independent)"
echo ""
echo "Monitor progress with:"
echo "  squeue -u \$USER"
echo "  tail -f logs/*.out"
echo "=============================================="
