#!/bin/bash
#SBATCH --job-name=05_func_enrich
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --output=logs/05_func_enrich.out
#SBATCH --error=logs/05_func_enrich.err
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it

# Module 5: Functional Enrichment
# FIXED: Absolute paths, proper SLURM params, correct conda env

set -e

# Navigate to script directory
SCRIPT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_2"
cd "$SCRIPT_DIR"

# Create directories
mkdir -p logs
mkdir -p results/05_functional_enrichment

echo "=== Module 5: Functional Enrichment ==="
echo "Working directory: $(pwd)"
echo "Date: $(date)"

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate seurat_full2

# Run script
Rscript 05_functional_enrichment.R

echo ""
echo "Module 5 completed: $(date)"
