#!/bin/bash
#SBATCH --job-name=03_medip_int
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --output=logs/03_medip_int.out
#SBATCH --error=logs/03_medip_int.err
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it

# Module 3: meDIP Integration
# FIXED: Absolute paths, proper SLURM params, correct conda env

set -e

# Navigate to script directory
SCRIPT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_2"
cd "$SCRIPT_DIR"

# Create directories
mkdir -p logs
mkdir -p results/03_medip_integration

echo "=== Module 3: meDIP Integration ==="
echo "Working directory: $(pwd)"
echo "Date: $(date)"

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate seurat_full2

# Run script
Rscript 03_medip_integration.R

echo ""
echo "Module 3 completed: $(date)"
