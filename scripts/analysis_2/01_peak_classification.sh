#!/bin/bash
#SBATCH --job-name=01_peak_class
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --output=logs/01_peak_class.out
#SBATCH --error=logs/01_peak_class.err
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it

# Module 1: Peak Classification
# FIXED: Absolute paths, proper SLURM params, correct conda env

set -e

# Navigate to script directory
SCRIPT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_2"
cd "$SCRIPT_DIR"

# Create logs and results directories
mkdir -p logs
mkdir -p results/01_peak_classification

echo "=== Module 1: Peak Classification ==="
echo "Working directory: $(pwd)"
echo "Date: $(date)"

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate chipseq_env

# Run script
Rscript 01_peak_classification.R

echo ""
echo "Module 1 completed: $(date)"
