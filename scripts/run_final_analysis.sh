#!/bin/bash
#SBATCH --job-name=final_integrative
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --time=6:00:00
#SBATCH --mem=82G
#SBATCH --cpus-per-task=16
#SBATCH --output=logs/final_integrative.out
#SBATCH --error=logs/final_integrative.err

echo "=========================================="
echo "FINAL INTEGRATIVE ANALYSIS"
echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Running on node: $SLURM_NODELIST"
echo "Started: $(date)"
echo ""

# Navigate to project directory
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top

# Load R environment
echo "Loading R environment..."
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate r_chipseq_env

# Run the final integrative analysis
echo "Running final integrative analysis..."
echo "Script: integrative_analysis/scripts/final_integrative_analysis.R"
echo ""

Rscript integrative_analysis/scripts/final_integrative_analysis.R

echo ""
echo "=========================================="
echo "Analysis completed: $(date)"
echo "Job ID: $SLURM_JOB_ID finished"
echo "=========================================="