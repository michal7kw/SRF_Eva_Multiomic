#!/bin/bash
#SBATCH --job-name=deg_only_integrative
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --time=4:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --output=logs/deg_only_integrative.out
#SBATCH --error=logs/deg_only_integrative.err

echo "=========================================="
echo "INTEGRATIVE ANALYSIS - DEGs ONLY"
echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Running on node: $SLURM_NODELIST"
echo "Started: $(date)"
echo ""
echo "Analysis mode: Using only significant DEGs"
echo "Approach: Direct targets = TF-bound significant DEGs"
echo ""

# Navigate to project directory
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top

# Load R environment
echo "Loading R environment..."
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate r_chipseq_env

# Run the DEG-only integrative analysis
echo "Running DEG-only integrative analysis..."
echo "Script: integrative_analysis/scripts/final_integrative_analysis.R"
echo ""

Rscript integrative_analysis/scripts/final_integrative_analysis.R

echo ""
echo "=========================================="
echo "Analysis completed: $(date)"
echo "Job ID: $SLURM_JOB_ID finished"
echo "=========================================="