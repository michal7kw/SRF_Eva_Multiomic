#!/bin/bash
#SBATCH --job-name=a1_10_final_integrative_analysis
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --time=4:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --output=logs/10_final_integrative_analysis.out
#SBATCH --error=logs/10_final_integrative_analysis.err

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
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_1

# Load R environment
echo "Loading R environment..."
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate r_chipseq_env

# Run the DEG-only integrative analysis
echo "Running DEG-only integrative analysis..."
echo "Script: final_integrative_analysis.R"
echo ""

Rscript 10_final_integrative_analysis.R

echo ""
echo "=========================================="
echo "Analysis completed: $(date)"
echo "Job ID: $SLURM_JOB_ID finished"
echo "=========================================="