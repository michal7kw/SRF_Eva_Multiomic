#!/bin/bash
#SBATCH --job-name=integrative_analysis
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --time=4:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --output=logs/integrative_analysis.out
#SBATCH --error=logs/integrative_analysis.err

echo "Starting integrative analysis: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Running on node: $SLURM_NODELIST"

# Load R environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate r_chipseq_env

# Run the integrative analysis
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top

echo "Running integrative analysis script..."
Rscript integrative_analysis/scripts/integrative_analysis.R

echo "Integrative analysis completed: $(date)"