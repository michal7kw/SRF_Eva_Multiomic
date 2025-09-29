#!/bin/bash
#SBATCH --job-name=all_genes_integrative
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --time=6:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --output=logs/all_genes_integrative.out
#SBATCH --error=logs/all_genes_integrative.err

echo "=========================================="
echo "INTEGRATIVE ANALYSIS - ALL GENES"
echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Running on node: $SLURM_NODELIST"
echo "Started: $(date)"
echo ""
echo "Analysis mode: Using ALL genes in the dataset"
echo "Approach: Direct targets = TF-bound + DE genes"
echo "Background: Complete transcriptome for pathway analysis"
echo ""

# Navigate to project directory
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top

# Load R environment
echo "Loading R environment..."
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate r_chipseq_env

# Run the all-genes integrative analysis
echo "Running all-genes integrative analysis..."
echo "Script: integrative_analysis/scripts/final_integrative_analysis_all_genes.R"
echo ""

Rscript integrative_analysis/scripts/final_integrative_analysis_all_genes.R

echo ""
echo "=========================================="
echo "Analysis completed: $(date)"
echo "Job ID: $SLURM_JOB_ID finished"
echo "=========================================="