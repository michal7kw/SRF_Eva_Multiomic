#!/bin/bash
#SBATCH --job-name=a1_11_final_integrative_analysis
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --time=6:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --output=logs/11_final_integrative_analysis.out
#SBATCH --error=logs/11_final_integrative_analysis.err

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
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_1

# Load R environment
echo "Loading R environment..."
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate r_chipseq_env

# Run the all-genes integrative analysis
echo "Running all-genes integrative analysis..."
echo "Script: final_integrative_analysis_all_genes.R"
echo ""

Rscript 11_final_integrative_analysis_all_genes.R

echo ""
echo "=========================================="
echo "Analysis completed: $(date)"
echo "Job ID: $SLURM_JOB_ID finished"
echo "=========================================="