#!/bin/bash
#SBATCH --job-name=phase5_3_crosstalk_a3_
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --cpus-per-task=8
#SBATCH --time=03:00:00
#SBATCH --partition=workq
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --output=logs/advanced_phase5_3_crosstalk.out
#SBATCH --error=logs/advanced_phase5_3_crosstalk.err

echo "=================================================="
echo "Phase 5.3: Pathway Crosstalk Analysis"
echo "Start time: $(date)"
echo "=================================================="
# Change to project directory
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis


# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate analysis3_env

# Run R script
Rscript scripts/analysis_3/advanced_phase5_3_pathway_crosstalk.R

echo "=================================================="
echo "Analysis complete: $(date)"
echo "=================================================="
