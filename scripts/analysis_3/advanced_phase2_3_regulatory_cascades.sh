#!/bin/bash
#SBATCH --job-name=phase2_3_cascades_a3_
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --cpus-per-task=8
#SBATCH --time=03:00:00
#SBATCH --partition=workq
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --output=logs/advanced_phase2_3_cascades.out
#SBATCH --error=logs/advanced_phase2_3_cascades.err

echo "=================================================="
echo "Phase 2.3: Regulatory Cascades and Secondary Targets"
echo "Start time: $(date)"
echo "=================================================="
# Change to project directory
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis


# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate analysis3_env

# Run R script
Rscript scripts/analysis_3/advanced_phase2_3_regulatory_cascades.R

echo "=================================================="
echo "Analysis complete: $(date)"
echo "=================================================="
