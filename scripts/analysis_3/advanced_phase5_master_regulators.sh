#!/bin/bash
#SBATCH --job-name=phase5_2_master_reg_a3_
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --cpus-per-task=8
#SBATCH --time=03:00:00
#SBATCH --partition=workq
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --output=logs/advanced_phase5_2_master_regulators.out
#SBATCH --error=logs/advanced_phase5_2_master_regulators.err

echo "=================================================="
echo "Phase 5.2: Master Regulator Analysis"
echo "Start time: $(date)"
echo "=================================================="

# Change to project directory
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate analysis3_env

# Run R script
Rscript scripts/analysis_3/advanced_phase5_master_regulators.R

echo "=================================================="
echo "Analysis complete: $(date)"
echo "=================================================="
