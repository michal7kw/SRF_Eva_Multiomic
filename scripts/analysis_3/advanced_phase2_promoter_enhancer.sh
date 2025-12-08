#!/bin/bash
#SBATCH --job-name=phase2_2_promoter_enh_a3_
#SBATCH --account=kubacki.michal
#SBATCH --mem=16GB
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --partition=workq
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --output=logs/advanced_phase2_2_promoter_enhancer.out
#SBATCH --error=logs/advanced_phase2_2_promoter_enhancer.err

echo "=================================================="
echo "Phase 2.2: Promoter vs Enhancer Binding Effects"
echo "Start time: $(date)"
echo "=================================================="
# Change to project directory
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis


# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate analysis3_env 

# Run R script
Rscript scripts/analysis_3/advanced_phase2_promoter_enhancer.R

echo "=================================================="
echo "Analysis complete: $(date)"
echo "=================================================="
