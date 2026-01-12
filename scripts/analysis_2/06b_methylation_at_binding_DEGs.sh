#!/bin/bash
#SBATCH --job-name=06b_meth_binding_DEGs
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --time=2:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/06b_methylation_at_binding_DEGs_%j.out
#SBATCH --error=logs/06b_methylation_at_binding_DEGs_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it

# 06b_methylation_at_binding_DEGs.sh
# Modified analysis: methylation at DEGs with TES/TEAD1 binding only

# Change to script directory
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_2

# Create logs directory if needed
mkdir -p logs

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate seurat_full2

echo "============================================"
echo "Module 6b: Methylation at Binding DEGs"
echo "============================================"
echo "Date: $(date)"
echo "Host: $(hostname)"
echo ""

# Run the R script
Rscript 06b_methylation_at_binding_DEGs.R

echo ""
echo "============================================"
echo "Completed: $(date)"
echo "============================================"
