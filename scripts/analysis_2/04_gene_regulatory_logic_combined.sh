#!/bin/bash
#SBATCH --job-name=04_gene_reg_comb
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --output=logs/04_gene_reg_combined.out
#SBATCH --error=logs/04_gene_reg_combined.err
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it

# Module 4 (Combined): Gene Regulatory Logic - All Peaks Combined
# Compares Promoter vs Enhancer/Distal expression changes with all peak categories merged

set -e

# Navigate to script directory
SCRIPT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_2"
cd "$SCRIPT_DIR"

# Create directories
mkdir -p logs
mkdir -p results/04_gene_regulatory_logic

echo "=== Module 4 (Combined): Gene Regulatory Logic - All Peaks ==="
echo "Working directory: $(pwd)"
echo "Date: $(date)"

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate seurat_full2

# Run script
Rscript 04_gene_regulatory_logic_combined.R

echo ""
echo "Module 4 (Combined) completed: $(date)"
