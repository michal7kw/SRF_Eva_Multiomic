#!/bin/bash
#SBATCH --job-name=paper_go_plot
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --mem=8GB
#SBATCH --cpus-per-task=2
#SBATCH --time=00:30:00
#SBATCH --output=logs/02b_paper_go_plot_%j.out
#SBATCH --error=logs/02b_paper_go_plot_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it

# Publication-ready GO enrichment plot
# Reuses precomputed results - fast execution

cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_1

# Create logs directory if needed
mkdir -p logs

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate r_chipseq_env

echo "Starting publication GO plot generation..."
echo "Date: $(date)"
echo ""

Rscript 02b_paper_go_plot.R

echo ""
echo "Completed: $(date)"
