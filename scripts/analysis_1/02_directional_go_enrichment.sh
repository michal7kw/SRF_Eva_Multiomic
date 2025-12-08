#!/bin/bash
#SBATCH --job-name=a1_02_dir_go
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --time=2:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/02_directional_go.out
#SBATCH --error=logs/02_directional_go.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it

echo "=========================================="
echo "DIRECTIONAL GO ENRICHMENT (UP vs DOWN)"
echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Started: $(date)"
echo ""

cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_1

# Create logs directory if missing
mkdir -p logs

echo "Loading R environment..."
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate r_chipseq_env

echo "Running directional GO enrichment..."
Rscript 02_directional_go_enrichment.R

if [ $? -eq 0 ]; then
    echo ""
    echo "=========================================="
    echo "Analysis completed successfully!"
    echo "=========================================="
    echo "Results: output/directional_go_enrichment/"
else
    echo "ERROR: Analysis failed!"
    exit 1
fi
