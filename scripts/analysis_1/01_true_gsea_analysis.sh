#!/bin/bash
#SBATCH --job-name=a1_01_true_gsea
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --time=2:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --output=logs/01_true_gsea.out
#SBATCH --error=logs/01_true_gsea.err

echo "=========================================="
echo "TRUE GSEA ANALYSIS (fgsea)"
echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Started: $(date)"
echo ""
echo "Method: Rank-based enrichment (proper GSEA)"
echo "All genes ranked by log2FoldChange"
echo ""

cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_1

echo "Loading R environment..."
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate r_chipseq_env

echo "Running TRUE GSEA analysis..."
Rscript 01_true_gsea_analysis.R

if [ $? -eq 0 ]; then
    echo ""
    echo "=========================================="
    echo "Analysis completed successfully!"
    echo "=========================================="
    echo "Results: output/true_gsea_analysis/"
else
    echo "ERROR: Analysis failed!"
    exit 1
fi
