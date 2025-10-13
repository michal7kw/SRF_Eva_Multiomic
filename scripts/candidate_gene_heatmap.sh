#!/bin/bash
#SBATCH --job-name=cand_hm
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --time=1:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/candidate_heatmap.out
#SBATCH --error=logs/candidate_heatmap.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it

echo "=========================================="
echo "CANDIDATE GENE HEATMAP ANALYSIS"
echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Started: $(date)"
echo ""
echo "Analyzing 223 genes from TES_degs.txt"
echo ""

cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis

echo "Loading R environment..."
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate r_chipseq_env

echo "Running candidate gene analysis..."
Rscript scripts/candidate_gene_heatmap.R

if [ $? -eq 0 ]; then
    echo ""
    echo "=========================================="
    echo "Analysis completed successfully!"
    echo "=========================================="
    echo "Results: output/candidate_gene_heatmap/"
else
    echo "ERROR: Analysis failed!"
    exit 1
fi
