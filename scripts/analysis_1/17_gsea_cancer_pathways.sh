#!/bin/bash
#SBATCH --job-name=a1_17_gsea_cancer_pathways
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --time=2:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --output=logs/17_gsea_cancer_pathways.out
#SBATCH --error=logs/17_gsea_cancer_pathways.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it

echo "=========================================="
echo "GSEA CANCER PATHWAYS ANALYSIS"
echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Running on node: $SLURM_NODELIST"
echo "Started: $(date)"
echo ""
echo "Analysis focus: Cell death, apoptosis, migration, proliferation"
echo "Target groups: TES-only, TEAD1-only, Shared TES/TEAD1"
echo ""

# Navigate to project directory
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_1

# Load R environment
echo "Loading R environment..."
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate r_chipseq_env

# Run the GSEA cancer pathways analysis
echo "Running GSEA cancer pathways analysis..."
echo "Script: 17_gsea_cancer_pathways.R"
echo ""

Rscript 17_gsea_cancer_pathways.R

# Check exit status
if [ $? -eq 0 ]; then
    echo ""
    echo "=========================================="
    echo "Analysis completed successfully!"
    echo "=========================================="
    echo "Completed: $(date)"
    echo "Job ID: $SLURM_JOB_ID"
    echo ""
    echo "Results available in:"
    echo "  output/gsea_cancer_pathways/"
    echo ""
    echo "Key outputs:"
    echo "  - TES_only_cancer_pathways.csv"
    echo "  - TEAD1_only_cancer_pathways.csv"
    echo "  - Shared_cancer_pathways.csv"
    echo "  - Multiple visualization PDFs"
    echo "  - GSEA_CANCER_PATHWAYS_SUMMARY.txt"
else
    echo ""
    echo "=========================================="
    echo "ERROR: Analysis failed!"
    echo "=========================================="
    echo "Check error log: logs/gsea_cancer_pathways.err"
    echo "Job ID: $SLURM_JOB_ID"
    exit 1
fi
