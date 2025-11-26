#!/bin/bash
#SBATCH --job-name=a1_12_msigdb_gsea_by_collection
#SBATCH --account=kubacki.michal
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=logs/12_msigdb_by_collection.out
#SBATCH --error=logs/12_msigdb_by_collection.err

# ============================================================================
# MSIGDB GSEA ANALYSIS - ORGANIZED BY COLLECTION
# ============================================================================
# Description: Run fgsea with MSigDB gene sets, save results per collection
# Author: Michal Kubacki
# Date: 2025-10-13
# ============================================================================

echo "============================================================================"
echo "  MSIGDB GSEA ANALYSIS - BY COLLECTION"
echo "============================================================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Job name: $SLURM_JOB_NAME"
echo "Node: $SLURM_NODELIST"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Memory: 64GB"
echo "Time limit: 8 hours"
echo "Started: $(date)"
echo "============================================================================"
echo ""

# Change to working directory
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_1

# Activate conda environment with R packages
echo "Loading R environment..."
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate r_chipseq_env

# Check packages
echo ""
echo "Checking R package availability..."
Rscript -e "packageVersion('fgsea')"
Rscript -e "packageVersion('dplyr')"
Rscript -e "packageVersion('ggplot2')"
echo ""

# Run the analysis
echo "Starting GSEA analysis with collection-based organization..."
echo "============================================================================"
echo ""

Rscript 12_msigdb_gsea_by_collection.R

EXIT_CODE=$?

echo ""
echo "============================================================================"
if [ $EXIT_CODE -eq 0 ]; then
    echo "  ANALYSIS COMPLETED SUCCESSFULLY"
    echo "  Results saved in: output/12_msigdb_by_collection/"
    echo "  Each MSigDB collection has its own subdirectory"
else
    echo "  ANALYSIS FAILED (exit code: $EXIT_CODE)"
fi
echo "============================================================================"
echo "Ended: $(date)"
echo ""

exit $EXIT_CODE
