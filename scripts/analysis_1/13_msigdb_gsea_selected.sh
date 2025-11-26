#!/bin/bash
#SBATCH --job-name=a1_13_msigdb_gsea_selected
#SBATCH --account=kubacki.michal
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/13_msigdb_gsea_selected.out
#SBATCH --error=logs/13_msigdb_gsea_selected.err

# ============================================================================
# MSIGDB GSEA ANALYSIS - SELECTED GENE SETS (Hallmark + C6)
# ============================================================================
# Description: Run fgsea with selected MSigDB gene sets (GENE_SETS_selected)
# Author: Michal Kubacki
# Date: 2025-10-13
# ============================================================================

echo "============================================================================"
echo "  MSIGDB GSEA ANALYSIS - SELECTED GENE SETS"
echo "============================================================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Job name: $SLURM_JOB_NAME"
echo "Node: $SLURM_NODELIST"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Memory: 32GB"
echo "Time limit: 2 hours"
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
echo "Starting GSEA analysis with selected gene sets (Hallmark + C6)..."
echo "============================================================================"
echo ""

Rscript 13_msigdb_gsea_selected.R

EXIT_CODE=$?

echo ""
echo "============================================================================"
if [ $EXIT_CODE -eq 0 ]; then
    echo "  ANALYSIS COMPLETED SUCCESSFULLY"
    echo "  Results saved in: output/13_msigdb_gsea_selected/"
    echo "  Gene sets used: GENE_SETS_selected/ (Hallmark + C6 Oncogenic)"
else
    echo "  ANALYSIS FAILED (exit code: $EXIT_CODE)"
fi
echo "============================================================================"
echo "Ended: $(date)"
echo ""

exit $EXIT_CODE
