#!/bin/bash
#SBATCH --job-name=a1_18_gsea_cancer_pathways_IMPROVED
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --time=2:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --output=logs/18_gsea_cancer_pathways_IMPROVED.out
#SBATCH --error=logs/18_gsea_cancer_pathways_IMPROVED.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it

echo "=========================================="
echo "IMPROVED GSEA CANCER PATHWAYS ANALYSIS"
echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Running on node: $SLURM_NODELIST"
echo "Started: $(date)"
echo ""
echo "Improvements over original:"
echo "  - Expression direction analysis (UP/DOWN)"
echo "  - Angiogenesis & metabolism pathways"
echo "  - Direction-aware visualizations"
echo "  - Biological interpretation"
echo ""

# Navigate to project directory
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_1

# Load R environment
echo "Loading R environment..."
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate r_chipseq_env

# Run the improved GSEA analysis
echo "Running improved GSEA cancer pathways analysis..."
echo ""

Rscript 18_gsea_cancer_pathways_IMPROVED.R

# Check exit status
if [ $? -eq 0 ]; then
    echo ""
    echo "=========================================="
    echo "Analysis completed successfully!"
    echo "=========================================="
    echo "Completed: $(date)"
    echo ""
    echo "Results available in:"
    echo "  output/gsea_cancer_pathways_improved/"
    echo ""
    echo "Key outputs:"
    echo "  - *_directional.csv (with UP/DOWN annotation)"
    echo "  - 01_pathway_direction_bias.pdf"
    echo "  - 02_key_pathways_with_direction.pdf"
    echo "  - BIOLOGICAL_INTERPRETATION.txt"
else
    echo ""
    echo "ERROR: Analysis failed!"
    echo "Check error log: logs/gsea_improved.err"
    exit 1
fi
