#!/bin/bash
#SBATCH --job-name=06_meth_at_degs
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --time=1:00:00
#SBATCH --output=logs/06_methylation_at_regulated_genes.out
#SBATCH --error=logs/06_methylation_at_regulated_genes.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it

# Module 6: Methylation at Regulated Genes
# Tests whether DMRs occur at gene bodies of DEGs (indirect methylation hypothesis)

set -e

# Navigate to script directory
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_2

# Create logs directory
mkdir -p logs

echo "=== Starting Module 6: Methylation at Regulated Genes ==="
echo "Date: $(date)"
echo "Host: $(hostname)"
echo ""

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate seurat_full2

# Check input files
echo "Checking input files..."
RNA_FILE="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_RNA/results/05_deseq2/deseq2_results_TES_vs_GFP.txt"
DMR_FILE="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP/results/07_differential_MEDIPS/TES_vs_GFP_DMRs_FDR05.csv"

if [ ! -f "$RNA_FILE" ]; then
    echo "ERROR: RNA-seq results not found: $RNA_FILE"
    exit 1
fi

if [ ! -f "$DMR_FILE" ]; then
    echo "ERROR: DMR file not found: $DMR_FILE"
    exit 1
fi

echo "  RNA-seq: OK"
echo "  DMRs: OK"
echo ""

# Run analysis
echo "Running R script..."
Rscript 06_methylation_at_regulated_genes.R

echo ""
echo "=== Module 6 Complete ==="
echo "End time: $(date)"
