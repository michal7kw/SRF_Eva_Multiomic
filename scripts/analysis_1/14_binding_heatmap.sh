#!/bin/bash
#SBATCH --job-name=a1_14_binding_heatmap
#SBATCH --account=kubacki.michal
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/14_binding_heatmap.out
#SBATCH --error=logs/14_binding_heatmap.err

# Cut&Tag Binding Heatmap Analysis
# Creates a clustered heatmap showing promoter binding signal for TES DEGs

set -e
set -u
set -o pipefail

echo "=================================================="
echo "Cut&Tag Binding Heatmap Analysis"
echo "=================================================="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Start time: $(date)"
echo "Running on node: $(hostname)"
echo "Working directory: $(pwd)"
echo "=================================================="

# Set working directory
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_1

# Activate R environment with ComplexHeatmap and genomics tools
echo "Activating conda environment..."
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate r_chipseq_env

# Check required R packages
echo "Checking R environment..."
Rscript -e "library(GenomicRanges); library(rtracklayer); library(ComplexHeatmap); library(circlize)" || {
    echo "ERROR: Required R packages not found"
    echo "Please ensure the r_chipseq_env has the following packages:"
    echo "  - GenomicRanges"
    echo "  - rtracklayer"
    echo "  - ComplexHeatmap"
    echo "  - circlize"
    echo "  - TxDb.Hsapiens.UCSC.hg38.knownGene"
    echo "  - org.Hs.eg.db"
    echo "  - dplyr"
    exit 1
}

# Run the analysis
echo ""
echo "Running Cut&Tag binding heatmap analysis..."
Rscript 14_binding_heatmap.R

# Check if analysis completed successfully
if [ $? -eq 0 ]; then
    echo ""
    echo "=================================================="
    echo "Analysis completed successfully!"
    echo "End time: $(date)"
    echo "Output directory: output/binding_heatmap/"
    echo ""
    echo "Generated files:"
    echo "  - promoter_binding_heatmap.pdf"
    echo "  - promoter_binding_heatmap_with_fc.pdf"
    echo "  - promoter_binding_signals.txt"
    echo "=================================================="
else
    echo ""
    echo "ERROR: Analysis failed!"
    exit 1
fi
