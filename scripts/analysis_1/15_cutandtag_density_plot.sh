#!/bin/bash
#SBATCH --job-name=a1_15_cutandtag_density
#SBATCH --account=kubacki.michal
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/15_cutandtag_density.out
#SBATCH --error=logs/15_cutandtag_density.err

# Cut&Tag Density Plot Analysis
# This script generates heatmaps showing Cut&Tag signal enrichment around TES DEGs

set -e
set -u
set -o pipefail

echo "=================================================="
echo "Cut&Tag Density Plot Analysis"
echo "=================================================="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Start time: $(date)"
echo "Running on node: $(hostname)"
echo "Working directory: $(pwd)"
echo "=================================================="

# Set working directory
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_1

# Activate R environment with ChIPseeker and visualization tools
echo "Activating conda environment..."
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate r_chipseq_env

# Check required R packages
echo "Checking R environment..."
Rscript -e "library(GenomicRanges); library(rtracklayer); library(EnrichedHeatmap); library(ComplexHeatmap)" || {
    echo "ERROR: Required R packages not found"
    echo "Please ensure the r_chipseq_env has the following packages:"
    echo "  - GenomicRanges"
    echo "  - rtracklayer"
    echo "  - EnrichedHeatmap"
    echo "  - ComplexHeatmap"
    echo "  - TxDb.Hsapiens.UCSC.hg38.knownGene"
    echo "  - org.Hs.eg.db"
    exit 1
}

# Run the analysis
echo ""
echo "Running Cut&Tag density plot analysis..."
Rscript 15_cutandtag_density_plot.R

# Check if analysis completed successfully
if [ $? -eq 0 ]; then
    echo ""
    echo "=================================================="
    echo "Analysis completed successfully!"
    echo "End time: $(date)"
    echo "Output directory: output/density_plots/"
    echo "=================================================="
else
    echo ""
    echo "ERROR: Analysis failed!"
    exit 1
fi
