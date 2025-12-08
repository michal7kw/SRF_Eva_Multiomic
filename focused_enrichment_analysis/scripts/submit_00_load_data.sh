#!/bin/bash
#SBATCH --job-name=load_data
#SBATCH --account=kubacki.michal
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=logs/00_load_data.out
#SBATCH --error=logs/00_load_data.err

echo "========================================================================"
echo "Loading Data for Enrichment Analysis"
echo "========================================================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Start time: $(date)"
echo "Running on node: $(hostname)"
echo "Working directory: $(pwd)"
echo "========================================================================"
echo ""

# Load required modules and activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate r-bio

# Set working directory
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top

# Run the R script
echo "Starting data loading..."
echo ""
Rscript SRF_Eva_integrated_analysis/focused_enrichment_analysis/scripts/00_load_data.R

# Check exit status
if [ $? -eq 0 ]; then
    echo ""
    echo "========================================================================"
    echo "Data loading completed successfully!"
    echo "End time: $(date)"
    echo "========================================================================"
else
    echo ""
    echo "========================================================================"
    echo "ERROR: Data loading failed!"
    echo "End time: $(date)"
    echo "========================================================================"
    exit 1
fi
