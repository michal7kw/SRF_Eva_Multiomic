#!/bin/bash
#SBATCH --job-name=approach4_DOWN
#SBATCH --account=kubacki.michal
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/approach4_DOWN.out
#SBATCH --error=logs/approach4_DOWN.err

echo "========================================================================"
echo "Approach 4-DOWN: High-Confidence Peaks + Downregulated DEGs Only"
echo "========================================================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Start time: $(date)"
echo "Running on node: $(hostname)"
echo "========================================================================"
echo ""

# Load conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate r-bio

# Set working directory
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top

# Run analysis
echo "Running Approach 4-DOWN analysis..."
Rscript SRF_Eva_integrated_analysis/focused_enrichment_analysis/scripts/approach4_high_confidence_DOWN.R

# Check exit status
if [ $? -eq 0 ]; then
    echo ""
    echo "========================================================================"
    echo "Approach 4-DOWN completed successfully!"
    echo "End time: $(date)"
    echo "========================================================================"
else
    echo ""
    echo "========================================================================"
    echo "ERROR: Approach 4-DOWN failed!"
    echo "End time: $(date)"
    echo "========================================================================"
    exit 1
fi
