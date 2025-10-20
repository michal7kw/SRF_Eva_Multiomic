#!/bin/bash
#SBATCH --job-name=comprehensive_enrichment
#SBATCH --account=kubacki.michal
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=06:00:00
#SBATCH --output=logs/comprehensive_analysis.out
#SBATCH --error=logs/comprehensive_analysis.err

echo "========================================================================"
echo "Comprehensive Enrichment Analysis"
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

# Check if R script exists
if [ ! -f "SRF_Eva_integrated_analysis/focused_enrichment_analysis/scripts/comprehensive_enrichment_analysis.R" ]; then
    echo "ERROR: R script not found!"
    exit 1
fi

echo "Starting R analysis..."
echo ""

# Run the comprehensive R script
Rscript SRF_Eva_integrated_analysis/focused_enrichment_analysis/scripts/comprehensive_enrichment_analysis.R

# Check exit status
if [ $? -eq 0 ]; then
    echo ""
    echo "========================================================================"
    echo "Analysis completed successfully!"
    echo "End time: $(date)"
    echo "========================================================================"
    echo ""
    echo "Results are available in:"
    echo "  SRF_Eva_integrated_analysis/focused_enrichment_analysis/"
    echo ""
    echo "Individual approach folders:"
    echo "  - approach1_direct_targets/"
    echo "  - approach2_downregulated/"
    echo "  - approach3_promoter_peaks/"
    echo "  - approach4_high_confidence/"
    echo "  - approach5_diffbind/"
    echo "  - approach6_migration_focused/"
    echo "  - tier1_progression/"
    echo "  - tier2_validation/"
    echo ""
else
    echo ""
    echo "========================================================================"
    echo "ERROR: Analysis failed!"
    echo "End time: $(date)"
    echo "Check the error log for details."
    echo "========================================================================"
    exit 1
fi
