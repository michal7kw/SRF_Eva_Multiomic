#!/bin/bash
#SBATCH --job-name=a1_22_peak_DMR_mapping
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --time=2:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --output=logs/22_peak_DMR_mapping.out
#SBATCH --error=logs/22_peak_DMR_mapping.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it

echo "=========================================="
echo "PEAK-DMR MAPPING ANALYSIS"
echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Started: $(date)"
echo ""
echo "Task: Map TES/TEAD1 peaks to DMRs"
echo "      - Find overlaps and distances"
echo "      - Annotate to genes"
echo "      - Create tables and visualizations"
echo ""

# Set working directory
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_1

# Create output directories
mkdir -p output/22_peak_DMR_mapping/tables
mkdir -p output/22_peak_DMR_mapping/plots
mkdir -p logs

# Activate environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate r_chipseq_env

# Run analysis
echo "Running Peak-DMR mapping analysis..."
Rscript 22_peak_DMR_mapping.R

if [ $? -eq 0 ]; then
    echo ""
    echo "=========================================="
    echo "ANALYSIS COMPLETED SUCCESSFULLY"
    echo "=========================================="
    echo "Completed: $(date)"
    echo ""
    echo "Output directory: output/22_peak_DMR_mapping/"
    echo ""
    echo "Generated tables:"
    ls -lh output/22_peak_DMR_mapping/tables/
    echo ""
    echo "Generated plots:"
    ls -lh output/22_peak_DMR_mapping/plots/
else
    echo "ERROR: Analysis failed!"
    exit 1
fi
