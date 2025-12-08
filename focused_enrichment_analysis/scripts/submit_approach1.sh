#!/bin/bash
#SBATCH --job-name=approach1_direct
#SBATCH --account=kubacki.michal
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/approach1.out
#SBATCH --error=logs/approach1.err

echo "========================================================================"
echo "APPROACH 1: Direct Targets Only (Baseline)"
echo "========================================================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Start time: $(date)"
echo "Running on node: $(hostname)"
echo "========================================================================"
echo ""

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate r-bio

cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top

Rscript SRF_Eva_integrated_analysis/focused_enrichment_analysis/scripts/approach1_direct_targets.R

if [ $? -eq 0 ]; then
    echo ""
    echo "========================================================================"
    echo "APPROACH 1 completed successfully!"
    echo "End time: $(date)"
    echo "========================================================================"
else
    echo ""
    echo "========================================================================"
    echo "ERROR: APPROACH 1 failed!"
    echo "End time: $(date)"
    echo "========================================================================"
    exit 1
fi
