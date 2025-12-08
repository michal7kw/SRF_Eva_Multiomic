#!/bin/bash
#SBATCH --job-name=tier1_progression
#SBATCH --account=kubacki.michal
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/tier1.out
#SBATCH --error=logs/tier1.err

echo "========================================================================"
echo "TIER 1: Show Progression of Specificity"
echo "========================================================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Start time: $(date)"
echo "========================================================================"

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate r-bio

cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top

Rscript SRF_Eva_integrated_analysis/focused_enrichment_analysis/scripts/tier1_progression.R

if [ $? -eq 0 ]; then
    echo ""
    echo "TIER 1 completed successfully!"
    echo "End time: $(date)"
else
    echo "ERROR: TIER 1 failed!"
    exit 1
fi
