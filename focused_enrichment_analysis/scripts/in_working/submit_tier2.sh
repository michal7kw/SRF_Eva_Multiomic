#!/bin/bash
#SBATCH --job-name=tier2_validation
#SBATCH --account=kubacki.michal
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/tier2.out
#SBATCH --error=logs/tier2.err

echo "========================================================================"
echo "TIER 2: Validation with Orthogonal Approach"
echo "========================================================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Start time: $(date)"
echo "========================================================================"

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate r-bio

cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top

Rscript SRF_Eva_integrated_analysis/focused_enrichment_analysis/scripts/tier2_validation.R

if [ $? -eq 0 ]; then
    echo ""
    echo "TIER 2 completed successfully!"
    echo "End time: $(date)"
else
    echo "ERROR: TIER 2 failed!"
    exit 1
fi
