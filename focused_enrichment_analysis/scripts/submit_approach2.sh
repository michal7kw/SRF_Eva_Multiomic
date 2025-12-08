#!/bin/bash
#SBATCH --job-name=approach2_down
#SBATCH --account=kubacki.michal
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/approach2.out
#SBATCH --error=logs/approach2.err

echo "========================================================================"
echo "APPROACH 2: Downregulated Direct Targets"
echo "========================================================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Start time: $(date)"
echo "========================================================================"

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate r-bio

cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top

Rscript SRF_Eva_integrated_analysis/focused_enrichment_analysis/scripts/approach2_downregulated.R

if [ $? -eq 0 ]; then
    echo ""
    echo "APPROACH 2 completed successfully!"
    echo "End time: $(date)"
else
    echo "ERROR: APPROACH 2 failed!"
    exit 1
fi
