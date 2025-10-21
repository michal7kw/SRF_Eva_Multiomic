#!/bin/bash
#SBATCH --job-name=yaptaz_enrichment
#SBATCH --account=kubacki.michal
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/yaptaz_analysis.out
#SBATCH --error=logs/yaptaz_analysis.err

echo "========================================================================"
echo "YAP/TAZ TARGET ENRICHMENT ANALYSIS"
echo "Using known YAP/TAZ targets from TES_degs.txt"
echo "========================================================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Start time: $(date)"
echo "========================================================================"

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate r-bio

cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top

Rscript SRF_Eva_integrated_analysis/focused_enrichment_analysis_list_targets/scripts/yaptaz_enrichment_analysis.R

if [ $? -eq 0 ]; then
    echo ""
    echo "YAP/TAZ target analysis completed successfully!"
    echo "End time: $(date)"
else
    echo "ERROR: YAP/TAZ target analysis failed!"
    exit 1
fi
