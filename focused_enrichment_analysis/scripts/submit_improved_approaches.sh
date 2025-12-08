#!/bin/bash
#SBATCH --job-name=improved_approaches
#SBATCH --account=kubacki.michal
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=SRF_Eva_integrated_analysis/focused_enrichment_analysis/logs/improved_approaches.out
#SBATCH --error=SRF_Eva_integrated_analysis/focused_enrichment_analysis/logs/improved_approaches.err

echo "========================================================================"
echo "IMPROVED ENRICHMENT ANALYSIS APPROACHES"
echo "Running Approaches 6 (Corrected), 7, and 8"
echo "========================================================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Start time: $(date)"
echo "========================================================================"

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate r-bio

cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top

echo ""
echo "========================================================================"
echo "APPROACH 6 CORRECTED: Migration Gene-Focused (Hypothesis-Driven)"
echo "========================================================================"
Rscript SRF_Eva_integrated_analysis/focused_enrichment_analysis/scripts/approach6_migration_focused_CORRECTED.R
APPROACH6_STATUS=$?

if [ $APPROACH6_STATUS -eq 0 ]; then
    echo "✓ APPROACH 6 CORRECTED completed successfully!"
else
    echo "✗ APPROACH 6 CORRECTED failed!"
fi

echo ""
echo "========================================================================"
echo "APPROACH 7: Glioblastoma & Cancer Marker Analysis"
echo "========================================================================"
Rscript SRF_Eva_integrated_analysis/focused_enrichment_analysis/scripts/approach7_glioblastoma_markers.R
APPROACH7_STATUS=$?

if [ $APPROACH7_STATUS -eq 0 ]; then
    echo "✓ APPROACH 7 completed successfully!"
else
    echo "✗ APPROACH 7 failed!"
fi

echo ""
echo "========================================================================"
echo "APPROACH 8: Comparative Enrichment Analysis"
echo "========================================================================"
Rscript SRF_Eva_integrated_analysis/focused_enrichment_analysis/scripts/approach8_comparative_enrichment.R
APPROACH8_STATUS=$?

if [ $APPROACH8_STATUS -eq 0 ]; then
    echo "✓ APPROACH 8 completed successfully!"
else
    echo "✗ APPROACH 8 failed!"
fi

echo ""
echo "========================================================================"
echo "SUMMARY"
echo "========================================================================"
echo "Approach 6 Corrected: $([ $APPROACH6_STATUS -eq 0 ] && echo 'SUCCESS' || echo 'FAILED')"
echo "Approach 7 (Glioblastoma): $([ $APPROACH7_STATUS -eq 0 ] && echo 'SUCCESS' || echo 'FAILED')"
echo "Approach 8 (Comparative): $([ $APPROACH8_STATUS -eq 0 ] && echo 'SUCCESS' || echo 'FAILED')"
echo ""
echo "End time: $(date)"
echo "========================================================================"

# Exit with error if any approach failed
if [ $APPROACH6_STATUS -ne 0 ] || [ $APPROACH7_STATUS -ne 0 ] || [ $APPROACH8_STATUS -ne 0 ]; then
    echo "ERROR: One or more approaches failed!"
    exit 1
else
    echo "All approaches completed successfully!"
    exit 0
fi
