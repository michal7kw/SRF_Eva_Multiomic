#!/bin/bash
#SBATCH --job-name=a1_23_genomicplot
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --time=4:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --output=logs/23_genomicplot_metagene.out
#SBATCH --error=logs/23_genomicplot_metagene.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it

echo "=========================================="
echo "GENOMICPLOT METAGENE PROFILES"
echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Started: $(date)"
echo ""
echo "Task: Create metagene profiles with gene body structure"
echo "      (Upstream -> 5'UTR -> CDS -> 3'UTR -> Downstream)"
echo ""

# Set working directory
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_1

# Create output directory
mkdir -p output/23_genomicplot_metagene
mkdir -p logs

# Activate environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate r_chipseq_env

# Run analysis
echo "Running GenomicPlot metagene analysis..."
Rscript 23_genomicplot_metagene.R

if [ $? -eq 0 ]; then
    echo ""
    echo "=========================================="
    echo "ANALYSIS COMPLETED SUCCESSFULLY"
    echo "=========================================="
    echo "Completed: $(date)"
    echo ""
    echo "Output directory: output/23_genomicplot_metagene/"
    ls -lh output/23_genomicplot_metagene/*.pdf 2>/dev/null || echo "Check output directory for results"
else
    echo "ERROR: Analysis failed!"
    exit 1
fi
