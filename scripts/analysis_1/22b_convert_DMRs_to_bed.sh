#!/bin/bash
#SBATCH --job-name=a1_22b_convert_DMRs_to_bed
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --time=2:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --output=logs/22b_convert_DMRs_to_bed.out
#SBATCH --error=logs/22b_convert_DMRs_to_bed.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it

# Set working directory
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_1

# Activate environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate r_chipseq_env

# Run analysis
echo "Running Peak-DMR mapping analysis..."
Rscript convert_DMRs_to_bed.R