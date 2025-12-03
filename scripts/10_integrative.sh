#!/bin/bash
#SBATCH --job-name=medip_peaks
#SBATCH --account=kubacki.michal
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/10_integrative.out
#SBATCH --error=logs/10_integrative.err


set -euo pipefail

echo "Start: $(date)"

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate r_chipseq_env


BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP"
cd $BASE_DIR


Rscript "./scripts/10_integrative.R"