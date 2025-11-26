#!/bin/bash
#SBATCH --job-name=binding_classification_a3_
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --cpus-per-task=8
#SBATCH --time=04:00:00
#SBATCH --partition=workq
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --output=logs/advanced_01_binding_classification.out
#SBATCH --error=logs/advanced_01_binding_classification.err

################################################################################
# Phase 1.1: Advanced Binding Site Classification
# TES vs TEAD1 Comparative Study (excluding TESmut)
#
# This script classifies peaks based on binding patterns:
#   - TES_unique: Peaks only found in TES
#   - TEAD1_unique: Peaks only found in TEAD1
#   - Shared_TES_dominant: Overlap with TES signal > TEAD1
#   - Shared_TEAD1_dominant: Overlap with TEAD1 signal > TES
#   - Shared_equivalent: Overlap with similar signal intensity
#   - Shared_high: Both TF signals in top 25%
#
# Prerequisites:
#   - Completed peak calling (SRF_Eva_CUTandTAG/results/05_peaks_narrow/)
#   - BigWig files generated (SRF_Eva_CUTandTAG/results/06_bigwig/)
#
# Output:
#   - results/01_binding_classification/binding_site_classification.csv
#   - results/01_binding_classification/<category>.bed (for motif analysis)
#   - Visualization PDFs
################################################################################

echo "=================================================="
echo "Phase 1.1: Advanced Binding Site Classification"
echo "Start time: $(date)"
echo "=================================================="

# Change to project directory
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate analysis3_env

# Run R script
Rscript scripts/analysis_3/advanced_01_binding_classification.R

echo "=================================================="
echo "Analysis complete: $(date)"
echo "=================================================="
