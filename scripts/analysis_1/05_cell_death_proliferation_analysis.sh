#!/bin/bash

#SBATCH --job-name=a1_05_cell_death_analysis
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --mem=16GB
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/05_cell_death_proliferation.err"
#SBATCH --output="logs/05_cell_death_proliferation.out"

# SCRIPT: 05_cell_death_proliferation_analysis.sh
# PURPOSE: Execute cell death and proliferation pathway analysis
#
# DESCRIPTION:
# This script runs specialized pathway analysis focused on cell death and
# proliferation mechanisms. This is particularly relevant for TES/TEAD1
# research in glioblastoma, where these transcription factors play crucial
# roles in cell survival, apoptosis, and proliferation control.
#
# KEY ANALYSES:
# 1. Focused GO enrichment for cell death and proliferation
# 2. KEGG pathway analysis for cancer-related pathways
# 3. Reactome pathway analysis
# 4. Comparative analysis between TES, TESmut, and TEAD1
# 5. Gene overlap and functional similarity analysis

# Set working directory
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_1

# Create directories
mkdir -p logs
mkdir -p output/05_cell_death_proliferation

echo "=========================================="
echo "Cell Death & Proliferation Analysis"
echo "Focus: TES/TEAD1 in Glioblastoma"
echo "Date: $(date)"
echo "=========================================="
echo ""
echo "NOTE: This analysis uses annotated peaks from Cut&Tag pipeline."
echo ""

# Set up conda environment with required tools
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate cell_death_prolif

# Input files (from CUTandTAG pipeline)
CUTANDTAG_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_CUTandTAG"
ANNOTATION_DIR="${CUTANDTAG_DIR}/results/07_analysis_narrow"

# Check if required input files exist
echo "Checking input files..."
echo "----------------------"

required_files=(
    "${ANNOTATION_DIR}/TES_peaks_annotated.csv"
    "${ANNOTATION_DIR}/TESmut_peaks_annotated.csv"
    "${ANNOTATION_DIR}/TEAD1_peaks_annotated.csv"
)

missing_files=0
for file in "${required_files[@]}"; do
    if [ -f "$file" ]; then
        echo "Found: $file"
    else
        echo "Missing: $file"
        missing_files=$((missing_files + 1))
    fi
done

if [ $missing_files -gt 0 ]; then
    echo ""
    echo "Warning: $missing_files required files are missing"
    echo "Please run the peak annotation script first."
    echo ""
    echo "Attempting to run analysis with available files..."
fi

# Run the cell death and proliferation analysis
echo ""
echo "Running cell death and proliferation analysis..."
echo "------------------------------------------------"

if [ -f "05_cell_death_proliferation_analysis.R" ]; then
    echo "Executing specialized pathway analysis..."

    Rscript 05_cell_death_proliferation_analysis.R

    if [ $? -eq 0 ]; then
        echo "Cell death and proliferation analysis completed successfully"
    else
        echo "Analysis encountered errors - check log files"
        exit 1
    fi
else
    echo "Error: 05_cell_death_proliferation_analysis.R script not found"
    exit 1
fi

echo ""
echo "=========================================="
echo "Analysis Summary"
echo "=========================================="

# Check output directory
OUTPUT_DIR="output/cell_death_proliferation"
if [ -d "$OUTPUT_DIR" ]; then
    echo ""
    echo "Generated files:"
    echo "---------------"

    # Count different types of outputs
    pdf_count=$(find "$OUTPUT_DIR" -name "*.pdf" 2>/dev/null | wc -l)
    png_count=$(find "$OUTPUT_DIR" -name "*.png" 2>/dev/null | wc -l)
    csv_count=$(find "$OUTPUT_DIR" -name "*.csv" 2>/dev/null | wc -l)

    echo "PDF plots: $pdf_count"
    echo "PNG plots: $png_count"
    echo "CSV results: $csv_count"

    echo ""
    echo "Key output files:"
    echo "----------------"

    # List important files
    if [ -f "$OUTPUT_DIR/analysis_summary.csv" ]; then
        echo "Analysis summary: $OUTPUT_DIR/analysis_summary.csv"
    fi

    if [ -f "$OUTPUT_DIR/gene_overlap_stats.csv" ]; then
        echo "Gene overlap statistics: $OUTPUT_DIR/gene_overlap_stats.csv"
    fi

    echo ""
    echo "Biological insights to examine:"
    echo "------------------------------"
    echo "- Cell death pathways regulated by TES vs TEAD1"
    echo "- Proliferation mechanisms specific to each factor"
    echo "- Apoptosis regulation differences between TES and TESmut"
    echo "- Cancer-related pathway enrichments"
    echo "- Gene overlap patterns between conditions"

else
    echo "Output directory not created - analysis may have failed"
fi

echo ""
echo "=========================================="
echo "Analysis completed: $(date)"
echo "Results location: $OUTPUT_DIR"
echo "=========================================="

# Quick summary if analysis summary exists
if [ -f "$OUTPUT_DIR/analysis_summary.csv" ]; then
    echo ""
    echo "Quick Results Summary:"
    echo "---------------------"
    cat "$OUTPUT_DIR/analysis_summary.csv"
fi
