#!/bin/bash
#SBATCH --job-name=a1_16_enhanced_integrative_plots
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --time=6:00:00
#SBATCH --mem=84G
#SBATCH --cpus-per-task=32
#SBATCH --output=logs/16_enhanced_integrative_plots.out
#SBATCH --error=logs/16_enhanced_integrative_plots.err

set -e  # Exit on error

# Working directory
WORK_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_1"
cd ${WORK_DIR}

echo "=========================================="
echo "ENHANCED VISUALIZATIONS MASTER SCRIPT"
echo "=========================================="
echo "Started: $(date)"
echo "Working directory: ${WORK_DIR}"
echo ""

# Create directories
mkdir -p logs
mkdir -p output/16_enhanced_gsea_visualizations
mkdir -p output/16_enhanced_integrative_plots
mkdir -p output/16_enhanced_directional_go

# Activate R environment
echo "Loading R environment..."
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate r-gsea-vis

# =============================================================================
# CHECK PREREQUISITES
# =============================================================================

echo ""
echo "=== Checking Prerequisites ==="
echo ""

# Check if standard analyses are complete
required_files=(
  "output/01_true_gsea_analysis/fgsea_all_pathways.csv"
  "output/02_directional_go_enrichment/upregulated_GO_BP.csv"
  "output/10_final_integrative_analysis/direct_targets/TES_direct_targets.csv"
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
  echo "WARNING: $missing_files required files are missing!"
  echo ""
  echo "Please run standard analyses first:"
  echo "  sbatch 01_true_gsea_analysis.sh"
  echo "  sbatch 02_directional_go_enrichment.sh"
  echo "  sbatch 10_final_integrative_analysis.sh"
  echo ""
  echo "Continuing anyway - scripts will skip missing data..."
fi

echo ""
echo "Prerequisites check complete!"

# =============================================================================
# RUN ENHANCED VISUALIZATIONS
# =============================================================================

echo ""
echo "=== Running Enhanced Visualization Scripts ==="
echo ""

# 1. Enhanced GSEA
echo "1/3 Running enhanced GSEA visualizations..."
if Rscript 16_enhanced_gsea_visualizations.R; then
  echo "   GSEA enhanced plots complete"
else
  echo "   GSEA visualization failed (continuing...)"
fi

echo ""

# 2. Enhanced Directional GO
echo "2/3 Running enhanced directional GO visualizations..."
if Rscript 16_enhanced_directional_go.R; then
  echo "   Directional GO enhanced plots complete"
else
  echo "   Directional GO visualization failed (continuing...)"
fi

echo ""

# 3. Enhanced Integrative Analysis
echo "3/3 Running enhanced integrative visualizations..."
if Rscript 16_enhanced_integrative_plots.R; then
  echo "   Integrative enhanced plots complete"
else
  echo "   Integrative visualization failed (continuing...)"
fi

echo ""

# =============================================================================
# SUMMARY
# =============================================================================

echo "=========================================="
echo "ENHANCED VISUALIZATIONS COMPLETE"
echo "=========================================="
echo "Completed: $(date)"
echo ""

echo "Output locations:"
echo "-----------------"
echo "GSEA plots:"
echo "  output/16_enhanced_gsea_visualizations/"
echo ""
echo "Directional GO plots:"
echo "  output/16_enhanced_directional_go/"
echo ""
echo "Integrative plots:"
echo "  output/16_enhanced_integrative_plots/"
echo ""

echo "Generated files summary:"
echo "------------------------"

# Count generated files
gsea_count=$(find output/16_enhanced_gsea_visualizations/ -type f 2>/dev/null | wc -l)
go_count=$(find output/16_enhanced_directional_go/ -type f 2>/dev/null | wc -l)
int_count=$(find output/16_enhanced_integrative_plots/ -type f 2>/dev/null | wc -l)

echo "  GSEA enhanced plots: $gsea_count files"
echo "  Directional GO enhanced plots: $go_count files"
echo "  Integrative enhanced plots: $int_count files"
echo ""

total_count=$((gsea_count + go_count + int_count))
echo "  Total new files: $total_count"
echo ""

echo "Next steps:"
echo "-----------"
echo "1. Review plots in output directories"
echo "2. Select figures for publication"
echo ""
