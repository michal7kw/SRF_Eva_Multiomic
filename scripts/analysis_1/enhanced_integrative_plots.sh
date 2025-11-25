#!/bin/bash
#SBATCH --job-name=run_enhanced_visualizations
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --time=6:00:00
#SBATCH --mem=84G
#SBATCH --cpus-per-task=32
#SBATCH --output=logs/run_enhanced_visualizations.out
#SBATCH --error=logs/run_enhanced_visualizations.err

set -e  # Exit on error

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top"
cd ${BASE_DIR}/SRF_Eva_integrated_analysis/scripts/analysis_1

echo "=========================================="
echo "ENHANCED VISUALIZATIONS MASTER SCRIPT"
echo "=========================================="
echo "Started: $(date)"
echo ""

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
  "${BASE_DIR}/SRF_Eva_integrated_analysis/output/true_gsea_analysis/fgsea_all_pathways.csv"
  "${BASE_DIR}/SRF_Eva_RNA/results/06_directional_go/upregulated_GO_BP.csv"
  "output/results/direct_targets/TES_direct_targets_all_genes.csv"
)

missing_files=0
for file in "${required_files[@]}"; do
  if [ -f "$file" ]; then
    echo "✓ Found: $file"
  else
    echo "✗ Missing: $file"
    missing_files=$((missing_files + 1))
  fi
done

if [ $missing_files -gt 0 ]; then
  echo ""
  echo "ERROR: $missing_files required files are missing!"
  echo ""
  echo "Please run standard analyses first:"
  echo "  sbatch true_gsea_analysis.sh"
  echo "  sbatch ${BASE_DIR}/SRF_Eva_RNA/scripts/run_directional_go.sh"
  echo "  sbatch final_integrative_analysis.sh"
  echo ""
  exit 1
fi

echo ""
echo "✓ All prerequisites met!"

# =============================================================================
# RUN ENHANCED VISUALIZATIONS
# =============================================================================

echo ""
echo "=== Running Enhanced Visualization Scripts ==="
echo ""

# 1. Enhanced GSEA
echo "1/3 Running enhanced GSEA visualizations..."
if Rscript enhanced_gsea_visualizations.R; then
  echo "   ✓ GSEA enhanced plots complete"
else
  echo "   ✗ GSEA visualization failed (continuing...)"
fi

echo ""

# 2. Enhanced Directional GO
echo "2/3 Running enhanced directional GO visualizations..."
if Rscript ${BASE_DIR}/SRF_Eva_RNA/scripts/enhanced_directional_go.R; then
  echo "   ✓ Directional GO enhanced plots complete"
else
  echo "   ✗ Directional GO visualization failed (continuing...)"
fi

echo ""

# 3. Enhanced Integrative Analysis
echo "3/3 Running enhanced integrative visualizations..."
if Rscript enhanced_integrative_plots.R; then
  echo "   ✓ Integrative enhanced plots complete"
else
  echo "   ✗ Integrative visualization failed (continuing...)"
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
echo "  output/true_gsea_analysis/enhanced_plots/"
echo ""
echo "Directional GO plots:"
echo "  ${BASE_DIR}/SRF_Eva_RNA/results/06_directional_go/enhanced_plots/"
echo "  ${BASE_DIR}/SRF_Eva_RNA/results/06_directional_go/interactive/"
echo ""
echo "Integrative plots:"
echo "  output/enhanced_integrative_plots/"
echo ""

echo "Generated files summary:"
echo "------------------------"

# Count generated files
gsea_count=$(find output/true_gsea_analysis/enhanced_plots/ -type f 2>/dev/null | wc -l)
go_count=$(find ${BASE_DIR}/SRF_Eva_RNA/results/06_directional_go/enhanced_plots/ -type f 2>/dev/null | wc -l)
int_count=$(find output/enhanced_integrative_plots/ -type f 2>/dev/null | wc -l)
interactive_count=$(find ${BASE_DIR}/SRF_Eva_RNA/results/06_directional_go/interactive/ -name "*.html" 2>/dev/null | wc -l)

echo "  GSEA enhanced plots: $gsea_count files"
echo "  Directional GO enhanced plots: $go_count files"
echo "  Interactive HTML plots: $interactive_count files"
echo "  Integrative enhanced plots: $int_count files"
echo ""

total_count=$((gsea_count + go_count + int_count + interactive_count))
echo "  Total new files: $total_count"
echo ""

echo "Next steps:"
echo "-----------"
echo "1. Review plots in output directories"
echo "2. Open interactive HTML in browser:"
echo "   firefox SRF_Eva_RNA/results/06_directional_go/interactive/interactive_volcano.html"
echo "3. Select figures for publication"
echo "4. Refer to ENHANCED_ANALYSIS_GUIDE.md for details"
echo ""
