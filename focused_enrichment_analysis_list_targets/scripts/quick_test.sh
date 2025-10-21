#!/bin/bash

echo "Quick test of plotting function..."
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate r-bio

# Run just the first part to test
Rscript -e "
source('SRF_Eva_integrated_analysis/focused_enrichment_analysis/scripts/00_helper_functions.R')
cat('Function loaded with', length(deparse(save_enrichment_results)), 'lines\n')

# Check for PNG commands
func_text <- paste(deparse(save_enrichment_results), collapse=' ')
has_png <- grepl('png', func_text, ignore.case=TRUE)
cat('Has PNG output:', has_png, '\n')

# Check for debug output
has_debug <- grepl('Generating plots for', func_text)
cat('Has debug output:', has_debug, '\n')
"

echo "Test complete"
