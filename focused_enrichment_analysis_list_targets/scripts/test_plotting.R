#!/usr/bin/env Rscript

# Test script to verify plotting updates
cat("Testing updated plotting function...\n")

# Load helper functions
source("SRF_Eva_integrated_analysis/focused_enrichment_analysis/scripts/00_helper_functions.R")

# Check if the function contains PNG saving code
func_text <- deparse(save_enrichment_results)
has_png <- any(grepl("png", func_text, ignore.case = TRUE))
has_height_calc <- any(grepl("plot_height.*max.*14", func_text))

cat("\nFunction checks:\n")
cat("Contains PNG saving code:", has_png, "\n")
cat("Contains updated height calculation:", has_height_calc, "\n")

# Print relevant lines
cat("\nSearching for ggsave calls in function:\n")
png_lines <- grep("ggsave.*png|Save as PNG", func_text, ignore.case = TRUE, value = TRUE)
for (line in png_lines) {
  cat(line, "\n")
}

cat("\nTest complete.\n")
