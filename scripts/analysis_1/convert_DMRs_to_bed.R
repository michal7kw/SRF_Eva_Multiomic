#!/usr/bin/env Rscript
# Convert DMRs_with_peak_info.csv to BED format for IGV visualization
# Creates multiple BED files with different coloring schemes

library(dplyr)
library(readr)

# CRITICAL: Disable scientific notation for genomic coordinates
options(scipen = 999)

# Input/output paths
input_file <- "output/22_peak_DMR_mapping/tables/DMRs_with_peak_info.csv"
output_dir <- "output/22_peak_DMR_mapping/igv_tracks"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Read DMR data
dmrs <- read_csv(input_file, show_col_types = FALSE)
cat(sprintf("Loaded %d DMRs\n", nrow(dmrs)))

# Ensure coordinates are integers (not scientific notation)
dmrs <- dmrs %>%
  mutate(
    start = as.integer(start),
    end = as.integer(end)
  )

# ============================================================================
# BED9 format: chr, start, end, name, score, strand, thickStart, thickEnd, itemRgb
# ============================================================================

# Helper function to write BED with proper formatting
write_bed <- function(df, filepath) {
  # Format coordinates as plain integers
  df <- df %>%
    mutate(
      start = format(start, scientific = FALSE, trim = TRUE),
      end = format(end, scientific = FALSE, trim = TRUE),
      thickStart = format(thickStart, scientific = FALSE, trim = TRUE),
      thickEnd = format(thickEnd, scientific = FALSE, trim = TRUE)
    )
  write.table(df, filepath, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

# 1. All DMRs colored by direction (Hyper=red, Hypo=blue)
bed_direction <- dmrs %>%
  mutate(
    name = paste0(dmr_id, "_", gene, "_", direction),
    score = as.integer(pmin(1000, abs(logFC) * 100)),  # scale logFC to 0-1000
    strand = ".",
    thickStart = start,
    thickEnd = end,
    color = ifelse(direction == "Hyper", "255,0,0", "0,0,255")  # Red/Blue
  ) %>%
  select(chr, start, end, name, score, strand, thickStart, thickEnd, color)

write_bed(bed_direction, file.path(output_dir, "DMRs_by_direction.bed"))

# 2. DMRs colored by TES peak presence (with peak=green, no peak=gray)
bed_TES <- dmrs %>%
  mutate(
    name = paste0(dmr_id, "_", gene, ifelse(has_TES_peak, "_TES+", "_TES-")),
    score = as.integer(pmin(1000, abs(logFC) * 100)),
    strand = ".",
    thickStart = start,
    thickEnd = end,
    color = ifelse(has_TES_peak, "0,180,0", "128,128,128")  # Green/Gray
  ) %>%
  select(chr, start, end, name, score, strand, thickStart, thickEnd, color)

write_bed(bed_TES, file.path(output_dir, "DMRs_by_TES_peak.bed"))

# 3. DMRs colored by TEAD1 peak presence
bed_TEAD1 <- dmrs %>%
  mutate(
    name = paste0(dmr_id, "_", gene, ifelse(has_TEAD1_peak, "_TEAD1+", "_TEAD1-")),
    score = as.integer(pmin(1000, abs(logFC) * 100)),
    strand = ".",
    thickStart = start,
    thickEnd = end,
    color = ifelse(has_TEAD1_peak, "180,0,180", "128,128,128")  # Purple/Gray
  ) %>%
  select(chr, start, end, name, score, strand, thickStart, thickEnd, color)

write_bed(bed_TEAD1, file.path(output_dir, "DMRs_by_TEAD1_peak.bed"))

# 4. Combined: 4-category coloring (both peaks, TES only, TEAD1 only, neither)
bed_combined <- dmrs %>%
  mutate(
    category = case_when(
      has_TES_peak & has_TEAD1_peak ~ "Both",
      has_TES_peak ~ "TES_only",
      has_TEAD1_peak ~ "TEAD1_only",
      TRUE ~ "Neither"
    ),
    name = paste0(dmr_id, "_", gene, "_", category),
    score = as.integer(pmin(1000, abs(logFC) * 100)),
    strand = ".",
    thickStart = start,
    thickEnd = end,
    color = case_when(
      category == "Both" ~ "255,165,0",      # Orange
      category == "TES_only" ~ "0,180,0",    # Green
      category == "TEAD1_only" ~ "180,0,180", # Purple
      TRUE ~ "128,128,128"                    # Gray
    )
  ) %>%
  select(chr, start, end, name, score, strand, thickStart, thickEnd, color)

write_bed(bed_combined, file.path(output_dir, "DMRs_by_peak_category.bed"))

# 5. Simple BED4 format (for basic viewing)
bed_simple <- dmrs %>%
  mutate(
    start = format(start, scientific = FALSE, trim = TRUE),
    end = format(end, scientific = FALSE, trim = TRUE)
  ) %>%
  select(chr, start, end, dmr_id)

write.table(bed_simple, file.path(output_dir, "DMRs_simple.bed"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Summary
cat("\n=== BED files created in", output_dir, "===\n")
cat("1. DMRs_by_direction.bed     - Hyper (red) vs Hypo (blue)\n")
cat("2. DMRs_by_TES_peak.bed      - With TES peak (green) vs without (gray)\n")
cat("3. DMRs_by_TEAD1_peak.bed    - With TEAD1 peak (purple) vs without (gray)\n")
cat("4. DMRs_by_peak_category.bed - Both (orange), TES only (green), TEAD1 only (purple), Neither (gray)\n")
cat("5. DMRs_simple.bed           - Basic BED4 format\n")
cat("\nTo view in IGV:\n")
cat("1. Open IGV and load hg38 genome\n")
cat("2. File > Load from File > select .bed file\n")
cat("3. Right-click track > Set Track Color > choose 'By strand or item RGB'\n")
