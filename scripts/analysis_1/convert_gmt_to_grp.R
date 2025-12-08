#!/usr/bin/env Rscript
#
# CONVERT GMT TO GRP FORMAT
# Converts MSigDB .gmt files (all pathways in one file) to individual .grp files
#
# Usage: Rscript convert_gmt_to_grp.R <input.gmt> [output_directory]
#

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  cat("CONVERT MSigDB GMT TO GRP FORMAT\n")
  cat("================================\n\n")
  cat("Usage: Rscript convert_gmt_to_grp.R <input.gmt> [output_directory]\n\n")
  cat("Arguments:\n")
  cat("  input.gmt          Path to MSigDB .gmt file\n")
  cat("  output_directory   Output directory (default: GENE_SETS/)\n\n")
  cat("Example:\n")
  cat("  Rscript convert_gmt_to_grp.R h.all.v2024.1.Hs.symbols.gmt GENE_SETS/\n\n")
  quit(status = 1)
}

gmt_file <- args[1]
output_dir <- ifelse(length(args) >= 2, args[2], "GENE_SETS_selected") # "GENE_SETS")

# Check if GMT file exists
if (!file.exists(gmt_file)) {
  stop("ERROR: GMT file not found: ", gmt_file)
}

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("=============================================================================\n")
cat("  CONVERT MSigDB GMT TO GRP FORMAT\n")
cat("=============================================================================\n")
cat("Input file:", gmt_file, "\n")
cat("Output directory:", output_dir, "\n\n")

# Read GMT file
gmt_lines <- readLines(gmt_file)
cat(sprintf("Read %d gene sets from GMT file\n\n", length(gmt_lines)))

# Extract version from filename
version <- "v2024.1"
species <- "Hs"
if (grepl("v[0-9]+\\.[0-9]+", basename(gmt_file))) {
  version_match <- regmatches(
    basename(gmt_file),
    regexpr("v[0-9]+\\.[0-9]+", basename(gmt_file))
  )
  if (length(version_match) > 0) version <- version_match[1]
}
if (grepl("\\.(Hs|Mm)\\.", basename(gmt_file))) {
  species_match <- regmatches(
    basename(gmt_file),
    regexpr("(Hs|Mm)", basename(gmt_file))
  )
  if (length(species_match) > 0) species <- species_match[1]
}

# Process each line
cat("Converting gene sets...\n")
success_count <- 0
error_count <- 0

for (i in seq_along(gmt_lines)) {
  tryCatch(
    {
      # Parse GMT line
      parts <- strsplit(gmt_lines[i], "\t")[[1]]

      if (length(parts) < 3) {
        cat(sprintf(
          "  [%d/%d] SKIP: %s (too few columns)\n",
          i, length(gmt_lines), parts[1]
        ))
        error_count <- error_count + 1
        next
      }

      pathway_name <- parts[1]
      description <- parts[2]
      genes <- parts[-(1:2)]

      # Remove empty gene entries
      genes <- genes[nchar(genes) > 0]

      if (length(genes) == 0) {
        cat(sprintf(
          "  [%d/%d] SKIP: %s (no genes)\n",
          i, length(gmt_lines), pathway_name
        ))
        error_count <- error_count + 1
        next
      }

      # Create .grp filename
      grp_filename <- paste0(pathway_name, ".", version, ".", species, ".grp")
      grp_filepath <- file.path(output_dir, grp_filename)

      # Write .grp file
      writeLines(c(
        pathway_name,
        paste("#", description),
        genes
      ), grp_filepath)

      if (i %% 100 == 0) {
        cat(sprintf(
          "  [%d/%d] %s (%d genes)\n",
          i, length(gmt_lines), pathway_name, length(genes)
        ))
      }

      success_count <- success_count + 1
    },
    error = function(e) {
      cat(sprintf("  [%d/%d] ERROR: %s\n", i, length(gmt_lines), e$message))
      error_count <- error_count + 1
    }
  )
}

cat("\n=============================================================================\n")
cat("  CONVERSION COMPLETE\n")
cat("=============================================================================\n")
cat(sprintf("Successfully converted: %d gene sets\n", success_count))
cat(sprintf("Errors: %d\n", error_count))
cat(sprintf("Output directory: %s\n", output_dir))
cat("\n")

# List a few output files
output_files <- list.files(output_dir, pattern = "\\.grp$", full.names = FALSE)
cat(sprintf("Total .grp files in output directory: %d\n", length(output_files)))

if (length(output_files) > 0) {
  cat("\nFirst 5 output files:\n")
  for (f in head(output_files, 5)) {
    cat(sprintf("  - %s\n", f))
  }
}

cat("\nYou can now run GSEA with:\n")
cat("  sbatch scripts/analysis_1/msigdb_gsea_by_collection.sh\n\n")
