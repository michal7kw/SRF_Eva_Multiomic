#!/bin/bash
#SBATCH --job-name=motif_analysis_a3_
#SBATCH --account=kubacki.michal
#SBATCH --mem=16GB
#SBATCH --cpus-per-task=8
#SBATCH --time=06:00:00
#SBATCH --partition=workq
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --output=logs/advanced_02_motif_analysis.out
#SBATCH --error=logs/advanced_02_motif_analysis.err

################################################################################
# Phase 1.2: Motif Analysis for Binding Categories
# TES vs TEAD1 Comparative Study
#
# Prerequisites:
#   - HOMER installed in genomics_env conda environment
#   - hg38 genome installed for HOMER via:
#     perl /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/genomics_env/share/homer/configureHomer.pl -install hg38
#   - Phase 1.1 binding classification completed
#
# This script performs de novo motif discovery on each binding category
# using HOMER's findMotifsGenome.pl
#
# Output:
#   - results/02_motif_analysis/<category>_motifs/ - HOMER output per category
#   - results/02_motif_analysis/motif_enrichment_heatmap.pdf
#   - results/02_motif_analysis/top_motifs_by_category.pdf
################################################################################

echo "=================================================="
echo "Phase 1.2: Motif Analysis for Binding Categories"
echo "Start time: $(date)"
echo "=================================================="

# Change to project directory
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top

# Set paths
INPUT_DIR="SRF_Eva_integrated_analysis/scripts/analysis_3/results/01_binding_classification"
OUTPUT_DIR="SRF_Eva_integrated_analysis/scripts/analysis_3/results/02_motif_analysis"
# Use HOMER's built-in hg38 genome (will be downloaded if not present)
GENOME="hg38"

# Create output directory
mkdir -p ${OUTPUT_DIR}

# Activate conda environment with HOMER
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate genomics_env  # This environment has HOMER installed

# Check if classification results exist
if [ ! -f "${INPUT_DIR}/binding_site_classification.csv" ]; then
    echo "ERROR: Binding classification results not found!"
    echo "Please run Phase 1.1 first: sbatch scripts/advanced_01_binding_classification.sh"
    exit 1
fi

echo ""
echo "Running HOMER motif analysis for each binding category..."
echo ""

# Categories to analyze
CATEGORIES=(
    "TES_unique"
    "TEAD1_unique"
    "Shared_high"
    "Shared_TES_dominant"
    "Shared_TEAD1_dominant"
    "Shared_equivalent"
)

# Run HOMER for each category
for CATEGORY in "${CATEGORIES[@]}"; do
    echo "----------------------------------------"
    echo "Processing category: ${CATEGORY}"
    echo "----------------------------------------"

    BED_FILE="${INPUT_DIR}/${CATEGORY}.bed"
    MOTIF_OUTPUT="${OUTPUT_DIR}/${CATEGORY}_motifs"

    if [ ! -f "${BED_FILE}" ]; then
        echo "WARNING: BED file not found for ${CATEGORY}, skipping..."
        continue
    fi

    # Count peaks
    N_PEAKS=$(wc -l < ${BED_FILE})
    echo "  Number of peaks: ${N_PEAKS}"

    if [ ${N_PEAKS} -lt 50 ]; then
        echo "  WARNING: Too few peaks for motif analysis, skipping..."
        continue
    fi

    # Run HOMER findMotifsGenome.pl
    echo "  Running HOMER motif discovery..."
    findMotifsGenome.pl \
        ${BED_FILE} \
        ${GENOME} \
        ${MOTIF_OUTPUT} \
        -size 200 \
        -mask \
        -p 8 \
        -len 8,10,12 \
        -S 25 \
        -noweight

    echo "  HOMER analysis complete for ${CATEGORY}"
    echo ""
done

echo ""
echo "=================================================="
echo "Creating motif comparison summary..."
echo "=================================================="

# Activate R environment for comparative analysis
conda activate analysis3_env

# Create R script for motif comparison
cat > ${OUTPUT_DIR}/compare_motifs.R << 'EOF'
#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)
library(RColorBrewer)

setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top")

OUTPUT_DIR <- "SRF_Eva_integrated_analysis/scripts/analysis_3/results/02_motif_analysis"

# Dynamically find categories that have results
all_dirs <- list.dirs(OUTPUT_DIR, full.names = FALSE, recursive = FALSE)
categories <- gsub("_motifs$", "", all_dirs[grepl("_motifs$", all_dirs)])

if (length(categories) == 0) {
  message("ERROR: No HOMER results found. Exiting.")
  quit(status = 1)
}

message("Found results for categories: ", paste(categories, collapse = ", "))

# Function to read HOMER known motifs
read_homer_motifs <- function(category) {
  motif_file <- file.path(OUTPUT_DIR, paste0(category, "_motifs"), "knownResults.txt")

  if (!file.exists(motif_file)) {
    message("  WARNING: ", motif_file, " not found")
    return(NULL)
  }

  # Read HOMER results
  motifs <- tryCatch({
    read.delim(motif_file, stringsAsFactors = FALSE)
  }, error = function(e) {
    message("  ERROR reading: ", motif_file)
    return(NULL)
  })

  if (is.null(motifs) || nrow(motifs) == 0) return(NULL)

  # Extract top 10 motifs
  motifs_top <- motifs[1:min(10, nrow(motifs)), ]
  motifs_top$category <- category
  motifs_top$rank <- 1:nrow(motifs_top)

  return(motifs_top)
}

message("Reading HOMER results for each category...")

# Read all categories
all_motifs <- lapply(categories, read_homer_motifs)
all_motifs <- do.call(rbind, all_motifs[!sapply(all_motifs, is.null)])

if (is.null(all_motifs) || nrow(all_motifs) == 0) {
  message("ERROR: No motif data found. Exiting.")
  quit(status = 1)
}

# Save combined results
write.csv(all_motifs, file.path(OUTPUT_DIR, "all_categories_top_motifs.csv"), row.names = FALSE)
message("Saved: all_categories_top_motifs.csv")

# Create enrichment heatmap for key TF families
message("Creating motif enrichment heatmap...")

# Extract key TF families
tf_families <- c("TEAD", "AP-1", "RUNX", "ETS", "NF-kB", "YAP", "CTCF", "E2F", "MYC", "TP53")

# Get unique categories that have data
valid_categories <- unique(all_motifs$category)

# Create enrichment matrix
enrichment_matrix <- matrix(0,
  nrow = length(tf_families),
  ncol = length(valid_categories)
)
rownames(enrichment_matrix) <- tf_families
colnames(enrichment_matrix) <- valid_categories

for (i in seq_along(valid_categories)) {
  cat_motifs <- all_motifs[all_motifs$category == valid_categories[i], ]

  for (j in seq_along(tf_families)) {
    # Find best match for TF family
    matches <- grep(tf_families[j], cat_motifs$Motif.Name, ignore.case = TRUE)

    if (length(matches) > 0) {
      # Use best p-value (Log.P.value is already negative log)
      best_match <- matches[which.min(cat_motifs$Log.P.value[matches])]
      enrichment_matrix[j, i] <- -cat_motifs$Log.P.value[best_match]
    }
  }
}

# Only create heatmap if we have variation in the data
if (length(unique(as.vector(enrichment_matrix))) > 1) {
  tryCatch({
    library(pheatmap)
    pdf(file.path(OUTPUT_DIR, "motif_enrichment_heatmap.pdf"), width = 10, height = 8)
    pheatmap(
      enrichment_matrix,
      color = colorRampPalette(c("white", "red", "darkred"))(100),
      cluster_rows = TRUE,
      cluster_cols = ifelse(ncol(enrichment_matrix) > 1, TRUE, FALSE),
      main = "TF Motif Enrichment Across Binding Categories",
      display_numbers = FALSE,
      fontsize = 10,
      angle_col = 45
    )
    dev.off()
    message("Saved: motif_enrichment_heatmap.pdf")
  }, error = function(e) {
    message("WARNING: Could not create heatmap: ", e$message)
  })
} else {
  message("WARNING: No variation in enrichment data, skipping heatmap")
}

# Create bar plot of top motifs per category
message("Creating bar plot...")

tryCatch({
  pdf(file.path(OUTPUT_DIR, "top_motifs_by_category.pdf"), width = 14, height = 10)

  # Select top 5 motifs per category
  top_motifs <- all_motifs %>%
    group_by(category) %>%
    slice_min(order_by = P.value, n = 5) %>%
    ungroup()

  p <- ggplot(top_motifs, aes(x = reorder(Motif.Name, -log10(P.value)),
                         y = -log10(P.value),
                         fill = category)) +
    geom_bar(stat = "identity") +
    facet_wrap(~ category, scales = "free", ncol = 2) +
    coord_flip() +
    scale_fill_brewer(palette = "Set2") +
    theme_classic() +
    theme(
      legend.position = "none",
      strip.text = element_text(face = "bold", size = 10)
    ) +
    labs(
      title = "Top 5 Enriched Motifs per Binding Category",
      x = "Motif",
      y = "-log10(P-value)"
    )
  print(p)
  dev.off()
  message("Saved: top_motifs_by_category.pdf")
}, error = function(e) {
  message("WARNING: Could not create bar plot: ", e$message)
})

message("Motif comparison analysis complete!")
EOF

# Run R comparison script
Rscript ${OUTPUT_DIR}/compare_motifs.R

echo ""
echo "=================================================="
echo "Motif analysis complete: $(date)"
echo "=================================================="
echo ""
echo "Results location: ${OUTPUT_DIR}"
echo ""
echo "Next steps:"
echo "  1. Review HOMER results in ${OUTPUT_DIR}/*_motifs/"
echo "  2. Check motif_enrichment_heatmap.pdf"
echo "  3. Proceed to Phase 1.3 (Signal Dynamics)"
