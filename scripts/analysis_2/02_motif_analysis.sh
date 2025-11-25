#!/bin/bash
#SBATCH --job-name=02_motif_analysis
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --output=logs/02_motif_analysis.out
#SBATCH --error=logs/02_motif_analysis.err
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it

# Module 2: Motif Enrichment Analysis
# Requires HOMER (findMotifsGenome.pl)
# FIXED: Absolute paths, proper input checking

set -e

# Navigate to script directory
SCRIPT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_2"
cd "$SCRIPT_DIR"

echo "=== Module 2: Motif Analysis ==="
echo "Working directory: $(pwd)"
echo "Date: $(date)"

# Activate environment containing HOMER
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate genomics_env

# Configuration
GENOME="hg38"
INPUT_DIR="${SCRIPT_DIR}/results/01_peak_classification"
OUTPUT_DIR="${SCRIPT_DIR}/results/02_motif_enrichment"

mkdir -p "$OUTPUT_DIR"

# Check if input files exist
echo ""
echo "Checking input files..."
for bed_file in "${INPUT_DIR}/Shared_TES_TEAD1.bed" "${INPUT_DIR}/TES_Unique.bed" "${INPUT_DIR}/TEAD1_Unique.bed"; do
    if [[ ! -f "$bed_file" ]]; then
        echo "ERROR: Input file not found: $bed_file"
        echo "Please run Module 1 (01_peak_classification.sh) first."
        exit 1
    fi
    echo "  Found: $bed_file ($(wc -l < "$bed_file") lines)"
done

# Check if HOMER is available
if ! command -v findMotifsGenome.pl &> /dev/null; then
    echo ""
    echo "Error: findMotifsGenome.pl not found."
    echo "Please install HOMER. See README.md for instructions."
    exit 1
fi

# Function to run HOMER
run_homer() {
    local bed_file=$1
    local name=$2

    echo ""
    echo "Running HOMER for $name..."
    mkdir -p "${OUTPUT_DIR}/${name}"

    # findMotifsGenome.pl <peak/BED file> <genome> <output directory> -size given -mask
    findMotifsGenome.pl "$bed_file" "$GENOME" "${OUTPUT_DIR}/${name}" \
        -size given \
        -mask \
        -p 8
}

# Run for each category
run_homer "${INPUT_DIR}/Shared_TES_TEAD1.bed" "Shared"
run_homer "${INPUT_DIR}/TES_Unique.bed" "TES_Unique"
run_homer "${INPUT_DIR}/TEAD1_Unique.bed" "TEAD1_Unique"

echo ""
echo "Module 2 completed: $(date)"
