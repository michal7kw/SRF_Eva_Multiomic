#!/bin/bash
# SBATCH --job-name=convert_gmt_to_grp
# SBATCH --account=kubacki.michal
# SBATCH --partition=workq
# SBATCH --time=1:00:00
# SBATCH --mem=16G
# SBATCH --cpus-per-task=4
# SBATCH --output=logs/convert_gmt_to_grp.out
# SBATCH --error=logs/convert_gmt_to_grp.err
# SBATCH --mail-type=END,FAIL
# SBATCH --mail-user=kubacki.michal@hsr.it
# CHECK GENE SETS - Quick diagnostic script
# Shows what MSigDB gene sets are available and ready for analysis

GENE_SET_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/GENE_SETS"
SCRIPT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_1"

echo "============================================================================="
echo "  MSIGDB GENE SET DIAGNOSTIC"
echo "============================================================================="
echo ""
echo "Gene set directory: $GENE_SET_DIR"
echo ""

# Check if directory exists
if [ ! -d "$GENE_SET_DIR" ]; then
    echo "ERROR: Gene set directory not found!"
    exit 1
fi

cd "$GENE_SET_DIR"

# Count files
n_grp=$(find . -name "*.grp" | wc -l)
n_gmt=$(find . -name "*.gmt" | wc -l)

echo "File counts:"
echo "  .grp files (ready for analysis): $n_grp"
echo "  .gmt files (need conversion): $n_gmt"
echo ""

if [ $n_grp -eq 0 ] && [ $n_gmt -eq 0 ]; then
    echo "⚠️  NO GENE SETS FOUND"
    echo ""
    echo "You need to download MSigDB gene sets."
    echo "See: GENE_SETS/README.md for instructions"
    echo ""
    exit 1
fi

if [ $n_grp -gt 0 ]; then
    echo "✓ Found $n_grp .grp gene sets (ready for analysis)"
    echo ""

    # Check species
    n_human=$(find . -name "*.Hs.grp" | wc -l)
    n_mouse=$(find . -name "*.Mm.grp" | wc -l)

    echo "Species breakdown:"
    echo "  Human (Hs): $n_human"
    echo "  Mouse (Mm): $n_mouse"
    echo ""

    if [ $n_human -gt 0 ]; then
        echo "✓ Using human gene sets for analysis"
        echo ""
    fi

    if [ $n_mouse -gt 0 ] && [ $n_human -eq 0 ]; then
        echo "⚠️  WARNING: Only mouse gene sets found"
        echo "   Your data is HUMAN (SNB19 cells)"
        echo "   Action: Download human gene sets from MSigDB"
        echo ""
    fi

    # Check collections
    echo "Collections detected:"
    for collection in HALLMARK GOBP GOCC GOMF KEGG REACTOME BIOCARTA WP PID
    do
        count=$(find . -name "${collection}_*.grp" | wc -l)
        if [ $count -gt 0 ]; then
            echo "  $collection: $count gene sets"
        fi
    done
    echo ""

    # Show sample gene sets
    echo "Sample gene sets (first 10):"
    find . -name "*.grp" | head -10 | while read file; do
        basename "$file"
    done
    echo ""

    # Check gene set sizes
    echo "Checking gene set sizes..."
    total_genes=0
    count=0
    for file in $(find . -name "*.grp" | head -100); do
        # Count non-comment lines excluding first line (pathway name)
        n_genes=$(grep -v "^#" "$file" | tail -n +2 | wc -l)
        total_genes=$((total_genes + n_genes))
        count=$((count + 1))
    done

    if [ $count -gt 0 ]; then
        avg_genes=$((total_genes / count))
        echo "  Average genes per set (sample of $count): $avg_genes"
        echo ""
    fi
fi

if [ $n_gmt -gt 0 ]; then
    echo "⚠️  Found $n_gmt .gmt files (need conversion)"
    echo ""
    echo "GMT files need to be converted to .grp format:"
    find . -name "*.gmt" | while read file; do
        echo "  - $(basename "$file")"
    done
    echo ""
    echo "Convert with:"
    echo "  Rscript ${SCRIPT_DIR}/convert_gmt_to_grp.R $GENE_SET_DIR/<file.gmt> $GENE_SET_DIR/"
    echo ""
fi

echo "============================================================================="
echo "  STATUS SUMMARY"
echo "============================================================================="

if [ $n_grp -eq 0 ]; then
    echo "❌ NOT READY - No .grp files found"
    echo "   Action: Download and convert MSigDB gene sets"
elif [ $n_grp -lt 50 ]; then
    echo "⚠️  LIMITED - Only $n_grp gene sets"
    echo "   Action: Consider downloading more collections for comprehensive analysis"
else
    echo "✅ READY - $n_grp gene sets available"
    echo "   Action: Run analysis with 'sbatch msigdb_gsea_by_collection.sh'"
fi

echo ""
echo "For detailed instructions, see:"
echo "  - GENE_SETS/README.md"
echo "  - MSIGDB_GSEA_GUIDE.md"
echo ""
