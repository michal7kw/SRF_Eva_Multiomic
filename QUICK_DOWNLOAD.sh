#!/bin/bash
#
# QUICK DOWNLOAD SCRIPT FOR MSIGDB COLLECTIONS
# Downloads and converts the collections you requested
#

cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis

echo "============================================================================="
echo "  MSigDB Gene Set Collections - Quick Download"
echo "============================================================================="
echo ""

GENE_SET_DIR="GENE_SETS"
cd "$GENE_SET_DIR"

echo "Downloading collections to: $(pwd)"
echo ""

# ============================================================================
# C3CA: Curated Cancer Cell Atlas (148 gene sets)
# ============================================================================
echo "1/6: Downloading C3CA - Curated Cancer Cell Atlas (148 sets)..."
wget -q --show-progress https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/c3.ca.v2024.1.Hs.symbols.gmt
echo "  ✓ Downloaded c3.ca.v2024.1.Hs.symbols.gmt"
echo ""

# ============================================================================
# C4.CGN: Cancer Gene Neighborhoods (427 gene sets)
# ============================================================================
echo "2/6: Downloading C4.CGN - Cancer Gene Neighborhoods (427 sets)..."
wget -q --show-progress https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/c4.cgn.v2024.1.Hs.symbols.gmt
echo "  ✓ Downloaded c4.cgn.v2024.1.Hs.symbols.gmt"
echo ""

# ============================================================================
# C4.CM: Cancer Modules (431 gene sets)
# ============================================================================
echo "3/6: Downloading C4.CM - Cancer Modules (431 sets)..."
wget -q --show-progress https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/c4.cm.v2024.1.Hs.symbols.gmt
echo "  ✓ Downloaded c4.cm.v2024.1.Hs.symbols.gmt"
echo ""

# ============================================================================
# C5: GO Biological Process (~7,800 gene sets)
# ============================================================================
echo "4/6: Downloading C5.GO.BP - GO Biological Process (~7800 sets)..."
echo "  (This may take a few minutes - large file)"
wget -q --show-progress https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/c5.go.bp.v2024.1.Hs.symbols.gmt
echo "  ✓ Downloaded c5.go.bp.v2024.1.Hs.symbols.gmt"
echo ""

# ============================================================================
# C2.CP: KEGG Pathways (~190 gene sets)
# ============================================================================
echo "5/6: Downloading C2.CP.KEGG - KEGG Pathways (~190 sets)..."
wget -q --show-progress https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/c2.cp.kegg.v2024.1.Hs.symbols.gmt
echo "  ✓ Downloaded c2.cp.kegg.v2024.1.Hs.symbols.gmt"
echo ""

# ============================================================================
# C2.CP: Reactome Pathways (~1,900 gene sets)
# ============================================================================
echo "6/6: Downloading C2.CP.Reactome - Reactome Pathways (~1900 sets)..."
echo "  (This may take a few minutes - large file)"
wget -q --show-progress https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/c2.cp.reactome.v2024.1.Hs.symbols.gmt
echo "  ✓ Downloaded c2.cp.reactome.v2024.1.Hs.symbols.gmt"
echo ""

echo "============================================================================="
echo "  All downloads complete!"
echo "============================================================================="
echo ""

cd ..

# ============================================================================
# Convert GMT to GRP format
# ============================================================================
echo "Converting GMT files to GRP format..."
echo ""

for gmt_file in "$GENE_SET_DIR"/*.gmt; do
    if [ -f "$gmt_file" ]; then
        basename=$(basename "$gmt_file")
        echo "Converting: $basename"
        Rscript scripts/convert_gmt_to_grp.R "$gmt_file" "$GENE_SET_DIR/"
    fi
done

echo ""
echo "============================================================================="
echo "  Conversion complete!"
echo "============================================================================="
echo ""

# Check what we have
echo "Checking gene set status..."
./scripts/check_gene_sets.sh

echo ""
echo "============================================================================="
echo "  READY TO RUN ANALYSIS"
echo "============================================================================="
echo ""
echo "Run GSEA with collection-based organization:"
echo "  sbatch scripts/run_msigdb_by_collection.sh"
echo ""
echo "Results will be saved in:"
echo "  output/12_msigdb_by_collection/"
echo ""
echo "Each collection will have its own subdirectory with:"
echo "  - all_pathways.csv"
echo "  - significant_pathways.csv"
echo "  - plot_01_top_pathways.pdf"
echo "  - plot_02_enrichment_scatter.pdf"
echo "  - plot_03_enrichment_curves.pdf"
echo "  - SUMMARY.txt"
echo ""
