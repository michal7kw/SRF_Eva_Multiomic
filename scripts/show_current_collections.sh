#!/bin/bash
#
# Shows what MSigDB collections are currently in GENE_SETS/
#

echo "============================================================================="
echo "  CURRENT GENE SET COLLECTIONS"
echo "============================================================================="
echo ""

cd GENE_SETS

total=$(ls -1 *.grp 2>/dev/null | wc -l)

if [ $total -eq 0 ]; then
    echo "ERROR: No .grp files found in GENE_SETS/"
    exit 1
fi

echo "Total gene sets: $total"
echo ""
echo "Collections by prefix:"
echo "-----------------------------------------------------------------------------"

# Count by prefix
ls -1 *.grp | sed 's/_.*//' | sort | uniq -c | sort -rn | while read count prefix; do
    case $prefix in
        HALLMARK)
            echo "  ✓ $count sets - HALLMARK (Hallmark gene sets)"
            ;;
        GOBP)
            echo "  ✓ $count sets - GOBP (GO Biological Process)"
            ;;
        GOCC)
            echo "  ✓ $count sets - GOCC (GO Cellular Component)"
            ;;
        GOMF)
            echo "  ✓ $count sets - GOMF (GO Molecular Function)"
            ;;
        KEGG)
            echo "  ✓ $count sets - KEGG (KEGG Pathways)"
            ;;
        REACTOME)
            echo "  ✓ $count sets - REACTOME (Reactome Pathways)"
            ;;
        WP)
            echo "  ✓ $count sets - WP (WikiPathways)"
            ;;
        BIOCARTA)
            echo "  ✓ $count sets - BIOCARTA (BioCarta Pathways)"
            ;;
        PID)
            echo "  ✓ $count sets - PID (PID Pathways)"
            ;;
        CAR)
            echo "  ✓ $count sets - CAR (C3 Curated Cancer Cell Atlas)"
            ;;
        MODULE)
            echo "  ✓ $count sets - MODULE (C4 Cancer Modules)"
            ;;
        GCM)
            echo "  ✓ $count sets - GCM (C4 Cancer Gene Neighborhoods)"
            ;;
        MORF)
            echo "  ✓ $count sets - MORF (C2 Curated - Morpheus)"
            ;;
        GNF2)
            echo "  ✓ $count sets - GNF2 (C2 Curated - GNF2 Atlas)"
            ;;
        *)
            echo "  ? $count sets - $prefix (Other/Unknown)"
            ;;
    esac
done

echo ""
echo "============================================================================="
echo "  COLLECTION MAPPING (for GSEA output folders)"
echo "============================================================================="
echo ""

# Show which output folders will be created
has_hallmark=$(ls -1 *.grp 2>/dev/null | grep -c "^HALLMARK_")
has_gobp=$(ls -1 *.grp 2>/dev/null | grep -c "^GOBP_")
has_gocc=$(ls -1 *.grp 2>/dev/null | grep -c "^GOCC_")
has_gomf=$(ls -1 *.grp 2>/dev/null | grep -c "^GOMF_")
has_kegg=$(ls -1 *.grp 2>/dev/null | grep -c "^KEGG_")
has_reactome=$(ls -1 *.grp 2>/dev/null | grep -c "^REACTOME_")
has_wp=$(ls -1 *.grp 2>/dev/null | grep -c "^WP_")
has_biocarta=$(ls -1 *.grp 2>/dev/null | grep -c "^BIOCARTA_")
has_pid=$(ls -1 *.grp 2>/dev/null | grep -c "^PID_")
has_car=$(ls -1 *.grp 2>/dev/null | grep -c "^CAR_")
has_module=$(ls -1 *.grp 2>/dev/null | grep -c "^MODULE_")
has_gcm=$(ls -1 *.grp 2>/dev/null | grep -c "^GCM_")
has_morf=$(ls -1 *.grp 2>/dev/null | grep -c "^MORF_")
has_gnf2=$(ls -1 *.grp 2>/dev/null | grep -c "^GNF2_")

if [ $has_hallmark -gt 0 ]; then
    echo "  01_Hallmark/ ($has_hallmark gene sets)"
fi

if [ $has_gobp -gt 0 ]; then
    echo "  02_GO_Biological_Process/ ($has_gobp gene sets)"
fi

if [ $has_gocc -gt 0 ]; then
    echo "  03_GO_Cellular_Component/ ($has_gocc gene sets)"
fi

if [ $has_gomf -gt 0 ]; then
    echo "  04_GO_Molecular_Function/ ($has_gomf gene sets)"
fi

if [ $has_kegg -gt 0 ]; then
    echo "  05_KEGG_Pathways/ ($has_kegg gene sets)"
fi

if [ $has_reactome -gt 0 ]; then
    echo "  06_Reactome_Pathways/ ($has_reactome gene sets)"
fi

if [ $has_wp -gt 0 ]; then
    echo "  07_WikiPathways/ ($has_wp gene sets)"
fi

if [ $has_biocarta -gt 0 ]; then
    echo "  08_BioCarta/ ($has_biocarta gene sets)"
fi

if [ $has_pid -gt 0 ]; then
    echo "  09_PID/ ($has_pid gene sets)"
fi

if [ $has_car -gt 0 ]; then
    echo "  11_C3_Cancer_Atlas/ ($has_car gene sets - CAR prefix)"
fi

c4_total=$((has_module + has_gcm))
if [ $c4_total -gt 0 ]; then
    echo "  12_C4_Computational/ ($c4_total gene sets - MODULE + GCM)"
fi

c2_curated=$((has_morf + has_gnf2))
if [ $c2_curated -gt 0 ]; then
    echo "  14_C2_Curated/ ($c2_curated gene sets - MORF + GNF2)"
fi

echo ""
echo "============================================================================="
echo "  MISSING COLLECTIONS (can be downloaded)"
echo "============================================================================="
echo ""

if [ $has_hallmark -eq 0 ]; then
    echo "  ✗ HALLMARK - Hallmark gene sets (50 sets)"
    echo "    wget https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/h.all.v2024.1.Hs.symbols.gmt"
    echo ""
fi

if [ $has_kegg -eq 0 ]; then
    echo "  ✗ KEGG - KEGG Pathways (~190 sets)"
    echo "    wget https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/c2.cp.kegg.v2024.1.Hs.symbols.gmt"
    echo ""
fi

if [ $has_gocc -eq 0 ]; then
    echo "  ✗ GO Cellular Component (~1,000 sets)"
    echo "    wget https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/c5.go.cc.v2024.1.Hs.symbols.gmt"
    echo ""
fi

if [ $has_gomf -eq 0 ]; then
    echo "  ✗ GO Molecular Function (~1,700 sets)"
    echo "    wget https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/c5.go.mf.v2024.1.Hs.symbols.gmt"
    echo ""
fi

echo "============================================================================="
echo ""
