#!/bin/bash
#SBATCH --job-name=a1_26_combined
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --time=2:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --output=logs/26_combined_binding_methylation.out
#SBATCH --error=logs/26_combined_binding_methylation.err

# =============================================================================
# COMBINED BINDING + METHYLATION GENE BODY PROFILES
# =============================================================================
#
# This script generates gene body profiles showing both TES/TEAD1 binding
# (from Cut&Tag) and methylation (from meDIP) for genes stratified by DMR status.
#
# Gene sets (from 25_genebody_methylation_DMR_genes.R):
#   - Hypermethylated: Genes with hypermethylated DMRs in gene body
#   - Hypomethylated: Genes with hypomethylated DMRs in gene body
#   - No DMR: Expressed genes without any DMRs (sampled control)
#   - DEGs DOWN + Hypermethylated: Downregulated genes with hypermethylation
#
# Output: PNG profile plots with 4 signals:
#   - TES binding (Cut&Tag)
#   - TEAD1 binding (Cut&Tag)
#   - TES methylation (meDIP)
#   - GFP methylation (meDIP control)
#
# =============================================================================

echo "=========================================="
echo "COMBINED BINDING + METHYLATION PROFILES"
echo "=========================================="
echo "Started: $(date)"
echo ""

cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_1

# -----------------------------------------------------------------------------
# CONFIGURATION
# -----------------------------------------------------------------------------

OUTDIR="output/26_combined_binding_methylation"
mkdir -p ${OUTDIR}

BED_DIR="output/25_genebody_methylation_DMR_genes"

# BigWig files - Binding (Cut&Tag)
TES_BIND="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_CUTandTAG/results/06_bigwig/TES_comb.bw"
TEAD1_BIND="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_CUTandTAG/results/06_bigwig/TEAD1_comb.bw"

# BigWig files - Methylation (meDIP)
TES_METH="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP/results/05_bigwig/TES_combined_RPKM.bw"
GFP_METH="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP/results/05_bigwig/GFP_combined_RPKM.bw"

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate tg

# Count genes in each category
N_HYPER=$(wc -l < ${BED_DIR}/genes_with_hypermethylated_DMRs.bed)
N_HYPO=$(wc -l < ${BED_DIR}/genes_with_hypomethylated_DMRs.bed)
N_NODMR=$(wc -l < ${BED_DIR}/genes_without_DMRs.bed)
N_DEGS_DOWN_HYPER=$(wc -l < ${BED_DIR}/DEGs_DOWN_hypermethylated.bed 2>/dev/null || echo "0")

echo "Gene counts:"
echo "  Hypermethylated DMRs: ${N_HYPER}"
echo "  Hypomethylated DMRs: ${N_HYPO}"
echo "  No DMR (control): ${N_NODMR}"
echo "  DEGs DOWN + Hypermethylated: ${N_DEGS_DOWN_HYPER}"
echo ""

# -----------------------------------------------------------------------------
# 1. HYPERMETHYLATED GENES
# -----------------------------------------------------------------------------
echo "=== Processing hypermethylated genes (n=${N_HYPER}) ==="

computeMatrix scale-regions \
    -S $TES_BIND $TEAD1_BIND $TES_METH $GFP_METH \
    -R ${BED_DIR}/genes_with_hypermethylated_DMRs.bed \
    --beforeRegionStartLength 3000 \
    --afterRegionStartLength 3000 \
    --regionBodyLength 5000 \
    --binSize 50 \
    --skipZeros \
    --missingDataAsZero \
    -o ${OUTDIR}/hypermethylated_matrix.gz \
    -p 16 \
    2>&1 | grep -v "Skipping\|did not match"

plotProfile -m ${OUTDIR}/hypermethylated_matrix.gz \
    -out ${OUTDIR}/hypermethylated_genes_profile.png \
    --perGroup \
    --colors "#E31A1C" "#377EB8" "#7B3294" "#636363" \
    --startLabel "TSS" \
    --endLabel "TES" \
    --samplesLabel "TES bind" "TEAD1 bind" "TES meth" "GFP meth" \
    --plotTitle "Binding + Methylation at Hypermethylated Genes (n=${N_HYPER})" \
    --plotHeight 10 \
    --plotWidth 16 \
    --legendLocation "lower-right" \
    --yMin 0 \
    --dpi 300

echo "  Done: hypermethylated_genes_profile.png"

# -----------------------------------------------------------------------------
# 2. HYPOMETHYLATED GENES
# -----------------------------------------------------------------------------
echo ""
echo "=== Processing hypomethylated genes (n=${N_HYPO}) ==="

computeMatrix scale-regions \
    -S $TES_BIND $TEAD1_BIND $TES_METH $GFP_METH \
    -R ${BED_DIR}/genes_with_hypomethylated_DMRs.bed \
    --beforeRegionStartLength 3000 \
    --afterRegionStartLength 3000 \
    --regionBodyLength 5000 \
    --binSize 50 \
    --skipZeros \
    --missingDataAsZero \
    -o ${OUTDIR}/hypomethylated_matrix.gz \
    -p 16 \
    2>&1 | grep -v "Skipping\|did not match"

plotProfile -m ${OUTDIR}/hypomethylated_matrix.gz \
    -out ${OUTDIR}/hypomethylated_genes_profile.png \
    --perGroup \
    --colors "#E31A1C" "#377EB8" "#7B3294" "#636363" \
    --startLabel "TSS" \
    --endLabel "TES" \
    --samplesLabel "TES bind" "TEAD1 bind" "TES meth" "GFP meth" \
    --plotTitle "Binding + Methylation at Hypomethylated Genes (n=${N_HYPO})" \
    --plotHeight 10 \
    --plotWidth 16 \
    --legendLocation "lower-right" \
    --yMin 0 \
    --dpi 300

echo "  Done: hypomethylated_genes_profile.png"

# -----------------------------------------------------------------------------
# 3. DEGs DOWN with HYPERMETHYLATED DMRs
# -----------------------------------------------------------------------------
if [ -f "${BED_DIR}/DEGs_DOWN_hypermethylated.bed" ]; then
    echo ""
    echo "=== Processing DEGs DOWN + hypermethylated (n=${N_DEGS_DOWN_HYPER}) ==="

    computeMatrix scale-regions \
        -S $TES_BIND $TEAD1_BIND $TES_METH $GFP_METH \
        -R ${BED_DIR}/DEGs_DOWN_hypermethylated.bed \
        --beforeRegionStartLength 3000 \
        --afterRegionStartLength 3000 \
        --regionBodyLength 5000 \
        --binSize 50 \
        --skipZeros \
        --missingDataAsZero \
        -o ${OUTDIR}/DEGs_DOWN_hyper_matrix.gz \
        -p 16 \
        2>&1 | grep -v "Skipping\|did not match"

    plotProfile -m ${OUTDIR}/DEGs_DOWN_hyper_matrix.gz \
        -out ${OUTDIR}/DEGs_DOWN_hyper_profile.png \
        --perGroup \
        --colors "#E31A1C" "#377EB8" "#7B3294" "#636363" \
        --startLabel "TSS" \
        --endLabel "TES" \
        --samplesLabel "TES bind" "TEAD1 bind" "TES meth" "GFP meth" \
        --plotTitle "DEGs DOWN with Hypermethylation (n=${N_DEGS_DOWN_HYPER})" \
        --plotHeight 10 \
        --plotWidth 16 \
        --legendLocation "lower-right" \
        --yMin 0 \
        --dpi 300

    echo "  Done: DEGs_DOWN_hyper_profile.png"
fi

# -----------------------------------------------------------------------------
# 4. COMPARISON - All gene sets
# -----------------------------------------------------------------------------
echo ""
echo "=== Creating comparison plot ==="

computeMatrix scale-regions \
    -S $TES_BIND $TEAD1_BIND $TES_METH $GFP_METH \
    -R ${BED_DIR}/genes_with_hypermethylated_DMRs.bed \
       ${BED_DIR}/genes_with_hypomethylated_DMRs.bed \
       ${BED_DIR}/genes_without_DMRs.bed \
    --beforeRegionStartLength 3000 \
    --afterRegionStartLength 3000 \
    --regionBodyLength 5000 \
    --binSize 50 \
    --skipZeros \
    --missingDataAsZero \
    -o ${OUTDIR}/comparison_matrix.gz \
    -p 16 \
    2>&1 | grep -v "Skipping\|did not match"

plotProfile -m ${OUTDIR}/comparison_matrix.gz \
    -out ${OUTDIR}/comparison_profile.png \
    --perGroup \
    --colors "#E31A1C" "#377EB8" "#7B3294" "#636363" \
    --startLabel "TSS" \
    --endLabel "TES" \
    --samplesLabel "TES bind" "TEAD1 bind" "TES meth" "GFP meth" \
    --regionsLabel "Hypermeth (n=${N_HYPER})" "Hypometh (n=${N_HYPO})" "No DMR (n=${N_NODMR})" \
    --plotTitle "Gene Body Profiles by DMR Status" \
    --plotHeight 12 \
    --plotWidth 18 \
    --legendLocation "upper-left" \
    --yMin 0 \
    --dpi 300

echo "  Done: comparison_profile.png"

# -----------------------------------------------------------------------------
# HEATMAPS (commented out - uncomment if needed)
# -----------------------------------------------------------------------------
# echo ""
# echo "=== Creating heatmaps ==="
#
# plotHeatmap -m ${OUTDIR}/hypermethylated_matrix.gz \
#     -out ${OUTDIR}/hypermethylated_genes_heatmap.png \
#     --colorList "white,#2166AC" "white,#B2182B" "white,#762A83" "white,#636363" \
#     --zMin 0 0 0 0 \
#     --zMax 3 3 200 200 \
#     --sortRegions descend \
#     --sortUsing mean \
#     --sortUsingSamples 1 2 \
#     --heatmapHeight 15 \
#     --heatmapWidth 3 \
#     --startLabel "TSS" \
#     --endLabel "TES" \
#     --plotTitle "Hypermethylated Genes" \
#     --dpi 300
#
# plotHeatmap -m ${OUTDIR}/DEGs_DOWN_hyper_matrix.gz \
#     -out ${OUTDIR}/DEGs_DOWN_hyper_heatmap.png \
#     --colorList "white,#2166AC" "white,#B2182B" "white,#762A83" "white,#636363" \
#     --zMin 0 0 0 0 \
#     --zMax 3 3 200 200 \
#     --sortRegions descend \
#     --sortUsing mean \
#     --sortUsingSamples 1 2 \
#     --heatmapHeight 15 \
#     --heatmapWidth 3 \
#     --startLabel "TSS" \
#     --endLabel "TES" \
#     --plotTitle "DEGs DOWN + Hypermethylated" \
#     --dpi 300

# -----------------------------------------------------------------------------
# SUMMARY
# -----------------------------------------------------------------------------
echo ""
echo "=========================================="
echo "COMPLETE"
echo "=========================================="
echo "Finished: $(date)"
echo ""
echo "Output directory: ${OUTDIR}"
echo ""
echo "PNG files generated:"
ls -lh ${OUTDIR}/*.png
echo ""
echo "Matrix files (for re-plotting):"
ls -lh ${OUTDIR}/*.gz
echo ""
