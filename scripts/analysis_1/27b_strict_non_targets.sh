#!/bin/bash
#SBATCH --job-name=a1_27b_strict
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --time=2:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --output=logs/27b_strict_non_targets.out
#SBATCH --error=logs/27b_strict_non_targets.err

# =============================================================================
# STRICT NON-TARGET GENES - NO BINDING ANYWHERE
# =============================================================================
#
# More stringent definition: genes with NO TES/TEAD1 peaks anywhere in
# gene body OR within 10kb flanks (truly unbound genes)
#
# =============================================================================

echo "=========================================="
echo "STRICT NON-TARGET GENES"
echo "=========================================="
echo "Started: $(date)"
echo ""

cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_1

OUTDIR="output/27_non_targets_extended_flanks"
mkdir -p ${OUTDIR}

# -----------------------------------------------------------------------------
# STEP 1: CREATE STRICT NON-TARGET GENE SETS
# -----------------------------------------------------------------------------
echo "=== Step 1: Creating strict non-target gene sets ==="

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate r_chipseq_env

Rscript - << 'EOF'
suppressPackageStartupMessages({
    library(GenomicRanges)
    library(rtracklayer)
    library(dplyr)
})

setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/scripts/analysis_1")

cat("Creating STRICT non-target gene sets...\n\n")

OUTPUT_BASE <- "output/27_non_targets_extended_flanks"

# Load GTF
GTF_FILE <- "/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/annotation/gencode.v44.annotation.gtf"
gtf <- rtracklayer::import(GTF_FILE)
genes_gtf <- gtf[gtf$type == "gene"]
genes_df <- as.data.frame(genes_gtf)

genes_df <- genes_df %>%
    filter(gene_type == "protein_coding") %>%
    filter(seqnames %in% paste0("chr", c(1:22, "X", "Y")))

genes_df$gene_id_clean <- gsub("\\..*", "", genes_df$gene_id)

genes_gr <- GRanges(
    seqnames = genes_df$seqnames,
    ranges = IRanges(start = genes_df$start, end = genes_df$end),
    strand = genes_df$strand,
    gene_name = genes_df$gene_name,
    gene_id = genes_df$gene_id_clean
)

# Extend genes by 10kb each side for strict check
genes_extended <- genes_gr
start(genes_extended) <- pmax(1, start(genes_extended) - 10000)
end(genes_extended) <- end(genes_extended) + 10000

cat(sprintf("  Protein-coding genes: %d\n", length(genes_gr)))

# Load ALL TES and TEAD1 peaks
load_peaks <- function(peak_file, name) {
    peaks <- read.table(peak_file, header = FALSE,
                        col.names = c("chr", "start", "end", "name", "score",
                                      "strand", "signalValue", "pValue", "qValue", "peak"))
    peaks$chr <- paste0("chr", peaks$chr)  # Add chr prefix
    cat(sprintf("  %s peaks: %d\n", name, nrow(peaks)))
    GRanges(seqnames = peaks$chr, ranges = IRanges(start = peaks$start + 1, end = peaks$end))
}

TES_PEAKS <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_CUTandTAG/results/05_peaks_narrow/TES_peaks.narrowPeak"
TEAD1_PEAKS <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_CUTandTAG/results/05_peaks_narrow/TEAD1_peaks.narrowPeak"

tes_peaks <- load_peaks(TES_PEAKS, "TES")
tead1_peaks <- load_peaks(TEAD1_PEAKS, "TEAD1")
all_peaks <- c(tes_peaks, tead1_peaks)

# Find genes with ANY peak within gene body + 10kb flanks
overlaps <- findOverlaps(genes_extended, all_peaks)
genes_with_any_binding <- unique(queryHits(overlaps))
genes_no_binding_strict <- setdiff(1:length(genes_gr), genes_with_any_binding)

cat(sprintf("\n  Genes with TES/TEAD1 within gene+10kb: %d\n", length(genes_with_any_binding)))
cat(sprintf("  STRICT non-targets (no binding anywhere): %d\n", length(genes_no_binding_strict)))

# Load DMRs
DMR_FILE <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP/results/07_differential_MEDIPS/TES_vs_GFP_DMRs_FDR05_FC2.csv"
dmr_data <- read.csv(DMR_FILE)
dmr_hyper <- dmr_data %>% filter(logFC > 0)
dmr_hypo <- dmr_data %>% filter(logFC < 0)

create_dmr_gr <- function(dmr_df) {
    chr_names <- paste0("chr", dmr_df$chr)
    GRanges(seqnames = chr_names, ranges = IRanges(start = dmr_df$start, end = dmr_df$stop))
}

dmr_hyper_gr <- create_dmr_gr(dmr_hyper)
dmr_hypo_gr <- create_dmr_gr(dmr_hypo)

# Get strict non-target genes
strict_non_target_gr <- genes_gr[genes_no_binding_strict]

# Stratify by DMR status
hyper_overlaps <- findOverlaps(strict_non_target_gr, dmr_hyper_gr)
hypo_overlaps <- findOverlaps(strict_non_target_gr, dmr_hypo_gr)

strict_with_hyper <- unique(queryHits(hyper_overlaps))
strict_with_hypo <- unique(queryHits(hypo_overlaps))

strict_hyper_only <- setdiff(strict_with_hyper, strict_with_hypo)
strict_hypo_only <- setdiff(strict_with_hypo, strict_with_hyper)

all_dmr_overlaps <- findOverlaps(strict_non_target_gr, c(dmr_hyper_gr, dmr_hypo_gr))
strict_with_any_dmr <- unique(queryHits(all_dmr_overlaps))
strict_no_dmr <- setdiff(1:length(strict_non_target_gr), strict_with_any_dmr)

cat(sprintf("\n  Strict non-targets + hypermethylated: %d\n", length(strict_hyper_only)))
cat(sprintf("  Strict non-targets + hypomethylated: %d\n", length(strict_hypo_only)))
cat(sprintf("  Strict non-targets + no DMR: %d\n", length(strict_no_dmr)))

# Filter for expressed genes
DESEQ2_FILE <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_RNA/results/05_deseq2/deseq2_results_TES_vs_GFP.txt"
deseq2 <- read.delim(DESEQ2_FILE)
deseq2$gene_id_clean <- gsub("\\..*", "", deseq2$gene_id)
expressed_genes <- deseq2 %>% filter(!is.na(padj)) %>% pull(gene_id_clean)

strict_gene_ids <- strict_non_target_gr$gene_id

hyper_ids <- strict_gene_ids[strict_hyper_only]
hyper_expressed <- hyper_ids[hyper_ids %in% expressed_genes]

hypo_ids <- strict_gene_ids[strict_hypo_only]
hypo_expressed <- hypo_ids[hypo_ids %in% expressed_genes]

no_dmr_ids <- strict_gene_ids[strict_no_dmr]
no_dmr_expressed <- no_dmr_ids[no_dmr_ids %in% expressed_genes]
set.seed(42)
if (length(no_dmr_expressed) > 3000) {
    no_dmr_expressed <- sample(no_dmr_expressed, 3000)
}

cat(sprintf("\n  STRICT + Hypermethylated (expressed): %d\n", length(hyper_expressed)))
cat(sprintf("  STRICT + Hypomethylated (expressed): %d\n", length(hypo_expressed)))
cat(sprintf("  STRICT + No DMR (expressed, sampled): %d\n", length(no_dmr_expressed)))

# Create BED files
create_bed <- function(gene_ids, genes_df, output_file) {
    selected <- genes_df[genes_df$gene_id_clean %in% gene_ids, ]
    bed <- data.frame(
        chr = selected$seqnames,
        start = selected$start - 1,
        end = selected$end,
        name = selected$gene_name,
        score = 0,
        strand = selected$strand
    )
    bed <- bed[!is.na(bed$start) & bed$start >= 0, ]
    bed <- bed[order(bed$chr, bed$start), ]
    bed <- bed[!duplicated(paste(bed$chr, bed$start, bed$end)), ]
    write.table(bed, output_file, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
    cat(sprintf("  Created: %s (%d genes)\n", basename(output_file), nrow(bed)))
    return(nrow(bed))
}

cat("\nCreating STRICT BED files:\n")
n_hyper <- create_bed(hyper_expressed, genes_df, file.path(OUTPUT_BASE, "strict_non_targets_hyper.bed"))
n_hypo <- create_bed(hypo_expressed, genes_df, file.path(OUTPUT_BASE, "strict_non_targets_hypo.bed"))
n_nodmr <- create_bed(no_dmr_expressed, genes_df, file.path(OUTPUT_BASE, "strict_non_targets_nodmr.bed"))

# Save counts
write.csv(data.frame(
    category = c("strict_hyper", "strict_hypo", "strict_nodmr"),
    n_genes = c(n_hyper, n_hypo, n_nodmr)
), file.path(OUTPUT_BASE, "strict_gene_counts.csv"), row.names = FALSE)

cat("\nDone!\n")
EOF

# Check if files were created
if [ ! -f "${OUTDIR}/strict_non_targets_hyper.bed" ]; then
    echo "ERROR: Strict BED files not created"
    exit 1
fi

# -----------------------------------------------------------------------------
# STEP 2: VISUALIZATION
# -----------------------------------------------------------------------------
echo ""
echo "=== Step 2: Creating visualizations ==="

conda activate tg

TES_METH="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP/results/05_bigwig/TES_combined_RPKM.bw"
GFP_METH="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP/results/05_bigwig/GFP_combined_RPKM.bw"
TES_BIND="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_CUTandTAG/results/06_bigwig/TES_comb.bw"
TEAD1_BIND="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_CUTandTAG/results/06_bigwig/TEAD1_comb.bw"

N_HYPER=$(wc -l < ${OUTDIR}/strict_non_targets_hyper.bed)
N_HYPO=$(wc -l < ${OUTDIR}/strict_non_targets_hypo.bed)
N_NODMR=$(wc -l < ${OUTDIR}/strict_non_targets_nodmr.bed)

echo "STRICT non-target gene counts:"
echo "  Hypermethylated: ${N_HYPER}"
echo "  Hypomethylated: ${N_HYPO}"
echo "  No DMR: ${N_NODMR}"
echo ""

# Comparison plot with 10kb flanks
echo "Creating STRICT comparison plot (10kb flanks)..."
computeMatrix scale-regions \
    -S $TES_METH $GFP_METH $TES_BIND $TEAD1_BIND \
    -R ${OUTDIR}/strict_non_targets_hyper.bed \
       ${OUTDIR}/strict_non_targets_hypo.bed \
       ${OUTDIR}/strict_non_targets_nodmr.bed \
    --beforeRegionStartLength 10000 \
    --afterRegionStartLength 10000 \
    --regionBodyLength 5000 \
    --binSize 100 \
    --skipZeros \
    --missingDataAsZero \
    -o ${OUTDIR}/strict_comparison_matrix.gz \
    -p 16 \
    2>&1 | grep -v "Skipping\|did not match"

plotProfile -m ${OUTDIR}/strict_comparison_matrix.gz \
    -out ${OUTDIR}/STRICT_non_targets_comparison.png \
    --perGroup \
    --colors "#7B3294" "#636363" "#E31A1C" "#377EB8" \
    --startLabel "TSS" \
    --endLabel "TES" \
    --samplesLabel "TES meth" "GFP meth" "TES bind" "TEAD1 bind" \
    --regionsLabel "Hyper (n=${N_HYPER})" "Hypo (n=${N_HYPO})" "No DMR (n=${N_NODMR})" \
    --plotTitle "STRICT NON-TARGETS: No TES/TEAD1 binding within gene+10kb" \
    --plotHeight 12 \
    --plotWidth 20 \
    --legendLocation "lower-right" \
    --yMin 0 \
    --dpi 300

echo "  Done: STRICT_non_targets_comparison.png"

# Compare bound vs strictly unbound hypermethylated genes
echo ""
echo "Creating bound vs unbound comparison..."

# Get genes WITH binding
BOUND_BED="output/25_genebody_methylation_DMR_genes/genes_with_hypermethylated_DMRs.bed"

if [ -f "${BOUND_BED}" ]; then
    N_BOUND=$(wc -l < ${BOUND_BED})

    computeMatrix scale-regions \
        -S $TES_METH $GFP_METH $TES_BIND $TEAD1_BIND \
        -R ${BOUND_BED} \
           ${OUTDIR}/strict_non_targets_hyper.bed \
        --beforeRegionStartLength 10000 \
        --afterRegionStartLength 10000 \
        --regionBodyLength 5000 \
        --binSize 100 \
        --skipZeros \
        --missingDataAsZero \
        -o ${OUTDIR}/bound_vs_strict_unbound_matrix.gz \
        -p 16 \
        2>&1 | grep -v "Skipping\|did not match"

    plotProfile -m ${OUTDIR}/bound_vs_strict_unbound_matrix.gz \
        -out ${OUTDIR}/bound_vs_STRICT_unbound_hyper.png \
        --perGroup \
        --colors "#7B3294" "#636363" "#E31A1C" "#377EB8" \
        --startLabel "TSS" \
        --endLabel "TES" \
        --samplesLabel "TES meth" "GFP meth" "TES bind" "TEAD1 bind" \
        --regionsLabel "All Hyper (n=${N_BOUND})" "STRICT unbound Hyper (n=${N_HYPER})" \
        --plotTitle "Hypermethylated: All Genes vs STRICT Non-Targets" \
        --plotHeight 12 \
        --plotWidth 20 \
        --legendLocation "upper-left" \
        --yMin 0 \
        --dpi 300

    echo "  Done: bound_vs_STRICT_unbound_hyper.png"
fi

# -----------------------------------------------------------------------------
# SUMMARY
# -----------------------------------------------------------------------------
echo ""
echo "=========================================="
echo "COMPLETE"
echo "=========================================="
echo "Finished: $(date)"
echo ""
echo "STRICT non-target files:"
ls -lh ${OUTDIR}/STRICT*.png ${OUTDIR}/bound_vs_STRICT*.png 2>/dev/null
echo ""
