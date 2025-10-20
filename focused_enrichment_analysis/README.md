# Enrichment Analysis - Migration/Cancer Focus

## Directory Structure

```
focused_enrichment_analysis/
├── approach1_direct_targets/      # Baseline: Direct targets (peaks + DEGs)
├── approach2_downregulated/       # Mechanistic: Downregulated direct targets
├── approach3_promoter_peaks/      # Positional: Promoter peaks only
├── approach4_high_confidence/     # Signal-based: Top 50% peaks
├── approach5_diffbind/            # Differential: TES vs TESmut specific
```

## Approaches Explained

### **Approach 1: Direct Targets Only**

**Rationale**: genes that are BOTH bound and differentially expressed

**Method**:
- TES peaks ∩ significant DEGs (padj < 0.05)

**Files**:
- `gene_lists/TES_direct_targets.txt`
- `results/approach1_TES_direct_GO_enrichment.csv`
- `results/approach1_migration_terms.csv`

---

### **Approach 2: Downregulated Direct Targets**

**Rationale**: TES is a REPRESSOR → focus on downregulated genes only

**Method**:
- Direct targets + log2FC < 0

**Files**:
- `gene_lists/TES_direct_downregulated.txt`
- `gene_lists/TES_direct_downregulated_detailed.csv`
- `results/approach2_TES_direct_down_GO_enrichment.csv`
- `results/approach2_migration_terms.csv`

---

### **Approach 3: Promoter Peaks Only**

**Rationale**: Promoter binding = more direct regulation than distal peaks

**Method**:
- Filter ChIPseeker annotations for "Promoter", "5' UTR", "1st Exon"
- Intersect with DEGs

**Files**:
- `gene_lists/TES_promoter_peak_genes.txt`
- `gene_lists/TES_promoter_DEGs.txt`
- `results/approach3_TES_promoter_DEGs_GO_enrichment.csv`

---

### **Approach 4: High-Confidence Peaks**

**Rationale**: Not all peaks are equal - use top 50% by signal strength

**Method**:
- Filter peaks by median qValue/fold_enrichment
- Map to genes, intersect with DEGs

**Files**:
- `gene_lists/TES_highconf_peak_genes.txt`
- `gene_lists/TES_highconf_DEGs.txt`
- `results/approach4_TES_highconf_DEGs_GO_enrichment.csv`

---

### **Approach 5: Differential Binding + Differential Expression**

**Rationale**: Sites differential between TES and TESmut are functional TES targets

**Method**:
- DiffBind: TES > TESmut (FDR < 0.05, Fold > 1.5)
- Intersect with downregulated DEGs

**Files**:
- `gene_lists/TES_diffbind_genes.txt`
- `gene_lists/TES_diffbind_down.txt`
- `results/approach5_diffbind_down_GO_enrichment.csv`