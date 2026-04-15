# Biological Review: Dr. Han Data Portal Datasets & Single-Cell Universe Integration

**Reviewer:** Cell Biologist (CellBiologist agent)
**Date:** 2026-04-15
**Issue:** [OMI-36](/OMI/issues/OMI-36)
**Parent:** [OMI-34](/OMI/issues/OMI-34)

---

## 1. Dataset Summary

The 4 datasets selected by Chien-Wei Lin (highlighted in red in `Data Portal_source datasets.xlsx`) are:

| # | Study | Species | Tissue | #Cells/Nuclei | Platform | Accession |
|---|-------|---------|--------|---------------|----------|-----------|
| 1 | Wang et al. 2020, *Nat Cell Biol* | Human (adult) | Heart (healthy + failing + post-LVAD) | 21,422 cells | 10x scRNA-seq | PRJNA548584 |
| 2 | Tucker et al. 2020, *Circulation* | Human (adult) | Heart (healthy donors) | 287,269 nuclei | 10x snRNA-seq | GSE109816 |
| 3 | Zhu et al. 2024, *Nature* | Human (fetal) | Heart (PCW 9–16) | 142,946 cells | 10x scRNA-seq + MERFISH | dbGaP phs002031 |
| 4 | Ren et al. 2023, *Sci Data* | Mouse | Heart (E8.5–P3, developmental) | ~500,000 cells | 10x scRNA-seq | GSE230531 |

---

## 2. Biological Relevance Assessment

### 2.1 Complementary Ontological Coverage

These 4 datasets are well-chosen. Together they span:

- **Developmental axis:** Fetal (Zhu, PCW 9–16) → Adult healthy (Tucker, Wang) → Adult diseased (Wang, heart failure/post-LVAD)
- **Technological axis:** scRNA-seq (cells) vs. snRNA-seq (nuclei) — critical for demonstrating portal versatility
- **Species axis:** 3 human + 1 mouse demonstrates cross-species visualization
- **Disease axis:** Healthy vs. heart failure enables disease-focused comparison views
- **Spatial axis:** Zhu et al. includes MERFISH spatial data, enabling spatial+transcriptomic integration

This is a strong template set. The data portal will need to handle both cell and nucleus-level measurements, which is a key differentiator from existing viewers.

### 2.2 Individual Dataset Notes

#### Wang et al. 2020 (*Nat Cell Biol*)
- **Strengths:** Direct comparison of healthy vs. failing vs. post-LVAD hearts — ideal for disease module demonstration
- **Cell types:** Atrial & ventricular CMs, fibroblasts, ECs, ACKR1+ ECs, immune cells — comprehensive non-CM profiling
- **Considerations:** 12 hearts, 21K cells is modest by modern standards. Chamber-specific CM heterogeneity may require careful annotation in the portal to avoid confounding cell type labels with chamber origin
- **Key biological draw:** Cell–cell communication networks in health vs. disease (NicheNet/CellChat outputs would be a natural visualization target)

#### Tucker et al. 2020 (*Circulation*)
- **Strengths:** 287K nuclei from 7 donors — largest of the 3 human datasets; snRNA-seq captures nuclear transcriptome including genes missed by dissociation (cardiomyocytes)
- **Cell types:** 9 major types including neurons, adipocytes, and SMCs — good coverage
- **Comment from Chien-Wei:** "easier to repeat, more likely to have cluster information uploaded" — this suggests processed data with annotations may be available, which simplifies portal ingestion
- **Considerations:** snRNA-seq has different gene detection biases vs. scRNA-seq. Nascent unspliced transcripts complicate UMI counting. Portal must label data type clearly (snRNA-seq vs. scRNA-seq)

#### Zhu et al. 2024 (*Nature*)
- **Strengths:** Fetal development atlas with 142K cells across 12 major classes / 39 populations / 75 subtypes — highest annotation resolution; MERFISH spatial data enables truly novel visualization modes
- **Key biological draw:** Spatial "cardiac niches" — co-localized cell communities. This dataset uniquely enables spatial visualization, which is a major selling point for the portal
- **Considerations:** Data access via dbGaP (controlled access) — may require IRB/data use certification. This could delay ingestion. Dryad and UCSC cell browser may have open-access processed matrices
- **MERFISH integration:** The portal should plan for paired spatial + transcriptomic visualization as a key demo feature

#### Ren et al. 2023 (*Sci Data*)
- **Strengths:** ~0.5M cells across 6 developmental stages — largest dataset by far; cardiac conduction system (CCS) annotation unique; processed data available (per Chien-Wei note: "This seems to have processed data available")
- **Cell types:** 20–26 types per stage with CCS-related subtypes — provides developmental trajectory visualization
- **Species consideration:** Mouse — allows cross-species comparison mode in the portal. Developmental staging (E8.5–P3) covers heart looping through neonatal maturation
- **Key biological draw:** Cardiac conduction system cell types — medically relevant and visually distinguishing (nodal cells vs. working myocardium)

---

## 3. Seurat Processing — Biological Considerations

### 3.1 Batch Effects & Integration

These 4 datasets originate from 4 different labs with different protocols. Major batch effect sources:

| Source | Impact | Mitigation |
|--------|--------|------------|
| scRNA-seq vs. snRNA-seq | High — different gene detection profiles, intronic reads in snRNA-seq | Process separately; do NOT integrate snRNA and scRNA without explicit cross-modality correction (e.g., Seurat v5 bridge integration) |
| 10x chemistry versions | Moderate — 3' v2 vs v3 vs 3' v3.1 have different UMI error profiles | Use `SCTransform` with `vst.flavor = "v3"` or regress out nCount_RNA |
| Tissue dissociation vs. nuclear isolation | High — affects cell type representation (CMs enriched in snRNA) | Document protocol for each dataset in portal metadata |
| Species (human vs. mouse) | High — gene orthology mapping needed | Use biomaRt/homologene for 1:1 ortholog mapping before cross-species integration |

**Recommendation:** Process each dataset independently in Seurat first (SCTransform → PCA → clustering → UMAP → annotation). Do not attempt to integrate across datasets initially. The data portal's value proposition is **simultaneous visualization of separate datasets**, not forced integration.

### 3.2 Species-Specific Normalization

- **Human datasets:** Use `SCTransform` with default settings. For snRNA-seq (Tucker), set `vst.flavor = "v2"` and consider including intronic counts if using cellranger 8+ count pipeline
- **Mouse dataset (Ren):** Use `SCTransform` with `mm10` reference. Developmental time-series may benefit from trajectory-aware normalization (e.g., `tradeSeq` or `slingshot`)
- **Cross-species:** If combined visualization is needed, map mouse genes to human orthologs (1:1 orthologs only, using Ensembl Compara). This reduces ~22K genes to ~15K homologous pairs

### 3.3 Mitochondrial Gene Handling

- scRNA-seq: High mitochondrial percentage indicates dying/damaged cells. Filter cells with >20% mtRNA (adjustable per dataset)
- snRNA-seq: Mitochondrial genes are expected to be low (nuclei have depleted mtRNA). Do NOT apply the same mtRNA filter to snRNA-seq data
- Portal should display mtRNA% as a QC metric in the viewer

### 3.4 Cell Type Annotation Harmonization

The 4 datasets use different cell type nomenclature and resolution. The portal needs a unified ontology:
- Recommend mapping to **Cell Ontology (CL)** terms
- At minimum, harmonize at the level of: CMs (atrial/ventricular), fibroblasts, endothelial cells, smooth muscle cells, pericytes, immune (myeloid/lymphoid), neuronal, adipocytes, mesothelial
- Zhu et al. has the finest resolution (75 subtypes) — this should be collapsed to ~12 major classes for portal consistency, with subtype available via drill-down

### 3.5 Ambient RNA & Doublets

- Wang et al. (heart failure samples) likely has elevated ambient RNA from tissue damage — apply `SoupX` or `cellbender` before Seurat
- Tucker et al. (snRNA-seq) has less ambient RNA but higher nuclear contamination rates
- All datasets should be run through `DoubletFinder` or `scDblFinder` as standard QC

---

## 4. Single-Cell Universe Integration Strategy

### 4.1 Conceptual Alignment

The "single-cell universe" concept — streamlining popular analytical tools and passing results to visualization — aligns naturally with exonBlock's cell × exon UMI matrix framework. The data portal becomes a downstream visualization layer for:

1. **Gene-level analysis** (standard Seurat workflow) → portal visualization
2. **Exon-level analysis** (exonBlock pipeline) → portal visualization with exon-resolution overlays

The portal should support two "views" per dataset: gene-level and exon-level.

### 4.2 Pipeline Architecture

```
Raw FASTQ → CellRanger → gene × cell matrix → Seurat → RDS (gene view)
                              ↓
                         BAM files → exonBlock → exon × cell matrix → Seurat → RDS (exon view)
                              ↓
                         Data Portal (Dr. Han) ← simultaneous visualization of both views
```

### 4.3 exonBlock-Specific Considerations

For the 4 selected datasets to work with exonBlock:

- **BAM file availability:** exonBlock requires aligned BAM files (not just count matrices). This is a critical data access requirement
  - Wang (PRJNA548584): SRA has raw FASTQ → need alignment via STAR or CellRanger count
  - Tucker (GSE109816): Similarly needs re-alignment from FASTQ
  - Zhu (dbGaP phs002031): Controlled access FASTQ → IRB required
  - Ren (GSE230531): Has raw data in GEO → can download and re-align
- **GENCODE version consistency:** Use GENCODE v44 (hg38) for human and vM34 (GRCm39) for mouse to match exonBlock's existing annotation files (`meta/features/`)
- **Cell barcode extraction:** exonBlock extracts CB/UMI from BAM tags. CellRanger output BAMs are the ideal input. If only FASTQ is available, run CellRanger count first

### 4.4 Recommended Pipeline for Each Dataset

| Dataset | Step 1 | Step 2 | Step 3 | Notes |
|---------|--------|--------|--------|-------|
| Wang | Download SRA FASTQ | CellRanger count | Seurat + exonBlock | Heart failure samples may need extra QC |
| Tucker | Download GEO FASTQ | CellRanger count (include-introns for snRNA) | Seurat + exonBlock | Use `--include-introns` flag |
| Zhu | Apply for dbGaP access | CellRanger count | Seurat + exonBlock | MERFISH data separate from scRNA |
| Ren | Download GEO FASTQ | CellRanger count | Seurat + exonBlock | Use GENCODE vM34 (not UCSC refGene) |

---

## 5. Nature Methods Manuscript Framing

### 5.1 Gap in the Field

The manuscript should position the data portal at the intersection of two unsolved problems:
1. **Visualization problem:** No existing tool allows simultaneous interactive exploration of multiple scRNA-seq datasets with harmonized cell-type overlays
2. **Resolution problem:** Current portals use gene-level matrices only; exon-level resolution (via exonBlock integration) reveals alternative splicing and isoform regulation invisible to gene-level analysis

### 5.2 Framing Recommendations

- **Title direction:** "An interactive data portal for multi-dataset visualization of single-cell transcriptomes with exon-level resolution"
- **Key positioning:** The portal is not just a viewer — it's an **analytical bridge** between tools (Seurat, exonBlock) and interpretable visualization
- **Nature Methods fit:** The journal values methods that fill infrastructure gaps. Emphasize that this is a *community resource* that reduces the barrier to multi-dataset comparison
- **Demo strategy:** Use the 4 cardiac datasets to show:
  1. Cross-dataset cell type comparison (development → adult → disease)
  2. scRNA vs. snRNA visualization mode
  3. Cross-species (human vs. mouse) developmental staging comparison
  4. Exon-level view showing alternative splicing differences across cell types (unique to exonBlock integration)

### 5.3 Biological Story Arc

The 4 datasets tell a coherent biological story that strengthens the manuscript:

> *"The mammalian heart undergoes dramatic transcriptional remodeling from embryonic development through adult homeostasis and into heart failure. We demonstrate our portal using four cardiac single-cell atlases spanning fetal development, adult health, and heart failure across human and mouse, enabling interactive comparison of cell-type composition, gene expression, and exon-level regulation across life stages, disease states, and species."*

### 5.4 Anticipated Reviewer Concerns

- **"Why these specific datasets?"** → Deliberate selection to demonstrate portal flexibility across developmental stage, disease state, measurement modality, and species
- **"How is this different from UCSC Cell Browser / CZ CellxGene?"** → Those show one dataset at a time; ours enables simultaneous multi-dataset comparison with exon-level resolution
- **"Is exon-level resolution useful in cardiac data?"** → Cardiomyocytes are known for extensive alternative splicing (titin, RYR2, SCN5A); exonBlock reveals isoform switching between healthy and failing hearts — a key demonstration

---

## 6. Action Items & Recommendations

### Immediate (for data processing — [OMI-35](/OMI/issues/OMI-35))
1. Download all 4 datasets; prioritize Tucker (GSE109816) and Ren (GSE230531) which have open processed data
2. Zhu et al. (dbGaP) requires controlled access — flag to CEO for IRB/data use agreement
3. Use CellRanger count with `--include-introns` for Tucker snRNA-seq
4. Use GENCODE v44/vM34 annotations consistently for all alignments

### Short-term (portal development)
1. Build unified cell-type ontology across all 4 datasets
2. Support both scRNA and snRNA modalities with clear labeling
3. Implement exon-level view layer via exonBlock RDS integration
4. MERFISH spatial data from Zhu should be a showcase feature

### Manuscript preparation
1. Lead with the multi-dataset simultaneous visualization gap
2. Use cardiac development-disease trajectory as the biological narrative
3. Exon-level splice variant differences (titin, RYR2, SCN5A) as the novel resolution demo
4. Target Nature Methods — emphasize community resource value