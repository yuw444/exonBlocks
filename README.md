# exonBlock

**Cell × Exon UMI Matrix Framework for scRNA-seq**

## Overview

exonBlock converts the standard cell × gene UMI matrix used in single-cell RNA-seq studies into a **cell × exon UMI matrix**, enhancing resolution to reveal exon-level regulation, alternative splicing, and cell subtype heterogeneity that is invisible to gene-level analysis.

## Motivation

Current scRNA-seq analyses rely on cell × gene UMI matrices that collapse mRNA transcript diversity and obscure exon- and isoform-level regulation. Preliminary data from iPSC knockout experiments demonstrate that gene-level expression can appear unchanged despite disruption of specific exons — highlighting a critical limitation of existing pipelines.

## Key Features

- Efficient pipeline to generate **cell × exon UMI matrices** from short-read sequencing
- Modular, user-friendly software platform (analogous to Seurat for gene-level analysis)
- Supports normalization, differential exon usage, and isoform inference
- Designed for scalability to biobank-level datasets

## Status

- [ ] Pipeline development (in progress)
- [ ] Exon-level normalization methods
- [ ] Differential exon usage framework
- [ ] Isoform inference module
- [ ] Software documentation

## Related Projects

- **K99/R00** — NIH grant proposal using exonBlock as the computational framework for Aim 2
- **MoSAIC** — HMM-based mCA detection framework (provides mCA calls for integration in Aim 3)
## Test Data File Hierarchy

```
meta/
├── bam/
│   ├── human_pbmc_Smart-seq3xpressV3.bam    — Human PBMC Smart-seq3 test data (gitignored)
│   ├── human_pbmc_Smart-seq3xpressV3.bam.bai
│   ├── mice_brain_10X.bam                   — Mouse brain 10X test data (gitignored)
│   ├── mice_brain_10X.bam.bai
│   ├── test.bam                              — Small test BAM
│   └── test.bam.bai
├── features/
│   ├── homo_sapiens_GRCh38_gene.tsv          — Human gene annotations (GENCODE v44)
│   ├── homo_sapiens_GRCh38_exon.tsv          — Human exon annotations (GENCODE v44)
│   ├── mus_musculus_GRCm39_gene.tsv         — Mouse gene annotations (GENCODE vM34)
│   ├── mus_musculus_GRCm39_exon.tsv          — Mouse exon annotations (GENCODE vM34)
│   ├── human_barcodes_chr19.txt              — Human barcodes for chr19 region
│   └── mouse_barcodes_chr19.txt              — Mouse barcodes for chr19 region
└── output/                                    — Test output directory (gitignored)
```

## Test Paths

All test scripts have been updated to use local paths under `meta/` instead of the old HPC paths (`/scratch/g/chlin/Yu/exonBlocks/`). Tests use `skip_if_not(file.exists(...))` to gracefully handle missing data files.


typedef struct {
    size_t n;
    int32_t *start;
    int32_t *end;
    int32_t *gene_id;
    int32_t *exon_cluster_id;
} ExonClusterIndex;
typedef struct {
    size_t n;
    int32_t *start;
    int32_t *end;
    int32_t *gene_id;
} GeneIndex;
typedef struct {
    size_t n_contigs;
    ExonClusterIndex *exon_index;
    GeneIndex *gene_index;
} GenomeAnnotation; 

| overlap exon cluster | within gene | interpretation |
| ------------ | ------------ | -------------- |
| yes          | yes          | exonic         |
| no           | yes          | unspliced      |
| yes          | no           | impossible     |
| no           | no           | intergenic     | 


same time, let's update all the notes about classify rule of our package using bit field follow

NSEIG 

N: has N in CIGAR
S: all gene-mapped blocks belong to same gene
E: at least one exon block present
I: at least one intron block present
G: at least one intergenic block present

H: has N in CIGAR
U: uniqueness of read blocks mapped gene
G: all blocks within gene region
E: all blocks within exon cluster

1111 ->  SPLICED
1011 -> CHIMERIC
_110 -> UNSPLICED
Other -> UNASSIGNED


