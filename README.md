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