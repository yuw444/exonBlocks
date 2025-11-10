# exonBlocks

Map cell barcode + UMI-level read blocks to overlapping exons. This repository provides utilities to extract read blocks from BAMs and assign each block (expanded from semicolon-delimited fields) to annotated exons.

## Requirements

- R (>= 4.0)
- Packages: data.table, dplyr, tidyr, testthat
- Access to an indexed BAM file

## Installation

Install required packages in R:

```r
if(!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
if(!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if(!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr")
devtools::install_github("https://github.com/yuw444/exonBlocks")
```

## Example usage

The following example is adapted for the manuscript analyzing Tnfrsf9 (CD137) exons in single-cell RNA-seq data. It:
1. extracts read blocks overlapping a genomic interval into a TSV (extract_exon_reads_hts),
2. expands semicolon-delimited block fields and maps blocks to exons (cb_umi_exons),
3. summarizes exon assignments per (CB, UMI).

```r
library(exonBlocks)
library(dplyr)

chrom <- "4"
exon_start <- 150920155
exon_end <- 150946102
bam <- "/scratch/g/chlin/Yu/CD137/data/AI/possorted_genome_bam.bam"
block_tsv <- "AI.tsv"

# extract read blocks overlapping the region -> writes AI.tsv
temp <- extract_exon_reads_hts(
  bam = bam,
  chr = chrom,
  start = exon_start,
  end = exon_end,
  out_bam = NULL,
  tsv = block_tsv
)

# map expanded blocks to annotated exons
df <- cb_umi_exons(
  tsv = block_tsv,
  meta_exons = "/scratch/g/chlin/Yu/exonBlocks/Tnfrsf9_exons_annotated.csv"
)

# summarize exon ids per (CB, UMI)
df_sum <- df %>%
  filter(!is.na(exon)) %>%
  group_by(CB, UMI, .drop = TRUE) %>%
  summarise(exon_ids = paste(sort(unique(exon)), collapse = ";"))

table(df_sum$exon_ids)
```

## Notes

- cb_umi_exons expects `block_start`, `block_end`, `block_seq` (semicolon-delimited allowed) and `CB`, `UMI` columns in the blocks TSV.
- meta exons file must contain columns: `chr`, `start`, `end`, `exon`.
- The implementation uses tidyr::separate_rows to expand semicolon lists and data.table::foverlaps to detect overlaps.
