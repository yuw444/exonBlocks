test_that("multiplication works", {
  library(dplyr)
  library(exonBlocks)
  chrom = "4"
  exon_start = 150920155
  exon_end = 150946102
  # bam_1 = "/scratch/g/chlin/Yu/CD137/data/AI/possorted_genome_bam.bam"
  bam_1 = "/scratch/g/chlin/Yu/exonBlocks/sub/sh/aligned_Aligned.sortedByCoord.out.bam"
  block_tsv_1 = "/scratch/g/chlin/Yu/exonBlocks/meta/AI.tsv"

  temp <- extract_exon_reads_hts(
    bam = bam_1,
    chr = chrom,
    start = exon_start,
    end = exon_end,
    out_bam = NULL,
    tsv = block_tsv_1,
    # xf_values = c(25L, 17L)
    tech = "Smart-seq"
  )
  
  df_exon <- read.table(
    "/scratch/g/chlin/Yu/exonBlocks/meta/Tnfrsf9_exons_annotated.csv",
    header = TRUE,
    sep = ",",
    stringsAsFactors = FALSE
  ) %>%
    mutate(length = end - start + 1) %>%
    filter(!is.na(exon)) %>%
    arrange(exon)

  df <- cb_umi_exons(
    tsv = block_tsv_1,
    meta_exons = "/scratch/g/chlin/Yu/exonBlocks/meta/Tnfrsf9_exons_annotated.csv"
  ) %>%
    mutate(length = end - start + 1) %>%
    filter(!is.na(exon))

  df %>%
    group_by(exon) %>%
    summarise(
      exon = n(),
      length = unique(length)
    )

  df %>%
    filter(!is.na(exon)) %>%
    group_by(CB, .drop = TRUE) %>%
    summarise(
      exon_ids = paste(sort(unique(exon)), collapse = ";")
    ) -> df_sum

  table(df_sum$exon_ids)

  
})
