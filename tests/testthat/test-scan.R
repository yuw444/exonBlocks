test_that("cb_umi_exons maps blocks to overlapping exons", {
  skip_on_cran()
  library(dplyr)
  library(exonBlocks)

  block_tsv_1 <- test_path("..", "..", "meta", "features", "AI.tsv")
  exon_csv <- test_path("..", "..", "meta", "features", "Xkr4_exons_annotated.csv")

  skip_if_not(file.exists(block_tsv_1), "AI.tsv not found")
  skip_if_not(file.exists(exon_csv), "Xkr4_exons_annotated.csv not found")

  df_exon <- read.table(
    exon_csv,
    header = TRUE,
    sep = ",",
    stringsAsFactors = FALSE
  ) %>%
    mutate(length = end - start + 1) %>%
    filter(!is.na(exon)) %>%
    arrange(exon)

  df <- cb_umi_exons(
    tsv = block_tsv_1,
    meta_exons = exon_csv
  ) %>%
    mutate(length = end - start + 1) %>%
    filter(!is.na(exon))

  expect_s3_class(df, "data.table")
  expect_true("exon" %in% names(df))

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

  expect_true(nrow(df_sum) > 0)
})