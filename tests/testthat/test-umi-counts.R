library(data.table)

test_that("fread reads umicounts file", {
  skip_on_cran()
  umi_file <- test_path("..", "..", "meta", "features", "PBMCs.allruns.umicounts_intronexon.txt")
  skip_if_not(file.exists(umi_file), "UMI counts file not found")

  temp <- data.table::fread(
    umi_file,
    header = TRUE,
    nrows = 1000
  )

  expect_s3_class(temp, "data.table")
  expect_true(nrow(temp) > 0)
})