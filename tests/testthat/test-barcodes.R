test_that("specific barcode coverage calculation", {
  skip_on_cran()
  skip_if_not_installed("GenomicAlignments")
  skip_if_not_installed("Rsamtools")
  library(dplyr)
  library(GenomicAlignments)
  library(Rsamtools)
  library(ggplot2)

  bam_file <- test_path("..", "..", "meta", "bam", "test.bam")
  skip_if_not(file.exists(bam_file), "test.bam not found")

  chrom <- "1"
  region_start <- 3600000
  region_end <- 3750000
  barcode <- "GCGCAGTGTGATGTCT-1"

  region <- GRanges(chrom, IRanges(region_start, region_end))
  param <- ScanBamParam(which = region, tag = "CB")

  bam <- readGAlignments(bam_file, param = param)
  cb_vals <- mcols(bam)$CB
  filtered <- bam[!is.na(cb_vals) & cb_vals == barcode]
  skip_if_not(length(filtered) > 0, "No reads found for barcode in this region")

  coverage <- coverage(filtered)
  cov_df <- as.data.frame(IRanges::as.data.frame(coverage[[chrom]]))
  if (ncol(cov_df) == 1) {
    cov_df$group <- 1
  }
  colnames(cov_df) <- c("coverage", "group")
  cov_df$position <- as.numeric(rownames(cov_df))
  cov_df <- cov_df %>% filter(position >= region_start & position <= region_end)

  p <- ggplot(cov_df, aes(x = position, y = coverage)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    labs(x = "Position", y = "Coverage", title = "Coverage Plot for Specific Barcode") +
    theme_minimal()

  print(p)

  expect_gt(length(filtered), 0, label = "No reads found for the barcode in this region")
  expect_gt(nrow(cov_df), 0, label = "Coverage calculation returned no data")
})