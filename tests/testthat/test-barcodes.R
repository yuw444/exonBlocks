test_that("specific barcode coverage calculation", {
  library(dplyr)
  library(GenomicAlignments)
  library(Rsamtools)
  library(ggplot2)

  # Input parameters
  chrom <- "1"
  region_start <- 10003
  region_end <- 1000000
  bam_file <- "/home/yu-wang/Documents/exonBlocks/meta/bam/test.bam"
  barcode <- "TGTCAAGAGGCCAACCTATA"

  # Setup region and barcode
  region <- GRanges(chrom, IRanges(region_start, region_end))
  param <- ScanBamParam(which = region, tag = "CB")

  # Read and filter BAM
  bam <- readGAlignments(bam_file, param = param)
  filtered <- bam[mcols(bam)$CB == barcode]

  # Calculate coverage
  coverage <- coverage(filtered)
  cov_df <- as.data.frame(IRanges::as.data.frame(coverage[[chrom]]))
  colnames(cov_df) <- c("position", "coverage")
  cov_df <- cov_df %>% filter(as.numeric(position) >= region_start & as.numeric(position) <= region_end)

  # Plot coverage
  p <- ggplot(cov_df, aes(x = as.numeric(position), y = coverage)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    labs(x = "Position", y = "Coverage", title = "Coverage Plot for Specific Barcode") +
    theme_minimal()

  # Save or print plot
  print(p)

  # Test conditions
  expect_gt(nrow(filtered), 0, "No reads found for the barcode in this region")
  expect_gt(nrow(cov_df), 0, "Coverage calculation returned no data")
})