#' Convert BAM to CellRanger-style matrix via SQLite
#'
#' Processes a BAM file, extracts cell barcodes and UMIs, maps reads to
#' genomic regions (genes, exons, or custom BED), and outputs a sparse
#' count matrix in 10x format.
#'
#' @param bam_file Path to BAM file (must be indexed with .bai)
#' @param db_file Path to output SQLite database
#' @param path_out Path to output directory for TSV/MTX files
#' @param barcodes_file Path to gzipped barcodes file (one per line)
#' @param regions_file Path to 8-col TSV: id, symbol, chr, start, end, annotation1, annotation2
#' @param rate_cell Fraction of cells to sample (0, 1]
#' @param rate_depth Fraction of reads to sample (0, 1]
#' @param seed Random seed for reproducibility
#' @param umi_copies_flag If 1, also output per-UMI copy counts
#' @param protocol Protocol type: "10X" or "zUMIs"
#' @param region_chr Optional chromosome for ROI filter (e.g. "chr19")
#' @param region_start Optional start position (1-based)
#' @param region_end Optional end position (1-based)
#' @return NULL (invisibly). Output written to db_file and path_out.
#' @export
bam2db <- function(
    bam_file,
    db_file,
    path_out,
    barcodes_file,
    regions_file,
    rate_cell = 1.0,
    rate_depth = 1.0,
    seed = 42,
    umi_copies_flag = 0L,
    protocol = c("10X", "zUMIs"),
    region_chr = NULL,
    region_start = 0L,
    region_end = 0L) {

  stopifnot(file.exists(bam_file))
  stopifnot(file.exists(barcodes_file))
  stopifnot(file.exists(regions_file))
  stopifnot(rate_cell > 0 && rate_cell <= 1)
  stopifnot(rate_depth > 0 && rate_depth <= 1)

  protocol <- match.arg(protocol)
  proto_int <- if (protocol == "10X") 0L else 1L

  if (!is.null(region_chr)) {
    region_chr <- as.character(region_chr)
    region_start <- as.integer(region_start)
    region_end <- as.integer(region_end)
  } else {
    region_chr <- character(0)
    region_start <- 0L
    region_end <- 0L
  }

  invisible(.Call(
    `_exonBlocks_bam2db`,
    bam_file,
    db_file,
    path_out,
    barcodes_file,
    regions_file,
    as.numeric(rate_cell),
    as.numeric(rate_depth),
    as.integer(seed),
    as.integer(umi_copies_flag),
    proto_int,
    region_chr,
    region_start,
    region_end
  ))
}
