#' htslib-core extractor (fast path)
#' @useDynLib exonBlocks, .registration = TRUE
#' @param bam Path to input BAM/CRAM
#' @param chr Region chrom/contig (e.g., "chr1")
#' @param start 1-based inclusive
#' @param end 1-based inclusive
#' @param out_bam Path for filtered BAM to write (or NULL to skip)
#' @param tsv Path to TSV to write
#' @param xf_values integer vector; default c(25L,17L)
#' @export
extract_exon_reads_hts <- function(bam, chr, start, end, out_bam, tsv,
                                   xf_values = c(25L,17L)) {
  if (is.null(out_bam)) out_bam <- ""
  stopifnot(file.exists(bam), nzchar(chr), start <= end)
  .Call(`_exonBlocks_scan_bam_blocks_hts`,
        bam, chr, as.integer(start), as.integer(end),
        out_bam, tsv, as.integer(xf_values))
}
