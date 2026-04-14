#' Extract Exon Reads from BAM using htslib
#'
#' @description
#' Scans a BAM/CRAM file for a specified genomic region, extracts read blocks
#' with CB/UB tags, and writes a TSV of block coordinates per (CB, UMI) pair.
#' Optionally writes a filtered BAM containing only overlapping reads.
#'
#' @param bam Path to input BAM/CRAM
#' @param chr Region chrom/contig (e.g., "chr1")
#' @param start 1-based inclusive start position
#' @param end 1-based inclusive end position
#' @param out_bam Path for filtered BAM to write (or NULL to skip)
#' @param tsv Path to TSV output file
#' @param xf_values Integer vector of allowed xf values (10X only). Default: c(25L, 17L)
#' @param tech Character; sequencing technology, one of "10X" (default) or "Smart-seq"
#' @return Character vector of TSV output path (invisibly)
#' @useDynLib exonBlocks, .registration = TRUE
#' @export
extract_exon_reads_hts <- function(
      bam,
      chr,
      start,
      end,
      out_bam,
      tsv,
      xf_values = c(25L, 17L),
      tech = "10X"
) {
      if (is.null(out_bam)) {
            out_bam <- ""
      }
      stopifnot(file.exists(bam), nzchar(chr), start <= end)
      .Call(
            `_exonBlocks_scan_bam_blocks_hts`,
            bam,
            chr,
            as.integer(start),
            as.integer(end),
            out_bam,
            tsv,
            tech,
            as.integer(xf_values)
      )
}
