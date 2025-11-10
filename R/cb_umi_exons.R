#' @title Map Cell Barcode + UMI to Overlapping Exons
#' @description Expand semicolon-delimited block fields and identify exons that overlap
#'   each block for every (CB, UMI) pair. The function reads a blocks TSV (produced by
#'   extract_exon_reads_hts()) where block_start, block_end and block_seq may contain
#'   multiple entries separated by ";" and a meta exons TSV describing exon coordinates.
#' @details The function:
#'   1. reads the blocks TSV and expands rows with semicolon-separated fields using
#'      tidyr::separate_rows();
#'   2. coerces block coordinates to integer;
#'   3. reads the exon metadata and selects the required columns;
#'   4. uses data.table::foverlaps() to find any overlap between blocks and exons.
#' @param tsv Character. Path to the blocks TSV. Required columns: CB, UMI,
#'   block_start, block_end, block_seq (block_* may contain \";\"-separated lists).
#' @param meta_exons Character. Path to the exon metadata TSV. Required columns:
#'   chr, start, end, exon (or exon_id). start/end must be integer coordinates.
#' @return A data.table of overlapping records. Columns include at minimum:
#'   CB, UMI, block_start, block_end, start, end, exon (columns produced by foverlaps).
#'   Each row represents one block-entry (after expansion) overlapping one exon.
#' @import data.table dplyr tidyr
#' @export
#' @examples
#' \dontrun{
#'   blocks_tsv <- "AI.tsv"
#'   exons_tsv  <- "meta_exons.tsv"
#'   res <- cb_umi_exons(tsv = blocks_tsv, meta_exons = exons_tsv)
#'   head(res)
#' }
#' @seealso data.table::foverlaps, tidyr::separate_rows
cb_umi_exons <- function(
    tsv,
    meta_exons
) {
    df_blocks <- data.table::fread(
        tsv,
        header = TRUE
    ) %>%
        as.data.frame() %>%
        tidyr::separate_rows(block_start, block_end, block_seq, sep = ";") %>%
        dplyr::mutate(
            block_start = as.integer(block_start),
            block_end = as.integer(block_end)
        ) %>%
        as.data.table()
    
    df_exons <- data.table::fread(
        meta_exons,
        header = TRUE
    ) %>%
        as.data.frame() %>%
        dplyr::select(chr, start, end, exon) %>%
        as.data.table()
    
    setDT(df_blocks, key = c("block_start", "block_end"))
    setDT(df_exons, key = c("start", "end"))
    
    # Merge blocks with exons based on overlap
    df_merged <- foverlaps(
        df_blocks[, .(CB, UMI, block_start, block_end)],
        df_exons[, .(start, end, exon)],
        by.x = c("block_start", "block_end"),
        by.y = c("start", "end"),
        type = "any",
        nomatch = 0L
    )
    
    return(df_merged)
}