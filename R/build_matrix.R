#' Build Cell × Cluster Matrix from BAM (CellRanger-style output)
#'
#' @description
#' Streams BAM file and outputs CellRanger-compatible files:
#' - barcodes.tsv.gz: Cell barcodes
#' - features.tsv.gz: Exon cluster coordinates
#' - matrix.mtx.gz: Sparse matrix in Matrix Market format
#'
#' @param bam Path to input BAM/CRAM with CB/UB tags
#' @param exon_clusters GRanges object with exon cluster coordinates
#' @param out_dir Output directory for CellRanger-style files
#' @param xf_values Integer vector of allowed xf values (10X only)
#' @param tech Character: "10X" or "Smart-seq"
#' @return List with n_cells, n_features, nnz, and paths
#' @export
#' @aliases build_cell_cluster_matrix
#' @examples
#' \dontrun{
#' result <- bam_to_cellranger(
#'   bam = "possorted.bam",
#'   exon_clusters = exon_clusters,
#'   out_dir = "cellranger_output"
#' )
#' }
bam_to_cellranger <- function(
    bam,
    exon_clusters,
    out_dir,
    xf_values = c(25L, 17L),
    tech = "10X"
) {
    # Create output directory
    if (!dir.exists(out_dir)) {
        dir.create(out_dir, recursive = TRUE)
    }

    # Get cluster data
    clusters <- GenomicRanges::granges(exon_clusters)
    cluster_chr <- as.character(GenomicRanges::seqnames(clusters))
    cluster_start <- GenomicRanges::start(clusters)
    cluster_end <- GenomicRanges::end(clusters)
    n_clusters <- length(clusters)

    # Get unique chromosomes
    unique_chr <- unique(cluster_chr)

    # Process each chromosome
    for (chrom in unique_chr) {
        message("Processing chromosome: ", chrom)

        # Get clusters for this chromosome
        chr_mask <- cluster_chr == chrom
        chr_start <- cluster_start[chr_mask]
        chr_end <- cluster_end[chr_mask]

        if (length(chr_start) == 0) next

        chr_min <- min(chr_start)
        chr_max <- max(chr_end)

        # Call C function for this chromosome
        result <- .Call(
            `_exonBlocks_bam_to_cellranger`,
            bam,
            chr,  # single chromosome
            as.integer(chr_start),
            as.integer(chr_end),
            "",  # cluster_id (optional)
            tech,
            as.integer(xf_values),
            out_dir
        )

        message("  Processed: ", result$n_cells, " cells, ", result$n_features, " features")
    }

    invisible(list(
        out_dir = out_dir,
        barcodes_path = file.path(out_dir, "barcodes.tsv.gz"),
        features_path = file.path(out_dir, "features.tsv.gz"),
        matrix_path = file.path(out_dir, "matrix.mtx.gz")
    ))
}


#' Convert Cell × Cluster Matrix to Standard Formats
#'
#' @param mat dgCMatrix from bam_to_cellranger
#' @param format Output format: "10x" (HDF5), "loom", "arrow", "csv"
#' @param path Output file path
#' @param cell_metadata Optional data.frame with cell annotations (must have rowname = CB)
#' @return Invisible NULL (writes to file)
#' @export
#' @examples
#' \dontrun{
#' export_matrix(mat, format = "loom", path = "matrix.loom")
#' }
export_cell_cluster_matrix <- function(
    mat,
    format = c("10x", "loom", "arrow", "csv"),
    path,
    cell_metadata = NULL
) {
    format <- match.arg(format)

    switch(format,
           "10x" = {
               if (!requireNamespace("zellkonverter", quietly = TRUE)) {
                   stop("zellkonverter package required for 10x format. ",
                        "Install with: remotes::install_github('cellgeni/zellkonverter')")
               }
               # zellkonverter writes H5AD/SCE from matrices
               zellkonverter::writeH5AD(mat, path)
           },
           "loom" = {
               if (!requireNamespace("loomR", quietly = TRUE)) {
                   stop("loomR package required for loom format")
               }
               # Create loom from scratch
               loomR::create_loom(path, matrix = mat, overwrite = TRUE)
           },
           "arrow" = {
               if (!requireNamespace("arrow", quietly = TRUE)) {
                   stop("arrow package required for arrow format")
               }
               # Export as Arrow IPC/Feather
               arrow::write_feather(
                   as.data.frame(mat),
                   path
               )
           },
           "csv" = {
               # Write as dense CSV (warning for large matrices)
               if (ncol(mat) > 1000) {
                   warning("CSV export of large matrices may be slow")
               }
               write.csv(as.matrix(mat), file = path)
           })

    invisible(NULL)
}


#' Save Exon Cluster Index for C-level Processing
#'
#' Serializes exon clusters to a compact binary format for fast C lookups.
#' This function is optional - used for optimized batch processing.
#'
#' @param exon_clusters GRanges object with exon clusters
#' @param path Output file path (.rds format)
#' @param overwrite Overwrite existing file?
#' @return Invisible NULL
#' @export
save_cluster_index <- function(exon_clusters, path, overwrite = FALSE) {
    if (file.exists(path) && !overwrite) {
        stop("File exists. Use overwrite = TRUE to replace.")
    }

    # Convert to compact list
    idx <- list(
        chr = as.character(GenomicRanges::seqnames(exon_clusters)),
        start = GenomicRanges::start(exon_clusters),
        end = GenomicRanges::end(exon_clusters),
        n_clusters = length(exon_clusters)
    )

    saveRDS(idx, path)
    message("Cluster index saved to: ", path)
    invisible(NULL)
}
