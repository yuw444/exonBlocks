#' Create Exon Clusters for a Species
#'
#' @description
#' Fetches exon coordinates from TxDb for the specified species and
#' creates exon clusters by merging overlapping/adjacent exons using
#' \code{\link[GenomicRanges:reduce]{GenomicRanges::reduce()}}.
#'
#' @param species Species genome assembly. Options: "hg38" (human, default),
#'   "mm10" (mouse), "rn6" (rat), "dm6" (fly), "ce11" (worm).
#'   Must have corresponding TxDb package installed (e.g., TxDb.Hsapiens.UCSC.hg38.knownGene).
#' @param min_gap_width Integer. Exons separated by gap <= this distance are merged
#'   into the same cluster. Default: 10. Set to \code{NA} to only merge overlapping exons.
#' @param verbose Print progress messages? Default: TRUE.
#' @return A \code{GRanges} object with exon clusters. The metadata column
#'   \code{revmap} contains indices of original exons that were merged into each cluster.
#' @importFrom GenomicRanges GRanges reduce
#' @importFrom IRanges IRanges
#' @export
#' @examples
#' \dontrun{
#' # Human hg38 with default min_gap_width=10
#' clusters <- make_exon_clusters(species = "hg38")
#'
#' # Mouse mm10 with stricter clustering
#' clusters <- make_exon_clusters(species = "mm10", min_gap_width = 5)
#'
#' # Only merge overlapping exons
#' clusters <- make_exon_clusters("hg38", min_gap_width = NA)
#' }
make_exon_clusters <- function(
    species = "hg38",
    min_gap_width = 10,
    verbose = TRUE
) {
    # -------------------------------------------------------------------------
    # Species to TxDb package mapping
    # -------------------------------------------------------------------------
    txdb_map <- list(
        "hg38" = "TxDb.Hsapiens.UCSC.hg38.knownGene",
        "hg19" = "TxDb.Hsapiens.UCSC.hg19.knownGene",
        "mm10" = "TxDb.Mmusculus.UCSC.mm10.knownGene",
        "mm9"  = "TxDb.Mmusculus.UCSC.mm9.knownGene",
        "rn6"  = "TxDb.Rnorvegicus.UCSC.rn6.refGene",
        "dm6"  = "TxDb.Dmelanogaster.UCSC.dm6.ensGene",
        "ce11" = "TxDb.Celegans.UCSC.ce11.ensGene"
    )

    species <- tolower(species)

    if (!species %in% names(txdb_map)) {
        stop(
            "Unsupported species: ", species, "\n",
            "Available species: ", paste(names(txdb_map), collapse = ", "), "\n",
            "Install required package with: BiocManager::install('", txdb_map[[species]], "')"
        )
    }

    txdb_pkg <- txdb_map[[species]]

    # -------------------------------------------------------------------------
    # Load TxDb
    # -------------------------------------------------------------------------
    if (verbose) message("Loading TxDb for ", species, "...")

    if (!requireNamespace(txdb_pkg, quietly = TRUE)) {
        stop("Package ", txdb_pkg, " not installed. Install with: BiocManager::install('", txdb_pkg, "')")
    }
    txdb <- get(txdb_pkg, envir = loadNamespace(txdb_pkg))

    # -------------------------------------------------------------------------
    # Extract exons
    # -------------------------------------------------------------------------
    if (verbose) message("Extracting exons...")

    gr_exons <- GenomicFeatures::exons(
        txdb,
        columns = c("exon_id", "gene_id", "tx_id")
    )

    n_exons <- length(gr_exons)

    if (verbose) message("  Found ", n_exons, " exons")

    # -------------------------------------------------------------------------
    # Reduce to clusters
    # -------------------------------------------------------------------------
    if (verbose) {
        if (is.na(min_gap_width)) {
            message("Reducing to clusters (overlapping only)...")
        } else {
            message("Reducing to clusters (min_gap_width = ", min_gap_width, ")...")
        }
    }

    # For GenomicRanges::reduce, min.gapwidth must be an integer
    if (is.na(min_gap_width)) {
        exon_clusters <- GenomicRanges::reduce(
            gr_exons,
            with.revmap = TRUE
        )
    } else {
        exon_clusters <- GenomicRanges::reduce(
            gr_exons,
            with.revmap = TRUE,
            min.gapwidth = as.integer(min_gap_width)
        )
    }

    n_clusters <- length(exon_clusters)

    if (verbose) {
        message("Done! Created ", n_clusters, " clusters from ", n_exons, " exons")
    }

    return(exon_clusters)
}
