#' Prepare Exon and Intron Annotation for Read Classification
#'
#' @description
#' Extracts exon and intron coordinates with gene IDs from a TxDb package
#' and prepares exon clusters. Uses \code{GenomicFeatures::exons(columns="gene_id")}
#' for direct exon-to-gene mapping and \code{GenomicFeatures::intronsBy(by="gene")}
#' for intron coordinates — no per-gene loops.
#'
#' @param species Character. Genome assembly: "hg38", "hg19", "mm10", "mm9", "rn6", "dm6", "ce11".
#' @param min_gap_width Integer. Exons within this distance are merged into same cluster. Default: 10.
#' @param chromosomes Character vector. Restrict to these chromosomes (e.g. "chr19"). NULL = all.
#' @param verbose Logical. Print progress. Default: TRUE.
#' @return A list with components: exons (data.frame), introns (data.frame),
#'   exon_clusters (GRanges), gene_ids, gene_to_idx.
#' @importFrom GenomicRanges GRanges reduce seqnames start end
#' @importFrom IRanges IRanges
#' @importFrom GenomicFeatures exons
#' @export
prepare_classification_annotation <- function(
    species = "hg38",
    min_gap_width = 10,
    chromosomes = NULL,
    verbose = TRUE
) {
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
        stop("Unsupported species: ", species)
    }

    txdb_pkg <- txdb_map[[species]]
    if (verbose) message("Loading TxDb: ", txdb_pkg)

    if (!requireNamespace(txdb_pkg, quietly = TRUE)) {
        stop("Package ", txdb_pkg, " not installed. Install with: BiocManager::install('", txdb_pkg, "')")
    }
    txdb <- get(txdb_pkg, envir = loadNamespace(txdb_pkg))

    # Get exons with gene_id column — single GRanges, no loop
    if (verbose) message("Extracting exons with gene IDs...")
    gr_exons <- GenomicFeatures::exons(txdb, columns = "gene_id")

    # gene_id is an IntegerList ( exon can map to multiple genes )
    exon_chr <- as.character(GenomicRanges::seqnames(gr_exons))
    exon_start <- GenomicRanges::start(gr_exons)
    exon_end <- GenomicRanges::end(gr_exons)
    exon_strand <- as.character(BiocGenerics::strand(gr_exons))
    gene_id_list <- S4Vectors::mcols(gr_exons)$gene_id

    # Expand multi-gene exons: one row per (exon, gene) pair
    gene_ids_vec <- unlist(gene_id_list)
    exon_rep <- rep(seq_along(gene_id_list), lengths(gene_id_list))

    exon_df <- data.frame(
        chr = exon_chr[exon_rep],
        start = exon_start[exon_rep],
        end = exon_end[exon_rep],
        strand = exon_strand[exon_rep],
        gene_idx = gene_ids_vec,
        stringsAsFactors = FALSE
    )

    gene_ids <- as.character(sort(unique(gene_ids_vec)))

    # Filter exons to requested chromosomes
    if (!is.null(chromosomes)) {
        keep <- exon_df$chr %in% chromosomes
        exon_df <- exon_df[keep, ]
    }

    # Compute introns from exon gaps using data.table — no per-gene loop
    if (verbose) message("Computing introns from exon gaps...")
    dt <- data.table::as.data.table(exon_df)
    data.table::setorder(dt, chr, strand, gene_idx, start, end)

    # Reduce overlapping/adjacent exons per (chr, strand, gene_idx): group into blocks
    # A new block starts when start > running max of previous ends
    dt[, prev_runmax := data.table::shift(cummax(as.integer(end)), type = "lag", fill = 0L),
       by = .(chr, strand, gene_idx)]
    dt[, new_block := as.integer(start) > prev_runmax, by = .(chr, strand, gene_idx)]
    dt[, block_id := cumsum(new_block), by = .(chr, strand, gene_idx)]
    exons_reduced <- dt[, .(start = min(start), end = max(end)), by = .(chr, strand, gene_idx, block_id)]
    exons_reduced[, block_id := NULL]

    # Introns = gaps between consecutive reduced exons per (chr, strand, gene_idx)
    data.table::setorder(exons_reduced, chr, strand, gene_idx, start)
    exons_reduced[, next_start := data.table::shift(start, type = "lead"), by = .(chr, strand, gene_idx)]
    exons_reduced[, intron_start := as.integer(end) + 1L]
    exons_reduced[, intron_end := as.integer(next_start) - 1L]
    intron_dt <- exons_reduced[!is.na(next_start) & intron_start <= intron_end,
                               .(chr, strand, start = intron_start, end = intron_end, gene_idx)]

    intron_df <- as.data.frame(intron_dt)
    if (nrow(intron_df) == 0) {
        intron_df <- data.frame(chr = character(0), strand = character(0),
                                start = integer(0), end = integer(0),
                                gene_idx = integer(0),
                                stringsAsFactors = FALSE)
    }

    if (verbose) message("  ", nrow(exon_df), " exon rows, ", nrow(intron_df), " intron rows")

    # Build exon clusters
    if (verbose) message("Building exon clusters...")
    exon_clusters <- make_exon_clusters(species = species, min_gap_width = min_gap_width, verbose = verbose)

    # Filter clusters to requested chromosomes
    if (!is.null(chromosomes)) {
        cluster_chr <- as.character(GenomicRanges::seqnames(exon_clusters))
        exon_clusters <- exon_clusters[cluster_chr %in% chromosomes]
        if (verbose) message("  Filtered to ", length(exon_clusters), " clusters on requested chromosomes")
    }

    # Sort by (chr, start) for downstream C binary search
    exon_df <- exon_df[order(exon_df$chr, exon_df$start), ]
    rownames(exon_df) <- NULL
    intron_df <- intron_df[order(intron_df$chr, intron_df$start), ]
    rownames(intron_df) <- NULL

    result <- list(
        exons = exon_df,
        introns = intron_df,
        exon_clusters = exon_clusters,
        gene_ids = gene_ids,
        gene_to_idx = setNames(seq_along(gene_ids), gene_ids)
    )

    return(result)
}


#' Classify Reads and Build Per-Category Cell x Exon Matrices
#'
#' @description
#' Classifies aligned reads into spliced, unspliced, chimeric, ambiguous, and
#' unassigned categories based on the HUGE bit field scheme:
#' H: has N in CIGAR; U: unique gene assignment; G: all blocks in gene; E: all blocks exonic.
#'   1111 -> SPLICED, 1011 -> CHIMERIC, _110 -> UNSPLICED, 0111 -> AMBIGUOUS, other -> UNASSIGNED.
#' Outputs CellRanger-style sparse matrices for each non-unassigned category.
#'
#' @param bam Character. Path to input BAM/CRAM.
#' @param annotation List from \code{prepare_classification_annotation}.
#' @param out_dir Character. Output directory for classified matrices.
#' @param xf_values Integer vector of allowed xf values (10X only). Default: c(25L, 17L).
#' @param tech Character. "10X", "Smart-seq3Xpress", or "Smart-seq". Default: "10X".
#' @param gene_span_pad Integer. Padding (bp) added to read span when checking gene overlap. Default: 0L.
#' @return Named list with classification counts, HUGE flag totals, and output directory.
#' @export
classify_reads <- function(
    bam,
    annotation,
    out_dir,
    xf_values = c(25L, 17L),
    tech = "10X",
    gene_span_pad = 0L
) {
    stopifnot(file.exists(bam))
    tech <- match.arg(tech, c("10X", "Smart-seq3Xpress", "Smart-seq"))
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

    exon_df <- annotation$exons
    intron_df <- annotation$introns
    exon_clusters <- annotation$exon_clusters

    # Detect BAM chromosome naming: "chr19" vs "19"
    bam_info <- Rsamtools::scanBamHeader(bam)
    bam_seqs <- names(bam_info$targets)
    bam_has_chr <- length(bam_seqs) > 0 && any(grepl("^chr", bam_seqs))

    # Normalize annotation chromosome names to match BAM
    normalize_chr <- function(chr_vec) {
        if (bam_has_chr) {
            ifelse(grepl("^chr", chr_vec), chr_vec, paste0("chr", chr_vec))
        } else {
            sub("^chr", "", chr_vec)
        }
    }

    exon_chr <- normalize_chr(exon_df$chr)
    intron_chr <- if (nrow(intron_df) > 0) normalize_chr(intron_df$chr) else character(0)
    cluster_chr <- normalize_chr(as.character(GenomicRanges::seqnames(exon_clusters)))
    cluster_start <- GenomicRanges::start(exon_clusters)
    cluster_end <- GenomicRanges::end(exon_clusters)
    cluster_strand <- as.character(BiocGenerics::strand(exon_clusters))
    cluster_ids <- seq_along(cluster_start) - 1L

    message("Calling C classify_reads (single BAM pass)...")
    t0 <- Sys.time()

    result <- .Call(
        `_exonBlocks_classify_reads`,
        bam,
        exon_chr,
        as.integer(exon_df$start),
        as.integer(exon_df$end),
        as.integer(exon_df$gene_idx),
        exon_df$strand,
        intron_chr,
        as.integer(intron_df$start),
        as.integer(intron_df$end),
        as.integer(intron_df$gene_idx),
        intron_df$strand,
        cluster_chr,
        as.integer(cluster_start),
        as.integer(cluster_end),
        cluster_strand,
        cluster_ids,
        tech,
        as.integer(xf_values),
        out_dir,
        as.integer(gene_span_pad)
    )

    t1 <- Sys.time()

    total <- result$n_spliced + result$n_unspliced + result$n_chimeric + result$n_ambiguous + result$n_unassigned
    message(sprintf(
        "Done in %.1f sec | Cells: %d | Spliced: %d | Unspliced: %d | Chimeric: %d | Ambiguous: %d | Unassigned: %d",
        as.numeric(difftime(t1, t0, units = "secs")),
        result$n_cells, result$n_spliced, result$n_unspliced,
        result$n_chimeric, result$n_ambiguous, result$n_unassigned
    ))

    if (total > 0) {
        message(sprintf("  Spliced:    %5.1f%%", 100 * result$n_spliced / total))
        message(sprintf("  Unspliced:  %5.1f%%", 100 * result$n_unspliced / total))
        message(sprintf("  Chimeric:   %5.1f%%", 100 * result$n_chimeric / total))
        message(sprintf("  Ambiguous:  %5.1f%%", 100 * result$n_ambiguous / total))
        message(sprintf("  Unassigned: %5.1f%%", 100 * result$n_unassigned / total))
        message(sprintf("  HUGE flags: H=%d U=%d G=%d E=%d",
                        result$flag_H, result$flag_U, result$flag_G, result$flag_E))
    }

    invisible(result)
}