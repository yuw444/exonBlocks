#' Prepare Exon and Intron Annotation for Read Classification
#'
#' @description
#' Extracts exon and intron coordinates with gene IDs from a TxDb package
#' and prepares exon clusters. Introns are computed as the gaps within
#' each gene's exon span (union strategy).
#'
#' @param species Character. Genome assembly: "hg38", "hg19", "mm10", "mm9", "rn6", "dm6", "ce11".
#' @param min_gap_width Integer. Exons within this distance are merged into same cluster. Default: 10.
#' @param verbose Logical. Print progress. Default: TRUE.
#' @return A list with components: exons (GRanges with gene_id), introns (GRanges with gene_id),
#'   exon_clusters (GRanges from \code{make_exon_clusters}).
#' @importFrom GenomicRanges GRanges reduce seqnames start end split gaps
#' @importFrom IRanges IRanges
#' @importFrom GenomicFeatures exons
#' @export
prepare_classification_annotation <- function(
    species = "hg38",
    min_gap_width = 10,
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
    suppressPackageStartupMessages({ txdb <- get(txdb_pkg) })

    if (verbose) message("Extracting exons with gene_id...")
    gr_exons <- GenomicFeatures::exons(txdb, columns = "gene_id")

    # Unlist gene_id (list column) to first gene per exon for simplicity
    # Some exons map to multiple genes (overlapping loci) — keep all
    exon_gene_ids <- lapply(gr_exons$gene_id, function(x) if (is.na(x[1])) NA_integer_ else x)
    # Assign numeric gene IDs
    all_genes <- sort(unique(unlist(exon_gene_ids)))
    all_genes <- all_genes[!is.na(all_genes)]
    gene_to_idx <- setNames(seq_along(all_genes), all_genes)

    if (verbose) message("  ", length(all_genes), " unique genes, ", length(gr_exons), " exons")

    # Expand multi-gene exons to one row per (exon, gene) pair
    exon_rows <- lapply(seq_along(gr_exons), function(i) {
        gids <- gr_exons$gene_id[[i]]
        if (is.null(gids) || all(is.na(gids))) return(NULL)
        data.frame(
            chr = as.character(GenomicRanges::seqnames(gr_exons)[i]),
            start = GenomicRanges::start(gr_exons)[i],
            end = GenomicRanges::end(gr_exons)[i],
            gene_id = gids[!is.na(gids)],
            stringsAsFactors = FALSE
        )
    })
    exon_df <- do.call(rbind, exon_rows)
    exon_df$gene_idx <- as.integer(gene_to_idx[as.character(exon_df$gene_id)])

    # Compute introns: for each gene, take exon union span, then gaps within
    if (verbose) message("Computing introns per gene...")
    exon_gr <- GenomicRanges::GRanges(
        seqnames = exon_df$chr,
        IRanges::IRanges(start = exon_df$start, end = exon_df$end)
    )
    GenomicRanges::mcols(exon_gr)$gene_idx <- exon_df$gene_idx

    # Split by gene, reduce exons per gene, compute gaps
    exon_by_gene <- GenomicRanges::split(exon_gr, exon_gr$gene_idx)
    gene_exon_reduced <- GenomicRanges::reduce(exon_by_gene, min.gapwidth = 1L)
    gene_span <- GenomicRanges::reduce(exon_by_gene, min.gapwidth = 0L)
    intron_gr <- GenomicRanges::gaps(gene_exon_reduced)

    # Filter introns to only those within gene spans
    # gaps() returns all gaps including between chromosomes; subset to gene-span regions
    intron_with_gene <- S4Vectors::mcols(intron_gr)
    # Re-attach gene_idx from names
    intron_gene_idx <- as.integer(names(intron_gr))
    if (is.null(intron_gene_idx)) {
        # gaps() drops mcols; reconstruct from gene_exon_reduced spans
        # Use a different approach: compute introns directly
        intron_rows <- lapply(seq_along(gene_exon_reduced), function(i) {
            gex <- gene_exon_reduced[[i]]
            gidx <- as.integer(names(gene_exon_reduced)[i])
            if (length(gex) <= 1) return(NULL)
            # Find gaps between consecutive reduced exons within gene
            gene_sp <- gene_span[[i]]
            ir <- GenomicRanges::ranges(gex)
            gsp_ir <- GenomicRanges::ranges(gene_sp)
            gap_ir <- IRanges::gaps(ir, start = min(IRanges::start(gsp_ir)), end = max(IRanges::end(gsp_ir)))
            if (length(gap_ir) == 0) return(NULL)
            data.frame(
                chr = as.character(GenomicRanges::seqnames(gene_sp)[1]),
                start = IRanges::start(gap_ir),
                end = IRanges::end(gap_ir),
                gene_idx = gidx,
                stringsAsFactors = FALSE
            )
        })
        intron_df <- do.call(rbind, intron_rows)
    } else {
        intron_df <- data.frame(
            chr = as.character(GenomicRanges::seqnames(intron_gr)),
            start = GenomicRanges::start(intron_gr),
            end = GenomicRanges::end(intron_gr),
            gene_idx = intron_gene_idx,
            stringsAsFactors = FALSE
        )
    }

    if (is.null(intron_df) || nrow(intron_df) == 0) {
        if (verbose) message("  No introns computed (single-exon genes only?)")
        intron_df <- data.frame(chr = character(0), start = integer(0),
                                end = integer(0), gene_idx = integer(0),
                                stringsAsFactors = FALSE)
    } else {
        if (verbose) message("  ", nrow(intron_df), " intron intervals")
    }

    # Build exon clusters
    if (verbose) message("Building exon clusters...")
    exon_clusters <- make_exon_clusters(species = species, min_gap_width = min_gap_width, verbose = verbose)

    # Sort exon and intron dfs by (chr, start) for C binary search
    exon_df <- exon_df[order(exon_df$chr, exon_df$start), ]
    intron_df <- intron_df[order(intron_df$chr, intron_df$start), ]

    result <- list(
        exons = exon_df,
        introns = intron_df,
        exon_clusters = exon_clusters,
        gene_ids = all_genes,
        gene_to_idx = gene_to_idx
    )

    return(result)
}


#' Classify Reads and Build Per-Category Cell x Exon Matrices
#'
#' @description
#' Classifies aligned reads into spliced, unspliced, ambiguous, and unassigned
#' categories based on CIGAR and annotation overlap. Outputs CellRanger-style
#' sparse matrices for each non-unassigned category.
#'
#' @param bam Character. Path to input BAM/CRAM with CB/UB tags.
#' @param annotation List from \code{prepare_classification_annotation}.
#' @param out_dir Character. Output directory for classified matrices.
#' @param xf_values Integer vector of allowed xf values (10X only). Default: c(25L, 17L).
#' @param tech Character. "10X" or "Smart-seq". Default: "10X".
#' @return Named list with classification counts and output directory.
#' @export
classify_reads <- function(
    bam,
    annotation,
    out_dir,
    xf_values = c(25L, 17L),
    tech = "10X"
) {
    stopifnot(file.exists(bam))
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

    exon_df <- annotation$exons
    intron_df <- annotation$introns
    exon_clusters <- annotation$exon_clusters

    cluster_chr <- as.character(GenomicRanges::seqnames(exon_clusters))
    cluster_start <- GenomicRanges::start(exon_clusters)
    cluster_end <- GenomicRanges::end(exon_clusters)

    unique_chr <- unique(cluster_chr)

    total_result <- list(
        n_cells = 0L,
        n_features = length(exon_clusters),
        n_spliced = 0L,
        n_unspliced = 0L,
        n_ambiguous = 0L,
        n_unassigned = 0L
    )

    for (chrom in unique_chr) {
        message("Classifying reads on: ", chrom)

        # Subset annotation for this chromosome
        exon_mask <- exon_df$chr == chrom
        intron_mask <- intron_df$chr == chrom
        cluster_mask <- cluster_chr == chrom

        if (sum(exon_mask) == 0) next

        chr_exon_starts <- as.integer(exon_df$start[exon_mask])
        chr_exon_ends <- as.integer(exon_df$end[exon_mask])
        chr_exon_gene_ids <- as.integer(exon_df$gene_idx[exon_mask])

        if (sum(intron_mask) > 0) {
            chr_intron_starts <- as.integer(intron_df$start[intron_mask])
            chr_intron_ends <- as.integer(intron_df$end[intron_mask])
            chr_intron_gene_ids <- as.integer(intron_df$gene_idx[intron_mask])
        } else {
            chr_intron_starts <- integer(0)
            chr_intron_ends <- integer(0)
            chr_intron_gene_ids <- integer(0)
        }

        chr_cluster_starts <- as.integer(cluster_start[cluster_mask])
        chr_cluster_ends <- as.integer(cluster_end[cluster_mask])

        result <- .Call(
            `_exonBlocks_classify_reads`,
            bam,
            chrom,
            chr_exon_starts,
            chr_exon_ends,
            chr_exon_gene_ids,
            chr_intron_starts,
            chr_intron_ends,
            chr_intron_gene_ids,
            chr_cluster_starts,
            chr_cluster_ends,
            tech,
            as.integer(xf_values),
            out_dir
        )

        message("  Cells: ", result$n_cells,
                " | Spliced: ", result$n_spliced,
                " | Unspliced: ", result$n_unspliced,
                " | Ambiguous: ", result$n_ambiguous,
                " | Unassigned: ", result$n_unassigned)

        total_result$n_cells <- total_result$n_cells + result$n_cells
        total_result$n_spliced <- total_result$n_spliced + result$n_spliced
        total_result$n_unspliced <- total_result$n_unspliced + result$n_unspliced
        total_result$n_ambiguous <- total_result$n_ambiguous + result$n_ambiguous
        total_result$n_unassigned <- total_result$n_unassigned + result$n_unassigned
    }

    total_result$out_dir <- out_dir
    message("Total: ", total_result$n_cells, " cells, ",
            total_result$n_spliced, " spliced, ",
            total_result$n_unspliced, " unspliced, ",
            total_result$n_ambiguous, " ambiguous, ",
            total_result$n_unassigned, " unassigned")

    invisible(total_result)
}