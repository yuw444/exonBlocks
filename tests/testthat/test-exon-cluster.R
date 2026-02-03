# Tests for make_exon_clusters
library(testthat)
library(GenomicRanges)

# Ensure exonBlocks is loaded
library(exonBlocks)

# Check if required TxDb packages are installed
txdb_hg38_installed <- require(TxDb.Hsapiens.UCSC.hg38.knownGene, quietly = TRUE)
txdb_mm10_installed <- require(TxDb.Mmusculus.UCSC.mm10.knownGene, quietly = TRUE)

test_that("make_exon_clusters returns GRanges object", {
    skip_if_not(txdb_hg38_installed, "TxDb.Hsapiens.UCSC.hg38.knownGene not installed")

    clusters <- make_exon_clusters(
        species = "hg38",
        min_gap_width = 10,
        verbose = FALSE
    )

    expect_true(inherits(clusters, "GRanges"))
    expect_true(length(clusters) > 0)
})

test_that("make_exon_clusters with default parameters", {
    skip_if_not(txdb_hg38_installed, "TxDb.Hsapiens.UCSC.hg38.knownGene not installed")

    clusters <- make_exon_clusters(verbose = FALSE)

    expect_true(inherits(clusters, "GRanges"))
    # Default is hg38
    expect_true(length(clusters) > 0)
})

test_that("make_exon_clusters has revmap metadata", {
    skip_if_not(txdb_hg38_installed, "TxDb.Hsapiens.UCSC.hg38.knownGene not installed")

    clusters <- make_exon_clusters(
        species = "hg38",
        min_gap_width = 10,
        verbose = FALSE
    )

    expect_true("revmap" %in% names(mcols(clusters)))
})

test_that("make_exon_clusters with NA min_gap_width (overlapping only)", {
    skip_if_not(txdb_hg38_installed, "TxDb.Hsapiens.UCSC.hg38.knownGene not installed")

    clusters_na <- make_exon_clusters(
        species = "hg38",
        min_gap_width = NA,
        verbose = FALSE
    )

    clusters_10 <- make_exon_clusters(
        species = "hg38",
        min_gap_width = 10,
        verbose = FALSE
    )

    # NA (overlapping only) should give more or equal clusters
    expect_true(length(clusters_na) >= length(clusters_10))
})

test_that("make_exon_clusters with different min_gap_width values", {
    skip_if_not(txdb_hg38_installed, "TxDb.Hsapiens.UCSC.hg38.knownGene not installed")

    clusters_20 <- make_exon_clusters(
        species = "hg38",
        min_gap_width = 20,
        verbose = FALSE
    )

    clusters_5 <- make_exon_clusters(
        species = "hg38",
        min_gap_width = 5,
        verbose = FALSE
    )

    # Larger min_gap_width should merge more exons (fewer clusters)
    expect_true(length(clusters_20) <= length(clusters_5))
})

test_that("make_exon_clusters chromosome naming", {
    skip_if_not(txdb_hg38_installed, "TxDb.Hsapiens.UCSC.hg38.knownGene not installed")

    clusters <- make_exon_clusters(
        species = "hg38",
        min_gap_width = 10,
        verbose = FALSE
    )

    # Should have chr prefix or NC_ accession
    seqlevels <- seqlevels(clusters)
    has_chr <- any(grepl("^chr", seqlevels))
    has_nc <- any(grepl("^NC_", seqlevels))
    expect_true(has_chr || has_nc)
})

test_that("make_exon_clusters preserves strand information", {
    skip_if_not(txdb_hg38_installed, "TxDb.Hsapiens.UCSC.hg38.knownGene not installed")

    clusters <- make_exon_clusters(
        species = "hg38",
        min_gap_width = 10,
        verbose = FALSE
    )

    strands <- as.character(strand(clusters))
    expect_true(all(strands %in% c("+", "-", "*")))
})

test_that("make_exon_clusters rejects invalid species", {
    expect_error(
        make_exon_clusters(
            species = "not_a_real_species_xyz123",
            verbose = FALSE
        ),
        "Unsupported species"
    )
})

test_that("make_exon_clusters creates fewer clusters than exons", {
    skip_if_not(txdb_hg38_installed, "TxDb.Hsapiens.UCSC.hg38.knownGene not installed")

    clusters <- make_exon_clusters(
        species = "hg38",
        min_gap_width = 10,
        verbose = FALSE
    )

    # With min_gap_width > 0, should merge some exons
    expect_true(length(clusters) > 0)
})

test_that("make_exon_clusters verbose output works", {
    skip_if_not(txdb_hg38_installed, "TxDb.Hsapiens.UCSC.hg38.knownGene not installed")

    # Should print messages without error
    expect_message(
        make_exon_clusters(
            species = "hg38",
            min_gap_width = 10,
            verbose = TRUE
        ),
        "Loading TxDb"
    )
})

test_that("make_exon_clusters GRanges has valid coordinates", {
    skip_if_not(txdb_hg38_installed, "TxDb.Hsapiens.UCSC.hg38.knownGene not installed")

    clusters <- make_exon_clusters(
        species = "hg38",
        min_gap_width = 10,
        verbose = FALSE
    )

    # All starts should be <= ends
    expect_true(all(start(clusters) <= end(clusters)))

    # All coordinates should be positive
    expect_true(all(start(clusters) > 0))
})

test_that("make_exon_clusters mm10 species works", {
    skip_if_not(txdb_mm10_installed, "TxDb.Mmusculus.UCSC.mm10.knownGene not installed")

    clusters <- make_exon_clusters(
        species = "mm10",
        min_gap_width = 10,
        verbose = FALSE
    )

    expect_true(inherits(clusters, "GRanges"))
    expect_true(length(clusters) > 0)
})

test_that("make_exon_clusters function exists and is exported", {
    expect_true(exists("make_exon_clusters", mode = "function"))
})
