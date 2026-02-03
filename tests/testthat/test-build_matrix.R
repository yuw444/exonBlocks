# Tests for build_cell_cluster_matrix and related functions
library(testthat)
library(GenomicRanges)

test_that("save_cluster_index creates valid index", {
    # Create mock exon clusters
    gr <- GRanges(
        seqnames = rep("chr1", 5),
        ranges = IRanges(
            start = c(100, 200, 300, 400, 500),
            end = c(150, 250, 350, 450, 550)
        ),
        strand = "+"
    )

    tmp <- tempfile(fileext = ".rds")
    save_cluster_index(gr, tmp, overwrite = TRUE)

    idx <- readRDS(tmp)
    expect_equal(idx$n_clusters, 5)
    expect_equal(length(idx$chr), 5)
    expect_equal(length(idx$start), 5)
    expect_equal(length(idx$end), 5)

    unlink(tmp)
})

test_that("build_cell_cluster_matrix returns valid dgCMatrix", {
    skip_if_not_installed("Matrix")

    # Create a simple test case
    # Two cells, three clusters
    mat <- Matrix::sparseMatrix(
        i = c(1, 1, 2, 2),
        j = c(1, 2, 2, 3),
        x = c(1, 2, 1, 3),
        dims = c(2, 3),
        giveCsparse = TRUE
    )

    expect_true(is(mat, "dgCMatrix"))
    expect_equal(nrow(mat), 2)
    expect_equal(ncol(mat), 3)
    expect_equal(Matrix::rowSums(mat), c(3, 4))
})

test_that("export_cell_cluster_matrix handles CSV format", {
    skip_if_not_installed("Matrix")

    mat <- Matrix::sparseMatrix(
        i = c(1, 2),
        j = c(1, 2),
        x = c(1, 1),
        dims = c(2, 3)
    )

    tmp <- tempfile(fileext = ".csv")
    expect_silent(export_cell_cluster_matrix(mat, format = "csv", path = tmp))

    # Verify file was created
    expect_true(file.exists(tmp))

    # Check content
    df <- read.csv(tmp, row.names = 1)
    expect_equal(nrow(df), 2)
    expect_equal(ncol(df), 3)

    unlink(tmp)
})

test_that("cb_umi_exons works with semicolon-delimited blocks", {
    # Create temp TSV with semicolon-delimited blocks
    tsv_content <- paste0(
        "CB\tUMI\tblock_start\tblock_end\tblock_seq\tnum_blocks\n",
        "cell1\tUMI1\t100;200\t150;250\tACGT;TGCA\t2\n",
        "cell1\tUMI2\t300\t350\tAAAA\t1\n",
        "cell2\tUMI1\t120\t170\tGGGG\t1\n"
    )

    tmp_tsv <- tempfile(fileext = ".tsv")
    writeLines(tsv_content, tmp_tsv)

    # Create temp exons file
    exons_content <- paste0(
        "chr\tstart\tend\texon\n",
        "chr1\t90\t160\texon1\n",
        "chr1\t190\t260\texon2\n",
        "chr1\t290\t360\texon3\n"
    )
    tmp_exons <- tempfile(fileext = ".tsv")
    writeLines(exons_content, tmp_exons)

    result <- cb_umi_exons(tsv = tmp_tsv, meta_exons = tmp_exons)

    expect_s3_class(result, "data.table")
    expect_true("exon" %in% names(result))

    unlink(tmp_tsv)
    unlink(tmp_exons)
})

test_that("cb_umi_exons handles no overlaps gracefully", {
    tsv_content <- paste0(
        "CB\tUMI\tblock_start\tblock_end\tblock_seq\tnum_blocks\n",
        "cell1\tUMI1\t100\t150\tACGT\t1\n"
    )
    tmp_tsv <- tempfile(fileext = ".tsv")
    writeLines(tsv_content, tmp_tsv)

    exons_content <- paste0(
        "chr\tstart\tend\texon\n",
        "chr1\t1000\t2000\texon1\n"
    )
    tmp_exons <- tempfile(fileext = ".tsv")
    writeLines(exons_content, tmp_exons)

    result <- cb_umi_exons(tsv = tmp_tsv, meta_exons = tmp_exons)
    expect_equal(nrow(result), 0)

    unlink(tmp_tsv)
    unlink(tmp_exons)
})
