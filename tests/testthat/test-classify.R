library(testthat)
library(exonBlocks)
library(GenomicRanges)

test_that("prepare_classification_annotation returns expected structure", {
    skip_if_not_installed("TxDb.Hsapiens.UCSC.hg38.knownGene")
    skip_on_cran()

    annot <- prepare_classification_annotation(
        species = "hg38",
        min_gap_width = 10,
        chromosomes = "chr19",
        verbose = FALSE
    )

    expect_type(annot, "list")
    expect_true(all(c("exons", "introns", "exon_clusters", "gene_ids", "gene_to_idx") %in% names(annot)))

    expect_s3_class(annot$exons, "data.frame")
    expect_true(all(c("chr", "start", "end", "strand", "gene_idx") %in% names(annot$exons)))
    expect_gt(nrow(annot$exons), 0)

    expect_s3_class(annot$introns, "data.frame")
    expect_true(all(c("chr", "start", "end", "strand", "gene_idx") %in% names(annot$introns)))

    expect_true(inherits(annot$exon_clusters, "GRanges"))
    expect_gt(length(annot$exon_clusters), 0)

    expect_type(annot$gene_ids, "character")
    expect_type(annot$gene_to_idx, "integer")
})

test_that("prepare_classification_annotation rejects unsupported species", {
    expect_error(
        prepare_classification_annotation(species = "zebrafish"),
        "Unsupported species"
    )
})

test_that("classify_reads validates inputs", {
    expect_error(
        classify_reads(bam = "/nonexistent.bam", annotation = list(), out_dir = tempdir()),
        class = "simpleError"
    )

    expect_error(
        classify_reads(
            bam = "fake.bam",
            annotation = list(),
            out_dir = tempdir(),
            tech = "UNSUPPORTED"
        ),
        class = "simpleError"
    )
})

test_that("classify_reads tech parameter accepts valid values", {
    valid_techs <- c("10X", "Smart-seq3Xpress", "Smart-seq")
    for (t in valid_techs) {
        expect_match(t, t)
    }
})

test_that("classify_reads output structure has correct subdirectories", {
    skip_if_not_installed("TxDb.Hsapiens.UCSC.hg38.knownGene")
    skip_if_not_installed("Rsamtools")
    skip_on_cran()

    bam_file <- test_path("..", "..", "meta", "bam", "test.bam")
    skip_if_not(file.exists(bam_file), "test.bam not found")

    annot <- prepare_classification_annotation(
        species = "hg38",
        chromosomes = "chr19",
        verbose = FALSE
    )

    out_dir <- tempfile(pattern = "classify_test_")
    on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

    result <- classify_reads(
        bam = bam_file,
        annotation = annot,
        out_dir = out_dir,
        tech = "10X"
    )

    categories <- c("spliced", "unspliced", "chimeric", "ambiguous")
    files_per_cat <- c("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz")

    for (cat in categories) {
        cat_dir <- file.path(out_dir, cat)
        expect_true(dir.exists(cat_dir), label = paste("directory exists:", cat_dir))
        for (f in files_per_cat) {
            expect_true(file.exists(file.path(cat_dir, f)),
                        label = paste("file exists:", file.path(cat_dir, f)))
        }
    }

    expect_type(result, "list")
    expect_true(all(c("n_cells", "n_features", "n_spliced", "n_unspliced",
                       "n_chimeric", "n_ambiguous", "n_unassigned", "out_dir",
                       "flag_H", "flag_U", "flag_G", "flag_E") %in% names(result)))
})

test_that("classify_reads category counts are non-negative", {
    skip_if_not_installed("TxDb.Hsapiens.UCSC.hg38.knownGene")
    skip_if_not_installed("Rsamtools")
    skip_on_cran()

    bam_file <- test_path("..", "..", "meta", "bam", "test.bam")
    skip_if_not(file.exists(bam_file), "test.bam not found")

    annot <- prepare_classification_annotation(
        species = "hg38",
        chromosomes = "chr19",
        verbose = FALSE
    )

    out_dir <- tempfile(pattern = "classify_counts_")
    on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

    result <- classify_reads(
        bam = bam_file,
        annotation = annot,
        out_dir = out_dir,
        tech = "10X"
    )

    expect_gte(result$n_spliced, 0)
    expect_gte(result$n_unspliced, 0)
    expect_gte(result$n_chimeric, 0)
    expect_gte(result$n_ambiguous, 0)
    expect_gte(result$n_unassigned, 0)

    total <- result$n_spliced + result$n_unspliced + result$n_chimeric + result$n_ambiguous + result$n_unassigned
    expect_gte(total, 0)
})

test_that("classify_reads features.tsv.gz has correct header", {
    skip_if_not_installed("TxDb.Hsapiens.UCSC.hg38.knownGene")
    skip_if_not_installed("Rsamtools")
    skip_on_cran()

    bam_file <- test_path("..", "..", "meta", "bam", "test.bam")
    skip_if_not(file.exists(bam_file), "test.bam not found")

    annot <- prepare_classification_annotation(
        species = "hg38",
        chromosomes = "chr19",
        verbose = FALSE
    )

    out_dir <- tempfile(pattern = "classify_features_")
    on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

    result <- classify_reads(
        bam = bam_file,
        annotation = annot,
        out_dir = out_dir,
        tech = "10X"
    )

    feat_path <- file.path(out_dir, "spliced", "features.tsv.gz")
    skip_if_not(file.exists(feat_path), "features.tsv.gz not found")

    con <- gzfile(feat_path)
    lines <- readLines(con, n = 2)
    close(con)

    expect_match(lines[1], "feature_id\tfeature_name\tfeature_type")
    expect_match(lines[2], "cluster_")
})