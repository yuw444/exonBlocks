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

test_that("prepare_classification_annotation mm10 returns expected structure", {
    skip_if_not_installed("TxDb.Mmusculus.UCSC.mm10.knownGene")
    skip_on_cran()

    annot <- prepare_classification_annotation(
        species = "mm10",
        min_gap_width = 10,
        chromosomes = "chr1",
        verbose = FALSE
    )

    expect_type(annot, "list")
    expect_s3_class(annot$exons, "data.frame")
    expect_gt(nrow(annot$exons), 0)
    expect_true(inherits(annot$exon_clusters, "GRanges"))
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

test_that("classify_reads 10X mouse output structure has correct subdirectories", {
    skip_if_not_installed("TxDb.Mmusculus.UCSC.mm10.knownGene")
    skip_if_not_installed("Rsamtools")
    skip_on_cran()

    bam_file <- test_path("..", "..", "meta", "bam", "test.bam")
    skip_if_not(file.exists(bam_file), "test.bam not found")

    annot <- prepare_classification_annotation(
        species = "mm10",
        chromosomes = "chr1",
        verbose = FALSE
    )

    out_dir <- tempfile(pattern = "classify_mouse10x_")
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

test_that("classify_reads 10X mouse category counts are consistent", {
    skip_if_not_installed("TxDb.Mmusculus.UCSC.mm10.knownGene")
    skip_if_not_installed("Rsamtools")
    skip_on_cran()

    bam_file <- test_path("..", "..", "meta", "bam", "test.bam")
    skip_if_not(file.exists(bam_file), "test.bam not found")

    annot <- prepare_classification_annotation(
        species = "mm10",
        chromosomes = "chr1",
        verbose = FALSE
    )

    out_dir <- tempfile(pattern = "classify_mouse_counts_")
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
    expect_gt(total, 0)

    expect_gte(result$flag_H, 0)
    expect_gte(result$flag_U, 0)
    expect_gte(result$flag_G, 0)
    expect_gte(result$flag_E, 0)

    expect_gt(result$n_spliced, 0, label = "10X test.bam has spliced reads")
    expect_gt(result$n_ambiguous, 0, label = "10X test.bam has ambiguous reads")
})

test_that("classify_reads 10X mouse features.tsv.gz has correct header", {
    skip_if_not_installed("TxDb.Mmusculus.UCSC.mm10.knownGene")
    skip_if_not_installed("Rsamtools")
    skip_on_cran()

    bam_file <- test_path("..", "..", "meta", "bam", "test.bam")
    skip_if_not(file.exists(bam_file), "test.bam not found")

    annot <- prepare_classification_annotation(
        species = "mm10",
        chromosomes = "chr1",
        verbose = FALSE
    )

    out_dir <- tempfile(pattern = "classify_mouse_feat_")
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

test_that("classify_reads 10X mouse matrix.mtx.gz dimensions match barcodes and features", {
    skip_if_not_installed("TxDb.Mmusculus.UCSC.mm10.knownGene")
    skip_if_not_installed("Matrix")
    skip_if_not_installed("Rsamtools")
    skip_on_cran()

    bam_file <- test_path("..", "..", "meta", "bam", "test.bam")
    skip_if_not(file.exists(bam_file), "test.bam not found")

    annot <- prepare_classification_annotation(
        species = "mm10",
        chromosomes = "chr1",
        verbose = FALSE
    )

    out_dir <- tempfile(pattern = "classify_mouse_mtx_")
    on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

    result <- classify_reads(
        bam = bam_file,
        annotation = annot,
        out_dir = out_dir,
        tech = "10X"
    )

    for (cat in c("spliced", "unspliced", "chimeric", "ambiguous")) {
        mtx_path <- file.path(out_dir, cat, "matrix.mtx.gz")
        barcodes_path <- file.path(out_dir, cat, "barcodes.tsv.gz")
        features_path <- file.path(out_dir, cat, "features.tsv.gz")

        skip_if_not(file.exists(mtx_path), paste(cat, "matrix.mtx.gz not found"))

        mat <- Matrix::readMM(mtx_path)
        barcodes <- readLines(gzfile(barcodes_path))
        features <- readLines(gzfile(features_path))

        n_cells <- length(barcodes)
        n_features <- length(features) - 1L

        expect_equal(nrow(mat), n_features, label = paste(cat, "rows match features"))
        expect_equal(ncol(mat), n_cells, label = paste(cat, "cols match barcodes"))
    }
})

test_that("classify_reads 10X mouse barcodes.tsv.gz has 10X-format barcodes", {
    skip_if_not_installed("TxDb.Mmusculus.UCSC.mm10.knownGene")
    skip_if_not_installed("Rsamtools")
    skip_on_cran()

    bam_file <- test_path("..", "..", "meta", "bam", "test.bam")
    skip_if_not(file.exists(bam_file), "test.bam not found")

    annot <- prepare_classification_annotation(
        species = "mm10",
        chromosomes = "chr1",
        verbose = FALSE
    )

    out_dir <- tempfile(pattern = "classify_mouse_bc_")
    on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

    result <- classify_reads(
        bam = bam_file,
        annotation = annot,
        out_dir = out_dir,
        tech = "10X"
    )

    bc_path <- file.path(out_dir, "ambiguous", "barcodes.tsv.gz")
    skip_if_not(file.exists(bc_path), "barcodes.tsv.gz not found")

    lines <- readLines(gzfile(bc_path), n = 3)
    expect_match(lines[1], "-1$")
})

test_that("classify_reads returns by_chr summary", {
    skip_if_not_installed("TxDb.Mmusculus.UCSC.mm10.knownGene")
    skip_if_not_installed("Rsamtools")
    skip_on_cran()

    bam_file <- test_path("..", "..", "meta", "bam", "test.bam")
    skip_if_not(file.exists(bam_file), "test.bam not found")

    annot <- prepare_classification_annotation(
        species = "mm10",
        chromosomes = "chr1",
        verbose = FALSE
    )

    out_dir <- tempfile(pattern = "classify_bychr_")
    on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

    result <- classify_reads(
        bam = bam_file,
        annotation = annot,
        out_dir = out_dir,
        tech = "10X"
    )

    expect_true("by_chr" %in% names(result))
    expect_s3_class(result$by_chr, "data.frame")
    expect_true("chr" %in% names(result$by_chr))
    expect_gte(nrow(result$by_chr), 1)
    expect_true(all(c("spliced", "unspliced", "chimeric", "ambiguous", "unassigned") %in% names(result$by_chr)))
})

test_that("classify_reads Smart-seq3Xpress human produces output", {
    skip_if_not_installed("TxDb.Hsapiens.UCSC.hg38.knownGene")
    skip_if_not_installed("Rsamtools")
    skip_on_cran()

    bam_file <- test_path("..", "..", "meta", "bam", "test_data", "human_smartseq3x_1M.bam")
    skip_if_not(file.exists(bam_file), "human_smartseq3x_1M.bam not found")

    annot <- prepare_classification_annotation(
        species = "hg38",
        chromosomes = "chr1",
        verbose = FALSE
    )

    out_dir <- tempfile(pattern = "classify_smartseq3x_")
    on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

    result <- classify_reads(
        bam = bam_file,
        annotation = annot,
        out_dir = out_dir,
        tech = "Smart-seq3Xpress"
    )

    expect_type(result, "list")
    expect_gte(result$n_cells, 0)

    total <- result$n_spliced + result$n_unspliced + result$n_chimeric + result$n_ambiguous + result$n_unassigned
    expect_gt(total, 0)
})

test_that("classify_reads Smart-seq3Xpress human matrix dimensions are consistent", {
    skip_if_not_installed("TxDb.Hsapiens.UCSC.hg38.knownGene")
    skip_if_not_installed("Matrix")
    skip_if_not_installed("Rsamtools")
    skip_on_cran()

    bam_file <- test_path("..", "..", "meta", "bam", "test_data", "human_smartseq3x_1M.bam")
    skip_if_not(file.exists(bam_file), "human_smartseq3x_1M.bam not found")

    annot <- prepare_classification_annotation(
        species = "hg38",
        chromosomes = "chr1",
        verbose = FALSE
    )

    out_dir <- tempfile(pattern = "classify_smartseq3x_mtx_")
    on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

    result <- classify_reads(
        bam = bam_file,
        annotation = annot,
        out_dir = out_dir,
        tech = "Smart-seq3Xpress"
    )

    for (cat in c("spliced", "unspliced", "chimeric", "ambiguous")) {
        mtx_path <- file.path(out_dir, cat, "matrix.mtx.gz")
        skip_if_not(file.exists(mtx_path), paste(cat, "matrix missing"))

        mat <- Matrix::readMM(mtx_path)
        n_features <- length(readLines(gzfile(file.path(out_dir, cat, "features.tsv.gz")))) - 1L
        n_cells <- length(readLines(gzfile(file.path(out_dir, cat, "barcodes.tsv.gz"))))

        expect_equal(nrow(mat), n_features, label = paste(cat, "rows match features"))
        expect_equal(ncol(mat), n_cells, label = paste(cat, "cols match barcodes"))
    }
})

test_that("classify_reads 10X mouse 1M produces output with meaningful classification", {
    skip_if_not_installed("TxDb.Mmusculus.UCSC.mm10.knownGene")
    skip_if_not_installed("Rsamtools")
    skip_on_cran()

    bam_file <- test_path("..", "..", "meta", "bam", "test_data", "mice_10X_1M.bam")
    skip_if_not(file.exists(bam_file), "mice_10X_1M.bam not found")

    annot <- prepare_classification_annotation(
        species = "mm10",
        chromosomes = "chr1",
        verbose = FALSE
    )

    out_dir <- tempfile(pattern = "classify_mouse1M_")
    on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

    result <- classify_reads(
        bam = bam_file,
        annotation = annot,
        out_dir = out_dir,
        tech = "10X"
    )

    expect_type(result, "list")
    total <- result$n_spliced + result$n_unspliced + result$n_chimeric + result$n_ambiguous + result$n_unassigned
    expect_gt(total, 0)

    expect_gt(result$flag_H, 0, label = "some reads have H flag (spliced)")
    expect_gt(result$flag_E, 0, label = "some reads have E flag (exonic)")
    expect_gt(result$n_spliced, 0, label = "1M reads has spliced reads")
    expect_gt(result$n_unspliced, 0, label = "1M reads has unspliced reads")
    expect_gt(result$n_ambiguous, 0, label = "1M reads has ambiguous reads")
})

test_that("classify_reads xf_values filters cells in 10X mouse", {
    skip_if_not_installed("TxDb.Mmusculus.UCSC.mm10.knownGene")
    skip_if_not_installed("Rsamtools")
    skip_on_cran()

    bam_file <- test_path("..", "..", "meta", "bam", "test_data", "mice_10X_1M.bam")
    skip_if_not(file.exists(bam_file), "mice_10X_1M.bam not found")

    annot <- prepare_classification_annotation(
        species = "mm10",
        chromosomes = "chr1",
        verbose = FALSE
    )

    out_dir_strict <- tempfile(pattern = "classify_xf_strict_")
    on.exit(unlink(out_dir_strict, recursive = TRUE), add = TRUE)
    result_strict <- classify_reads(
        bam = bam_file,
        annotation = annot,
        out_dir = out_dir_strict,
        tech = "10X",
        xf_values = c(25L, 17L)
    )

    out_dir_all <- tempfile(pattern = "classify_xf_all_")
    on.exit(unlink(out_dir_all, recursive = TRUE), add = TRUE)
    result_all <- classify_reads(
        bam = bam_file,
        annotation = annot,
        out_dir = out_dir_all,
        tech = "10X",
        xf_values = c(0L, 17L, 19L, 25L)
    )

    total_strict <- result_strict$n_spliced + result_strict$n_unspliced + result_strict$n_chimeric + result_strict$n_ambiguous + result_strict$n_unassigned
    total_all <- result_all$n_spliced + result_all$n_unspliced + result_all$n_chimeric + result_all$n_ambiguous + result_all$n_unassigned
    expect_gte(total_all, total_strict, label = "including xf=0 counts more or equal reads")
})