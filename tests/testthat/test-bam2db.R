test_that("bam2db works with zUMIs protocol + gene features", {
  skip_on_cran()
  library(exonBlocks)
  human_bam <- test_path("/scratch/g/chlin/Yu/exonBlocks/meta/bam/human_pbmc_Smart-seq3xpressV3.bam")
  skip_if_not(file.exists(human_bam), "Human BAM not found")
  gene_tsv <- test_path("/scratch/g/chlin/Yu/exonBlocks/meta/features/homo_sapiens_GRCh38_gene.tsv")
  skip_if_not(file.exists(gene_tsv), "Gene features not found")
  bc_file <- test_path("/scratch/g/chlin/Yu/exonBlocks/meta/features/human_barcodes_chr19.txt")
  skip_if_not(file.exists(bc_file), "Human barcodes not found")

  out_dir <- "/scratch/g/chlin/Yu/exonBlocks/meta/output/bam2db_zumis_gene/"
  dir.create(out_dir, recursive = TRUE)

  db_file <- file.path(out_dir, "test.db")

  bam2db(
    bam_file = human_bam,
    db_file = db_file,
    path_out = out_dir,
    barcodes_file = bc_file,
    regions_file = gene_tsv,
    rate_cell = 0.1,
    rate_depth = 0.1,
    seed = 42,
    umi_copies_flag = 0L,
    protocol = "zUMIs",
    region_chr = "19",
    region_start = 1L,
    region_end = 10000000L
  )

  expect_true(file.exists(db_file))
  expect_true(file.exists(file.path(out_dir, "barcodes.tsv.gz")))
  expect_true(file.exists(file.path(out_dir, "matrix.mtx.gz")))
  expect_true(file.exists(file.path(out_dir, "features.tsv.gz")))
})

test_that("bam2db works with zUMIs protocol + exon features", {
  skip_on_cran()
  library(exonBlocks)
  human_bam <- test_path("/scratch/g/chlin/Yu/exonBlocks/meta/bam/human_pbmc_Smart-seq3xpressV3.bam")
  skip_if_not(file.exists(human_bam), "Human BAM not found")
  exon_tsv <- test_path("/scratch/g/chlin/Yu/exonBlocks/meta/features/homo_sapiens_GRCh38_exon.tsv")
  skip_if_not(file.exists(exon_tsv), "Exon features not found")
  bc_file <- test_path("/scratch/g/chlin/Yu/exonBlocks/meta/features/human_barcodes_chr19.txt")
  skip_if_not(file.exists(bc_file), "Human barcodes not found")

  out_dir <- "/scratch/g/chlin/Yu/exonBlocks/meta/output/bam2db_zumis_exon/"
  dir.create(out_dir, recursive = TRUE)

  db_file <- file.path(out_dir, "test.db")

  bam2db(
    bam_file = human_bam,
    db_file = db_file,
    path_out = out_dir,
    barcodes_file = bc_file,
    regions_file = exon_tsv,
    rate_cell = 0.1,
    rate_depth = 0.1,
    seed = 42,
    umi_copies_flag = 0L,
    protocol = "zUMIs",
    region_chr = "19",
    region_start = 1L,
    region_end = 10000000L
  )

  expect_true(file.exists(db_file))
  expect_true(file.exists(file.path(out_dir, "matrix.mtx.gz")))
})

test_that("bam2db works with 10X protocol + gene features", {
  skip_on_cran()
  library(exonBlocks)
  mouse_bam <- test_path("/scratch/g/chlin/Yu/exonBlocks/meta/bam/mice_brain_10X.bam")
  skip_if_not(file.exists(mouse_bam), "Mouse BAM not found")
  gene_tsv <- test_path("/scratch/g/chlin/Yu/exonBlocks/meta/features/mus_musculus_GRCm39_gene.tsv")
  skip_if_not(file.exists(gene_tsv), "Mouse gene features not found")
  bc_file <- test_path("/scratch/g/chlin/Yu/exonBlocks/meta/features/mouse_barcodes_chr19.txt")
  skip_if_not(file.exists(bc_file), "Mouse barcodes not found")

  out_dir <- "/scratch/g/chlin/Yu/exonBlocks/meta/output/bam2db_10x_gene/"
  dir.create(out_dir, recursive = TRUE)

  db_file <- file.path(out_dir, "test.db")

  bam2db(
    bam_file = mouse_bam,
    db_file = db_file,
    path_out = out_dir,
    barcodes_file = bc_file,
    regions_file = gene_tsv,
    rate_cell = 0.1,
    rate_depth = 0.1,
    seed = 42,
    umi_copies_flag = 0L,
    protocol = "10X",
    region_chr = "19",
    region_start = 1L,
    region_end = 10000000L
  )

  expect_true(file.exists(db_file))
  expect_true(file.exists(file.path(out_dir, "matrix.mtx.gz")))
})

test_that("bam2db works with 10X protocol + exon features", {
  skip_on_cran()
  library(exonBlocks)
  mouse_bam <- test_path("/scratch/g/chlin/Yu/exonBlocks/meta/bam/mice_brain_10X.bam")
  skip_if_not(file.exists(mouse_bam), "Mouse BAM not found")
  exon_tsv <- test_path("/scratch/g/chlin/Yu/exonBlocks/meta/features/mus_musculus_GRCm39_exon.tsv")
  skip_if_not(file.exists(exon_tsv), "Mouse exon features not found")
  bc_file <- test_path("/scratch/g/chlin/Yu/exonBlocks/meta/features/mouse_barcodes_chr19.txt")
  skip_if_not(file.exists(bc_file), "Mouse barcodes not found")

  out_dir <- "/scratch/g/chlin/Yu/exonBlocks/meta/output/bam2db_10x_exon/"
  dir.create(out_dir, recursive = TRUE)

  db_file <- file.path(out_dir, "test.db")

  bam2db(
    bam_file = mouse_bam,
    db_file = db_file,
    path_out = out_dir,
    barcodes_file = bc_file,
    regions_file = exon_tsv,
    rate_cell = 0.1,
    rate_depth = 0.1,
    seed = 42,
    umi_copies_flag = 0L,
    protocol = "10X",
    region_chr = "19",
    region_start = 1L,
    region_end = 10000000L
  )

  expect_true(file.exists(db_file))
  expect_true(file.exists(file.path(out_dir, "matrix.mtx.gz")))
})
