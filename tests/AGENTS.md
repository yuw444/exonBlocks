# tests/ — Test Suite

testthat Edition 3 tests for R functions and C bindings.

## STRUCTURE
```
tests/
├── testthat.R              # testthat runner
└── testthat/
    ├── test-scan.R         # BAM scanning, C binding tests
    ├── test-barcodes.R     # Barcode extraction tests
    ├── test-build_matrix.R # CellRanger matrix output tests
    ├── test-exon-cluster.R # Exon clustering tests
    └── test-umi-counts.R   # UMI counting tests
```

## WHERE TO LOOK
| Task | File | Notes |
|------|------|-------|
| BAM scan tests | `test-scan.R` | `extract_exon_reads_hts` C binding |
| Barcode tests | `test-barcodes.R` | CB/UB tag extraction |
| Matrix output tests | `test-build_matrix.R` | `bam_to_cellranger`, export formats |
| Cluster tests | `test-exon-cluster.R` | `make_exon_clusters` |
| UMI tests | `test-umi-counts.R` | `extract_unique_tags` |

## CONVENTIONS
- **Framework:** testthat Edition 3 (`Config/testthat/edition: 3` in DESCRIPTION)
- **Test structure:** `test_that("description", { ... })`
- **File naming:** `test-<feature>.R`
- **C code testing:** Package must be installed/loaded for DLL availability
- **Fixtures:** Use `meta/` directory sample data when possible

## COMMANDS
```bash
# Run all tests
R -e 'devtools::test()'

# Run single test file
R -e 'testthat::test_file("tests/testthat/test-scan.R")'

# Run filtered by description
R -e 'testthat::test_file("tests/testthat/test-scan.R", filter="scan_bam works")'
```

## NOTES
- Tests requiring BAM files need actual indexed BAM data — check for skip conditions
- C binding tests require package reinstall after `.c` changes
- testthat Edition 3 uses `testthat::test_file()` directly, not `source()`
