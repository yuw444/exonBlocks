# CRAN Submission Comments — exonBlocks

Package: exonBlocks  
Version submitted: 0.1.1  
Maintainer: Yu Wang <ywang@mcw.edu>

## Purpose
exonBlocks maps cell-barcode + UMI read blocks (from BAM-derived TSVs) to annotated exons. It expands semicolon-delimited block fields (block_start, block_end, block_seq) and uses data.table::foverlaps() to detect overlaps.

## Changes since last release
- Fix DESCRIPTION formatting (valid DCF).
- Remove non-portable debug flags (-g -O0) from src/Makevars.
- Add R/globals.R to declare NSE column names used by dplyr/data.table to silence R CMD check NOTES.
- Import magrittr pipe (%>%) in NAMESPACE via roxygen.
- Improve roxygen docs for cb_umi_exons() and add examples.
- Add NEWS.md and a README with usage examples.
- Add/repair basic testthat tests.

## Package contents / dependencies
- R package with C code that required for functionality.
- Imports: data.table, dplyr, tidyr, magrittr.
- Suggests: testthat.
- License: MIT (LICENSE included).

## Known check notes and how they were addressed
- cb_umi_exons: no visible binding for global variable ‘.data’
Undefined global functions or variables:
  .data
    - no treatement needed as .data are reserved variable within `dplyr` package

## Tests and examples
- testthat tests are included under tests/testthat and pass locally with devtools::check().
- Examples in Rd documentation are guarded with \dontrun{} where they require large BAM files and user data (so R CMD check on CRAN will not try to run long, data-dependent examples).

## System requirements / external resources
- An indexed BAM is required to run the example workflow that calls extract_exon_reads_hts(); however, no external data or network access is required for package installation or running unit tests.
- htslib(https://github.com/samtools/htslib) is required.

## Additional notes for reviewers
- If you see the git related file, it will be removed in the final submission.
- If CRAN requires changes to any packaging metadata, please advise and I will address them promptly.

For questions or follow-ups: Yu Wang <ywang@mcw.edu>