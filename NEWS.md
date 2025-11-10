# exonBlocks NEWS

All notable changes to this project will be documented in this file.


## 0.1.0 - Initial release
- Features
  - cb_umi_exons(): Expand semicolon-delimited block fields (block_start, block_end, block_seq) and map each expanded block to overlapping exons using `data.table::foverlaps()`.
  - extract_exon_reads_hts(): Helper to extract read blocks overlapping a genomic interval and write block TSVs (example workflow included).
- Expectations
  - Blocks TSV must contain columns: `CB`, `UMI`, `block_start`, `block_end`, `block_seq`.
  - Exon metadata TSV must contain columns: `chr`, `start`, `end`, `exon`.
- Testing
  - Basic testthat tests and example scripts included.

---

Reporting issues
- Please report bugs or feature requests via the repository issue tracker with a minimal reproducible example and the output of sessionInfo().

License
- MIT (see LICENSE file)