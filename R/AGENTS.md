# R/ — R Functions

R interface layer: C bindings, exon overlap mapping, CellRanger output, and exon clustering.

## STRUCTURE
```
R/
├── extract_exon_reads_hts.R  # C binding: BAM region scanning → TSV output
├── cb_umi_exons.R            # Pure R: block-to-exon overlap via data.table::foverlaps
├── build_matrix.R            # CellRanger output, matrix export, tag extraction wrapper
└── make_exon_clusters.R      # Exon clustering logic
```

## WHERE TO LOOK
| Task | File | Notes |
|------|------|-------|
| Add new C binding | `extract_exon_reads_hts.R` | Pattern: `.Call()` wrapper with input validation |
| Exon overlap logic | `cb_umi_exons.R` | `tidyr::separate_rows` + `data.table::foverlaps` |
| Classification wrapper | `classify_reads.R` | `classify_reads` calls C with HUGE bit field logic |
| CellRanger matrix | `build_matrix.R` | `bam_to_cellranger` calls C per chromosome |
| Matrix export formats | `build_matrix.R` | 10x (H5AD), loom, arrow, CSV |
| Tag extraction wrapper | `build_matrix.R` | `extract_unique_tags` — OpenMP C function |
| Cluster index serialization | `build_matrix.R` | `save_cluster_index` — RDS for C lookup |
| Exon clustering | `make_exon_clusters.R` | `make_exon_clusters` function |

## EXPORTED SYMBOLS
| Function | Type | Dependencies |
|----------|------|-------------|
| `extract_exon_reads_hts` | R → C `.Call` | htslib via `scan_core.c` |
| `cb_umi_exons` | Pure R | data.table, tidyr, dplyr |
| `bam_to_cellranger` | R → C `.Call` | GenomicRanges, htslib |
| `export_cell_cluster_matrix` | Pure R | zellkonverter/loomR/arrow (optional) |
| `save_cluster_index` | Pure R | GenomicRanges |
| `extract_unique_tags` | R → C `.Call` | OpenMP via `scan_core.c` |
| `make_exon_clusters` | Pure R | GenomicRanges, IRanges |

## CONVENTIONS
- **C bindings:** Validate inputs with `stopifnot()` before `.Call()`; convert NULL to `""` for C strings
- **Data flow:** TSV intermediate format between C output → R processing
- **Semicolon-delimited fields:** `block_start`, `block_end`, `block_seq` use `;` separator; expanded via `tidyr::separate_rows`
- **Overlap detection:** `data.table::foverlaps` with `type = "any"`, `nomatch = 0L`
- **GenomicRanges:** Used for cluster coordinates; extract via `seqnames()`, `start()`, `end()`

## ANTI-PATTERNS
- Manual `NAMESPACE` edits — use roxygen2 `@export` + `devtools::document()`
- Calling `.Call()` without input validation — always check file existence, types
- Mixing `data.table` and `dplyr` on same object — pick one per operation
- Forgetting `as.integer()` on numeric args passed to C

## NOTES
- `cb_umi_exons` expects TSV with columns: CB, UMI, block_start, block_end, block_seq
- Meta exons file must have: chr, start, end, exon columns
- Export formats (loom, arrow, 10x) require optional packages — check with `requireNamespace`
- C functions require package reinstall after `.c` file changes
