# src/ — C Extensions

C code for BAM/CRAM processing, parallel tag extraction, and CellRanger matrix generation.

## STRUCTURE
```
src/
├── scan_core.c          # Core BAM scanning, CIGAR parsing, OpenMP parallel tag extraction
├── scan_core_stream.c   # `_exonBlocks_bam_to_cellranger` — CellRanger matrix output
├── classify_reads.c     # Read classification (NSEIG), OpenMP per-chromosome BAM pass
├── bam2db_ds.c          # Legacy: SQLite + CellRanger (EXCLUDED from build — broken deps)
├── bam2db_ds.h          # Header for bam2db_ds.c
├── hashtable.c          # Generic hash table (chaining, DJB2 hash)
├── hashtable.h          # Hash table header
├── exonblocks_shared.h  # Shared definitions (ExonClusterIndex, GeneIndex, CBHash, etc.)
├── init.c               # R native routine registration (5 entries)
├── Makevars             # Explicit SOURCES + htslib autodetect chain
└── Makevars.win         # Windows build config
```

## WHERE TO LOOK
| Task | File | Notes |
|------|------|-------|
| Add new `.Call` entry | `init.c` | Must add to `CallEntries[]` + declare function |
| BAM region scanning | `scan_core.c` | `scan_bam_blocks_hts` — htslib iterator, CIGAR parsing |
| Read classification | `classify_reads.c` | HUGE bit field, `ExonClusterIndex`+`GeneIndex`, OpenMP per-chr, 5 categories |
| Parallel tag extraction | `scan_core.c` | `extract_unique_tags` — OpenMP round-robin contig scanning |
| Block-to-cluster mapping | `scan_core.c` | `map_blocks_to_clusters` — binary search overlap |
| Hash table operations | `hashtable.c` | Create/insert/lookup/delete, chaining collision |
| Legacy SQLite pipeline | `bam2db_ds.c` | NOT in init.c, uses `fprintf(stderr)` + `exit()` |
| Build flags | `Makevars` | htslib autodetect: env → pkg-config → HTSLIB_DIR → htsfile → Homebrew/conda |

## REGISTERED C FUNCTIONS
| R Symbol | C Function | Args | Source |
|----------|-----------|------|--------|
| `_exonBlocks_scan_bam_blocks_hts` | `scan_bam_blocks_hts` | 8 | `scan_core.c` |
| `_exonBlocks_bam_to_cellranger` | `bam_to_cellranger` | 8 | `scan_core_stream.c` |
| `_exonBlocks_extract_unique_tags` | `extract_unique_tags` | 3 | `scan_core.c` |
| `_exonBlocks_classify_reads` | `_exonBlocks_classify_reads` | 20 | `classify_reads.c` |

## KEY STRUCTURES (`exonblocks_shared.h`)
| Struct | Fields | Purpose |
|--------|--------|---------|
| `ExonClusterIndex` | `n`, `start`, `end`, `gene_id`, `exon_cluster_id` | Sorted exon intervals for overlap queries + cluster mapping |
| `GeneIndex` | `n`, `start`, `end`, `gene_id` | Sorted gene span intervals for intron/intergenic detection |
| `GenomeAnnotation` | `n_contigs`, `exon_index`, `gene_index` | Top-level annotation container (not yet used in C, planned) |
| `CBHash` / `CBEntry` | DJB2 hash table | Cell barcode → integer ID |
| `TripletCOO` | `cell_idx`, `cluster_idx`, `counts` | Sparse COO matrix accumulator |

## HUGE BIT FIELD (`exonblocks_shared.h`)
| Bit | Flag | Meaning |
|-----|------|---------|
| 3 | `H` | CIGAR has N (spliced) |
| 2 | `U` | Unique gene assignment across all blocks |
| 1 | `G` | All blocks within gene region |
| 0 | `E` | All blocks within exon cluster |

Classification: HUGE(1111)→SPLICED, HGE(1011)→CHIMERIC, _110(UGE?)→UNSPLICED (UG no E), 0111(UGE)→AMBIGUOUS, other→UNASSIGNED

## CONVENTIONS
- **Error handling:** Use `error()` / `warning()` from R API — NEVER `exit()` or `fprintf(stderr)` in registered entry points
- **Memory:** Use `PROTECT`/`UNPROTECT` for R objects; `malloc`/`free` for C memory; always cleanup on error paths
- **htslib cleanup pattern:** Destroy iterators → headers → close files (order matters)
- **OpenMP:** Round-robin contig distribution across threads; `#pragma omp critical` for hashtable merge

## ANTI-PATTERNS
- `exit()` in registered `.Call` entry points — use `error()` instead
- `fprintf(stderr)` in new code — `bam2db_ds.c` uses it (legacy), don't replicate
- Skipping cleanup on early returns — every error path must free all allocated resources
- Manual `NAMESPACE` edits — roxygen2 handles this
- Forgetting to register new C functions in `init.c` `CallEntries[]`

## NOTES
- `bam2db_ds.c` is legacy/excluded — missing Mersenne Twister RNG deps, NOT in Makevars SOURCES
- `Makevars` uses explicit `SOURCES = init.c scan_core.c scan_core_stream.c hashtable.c` — add new `.c` files here
- After editing ANY `.c` file: `R CMD INSTALL .` required to rebuild shared object
- `Makevars` has 6-step htslib detection — test with `make show-flags`
- OpenMP compilation requires `-fopenmp` flag (set in Makevars)
