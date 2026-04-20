# exonBlocks AGENTS.md

R package + C extensions for single-cell RNA-seq exon block analysis from BAM files. Maps cell barcodes + UMI-level read blocks to overlapping exons.

## STRUCTURE
```
exonBlocks/
├── R/                      # R functions (4 files)
├── src/                    # C extensions (htslib, OpenMP, SQLite)
├── tests/testthat/         # testthat Edition 3 tests
├── man/                    # Roxygen2-generated .Rd docs (auto)
├── meta/                   # Sample data (TSV, CSV)
└── .vscode/                # Pre-configured IDE settings
```

## WHERE TO LOOK
| Task | Location | Notes |
|------|----------|-------|
| BAM scanning, block extraction | `src/scan_core.c` | htslib iterator, CIGAR parsing |
| OpenMP parallel tag extraction | `src/scan_core.c` | `extract_unique_tags` |
| CellRanger matrix output | `src/bam2db_ds.c` | SQLite + Matrix Market format |
| Hash table utility | `src/hashtable.c` | Chaining, DJB2 hash |
| C→R registration | `src/init.c` | 3 registered `.Call` entries |
| Exon overlap mapping | `R/cb_umi_exons.R` | data.table::foverlaps |
| CellRanger output wrapper | `R/build_matrix.R` | `bam_to_cellranger`, export formats |
| Build config | `src/Makevars` | htslib autodetect chain |

## CONVENTIONS
- **R style:** Mixture of `data.table` (performance) + `dplyr/tidyr` (readability)
- **C error handling:** Use `error()` from R API, NEVER `exit()` or `fprintf(stderr)` + return in registered entry points
- **Documentation:** Roxygen2 `#'` comments → `devtools::document()` → auto-generates `NAMESPACE` + `man/`
- **NEVER edit `NAMESPACE` manually**
- **C registration:** All `.Call` functions must be declared in `src/init.c` `CallEntries[]`

## ANTI-PATTERNS
- `exit()` in C code — use `error()` (R API) instead
- `@ts-ignore` / `as any` — N/A but suppress type errors
- Manual `NAMESPACE` edits
- Using `fprintf(stderr)` in registered `.Call` entry points (use `error()`/`warning()`)
- `bam2db_ds.c` uses `fprintf(stderr)` — legacy pattern, don't replicate

## COMMANDS
```bash
# Install/rebuild (micromamba env: exonBlock)
R -e 'devtools::install()'

# Tests
R -e 'devtools::test()'

# Docs
R -e 'devtools::document()'

# C syntax check (micromamba paths)
gcc -fsyntax-only -c src/scan_core.c \
    -I/home/yu-wang/micromamba/envs/exonBlock/lib/R/include \
    -I/home/yu-wang/micromamba/envs/exonBlock/include \
    -Isrc

# C syntax check (HPC module paths, if using R/4.5.0)
gcc -fsyntax-only -c src/scan_core.c \
    -I/usr/lib/R/lib64/R/include \
    -I/usr/local/include \
    -Isrc
```

## KEY PATHS
| Component | Path (micromamba) | Path (HPC module) |
|-----------|-------------------|-------------------|
| R Runtime | `/home/yu-wang/micromamba/envs/exonBlock/` | `/usr/lib/R/` |
| R Headers | `.../lib/R/include` | `/usr/lib/R/lib64/R/include` |
| htslib | `.../include` | `/usr/local/include` |
| R Library | `.../lib/R/library` | `~/R/x86_64-pc-linux-gnu-library/4.5.0/` |

## EXPORTED SYMBOLS
| Symbol | Type | Source |
|--------|------|--------|
| `extract_exon_reads_hts` | R → C `.Call` | `R/extract_exon_reads_hts.R` |
| `cb_umi_exons` | Pure R | `R/cb_umi_exons.R` |
| `bam_to_cellranger` | R → C `.Call` | `R/build_matrix.R` |
| `export_cell_cluster_matrix` | Pure R | `R/build_matrix.R` |
| `save_cluster_index` | Pure R | `R/build_matrix.R` |
| `extract_unique_tags` | R → C `.Call` | `R/build_matrix.R` |
| `make_exon_clusters` | Pure R | `R/make_exon_clusters.R` |

## NOTES
- Modifying any `.c` file requires reinstall (`R CMD INSTALL .`) to rebuild shared object
- `src/bam2db_ds.c` is legacy — not registered in `init.c`, not exported. Contains SQLite-based pipeline
- `bam2db_ds.c` uses `fprintf(stderr)` pattern — do NOT replicate in new code
- See `COMPILE_ENV.md` for full HPC environment docs

## featureCounts ES/IS/GE/GI Tags (human_pbmc_Smart-seq3xpressV3.bam)

featureCounts assigns per-read tags for exon and intron classification:

| Tag | Meaning | Set when |
|-----|---------|----------|
| `ES` | Exon Status | Assigned / Ambiguity / NoFeatures |
| `IS` | Intron Status | Assigned / Ambiguity / NoFeatures |
| `GE` | Gene ID (exon-assigned) | Exon confidently maps to one gene |
| `GI` | Gene ID (intron-assigned) | Intron confidently maps to one gene |

`GE` and `GI` are only populated when the respective level gets a confident assignment. Unassigned reads have empty GE/GI.

### ES x IS cross-tabulation (chr1, first 10M reads)

| ES ↓ / IS → | Assigned3 | NoFeatures | Ambiguity | Total |
|---|---|---|---|---|
| **Assigned3** | 332,049 | 5,754,106 | — | 6,086,155 |
| **NoFeatures** | 351,316 | 3,521,114 | 16 | 3,872,446 |
| **Ambiguity** | 60 | 41,329 | — | 41,389 |

### Read categories

- **5.75M (57.5%)** ES=Assigned / IS=NoFeatures → **exon-only reads** (mature mRNA, no intron overlap)
- **3.52M (35.2%)** ES=NoFeatures / IS=NoFeatures → **intergenic/UTR/unannotated** (outside GTF features)
- **351K (3.5%)** ES=NoFeatures / IS=Assigned → **intron-only reads** (pre-mRNA/nascent transcripts with GI but no GE)
- **332K (3.3%)** ES=Assigned / IS=Assigned → **exon-intron spanning reads** (spliced junction reads; both GE and GI populated — key for RNA velocity)
- **41K** ES=Ambiguity / IS=NoFeatures → reads overlapping exons of multiple genes, no intron overlap to disambiguate
- **154 reads** in first 500K of chr19 had GE ≠ GI → reads bridging exons of one gene and introns of another (overlapping/antisense gene loci)

### Biological interpretation

- Exon-only reads = high-confidence mature mRNA expression
- Intron-only reads = transcriptional activity but not spliced; nuclear/pre-mRNA signal
- Exon-intron spanning = junction reads confirming splicing; most informative for isoform resolution and RNA velocity (spliced vs unspliced)
- GE ≠ GI reads = intergenic splicing ambiguity or overlapping gene architectures; should be flagged/excluded to avoid misattribution

## Git Commit Policy (Standing Rule)

No issue may be marked done without a corresponding git commit pushed to the `git-ai` remote. This is non-negotiable.

### Per-repo git identity

Before working in any repo, verify git config is set correctly:

- exonBlock: `user.name "Senior Bioinformatics Engineer"`, `user.email "sre@omilab"`

If misconfigured, fix it immediately before committing.

### Commit discipline

- Every completed task must result in at least one commit pushed to `git-ai`
- Commit messages must be descriptive
- Include `Co-Authored-By: Paperclip <noreply@paperclip.ing>` on every commit
