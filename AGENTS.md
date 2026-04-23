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
| Read classification (HUGE) | `src/classify_reads.c` | ExonClusterIndex + GeneIndex, HUGE bit field, OpenMP per-chr |
| OpenMP parallel tag extraction | `src/scan_core.c` | `extract_unique_tags` |
| CellRanger matrix output | `src/bam2db_ds.c` | SQLite + Matrix Market format |
| Hash table utility | `src/hashtable.c` | Chaining, DJB2 hash |
| C→R registration | `src/init.c` | 5 registered `.Call` entries |
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
- `system("mkdir -p")` in C code — use `mkdir_p()` from shared header instead

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
| `bam2db` | R → C `.Call` | `R/bam2db.R` |
| `classify_reads` | R → C `.Call` | `R/classify_reads.R` |
| `prepare_classification_annotation` | Pure R | `R/classify_reads.R` |

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

## Read Classification Logic (`src/classify_reads.c`)

OpenMP-parallel BAM iteration classifies each read into one of four categories using the **HUGE bit field** and an exon/gene overlap table. Each chromosome is processed independently by its own thread with a dedicated BAM file handle.

### HUGE bit field

Each read is annotated with a 4-bit flag:

| Bit | Flag | Meaning |
|-----|------|---------|
| 3   | `H`  | CIGAR contains `N` (BAM_CREF_SKIP) — spliced read |
| 2   | `U`  | Unique — all gene-mapped blocks belong to the same gene |
| 1   | `G`  | All blocks within gene region (no intergenic blocks) |
| 0   | `E`  | All blocks within exon cluster (every block hits an exon) |

### Block-level overlap classification

For each aligned block (CIGAR M/=/X segment):

| Exon overlap | Gene overlap | Effect on flags |
|-------------|-------------|----------------|
| yes | yes | Sets `E` (exonic block — if all blocks exonic, E=1) |
| no | yes | Clears `E`, sets `G` (intronic — within gene but not on exon) |
| no | no | Clears both `E` and `G` (intergenic — outside any gene) |

G and E are mutually exclusive: G is set only when all blocks are in gene regions **but not all are exonic** (some intronic). E is set when all blocks are exonic.

### Classification rules

| HUGE pattern | Category | Rationale |
|---------------|----------|-----------|
| `1111` (HUGE) | **SPLICED** | Junction read with confident unique gene + all blocks exonic + all in gene |
| `1011` (HGE, no U) | **CHIMERIC** | Junction read, all exonic and in gene, but gene assignment ambiguous across blocks |
| `_110` (UG, no E) | **UNSPLICED** | All blocks in gene, unique gene, not all exonic (some intronic) → pre-mRNA / nascent transcript |
| `0111` (UGE, no H) | **AMBIGUOUS** | Exon-only read with unique gene but no junction → mature mRNA signal |
| All other patterns | **UNASSIGNED** | Noise — mixed signals, no clear gene, or intergenic |

### Pre-filtering
- Skip unmapped, secondary, supplementary reads
- 10X: require `xf` tag in allowed values + `CB` barcode
- Smart-seq3Xpress: cell barcode via `BX` tag, UMI via `UX` tag
- Vanilla Smart-seq: use read name as cell ID (no barcode/UMI)

### Strand-aware annotation
- Per-chromosome `ChrAnnotNew` holds forward/reverse `ExonClusterIndex` and `GeneIndex` plus a `StrandGeneSpan` for gene span lookups
- Read routing: `BAM_FREVERSE` flag selects the correct strand's annotation for overlap queries
- `ExonClusterIndex`: sorted by end coordinate, stores `{start, end, gene_id, exon_cluster_id}` per exon — both overlap lookups and cluster counting use this single index
- `GeneIndex`: sorted by end coordinate, stores `{start, end, gene_id}` — used for intron overlap detection (intron intervals from exon gaps within genes)
- R side: `prepare_classification_annotation()` extracts strand from `BiocGenerics::strand()` for exons and introns; intron grouping key is `(chr, strand, gene_idx)`

### Classification (branching on CIGAR `N` skip)

**Spliced reads** (CIGAR contains `N` / `BAM_CREF_SKIP`):

Logic: intersect gene IDs across all aligned blocks via `ExonClusterIndex` overlap (strand-specific).
- All blocks share ≥1 common gene → `U` set → check `G` and `E` from block overlaps
- Gene intersection empty → `U` not set

Per-block check: exon overlap → contributes to `E`; if no exon overlap but gene overlap → sets `G`; if no exon and no gene overlap → clears both `E` and `G`.

- `HUGE` (1111) → **SPLICED** (confident junction read, all blocks exonic and in gene with unique gene)
- `HGE` (1011, no U) → **CHIMERIC** (junction read, all blocks exonic and in gene, but gene assignment ambiguous)
- `_110` (UG, no E) → **UNSPLICED** (all in gene, unique gene, not all exonic → intronic signal)
- `0111` (UGE, no H) → **AMBIGUOUS** (exon-only read with unique gene, no junction → mature mRNA)
- Other patterns → **UNASSIGNED**

### Key design point
Spliced vs unspliced is determined by CIGAR `N` operations (junction evidence), not by annotation overlap alone. A contiguous read overlapping an intron is unspliced (pre-mRNA). A spliced read whose blocks don't share a gene is ambiguous, not unspliced.

### OpenMP architecture
- **Per-chromosome parallelism**: `#pragma omp parallel for schedule(dynamic, 1)` over unique chromosomes
- Each thread opens its own `htsFile` (htslib not thread-safe for shared handles)
- Thread-local `CBHash` + `TripletCOO` per chromosome → merged after parallel region
- **CB remap**: build `local_id → global_id` array per chromosome, single-pass remap over COO triplets (O(n), not O(cells × triplets))
- **COO merge**: `triplet_add_count()` for O(1) merge across chromosomes
- Full chromosome range: `sam_itr_queryi(idx, tid, 0, INT_MAX)` (NOT `0, 0` which is empty)

### Output
Four CellRanger-format subdirectories: `spliced/`, `unspliced/`, `chimeric/`, `ambiguous/`, each with `barcodes.tsv.gz` + `features.tsv.gz` + `matrix.mtx.gz`. Unassigned reads are counted but not output.

Return value includes `by_chr` data.frame with per-chromosome counts (chr, spliced, unspliced, chimeric, ambiguous, unassigned, flag_H, flag_U, flag_G, flag_E) and global HUGE flag totals (flag_H, flag_U, flag_G, flag_E).

### `.Call` signature (20 args)
`_exonBlocks_classify_reads(bam_path, exon_chr, exon_starts, exon_ends, exon_gene_ids, exon_strand, intron_chr, intron_starts, intron_ends, intron_gene_ids, intron_strand, cluster_chr, cluster_starts, cluster_ends, cluster_strand, cluster_ids, tech, xf_values, out_dir, gene_span_pad)`

### Architecture
- **Shared utilities** (`src/exonblocks_shared.h`): CBHash, TripletCOO, CIGAR helpers, `find_first_overlap_end`, `mkdir_p`, `triplet_add_count`, `ExonClusterIndex`, `GeneIndex`, `GenomeAnnotation`, HUGE flag macros
- **Protocol-agnostic pre-filter** (`src/exonblocks_shared.h`): `TechProtocol` enum (10X / Smart-seq3Xpress / Smart-seq), `parse_tech()`, `read_passes_filter()`, `read_cell_id()` — reusable across all `.Call` entry points
- **`src/classify_reads.c`**: `ChrAnnotNew` (per-contig bundle with strand-split `ExonClusterIndex` + `GeneIndex` + `StrandGeneSpan`), HUGE bit field classification, OpenMP parallel loop, CB remap, COO merge
- **No `system()` calls**: `mkdir_p()` replaces `system("mkdir -p")`

### Performance (mice_brain_10X.bam, mm10, 16 OpenMP threads)
- 171.5M reads classified in **57 sec** → **3.0M reads/sec**
- 1,148,095 cells (with xf={0,17,19,25}); 784,294 cells (xf={17,25} only)
- 286,305 exon cluster features
- Breakdown: 5.1% spliced, 8.8% unspliced, 0.0% chimeric, 83.1% ambiguous, 3.0% unassigned

### Cell count caveat: xf=0 inflates cell barcodes
Including `xf=0` (CellRanger "unclassified") reads inflates the cell count because xf=0 barcodes include ambient RNA / empty droplets rejected by CellRanger's cell calling. In mice_brain_10X:
- `xf={17,25}` → 784K cells (high-confidence)
- `xf={0,17,19,25}` → 1.15M cells (includes ~370K likely empty droplets)
- The same cell barcode can appear across multiple xf values — this is correct (not double-counting); the inflation comes from xf=0-only barcodes that never appear with xf=17/25
- **Recommendation**: default `xf_values = c(25L, 17L)` for cell-level analysis; include xf=0 only for total-read quantification

### Per-chromosome finding: chr19 anomaly
chr19 shows 71.6% chimeric reads — much higher than the genome average of 29.9% — with only 9.6% unspliced vs the 30.1% average. This suggests chr19's exon annotation is dense with short exons that catch many contiguous reads as chimeric (exon-only, no intron overlap) rather than classifying them as unspliced. This pattern may reflect biological features (gene-dense chromosome with short introns) or annotation artifacts worth investigating.

### Per-chromosome summary (mice_brain_10X.bam, mm10, xf={0,17,19,25})

| Chr | Total reads | Spliced | Unspliced | Chimeric | Unassigned | %Spliced | %Unspliced | %Chimeric |
|-----|-----------:|--------:|----------:|---------:|-----------:|--------:|-----------:|----------:|
| 2 | 14,319,623 | 797,917 | 4,643,919 | 3,982,238 | 4,895,549 | 5.6 | 32.4 | 27.8 |
| 11 | 10,636,688 | 929,206 | 2,875,674 | 3,586,474 | 3,245,334 | 8.7 | 27.0 | 33.7 |
| 7 | 10,543,099 | 1,136,220 | 2,653,414 | 3,550,635 | 3,202,830 | 10.8 | 25.2 | 33.7 |
| 19 | 10,416,493 | 462,318 | 1,003,129 | 7,457,650 | 1,493,396 | 4.4 | 9.6 | 71.6 |
| 5 | 10,081,717 | 626,051 | 3,139,776 | 2,635,331 | 3,680,559 | 6.2 | 31.1 | 26.1 |
| 1 | 10,073,616 | 547,317 | 3,482,486 | 2,368,895 | 3,674,918 | 5.4 | 34.6 | 23.5 |
| 4 | 9,360,649 | 580,344 | 3,030,862 | 2,549,936 | 3,199,507 | 6.2 | 32.4 | 27.2 |
| 9 | 9,302,498 | 633,932 | 2,771,263 | 2,766,790 | 3,130,513 | 6.8 | 29.8 | 29.7 |
| 6 | 9,157,759 | 395,821 | 3,329,362 | 2,028,458 | 3,404,118 | 4.3 | 36.4 | 22.2 |
| 8 | 8,566,182 | 631,771 | 2,809,715 | 2,258,055 | 2,866,641 | 7.4 | 32.8 | 26.4 |
| 12 | 8,347,230 | 450,983 | 2,352,284 | 2,510,618 | 3,033,345 | 5.4 | 28.2 | 30.0 |
| 10 | 8,234,341 | 569,910 | 2,714,194 | 2,253,476 | 2,696,761 | 6.9 | 33.0 | 27.4 |
| 17 | 7,955,216 | 589,645 | 2,125,439 | 2,410,918 | 2,829,214 | 7.4 | 26.7 | 30.3 |
| 3 | 7,821,242 | 481,090 | 2,617,899 | 1,922,182 | 2,800,071 | 6.2 | 33.5 | 24.6 |
| 16 | 7,238,289 | 353,668 | 2,629,772 | 1,385,805 | 2,869,044 | 4.9 | 36.3 | 19.1 |
| 14 | 7,181,829 | 311,084 | 2,549,950 | 1,433,826 | 2,886,969 | 4.3 | 35.5 | 20.0 |
| 15 | 6,318,400 | 578,890 | 1,866,130 | 1,792,241 | 2,081,139 | 9.2 | 29.5 | 28.4 |
| 13 | 6,069,740 | 261,553 | 2,082,063 | 1,452,201 | 2,273,923 | 4.3 | 34.3 | 23.9 |
| X | 5,113,318 | 286,634 | 1,390,825 | 1,709,463 | 1,726,396 | 5.6 | 27.2 | 33.4 |
| 18 | 4,761,188 | 295,418 | 1,506,923 | 1,244,598 | 1,714,249 | 6.2 | 31.7 | 26.1 |
| Y | 30,298 | 2,867 | 5,705 | 3,960 | 17,766 | 9.5 | 18.8 | 13.1 |

**Overall**: 171.5M reads, 6.4% spliced, 30.1% unspliced, 29.9% chimeric, 33.7% unassigned.

Notable observations:
- chr19 is the clear outlier: 71.6% chimeric (2.4x genome average), only 9.6% unspliced (0.3x average) — dense short-exon annotation causes contiguous reads to be caught as chimeric rather than unspliced
- chr7 and chr11 have the highest spliced proportions (10.8%, 8.7%) — longer genes with more junction reads
- chr6, chr16, chr14 have the highest unspliced proportions (36.4%, 36.3%, 35.5%) — gene-sparse regions where reads often span introns
- chr15 and chr7 have the most balanced spliced/unspliced ratios — potential sweet spots for RNA velocity

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
