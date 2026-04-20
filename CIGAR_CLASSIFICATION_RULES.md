# CIGAR-Based Read Classification Rules

## CIGAR Operators

| Op | Name | Ref | Read | Meaning |
|----|------|-----|------|---------|
| M  | Match/mismatch | ✓ | ✓ | Aligned base (may include mismatches) |
| =  | Sequence match | ✓ | ✓ | Confirmed match |
| X  | Sequence mismatch | ✓ | ✓ | Confirmed mismatch |
| N  | Skip (intron) | ✓ | ✗ | **Splice junction — intron removed** |
| D  | Deletion | ✓ | ✗ | Deletion in reference (not a splice) |
| I  | Insertion | ✗ | ✓ | Insertion in read (not in reference) |
| S  | Soft clip | ✗ | ✓ | Clipped sequence present in read |
| H  | Hard clip | ✗ | ✗ | Clipped sequence absent from read |
| P  | Padding | ✗ | ✗ | Padding for multiple alignments |

## Classification Rules

### Rule 0: Extract aligned blocks from CIGAR

Split each read into aligned reference intervals (blocks) by iterating CIGAR:

- `M`, `=`, `X` → emit block `[pos, pos+len)`, advance pos
- `N` → advance pos (gap between blocks = intron)
- `D` → advance pos (deletion, NOT a splice)
- `I`, `S` → do not advance pos (read-only operations)
- `H`, `P` → no-op

```
CIGAR: 76M          → 1 block  [pos, pos+76)
CIGAR: 50M150N26M   → 2 blocks [pos, pos+50), [pos+50+150, pos+226)
CIGAR: 30M3I20M     → 1 block  [pos, pos+50)  (I does not advance ref)
CIGAR: 10S76M       → 1 block  [pos, pos+76)
CIGAR: 50M5D30M     → 1 block  [pos, pos+85)  (D treated as contiguous — deletion, not splice)
```

### Rule 1: Splice status from CIGAR

**Check: Does CIGAR contain `N` operator?**

```
has_N = any CIGAR operation == N
```

| Condition | `has_N` | `num_blocks` | Classification | Rationale |
|-----------|---------|--------------|----------------|-----------|
| N in CIGAR | true | ≥ 2 | **SPLICED** | Read spans intron(s); spliced junction read |
| No N, single block | false | 1 | **Requires exon annotation** (see Rule 2) | No splice junction; could be exon-only or intron-only |
| No N, D splits | false | 1* | **Requires exon annotation** | D is a deletion, not a splice; treated as single contiguous block |

\* D does not split blocks. A read with `50M5D30M` is one contiguous reference interval `[pos, pos+85)`.

### Rule 2: Single-block reads — exon annotation required

For reads with `has_N == false` and `num_blocks == 1`, classify using exon annotation overlap:

```
block = [block_start, block_end]
overlap = foverlaps(block, exons)
```

| Overlap result | Classification | Meaning |
|---------------|----------------|---------|
| Block overlaps ≥1 exon, no intron | **SPLICED** | Exon-only mature mRNA |
| Block overlaps ≥1 intron, no exon | **UNSPLICED** | Intronic pre-mRNA/nascent transcript |
| Block overlaps both exon and intron of same gene | **AMBIGUOUS** | Pre-mRNA with partial exon overlap |
| Block overlaps exon of gene A and intron of gene B | **AMBIGUOUS_CONFLICT** | Overlapping/antisense gene loci |
| Block overlaps no annotation | **UNASSIGNED** | Intergenic / UTR / unannotated |

### Rule 3: Multi-block spliced reads — per-block exon assignment

For reads with `has_N == true` and `num_blocks ≥ 2`, each block is independently overlapped against exons:

```
for each block in read:
    exon_hits = foverlaps(block, exons)
```

| Block position | Overlap result | Assigned to |
|---------------|---------------|-------------|
| Block overlaps exon E1 | foverlaps hit | Exon E1 (row in output) |
| Gap between blocks | No overlap | Intronic — no row emitted |
| Block overlaps no exon | No row | Unannotated region of spliced read |

A spliced read hitting exons E1 and E2 with an intron gap produces **2 rows** in the output, one per exon block. The intron gap produces no row.

### Rule 4: D (deletion) is NOT a splice

```
CIGAR: 50M3D30M → single block [pos, pos+83)
```

D represents a small deletion in the reference genome, not an intron. It does NOT split the read into multiple blocks and does NOT make the read spliced. Only N creates multiple blocks.

### Rule 5: Soft/hard clip handling

```
CIGAR: 10S76M → single block [pos, pos+76)   (S bases ignored for block)
CIGAR: 5H76M  → single block [pos, pos+76)   (H bases not in read at all)
```

Clipped bases are excluded from block coordinates. They do not affect splice classification.

### Rule 6: GE vs GI consistency check (when featureCounts tags available)

For reads where both GE and GI are populated:

| Condition | Action |
|-----------|--------|
| GE == GI | Normal — same gene, no conflict |
| GE != GI | Flag as CONFLICT — read bridges exon of gene A and intron of gene B. Exclude from gene-level counting or assign fractionally. |

## Decision Flowchart

```
Read
 │
 ├─ Parse CIGAR → extract blocks
 │
 ├─ CIGAR has N?
 │    │
 │    YES → SPLICED (junction read)
 │    │      Each block → foverlaps with exons
 │    │      Emit 1 row per exon-overlapping block
 │    │
 │    NO → Single contiguous block
 │           │
 │           ├─ foverlaps with exon annotation
 │           │
 │           ├─ Overlaps exon only    → SPLICED (mature mRNA)
 │           ├─ Overlaps intron only   → UNSPLICED (pre-mRNA)
 │           ├─ Overlaps both          → AMBIGUOUS (boundary read)
 │           └─ Overlaps nothing       → UNASSIGNED (intergenic/UTR)
 │
 └─ Output: CB  UMI  block_start  block_end  block_seq  num_blocks  spliced
            spliced = 1 if has_N, else 0
```

## Relationship to Existing Tools

### zUMIs (annotation-only, no CIGAR)

zUMIs uses two-pass featureCounts (exon SAF → intron SAF) to classify reads. It does **not** parse CIGAR. Categories:

| zUMIs Category | Rule |
|---|---|
| **Unused-BC** | Barcode not in kept cell whitelist |
| **Unmapped** | STAR could not map the read |
| **Ambiguity** | ES or IS = `Unassigned_Ambiguity` |
| **Intergenic** | ES=NoFeatures AND IS=NoFeatures |
| **Intron** | ES=NoFeatures AND IS=Assigned |
| **Exon** | ES=Assigned (regardless of GI) |

Key design choices:
- **Exon overrides intron**: If a read has both GE and GI, it is classified as **Exon**, not "both"
- Intron SAF constructed from gene-internal gaps (GTF exon gaps filtered to be within gene boundaries, width 10–100kb)
- Two-pass: first assign exons, then feed unassigned reads to introns
- Uses `largestOverlap=TRUE` and `allowMultiOverlap` for tie-breaking
- Offers intron probability score (truncated Poisson) to validate intronic signal above intergenic background

### Velocyto (CIGAR-based, annotation-assisted)

Velocyto uses **CIGAR N-operations** (splice junctions) plus GTF annotation overlap:

| Velocyto Category | Rule |
|---|---|
| **Spliced** | CIGAR contains N, or exon-overlapping with no intron overlap |
| **Unspliced** | Single block overlapping intron only |
| **Ambiguous** | Overlaps both exon and intron of same gene (boundary read) |

This is closest to our CIGAR-based approach. The key difference: velocyto uses `pysam` to parse CIGAR directly, while we use `write_blocks_one()` in `scan_core.c`.

### Cell Ranger / STARsolo

Cell Ranger and STARsolo use STAR's splice junction database (SJ.out.tab) plus GTF annotation. Classification is similar to zUMIs but integrated into the mapper. `--include-introns` mode adds intronic counting.

### Comparison Table

| Tool | Method | Spliced detection | Intron detection | Multi-exon reads |
|------|--------|-------------------|-----------------|------------------|
| **zUMIs** | featureCounts annotation | Annotation overlap | Annotation overlap (intron SAF) | Single assignment per read |
| **Velocyto** | CIGAR + GTF | CIGAR N-operator | GTF intron overlap | Per-exon from CIGAR blocks |
| **STARsolo** | STAR SJ + GTF | Splice junction database | GTF intron regions | Integrated in mapper |
| **alevin-fry** | Pseudoalignment | Extended reference index | Pre-mRNA transcripts in index | Spliced vs unspliced isoforms |
| **exonBlocks (ours)** | CIGAR + foverlaps | CIGAR N-operator | GTF exon gaps (intron SAF) | Per-block foverlaps assignment |

### Why our CIGAR-based approach is better for exon-level analysis

1. **Per-block exon assignment**: A spliced read spanning 3 exons produces 3 independent (CB, UMI, exon) rows. zUMIs assigns the whole read to one gene.
2. **No two-pass required**: Single CIGAR parse + single foverlaps. zUMIs requires two featureCounts passes.
3. **Junction precision**: CIGAR N-operations directly identify splice junctions. Annotation-only methods (zUMIs) cannot distinguish a junction read from a read that happens to overlap both exon and intron.
4. **Exon cluster resolution**: Our `make_exon_clusters()` merges overlapping exons, preventing double-counting from reads hitting multiple exons of the same gene.

### Where zUMIs annotation approach is still useful

- **ES/IS/GE/GI tags** provide a quality filter: reads with ES=Ambiguity should be flagged
- **Intron probability score** validates that intronic signal is genuine, not noise
- **Two-pass featureCounts** is a good cross-check against CIGAR-based classification for Smart-seq data where featureCounts tags are available

## Implementation Notes

- `spliced` column: `1` if CIGAR contains N, `0` otherwise
- Added to TSV after `num_blocks`: `CB\tUMI\tblock_start\tblock_end\tblock_seq\tnum_blocks\tspliced`
- R layer uses `spliced` + exon overlap to assign reads to spliced/unspliced matrices
- D operator does NOT produce multiple blocks and does NOT set `spliced=1`
- Multi-exon-mapping reads (same UMI, different exon blocks) require gene-level UMI deduplication before counting