# Velocyto Read Classification Logic — Detailed Notes

Source: [velocyto-team/velocyto.py](https://github.com/velocyto-team/velocyto.py)
Paper: La Manno et al., *Nature* 560:494–498, 2018

**Note:** This repo is no longer maintained. The authors recommend STARsolo instead.

---

## Architecture Overview

Velocyto classifies reads at **two levels**:

1. **CIGAR parsing** → splits each read into aligned reference segments (blocks)
2. **GTF annotation overlap** → each segment is checked against exon/intron features per transcript model

The classification is done per-UMI (per-molecule), not per-read. Multiple reads with the same (CB, UMI) are merged into a `Molitem`, and the union of their segment-mapping records determines the final category.

---

## Step 1: CIGAR Parsing → Segments

Velocyto's `parse_cigar_tuple()` (in `counter.py`) produces:

```python
segments, ref_skipped, clip5, clip3 = parse_cigar_tuple(read.cigartuples, pos)
```

### Segment extraction rules

| CIGAR Op | Action | Effect on segments |
|----------|--------|--------------------|
| `M`/`=`/`X` (0,7,8) | Emit segment `[p, p+len-1]`, advance pos | New aligned block |
| `N` (3) | Advance pos, set `ref_skipped=True` | **Splice junction detected** |
| `D` (2) | If `len <= PATCH_INDELS` and flanked by M ops → merge into previous segment; else advance pos | Treated as contiguous (NOT a splice) |
| `I` (1) | If `len <= PATCH_INDELS` and flanked by M ops → merge adjacent segments; else ignored | Read-only, no ref advance |
| `S` (4) | Track as `clip5`/`clip3` | Not part of alignment |
| `H` (5) | Log warning | Not part of alignment |

**Key difference from our implementation:**

- Velocyto **merges small indels**: if a `D` or `I` is ≤ `PATCH_INDELS` (default 5bp) and flanked by M operations, it bridges the adjacent segments into one. Our `write_blocks_one()` does NOT do this — we treat D as always advancing position but never splitting blocks.
- `ref_skipped` flag is the equivalent of our `num_blocks >= 2` spliced indicator.

---

## Step 2: Transcript Model Construction from GTF

Velocyto reads the GTF and builds `TranscriptModel` objects:

1. Each `exon` feature in the GTF → `Feature(kind='e', start, end, exon_number)`
2. **Introns** are derived as gaps between consecutive exons **within the same transcript**: `Feature(kind='i', start, end, is_validated=False)`
3. Introns longer than `LONGEST_INTRON_ALLOWED` (default 1Mbp) are split ("chopped") to avoid masking entire gene loci

This means velocyto uses **transcript-level annotation** (not gene-level), which is more precise for isoform resolution.

---

## Step 3: Intron Validation (mark_up_introns)

After building transcript models, velocyto scans the BAM a second time to **validate introns**:

- For each **unspliced read** (`ref_skipped == False`, i.e., single contiguous segment):
  - If the read overlaps an intron feature, the read's segment **spans** the intron boundary
  - This means the intron has read evidence that confirms it's real (not a genomic annotation artifact)
  - The intron's `is_validated` flag is set to `True`

**Validated introns** are those confirmed by at least one spanning read. This is a critical concept unique to velocyto.

---

## Step 4: Read Classification — The Logic Classes

Velocyto provides **5 logic modes** (strictness levels), all sharing the same underlying classification but with different rules for unvalidated introns and singleton reads:

### Per-segment classification

For each segment of a read, velocyto creates a `SegmentMatch` object with:

```python
class SegmentMatch:
    segment: Tuple[int, int]      # (start, end) of the aligned block
    feature: Feature               # The exon/intron feature it overlaps
    maps_to_intron: bool           # True if feature.kind == ord('i')
    maps_to_exon: bool             # True if feature.kind == ord('e')
    is_spliced: bool               # True if CIGAR contains N (ref_skipped flag)
```

### Per-molecule flags (computed from all segment mappings)

For a given (CB, UMI) molitem, velocyto evaluates all compatible transcript models and computes these flags:

| Flag | Meaning |
|------|---------|
| `has_introns` | At least one segment maps to an intron |
| `has_exons` | At least one segment maps to an exon |
| `has_exseg_with_spliced_flag` | At least one exon-overlapping segment has `is_spliced=True` (CIGAR contains N) |
| `has_validated_intron` | At least one intron segment is confirmed by spanning reads |
| `has_exin_intron_span` | A segment spanning the exon-intron boundary (segment overlaps both the exon AND the adjacent intron) |
| `has_non3prime` | At least one exon segment that is NOT the last 3' exon |

### Combining flags into transcript-level categories

For each compatible transcript model, velocyto derives:

| Composite Flag | Condition | Meaning |
|---------------|-----------|---------|
| `has_onlyexo_model` | `has_exons and not has_introns` | All segments hit exons only |
| `has_onlyintron_model` | `has_introns and not has_exons` | All segments hit introns only |
| `has_onlyintron_and_valid_model` | `has_validated_intron and not has_exons` | Introns only, and at least one is validated |
| `has_mixed_model` | `has_exons and has_introns and (not has_exin_intron_span)` | Mix of exons and introns, no exon-intron spanning boundary |
| `has_valid_mixed_model` | `has_exons and has_introns and has_validated_intron and (not has_exin_intron_span)` | Mix, with validated intron, no boundary spanning |
| `has_invalid_mixed_model` | `has_exons and has_introns and (not has_validated_intron) and (not has_exin_intron_span)` | Mix, no validated intron, no boundary spanning |
| `has_only_span_exin_model` | `not has_exin_intron_span` across ALL models | Every compatible model has exon-intron boundary spanning |

Wait — `has_only_span_exin_model` is initialized to `1` and set to `0` if any model has `has_exin_intron_span`. So `has_only_span_exin_model=True` means **no model has an exon-intron spanning boundary**.

Actually re-reading more carefully: `has_exin_intron_span` means a segment **spans the boundary** between an exon and its adjacent intron. This happens when a read's aligned segment partially overlaps an intron such that the segment also overlaps the downstream (or upstream) exon. This is the key signal for **unspliced pre-mRNA** — the read is physically present in both the exon and the intron at the splice junction.

### Final classification across all transcript models

After evaluating all transcript models for a molitem, velocyto checks:

1. **Multi-gene?** If segments map to >1 gene → **discard** (ambiguous across genes)
2. **Single gene, multiple transcript models** → combine flags across models:

| Condition | Category | Notes |
|-----------|-----------|-------|
| `has_onlyexo_model` and no `has_onlyintron_model` and no `has_mixed_model` | **SPLICED** | Exon-only, no intron overlap |
| `has_only_span_exin_model` | **UNSPLICED** | All models show exon-intron boundary spanning |
| `has_onlyintron_and_valid_model` and no mixed/exon | **UNSPLICED** | Validated intron only (singleton or non-singleton depends on logic) |
| `has_onlyintron_model` (no validated, no mixed, no exon) | **UNSPLICED** or **DISCARD** | Depends on logic strictness |
| `has_invalid_mixed_model` only | **UNSPLICED** or **DISCARD** | Unvalidated intron+exon, depends on logic |
| `has_valid_mixed_model` only | **UNSPLICED** | Validated intron+exon mix |
| `has_onlyintron_model` and `has_onlyexo_model` | **AMBIGUOUS** | Different transcript models disagree |
| Any combination with `has_mixed_model` | **AMBIGUOUS** or **UNSPLICED** | Depends on logic |

---

## Step 5: The 5 Logic Modes

| Logic | Singleton in non-validated intron | Non-singleton in non-validated intron | Singleton in validated intron | Non-singleton in validated intron | Invalid mixed | Validated mixed | Boundary spanning |
|-------|-----------------------------------|--------------------------------------|-------------------------------|------------------------------------|---------------|-----------------|-------------------|
| **Permissive10X** | COUNT unsPLICED | COUNT unsPLICED | COUNT unsPLICED | COUNT unsPLICED | COUNT unsPLICED | COUNT unsPLICED | COUNT unsPLICED |
| **Intermediate10X** | DISCARD | COUNT unsPLICED | COUNT unsPLICED | COUNT unsPLICED | DISCARD | COUNT unsPLICED | COUNT unsPLICED |
| **ValidatedIntrons10X** | DISCARD | DISCARD | COUNT unsPLICED | COUNT unsPLICED | DISCARD | COUNT unsPLICED | COUNT unsPLICED |
| **Stricter10X** | DISCARD | DISCARD | DISCARD | COUNT unsPLICED | DISCARD | COUNT unsPLICED | COUNT unsPLICED |
| **ObservedSpanning10X** | DISCARD | DISCARD | DISCARD | DISCARD | DISCARD | COUNT unsPLICED | COUNT unsPLICED |

### What each logic is for

- **Permissive**: Count everything. For 10X where UMI dedup provides noise filtering.
- **Intermediate**: Require ≥2 UMIs for non-validated introns. Good balance for 10X.
- **ValidatedIntrons**: Only count introns with read-through validation. Recommended for 10X.
- **Stricter**: Only count validated introns with >1 UMI. Most conservative 10X mode.
- **ObservedSpanning**: Only count introns where a segment physically spans the exon-intron boundary. The most stringent — only pre-mRNA molecules with direct splice-site evidence.

---

## Key Differences from exonBlocks

| Aspect | Velocyto | exonBlocks |
|--------|----------|-----------|
| **CIGAR parsing** | Merges small indels (D/I ≤ 5bp) into adjacent segments | Treats D as continuous block advancement, never splits |
| **Annotation** | Transcript-level (intron features derived from GTF exon gaps) | Gene-level exon clusters (`make_exon_clusters`) or raw exon TSV |
| **Intron validation** | Two-pass: scan BAM to mark introns with spanning evidence | No intron validation (not yet) |
| **Classification** | Per-molecule (UMI dedup first, then classify) | Per-read (each read independently via CIGAR + foverlaps) |
| **Output** | 3 matrices: spliced, unspliced, ambiguous | Single TSV per scan, then R-level classification |
| **Spliced detection** | `ref_skipped` flag from CIGAR N-operator | `num_blocks >= 2` (equivalent) |
| **Splice boundary detection** | `has_exin_intron_span` — segment physically bridges exon-intron boundary | Not explicitly detected (we check block-level exon overlap only) |
| **Stratification** | 5 strictness levels (logic classes) | Single mode (could add logic classes) |
| **Gene ambiguity** | Discards multi-gene molitems | GE ≠ GI flagged, not discarded |

### What exonBlocks should adopt from velocyto

1. **Intron validation** — two-pass scan where unspliced reads that span an intron boundary validate that intron. Prevents counting annotation artifacts as pre-mRNA.

2. **Exon-intron boundary spanning** (`has_exin_intron_span`) — detects when a segment overlaps both an exon and its adjacent intron. This is the strongest signal for unspliced pre-mRNA.

3. **Small indel patching** (`PATCH_INDELS`) — D/I operations ≤5bp bridging two M segments should be merged into one block. Our current implementation splits on D, which means a read like `50M3D30M` becomes one block, but `50M150N26M` becomes two blocks. This is correct, but we should document it.

4. **Logic classes** — offer Permissive/Intermediate/Validated/Strict/Observed as user-selectable modes rather than hardcoding one classification.

5. **Per-molecule aggregation** — classify after UMI dedup, not per-read. A UMI with 3 reads where 2 hit exons and 1 hits intron should be classified based on the combined evidence, not 3 independent classifications.