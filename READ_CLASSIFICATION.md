# Read Classification Logic

How each alignment read is processed from BAM through to spliced/unspliced matrix assignment.

## CIGAR-Based Classification

```mermaid
flowchart TD
    READ[BAM Read] --> CIGAR{CIGAR Parsing}
    CIGAR -->|Single block| SBLK[num_blocks = 1]
    CIGAR -->|Multiple blocks N operator| MBLK[num_blocks >= 2]

    SBLK --> OVLP{foverlaps with exon annotation?}
    OVLP -->|Yes| SPLICED1[Spliced exon-only read]
    OVLP -->|No| UNSPLICED1[Unspliced intronic read]

    MBLK --> SPLICED2[Spliced junction read]
    MBLK --> OVLP2{foverlaps each block with exons?}
    OVLP2 -->|Block hits exon E1| ROW1[Row: CB UMI block E1]
    OVLP2 -->|Block hits exon E2| ROW2[Row: CB UMI block E2]
    OVLP2 -->|Block hits intron gap| SKIP[No exon row for intron]

    SPLICED1 --> MAT_S[SPLICED matrix]
    SPLICED2 --> MAT_S
    UNSPLICED1 --> MAT_U[UNSPLICED matrix]

    READ --> TAGS{featureCounts Tags}
    TAGS --> ES[ES Exon Status]
    TAGS --> IS[IS Intron Status]
    TAGS --> GE[GE Gene ID exon]
    TAGS --> GI[GI Gene ID intron]

    ES -->|Assigned| IS_CHK{IS Intron?}
    IS_CHK -->|NoFeatures| EXON_ONLY[Exon-only mature mRNA]
    IS_CHK -->|Assigned| JUNCTION[Junction read GE eq GI?]
    IS_CHK -->|Ambiguity| AMBIG1[Flag ambiguous]

    ES -->|NoFeatures| IS_CHK2{IS Intron?}
    IS_CHK2 -->|Assigned| INTRON_ONLY[Intron-only pre-mRNA]
    IS_CHK2 -->|NoFeatures| INTERG[Discard or UTR]
    IS_CHK2 -->|Ambiguity| AMBIG2[Flag ambiguous]

    ES -->|Ambiguity| DISCARD1[Exclude multi-gene]

    JUNCTION -->|GE eq GI| CONFIRM[Same gene confirmed]
    JUNCTION -->|GE ne GI| CONFLICT[Flag conflicting gene loci]

    EXON_ONLY --> MAT_S
    INTRON_ONLY --> MAT_U
    CONFIRM --> MAT_S
    CONFLICT --> FLAG[Flag or exclude]
```

## Read Classification Mindmap

```mermaid
mindmap
  root((Alignment Read))
    CIGAR Parsing
      Single block
        1 M/X/= operation
        Contiguous alignment
      Multiple blocks
        N splice operator splits
        Each block = exon-aligned segment
        num_blocks >= 2 = spliced
    Tag-based Classification
      ES Exon Status
        Assigned
        Unassigned_NoFeatures
        Unassigned_Ambiguity
      IS Intron Status
        Assigned
        Unassigned_NoFeatures
        Unassigned_Ambiguity
      GE Gene ID exon-assigned
      GI Gene ID intron-assigned
    Splicing Status
      Spliced
        num_blocks >= 2 junction read
        OR single block overlapping exon only
        Assigned to spliced matrix
      Unspliced
        Single block
        No exon overlap
        Falls in intron or pre-mRNA region
        Assigned to unspliced matrix
    GE vs GI
      Same gene
        Standard assignment
        No conflict
      GE ne GI
        Exon from gene A intron from gene B
        Overlapping or antisense loci
        Flag or exclude to avoid misattribution
        ~154 reads in 10M sample
    Exon Overlap via foverlaps
      Each block independently matched
      Spliced read hitting 2 exons yields 2 rows
      Single-block read yields 1 row or 0
      Unspliced read in intron yields 0 rows
    Multi-Exon Mapping
      Read spans N exons
        N foverlaps rows for same CB-UMI
        Gene-level dedup via GE or grouping
        Exon-level counts may double-count UMIs
      Ambiguity
        ES Ambiguity read overlaps exons of different genes
        No GE assigned
        Exclude or distribute fractionally
    Matrix Output
      Spliced matrix
        CB x exon for spliced reads only
      Unspliced matrix
        CB x exon for unspliced reads only
      Gene-level collapse
        Sum exon counts per gene per CB
        Deduplicate UMIs per gene not per exon
    Discard or Flag Categories
      ES NoFeatures IS NoFeatures
        Intergenic or UTR
        No gene assignment
        3.52M reads 35.2%
      ES Ambiguity
        Multi-gene exon overlap
        No confident GE
        41K reads 0.4%
      GE ne GI
        Conflicting gene assignments
        154 reads in 10M
```

## ES x IS Cross-tabulation (chr1, first 10M reads)

| ES ↓ / IS → | Assigned3 | NoFeatures | Ambiguity | Total |
|---|---|---|---|---|
| **Assigned3** | 332,049 | 5,754,106 | — | 6,086,155 |
| **NoFeatures** | 351,316 | 3,521,114 | 16 | 3,872,446 |
| **Ambiguity** | 60 | 41,329 | — | 41,389 |

### Read categories

- **5.75M (57.5%)** ES=Assigned / IS=NoFeatures → exon-only reads (mature mRNA, no intron overlap)
- **3.52M (35.2%)** ES=NoFeatures / IS=NoFeatures → intergenic/UTR/unannotated (outside GTF features)
- **351K (3.5%)** ES=NoFeatures / IS=Assigned → intron-only reads (pre-mRNA/nascent transcripts with GI but no GE)
- **332K (3.3%)** ES=Assigned / IS=Assigned → exon-intron spanning reads (spliced junction reads; both GE and GI populated — key for RNA velocity)
- **41K** ES=Ambiguity / IS=NoFeatures → reads overlapping exons of multiple genes, no intron overlap to disambiguate
- **154 reads** in first 500K of chr19 had GE ≠ GI → reads bridging exons of one gene and introns of another (overlapping/antisense gene loci)

### Biological interpretation

- Exon-only reads = high-confidence mature mRNA expression
- Intron-only reads = transcriptional activity but not spliced; nuclear/pre-mRNA signal
- Exon-intron spanning = junction reads confirming splicing; most informative for isoform resolution and RNA velocity (spliced vs unspliced)
- GE ≠ GI reads = intergenic splicing ambiguity or overlapping gene architectures; should be flagged/excluded to avoid misattribution