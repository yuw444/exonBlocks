#ifndef EXONBLOCKS_H
#define EXONBLOCKS_H

#include <R.h>
#include <Rinternals.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/kstring.h>
#include <string.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <errno.h>

// ============================================================================
// BAM tag helpers
// ============================================================================

static inline char base_at(const uint8_t *s, int qpos) {
    return seq_nt16_str[bam_seqi(s, qpos)];
}

int write_blocks_one(const bam1_t *b, kstring_t *ks_start,
                    kstring_t *ks_end, kstring_t *ks_seq);

static inline const char *get_tagZ(const bam1_t *b, const char *tag) {
    uint8_t *p = bam_aux_get(b, tag);
    if (!p) return NULL;
    return bam_aux2Z(p);
}

static inline long get_tagi(const bam1_t *b, const char *tag, int *ok) {
    *ok = 0;
    uint8_t *p = bam_aux_get(b, tag);
    if (!p) return 0;
    *ok = 1;
    return bam_aux2i(p);
}

// ============================================================================
// CIGAR block extraction (1-based inclusive coordinates)
// ============================================================================

static inline int extract_blocks_int(const bam1_t *b, int *blk_start, int *blk_end, int max_blocks) {
    uint32_t *cigar = bam_get_cigar(b);
    int n_cigar = b->core.n_cigar;
    int n_blocks = 0;
    int pos = b->core.pos;
    for (int i = 0; i < n_cigar && n_blocks < max_blocks; i++) {
        int op = bam_cigar_op(cigar[i]);
        int len = bam_cigar_oplen(cigar[i]);
        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
            blk_start[n_blocks] = pos + 1;
            blk_end[n_blocks] = pos + len;
            n_blocks++;
            pos += len;
        } else if (op == BAM_CDEL || op == BAM_CREF_SKIP) {
            pos += len;
        } else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) {
        } else if (op == BAM_CHARD_CLIP) {
        } else if (op == BAM_CPAD) {
            pos += len;
        }
    }
    return n_blocks;
}

static inline int has_N_cigar(const bam1_t *b) {
    uint32_t *cigar = bam_get_cigar(b);
    int n_cigar = b->core.n_cigar;
    for (int i = 0; i < n_cigar; i++)
        if (bam_cigar_op(cigar[i]) == BAM_CREF_SKIP) return 1;
    return 0;
}

// ============================================================================
// Cell barcode hash (DJB2, chaining)
// ============================================================================

#define CB_HASH_SIZE 65536
#define CB_HASH_MASK (CB_HASH_SIZE - 1)

typedef struct CBEntry {
    char *key;
    int key_len;
    int cell_id;
    struct CBEntry *next;
} CBEntry;

typedef struct {
    CBEntry *buckets[CB_HASH_SIZE];
    int count;
    int next_id;
} CBHash;

static inline unsigned int hash_djb2(const char *str, int len) {
    unsigned int hash = 5381;
    for (int i = 0; i < len; i++)
        hash = ((hash << 5) + hash) + str[i];
    return hash;
}

static inline void cb_hash_init(CBHash *h) {
    memset(h, 0, sizeof(CBHash));
    h->next_id = 1;
}

static inline int cb_hash_get(CBHash *h, const char *key, int key_len) {
    unsigned int idx = hash_djb2(key, key_len) & CB_HASH_MASK;
    CBEntry *e = h->buckets[idx];
    while (e) {
        if (e->key_len == key_len && memcmp(e->key, key, key_len) == 0)
            return e->cell_id;
        e = e->next;
    }
    e = (CBEntry *)malloc(sizeof(CBEntry));
    e->key = (char *)malloc(key_len + 1);
    memcpy(e->key, key, key_len);
    e->key[key_len] = '\0';
    e->key_len = key_len;
    e->cell_id = h->next_id++;
    e->next = h->buckets[idx];
    h->buckets[idx] = e;
    h->count++;
    return e->cell_id;
}

static inline void cb_hash_free(CBHash *h) {
    for (int i = 0; i < CB_HASH_SIZE; i++) {
        CBEntry *e = h->buckets[i];
        while (e) {
            CBEntry *next = e->next;
            free(e->key);
            free(e);
            e = next;
        }
    }
}

// ============================================================================
// COO triplet storage (cell_idx, cluster_idx, counts)
// ============================================================================

typedef struct {
    int *cell_idx;
    int *cluster_idx;
    int *counts;
    int capacity;
    int size;
} TripletCOO;

static inline void triplet_init(TripletCOO *t, int initial_cap) {
    t->capacity = initial_cap;
    t->size = 0;
    t->cell_idx = (int *)malloc(initial_cap * sizeof(int));
    t->cluster_idx = (int *)malloc(initial_cap * sizeof(int));
    t->counts = (int *)malloc(initial_cap * sizeof(int));
}

static inline void triplet_free(TripletCOO *t) {
    free(t->cell_idx);
    free(t->cluster_idx);
    free(t->counts);
}

static inline void triplet_add(TripletCOO *t, int cell_id, int cluster_id) {
    if (t->size > 0) {
        int last = t->size - 1;
        if (t->cell_idx[last] == cell_id && t->cluster_idx[last] == cluster_id) {
            t->counts[last]++;
            return;
        }
    }
    for (int i = t->size - 2; i >= 0 && i >= t->size - 10; i--) {
        if (t->cell_idx[i] == cell_id && t->cluster_idx[i] == cluster_id) {
            t->counts[i]++;
            return;
        }
    }
    if (t->size >= t->capacity) {
        t->capacity *= 2;
        t->cell_idx = (int *)realloc(t->cell_idx, t->capacity * sizeof(int));
        t->cluster_idx = (int *)realloc(t->cluster_idx, t->capacity * sizeof(int));
        t->counts = (int *)realloc(t->counts, t->capacity * sizeof(int));
    }
    t->cell_idx[t->size] = cell_id;
    t->cluster_idx[t->size] = cluster_id;
    t->counts[t->size] = 1;
    t->size++;
}

static inline void triplet_add_count(TripletCOO *t, int cell_id, int cluster_id, int count) {
    if (t->size > 0) {
        int last = t->size - 1;
        if (t->cell_idx[last] == cell_id && t->cluster_idx[last] == cluster_id) {
            t->counts[last] += count;
            return;
        }
    }
    for (int i = t->size - 2; i >= 0 && i >= t->size - 10; i--) {
        if (t->cell_idx[i] == cell_id && t->cluster_idx[i] == cluster_id) {
            t->counts[i] += count;
            return;
        }
    }
    if (t->size >= t->capacity) {
        t->capacity *= 2;
        t->cell_idx = (int *)realloc(t->cell_idx, t->capacity * sizeof(int));
        t->cluster_idx = (int *)realloc(t->cluster_idx, t->capacity * sizeof(int));
        t->counts = (int *)realloc(t->counts, t->capacity * sizeof(int));
    }
    t->cell_idx[t->size] = cell_id;
    t->cluster_idx[t->size] = cluster_id;
    t->counts[t->size] = count;
    t->size++;
}

// ============================================================================
// Binary search: first index where ends[i] >= query_start
// ============================================================================

static inline int find_first_overlap_end(int *ends, int n, int query_start) {
    int lo = 0, hi = n;
    while (lo < hi) {
        int mid = (lo + hi) / 2;
        if (ends[mid] < query_start) lo = mid + 1;
        else hi = mid;
    }
    return lo;
}

// ============================================================================
// Protocol-agnostic read pre-filtering
// ============================================================================

typedef enum {
    TECH_10X,
    TECH_SMARTSEQ3XPRESS,
    TECH_SMARTSEQ
} TechProtocol;

static inline TechProtocol parse_tech(const char *tech) {
    if (strcmp(tech, "10X") == 0) return TECH_10X;
    if (strcmp(tech, "Smart-seq3Xpress") == 0) return TECH_SMARTSEQ3XPRESS;
    return TECH_SMARTSEQ;
}

static inline int read_passes_filter(const bam1_t *b, TechProtocol tech,
                                      long *xf_allow, int nxf) {
    uint16_t flag = b->core.flag;
    if (flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) return 0;

    if (tech == TECH_10X) {
        int tagok = 0;
        long xf = get_tagi(b, "xf", &tagok);
        if (!tagok) return 0;
        for (int i = 0; i < nxf; i++)
            if (xf == xf_allow[i]) return 1;
        return 0;
    }

    return 1;
}

static inline const char *read_cell_id(const bam1_t *b, TechProtocol tech) {
    if (tech == TECH_10X) {
        const char *cb = get_tagZ(b, "CB");
        if (cb) return cb;
        return NULL;
    }
    if (tech == TECH_SMARTSEQ3XPRESS) {
        const char *bx = get_tagZ(b, "BX");
        if (bx) return bx;
        const char *cb = get_tagZ(b, "CB");
        if (cb) return cb;
        return NULL;
    }
    return bam_get_qname(b);
}

// ============================================================================
// Genome annotation indices (per-contig, sorted for binary search)
// ============================================================================

typedef struct {
    size_t n;
    int32_t *start;
    int32_t *end;
    int32_t *gene_id;
    int32_t *exon_cluster_id;
} ExonClusterIndex;

typedef struct {
    size_t n;
    int32_t *start;
    int32_t *end;
    int32_t *gene_id;
} GeneIndex;

typedef struct {
    size_t n_contigs;
    ExonClusterIndex *exon_index;
    GeneIndex *gene_index;
} GenomeAnnotation;

// ============================================================================
// HUGE read classification flags (bit field)
// ============================================================================
// H: has N in CIGAR (spliced read)
// U: unique — all gene-mapped blocks belong to same gene
// G: all blocks within gene region (no intergenic blocks)
// E: all blocks within exon cluster (every block hits an exon)
//
// Classification:
//   1111 -> SPLICED     (junction read with confident gene + exon assignment)
//   1011 -> CHIMERIC    (junction read, ambiguous gene, but all exonic)
//   _110 -> UNSPLICED   (all in gene region, not all exonic — intronic signal)
//   Other -> UNASSIGNED

#define HUGE_H  0x08  /* has N in CIGAR (spliced read) */
#define HUGE_U  0x04  /* unique gene assignment across all blocks */
#define HUGE_G  0x02  /* all blocks within gene region */
#define HUGE_E  0x01  /* all blocks within exon cluster */

// ============================================================================
// mkdir -p helper (no shell injection)
// ============================================================================

static inline int mkdir_p(const char *path) {
    char tmp[2048];
    snprintf(tmp, sizeof(tmp), "%s", path);
    size_t len = strlen(tmp);
    if (len > 0 && tmp[len - 1] == '/') tmp[len - 1] = '\0';
    for (char *p = tmp + 1; *p; p++) {
        if (*p == '/') {
            *p = '\0';
            if (mkdir(tmp, 0755) != 0 && errno != EEXIST) return -1;
            *p = '/';
        }
    }
    if (mkdir(tmp, 0755) != 0 && errno != EEXIST) return -1;
    return 0;
}

#endif // EXONBLOCKS_H
