#include "exonblocks_shared.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <zlib.h>

// ============================================================================
// Read classification categories
// ============================================================================

typedef enum {
    READ_SPLICED,
    READ_UNSPLICED,
    READ_AMBIGUOUS,
    READ_UNASSIGNED
} ReadClass;

// ============================================================================
// Cell barcode hash (copied from scan_core_stream.c pattern)
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

static unsigned int hash_djb2(const char *str, int len) {
    unsigned int hash = 5381;
    for (int i = 0; i < len; i++) {
        hash = ((hash << 5) + hash) + str[i];
    }
    return hash;
}

static void cb_hash_init(CBHash *h) {
    memset(h, 0, sizeof(CBHash));
    h->next_id = 1;
}

static int cb_hash_get(CBHash *h, const char *key, int key_len) {
    unsigned int idx = hash_djb2(key, key_len) & CB_HASH_MASK;
    CBEntry *e = h->buckets[idx];
    while (e) {
        if (e->key_len == key_len && memcmp(e->key, key, key_len) == 0) {
            return e->cell_id;
        }
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

static void cb_hash_free(CBHash *h) {
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
// Per-category COO triplet storage
// ============================================================================

typedef struct {
    int *cell_idx;
    int *cluster_idx;
    int *counts;
    int capacity;
    int size;
} TripletCOO;

static void triplet_init(TripletCOO *t, int initial_cap) {
    t->capacity = initial_cap;
    t->size = 0;
    t->cell_idx = (int *)malloc(initial_cap * sizeof(int));
    t->cluster_idx = (int *)malloc(initial_cap * sizeof(int));
    t->counts = (int *)malloc(initial_cap * sizeof(int));
}

static void triplet_free(TripletCOO *t) {
    free(t->cell_idx);
    free(t->cluster_idx);
    free(t->counts);
}

static void triplet_add(TripletCOO *t, int cell_id, int cluster_id) {
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

// ============================================================================
// Annotation interval (exon or intron with gene_id)
// ============================================================================

typedef struct {
    int32_t start;
    int32_t end;
    int32_t gene_id;
} AnnotInterval;

// ============================================================================
// Extract aligned blocks as integer arrays (avoids kstring_t round-trip)
// ============================================================================

static int extract_blocks_int(const bam1_t *b, int *blk_start, int *blk_end, int max_blocks) {
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
            // no reference advance
        } else if (op == BAM_CHARD_CLIP) {
            // nothing
        } else if (op == BAM_CPAD) {
            pos += len;
        }
    }
    return n_blocks;
}

// ============================================================================
// Check if CIGAR contains N (intron skip) operation
// ============================================================================

static int has_N_cigar(const bam1_t *b) {
    uint32_t *cigar = bam_get_cigar(b);
    int n_cigar = b->core.n_cigar;
    for (int i = 0; i < n_cigar; i++) {
        if (bam_cigar_op(cigar[i]) == BAM_CREF_SKIP) return 1;
    }
    return 0;
}

// ============================================================================
// Binary search: first interval with end >= query_start (sorted by end)
// ============================================================================

static int find_first_overlap_end(int *ends, int n, int query_start) {
    int lo = 0, hi = n;
    while (lo < hi) {
        int mid = (lo + hi) / 2;
        if (ends[mid] < query_start) {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    return lo;
}

// Binary search: first interval with start > query_end (sorted by start)
// Returns index of first interval PAST the overlap zone
// ============================================================================
// Find gene IDs overlapping a single interval [qstart, qend]
// from a sorted AnnotInterval array (sorted by end)
// Returns gene IDs in out_ids (sorted, deduplicated), returns count
// ============================================================================

static int find_overlapping_gene_ids(AnnotInterval *ivs, int n_ivs,
                                     int *iv_ends, // precomputed sorted ends
                                     int qstart, int qend,
                                     int *out_ids, int max_out) {
    int n_out = 0;
    int start_idx = find_first_overlap_end(iv_ends, n_ivs, qstart);

    for (int i = start_idx; i < n_ivs; i++) {
        if (ivs[i].start > qend) break;
        if (ivs[i].end >= qstart && ivs[i].start <= qend) {
            int gid = ivs[i].gene_id;
            // insert sorted (dedup)
            int inserted = 0;
            for (int j = 0; j < n_out; j++) {
                if (out_ids[j] == gid) { inserted = 1; break; }
                if (out_ids[j] > gid) {
                    if (n_out < max_out) {
                        memmove(out_ids + j + 1, out_ids + j, (n_out - j) * sizeof(int));
                    }
                    out_ids[j] = gid;
                    n_out++;
                    inserted = 1;
                    break;
                }
            }
            if (!inserted && n_out < max_out) {
                out_ids[n_out++] = gid;
            }
            if (n_out >= max_out) break;
        }
    }
    return n_out;
}

// ============================================================================
// In-place sorted intersection of two integer arrays
// Returns result in `a`, returns new length
// ============================================================================

static int intersect_sorted_ints(int *a, int na, int *b, int nb) {
    int i = 0, j = 0, k = 0;
    while (i < na && j < nb) {
        if (a[i] < b[j]) { i++; }
        else if (a[i] > b[j]) { j++; }
        else { a[k++] = a[i]; i++; j++; }
    }
    return k;
}

// ============================================================================
// Check if any interval in sorted annotation overlaps [qstart, qend]
// Returns 1 if overlap found, 0 otherwise
// ============================================================================

static int any_overlap(AnnotInterval *ivs, int n_ivs, int *iv_ends,
                       int qstart, int qend) {
    int start_idx = find_first_overlap_end(iv_ends, n_ivs, qstart);
    for (int i = start_idx; i < n_ivs; i++) {
        if (ivs[i].start > qend) break;
        if (ivs[i].end >= qstart && ivs[i].start <= qend) return 1;
    }
    return 0;
}

// ============================================================================
// Write CellRanger-style output files for one category
// ============================================================================

static void write_category_output(const char *out_dir, const char *category,
                                  TripletCOO *triplets, CBHash *cb_hash,
                                  int *cluster_start, int *cluster_end,
                                  int n_clusters, const char *chr) {
    char path[2048];
    char subdir[1024];
    snprintf(subdir, 1024, "%s/%s", out_dir, category);

    // Create subdirectory (mkdir -p equivalent)
    {
        char cmd[2048];
        snprintf(cmd, 2048, "mkdir -p %s", subdir);
        (void)system(cmd);
    }

    // barcodes.tsv.gz
    snprintf(path, 2048, "%s/barcodes.tsv.gz", subdir);
    gzFile fp_bc = gzopen(path, "wb");
    if (!fp_bc) error("Cannot open %s", path);
    for (int i = 0; i < CB_HASH_SIZE; i++) {
        CBEntry *e = cb_hash->buckets[i];
        while (e) {
            if (e->cell_id >= 1 && e->cell_id <= cb_hash->count) {
                gzprintf(fp_bc, "%s\n", e->key);
            }
            e = e->next;
        }
    }
    gzclose(fp_bc);

    // features.tsv.gz
    snprintf(path, 2048, "%s/features.tsv.gz", subdir);
    gzFile fp_feat = gzopen(path, "wb");
    if (!fp_feat) error("Cannot open %s", path);
    gzprintf(fp_feat, "feature_id\tfeature_name\tfeature_type\n");
    for (int i = 0; i < n_clusters; i++) {
        gzprintf(fp_feat, "cluster_%d\t%s:%d-%d\texon_cluster\n",
                 i + 1, chr, cluster_start[i], cluster_end[i]);
    }
    gzclose(fp_feat);

    // matrix.mtx.gz
    snprintf(path, 2048, "%s/matrix.mtx.gz", subdir);
    gzFile fp_mtx = gzopen(path, "wb");
    if (!fp_mtx) error("Cannot open %s", path);
    gzprintf(fp_mtx, "%%%%MatrixMarket matrix coordinate integer general\n");
    gzprintf(fp_mtx, "%%%%metadata_json: {\"software_version\": \"exonBlocks\", \"category\": \"%s\"}\n", category);
    gzprintf(fp_mtx, "%d %d %d\n", n_clusters, cb_hash->count, triplets->size);

    for (int k = 0; k < triplets->size; k++) {
        gzprintf(fp_mtx, "%d %d %d\n",
                 triplets->cluster_idx[k] + 1,
                 triplets->cell_idx[k] + 1,
                 triplets->counts[k]);
    }
    gzclose(fp_mtx);
}

// ============================================================================
// Build AnnotInterval array from parallel R vectors
// Sort by end (for binary search pattern)
// ============================================================================

static AnnotInterval *build_annot_intervals(int *starts, int *ends, int *gene_ids, int n,
                                            int **out_ends) {
    AnnotInterval *ivs = (AnnotInterval *)malloc(n * sizeof(AnnotInterval));
    *out_ends = (int *)malloc(n * sizeof(int));

    for (int i = 0; i < n; i++) {
        ivs[i].start = starts[i];
        ivs[i].end = ends[i];
        ivs[i].gene_id = gene_ids[i];
    }

    // Sort by end, then by start, then by gene_id
    // Simple insertion sort for small arrays, or qsort
    // Using qsort with a comparison function
    // We need stable sort for the binary search on ends
    for (int i = 1; i < n; i++) {
        AnnotInterval tmp = ivs[i];
        int j = i - 1;
        while (j >= 0 && (ivs[j].end > tmp.end ||
                          (ivs[j].end == tmp.end && ivs[j].start > tmp.start))) {
            ivs[j + 1] = ivs[j];
            j--;
        }
        ivs[j + 1] = tmp;
    }

    for (int i = 0; i < n; i++) {
        (*out_ends)[i] = ivs[i].end;
    }

    return ivs;
}

// ============================================================================
// Main .Call entry: classify reads and build per-category matrices
// ============================================================================

SEXP _exonBlocks_classify_reads(
    SEXP bam_path_,
    SEXP chr_,
    SEXP exon_starts_,
    SEXP exon_ends_,
    SEXP exon_gene_ids_,
    SEXP intron_starts_,
    SEXP intron_ends_,
    SEXP intron_gene_ids_,
    SEXP cluster_starts_,
    SEXP cluster_ends_,
    SEXP tech_,
    SEXP xf_vals_,
    SEXP out_dir_)
{
    const char *bam_path = CHAR(STRING_ELT(bam_path_, 0));
    const char *chr = CHAR(STRING_ELT(chr_, 0));
    const char *tech = CHAR(STRING_ELT(tech_, 0));
    const char *out_dir = CHAR(STRING_ELT(out_dir_, 0));

    int n_exons = LENGTH(exon_starts_);
    int n_introns = LENGTH(intron_starts_);
    int n_clusters = LENGTH(cluster_starts_);

    int *exon_starts = INTEGER(exon_starts_);
    int *exon_ends = INTEGER(exon_ends_);
    int *exon_gene_ids = INTEGER(exon_gene_ids_);
    int *intron_starts = INTEGER(intron_starts_);
    int *intron_ends = INTEGER(intron_ends_);
    int *intron_gene_ids = INTEGER(intron_gene_ids_);
    int *cluster_starts = INTEGER(cluster_starts_);
    int *cluster_ends_ptr = INTEGER(cluster_ends_);

    int nxf = LENGTH(xf_vals_);
    long xf_allow[16];
    for (int i = 0; i < nxf && i < 16; ++i)
        xf_allow[i] = (long)INTEGER(xf_vals_)[i];

    // Build sorted annotation arrays
    int *exon_sort_ends = NULL;
    int *intron_sort_ends = NULL;
    AnnotInterval *exon_ivs = build_annot_intervals(exon_starts, exon_ends, exon_gene_ids, n_exons, &exon_sort_ends);
    AnnotInterval *intron_ivs = NULL;
    if (n_introns > 0) {
        intron_ivs = build_annot_intervals(intron_starts, intron_ends, intron_gene_ids, n_introns, &intron_sort_ends);
    }

    // Sort clusters by end (for binary search)
    // clusters are already expected sorted from R, but let's ensure
    int *cluster_end_sorted = (int *)malloc(n_clusters * sizeof(int));
    int *cluster_start_sorted = (int *)malloc(n_clusters * sizeof(int));
    int *cluster_orig_idx = (int *)malloc(n_clusters * sizeof(int));
    for (int i = 0; i < n_clusters; i++) {
        cluster_start_sorted[i] = cluster_starts[i];
        cluster_end_sorted[i] = cluster_ends_ptr[i];
        cluster_orig_idx[i] = i;
    }
    // Sort clusters by end (insertion sort, stable)
    for (int i = 1; i < n_clusters; i++) {
        int tmp_end = cluster_end_sorted[i];
        int tmp_start = cluster_start_sorted[i];
        int tmp_idx = cluster_orig_idx[i];
        int j = i - 1;
        while (j >= 0 && cluster_end_sorted[j] > tmp_end) {
            cluster_end_sorted[j + 1] = cluster_end_sorted[j];
            cluster_start_sorted[j + 1] = cluster_start_sorted[j];
            cluster_orig_idx[j + 1] = cluster_orig_idx[j];
            j--;
        }
        cluster_end_sorted[j + 1] = tmp_end;
        cluster_start_sorted[j + 1] = tmp_start;
        cluster_orig_idx[j + 1] = tmp_idx;
    }

    // Open BAM
    htsFile *in = sam_open(bam_path, "r");
    if (!in) error("Failed to open BAM: %s", bam_path);

    bam_hdr_t *hdr = sam_hdr_read(in);
    if (!hdr) { sam_close(in); error("Failed to read BAM header"); }

    hts_idx_t *idx = sam_index_load(in, bam_path);
    if (!idx) { bam_hdr_destroy(hdr); sam_close(in); error("No BAM index"); }

    int tid = bam_name2id(hdr, chr);
    if (tid < 0) {
        hts_idx_destroy(idx); bam_hdr_destroy(hdr); sam_close(in);
        error("Unknown contig: %s", chr);
    }
    hts_itr_t *itr = sam_itr_queryi(idx, tid, 0, 0);

    // Initialize data structures
    CBHash cb_hash;
    cb_hash_init(&cb_hash);

    TripletCOO spliced_triplets, unspliced_triplets, ambiguous_triplets;
    triplet_init(&spliced_triplets, 100000);
    triplet_init(&unspliced_triplets, 100000);
    triplet_init(&ambiguous_triplets, 100000);

    long n_spliced = 0, n_unspliced = 0, n_ambiguous = 0, n_unassigned = 0;

    bam1_t *b = bam_init1();
    int blk_start[64], blk_end[64];
    int gene_set[32], block_genes[32], intersect_buf[32];

    while (sam_itr_next(in, itr, b) >= 0) {
        uint16_t flag = b->core.flag;
        if ((flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) != 0) continue;

        // xf filter for 10X
        if (strcmp(tech, "10X") == 0) {
            int okxf = 0, tagok = 0;
            long xf = get_tagi(b, "xf", &tagok);
            if (!tagok) continue;
            for (int i = 0; i < nxf; ++i)
                if (xf == xf_allow[i]) { okxf = 1; break; }
            if (!okxf) continue;
        }

        // Get CB
        const char *CB = get_tagZ(b, "CB");
        if (!CB) continue;
        int cell_id = cb_hash_get(&cb_hash, CB, strlen(CB));

        // Get UMI
        const char *UMI = (strcmp(tech, "10X") == 0) ? get_tagZ(b, "UB") : get_tagZ(b, "RX");
        (void)UMI; // UMI available for future dedup

        // Extract blocks
        int nblk = extract_blocks_int(b, blk_start, blk_end, 64);
        if (nblk <= 0) continue;

        // Classify
        int is_spliced_cigar = has_N_cigar(b);
        ReadClass category;

        if (is_spliced_cigar) {
            // SPLICED path: all blocks must overlap exons of same gene
            int n_gene_set = 0;

            // Get genes for first block
            n_gene_set = find_overlapping_gene_ids(exon_ivs, n_exons, exon_sort_ends,
                                                    blk_start[0], blk_end[0],
                                                    gene_set, 32);

            if (n_gene_set == 0) {
                category = READ_AMBIGUOUS;
            } else {
                for (int bi = 1; bi < nblk && n_gene_set > 0; bi++) {
                    int n_bg = find_overlapping_gene_ids(exon_ivs, n_exons, exon_sort_ends,
                                                          blk_start[bi], blk_end[bi],
                                                          block_genes, 32);
                    if (n_bg == 0) {
                        n_gene_set = 0;
                        break;
                    }
                    // Intersect gene_set with block_genes
                    // Both are sorted already from find_overlapping_gene_ids
                    memcpy(intersect_buf, gene_set, n_gene_set * sizeof(int));
                    n_gene_set = intersect_sorted_ints(intersect_buf, n_gene_set, block_genes, n_bg);
                    memcpy(gene_set, intersect_buf, n_gene_set * sizeof(int));
                }

                if (n_gene_set > 0) {
                    category = READ_SPLICED;
                } else {
                    category = READ_AMBIGUOUS;
                }
            }
        } else {
            // UNSPLICED path: single block or contiguous alignment
            // Check intron overlap first, then exon
            int overlap_intron = 0;
            if (n_introns > 0) {
                // Read spans entire interval from first block start to last block end
                // Actually for no-N cigar, there's typically one block, but handle multi-M + D
                int read_start = blk_start[0];
                int read_end = blk_end[nblk - 1];
                overlap_intron = any_overlap(intron_ivs, n_introns, intron_sort_ends,
                                             read_start, read_end);
            }

            if (overlap_intron) {
                category = READ_UNSPLICED;
            } else {
                int overlap_exon = 0;
                for (int bi = 0; bi < nblk; bi++) {
                    if (any_overlap(exon_ivs, n_exons, exon_sort_ends,
                                    blk_start[bi], blk_end[bi])) {
                        overlap_exon = 1;
                        break;
                    }
                }
                if (overlap_exon) {
                    category = READ_AMBIGUOUS;
                } else {
                    category = READ_UNASSIGNED;
                }
            }
        }

        // Map blocks to clusters and add to category-specific COO
        if (category != READ_UNASSIGNED) {
            TripletCOO *coo;
            switch (category) {
                case READ_SPLICED:   coo = &spliced_triplets;   n_spliced++;   break;
                case READ_UNSPLICED: coo = &unspliced_triplets;  n_unspliced++; break;
                case READ_AMBIGUOUS: coo = &ambiguous_triplets;   n_ambiguous++; break;
                default: coo = NULL; break;
            }

            if (coo) {
                for (int bi = 0; bi < nblk; bi++) {
                    int start_idx = find_first_overlap_end(cluster_end_sorted, n_clusters, blk_start[bi]);
                    for (int ci = start_idx; ci < n_clusters; ci++) {
                        if (cluster_start_sorted[ci] > blk_end[bi]) break;
                        if (cluster_end_sorted[ci] >= blk_start[bi] && cluster_start_sorted[ci] <= blk_end[bi]) {
                            triplet_add(coo, cell_id, cluster_orig_idx[ci]);
                        }
                    }
                }
            }
        } else {
            n_unassigned++;
        }
    }

    // Cleanup BAM
    bam_destroy1(b);
    hts_itr_destroy(itr);
    hts_idx_destroy(idx);
    bam_hdr_destroy(hdr);
    sam_close(in);

    // Free annotation arrays
    free(exon_ivs);
    free(exon_sort_ends);
    if (intron_ivs) { free(intron_ivs); free(intron_sort_ends); }
    free(cluster_end_sorted);
    free(cluster_start_sorted);
    free(cluster_orig_idx);

    // Write output files for each category
    write_category_output(out_dir, "spliced", &spliced_triplets, &cb_hash,
                          cluster_starts, cluster_ends_ptr, n_clusters, chr);
    write_category_output(out_dir, "unspliced", &unspliced_triplets, &cb_hash,
                          cluster_starts, cluster_ends_ptr, n_clusters, chr);
    write_category_output(out_dir, "ambiguous", &ambiguous_triplets, &cb_hash,
                          cluster_starts, cluster_ends_ptr, n_clusters, chr);

    // Free COO triplets
    triplet_free(&spliced_triplets);
    triplet_free(&unspliced_triplets);
    triplet_free(&ambiguous_triplets);

    // Free CB hash
    cb_hash_free(&cb_hash);

    // Return summary
    SEXP result = PROTECT(allocVector(VECSXP, 7));
    SEXP names = PROTECT(allocVector(STRSXP, 7));

    SET_VECTOR_ELT(result, 0, ScalarInteger(cb_hash.count));
    SET_VECTOR_ELT(result, 1, ScalarInteger(n_clusters));
    SET_VECTOR_ELT(result, 2, ScalarInteger((int)n_spliced));
    SET_VECTOR_ELT(result, 3, ScalarInteger((int)n_unspliced));
    SET_VECTOR_ELT(result, 4, ScalarInteger((int)n_ambiguous));
    SET_VECTOR_ELT(result, 5, ScalarInteger((int)n_unassigned));
    SET_VECTOR_ELT(result, 6, mkString(out_dir));

    SET_STRING_ELT(names, 0, mkChar("n_cells"));
    SET_STRING_ELT(names, 1, mkChar("n_features"));
    SET_STRING_ELT(names, 2, mkChar("n_spliced"));
    SET_STRING_ELT(names, 3, mkChar("n_unspliced"));
    SET_STRING_ELT(names, 4, mkChar("n_ambiguous"));
    SET_STRING_ELT(names, 5, mkChar("n_unassigned"));
    SET_STRING_ELT(names, 6, mkChar("out_dir"));

    setAttrib(result, R_NamesSymbol, names);
    UNPROTECT(2);
    return result;
}