#include "exonblocks_shared.h"
#include <Rinternals.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/kstring.h>
#include <htslib/bgzf.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <zlib.h>

// ============================================================================
// Hash table for cell barcodes (simple DJB2 hash)
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
    // Not found, add new
    e = (CBEntry*)malloc(sizeof(CBEntry));
    e->key = (char*)malloc(key_len + 1);
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
// COO triplet storage (for in-memory aggregation before writing)
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
    t->cell_idx = (int*)malloc(initial_cap * sizeof(int));
    t->cluster_idx = (int*)malloc(initial_cap * sizeof(int));
    t->counts = (int*)malloc(initial_cap * sizeof(int));
}

static void triplet_free(TripletCOO *t) {
    free(t->cell_idx);
    free(t->cluster_idx);
    free(t->counts);
}

static void triplet_add(TripletCOO *t, int cell_id, int cluster_id) {
    // Check last entry first (common case)
    if (t->size > 0) {
        int last = t->size - 1;
        if (t->cell_idx[last] == cell_id && t->cluster_idx[last] == cluster_id) {
            t->counts[last]++;
            return;
        }
    }
    // Linear search backward (for nearby clusters from same read)
    for (int i = t->size - 2; i >= 0 && i >= t->size - 10; i--) {
        if (t->cell_idx[i] == cell_id && t->cluster_idx[i] == cluster_id) {
            t->counts[i]++;
            return;
        }
    }
    // Add new
    if (t->size >= t->capacity) {
        t->capacity *= 2;
        t->cell_idx = (int*)realloc(t->cell_idx, t->capacity * sizeof(int));
        t->cluster_idx = (int*)realloc(t->cluster_idx, t->capacity * sizeof(int));
        t->counts = (int*)realloc(t->counts, t->capacity * sizeof(int));
    }
    t->cell_idx[t->size] = cell_id;
    t->cluster_idx[t->size] = cluster_id;
    t->counts[t->size] = 1;
    t->size++;
}

// ============================================================================
// Binary search for cluster overlap
// ============================================================================

static int find_first_overlap(int *cluster_end, int n_clusters, int block_start) {
    int lo = 0, hi = n_clusters;
    while (lo < hi) {
        int mid = (lo + hi) / 2;
        if (cluster_end[mid] < block_start) {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    return lo;
}

// ============================================================================
// Parse comma/semicolon-separated coordinates
// ============================================================================

static int parse_coords(const char *str, int *coords, int max_coords) {
    int count = 0;
    const char *p = str;
    char *end;
    while (*p && count < max_coords) {
        while (*p == ',' || *p == ';' || *p == ' ') p++;
        if (!*p) break;
        coords[count++] = strtol(p, &end, 10);
        p = end;
    }
    return count;
}

// ============================================================================
// Write gzipped file helper
// ============================================================================

static void write_gzopen(const char *path, const char *data, size_t len) {
    gzFile fp = gzopen(path, "wb");
    if (!fp) error("Cannot open %s for writing", path);
    gzwrite(fp, data, len);
    gzclose(fp);
}

// ============================================================================
// Main streaming function: BAM -> CellRanger-style output
// Outputs: barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz
// ============================================================================

SEXP _exonBlocks_bam_to_cellranger(
    SEXP bam_,
    SEXP cluster_chr_,
    SEXP cluster_start_,
    SEXP cluster_end_,
    SEXP cluster_id_,
    SEXP tech_,
    SEXP xf_vals_,
    SEXP out_dir_)
{
    const char *bam_path = CHAR(STRING_ELT(bam_, 0));
    const char *out_dir = CHAR(STRING_ELT(out_dir_, 0));
    const char *tech = CHAR(STRING_ELT(tech_, 0));

    // Cluster data
    int n_clusters = LENGTH(cluster_start_);
    int *cluster_start = INTEGER(cluster_start_);
    int *cluster_end = INTEGER(cluster_end_);
    const char *cluster_chr = NULL;
    if (LENGTH(cluster_chr_) > 0) {
        cluster_chr = CHAR(STRING_ELT(cluster_chr_, 0));
    }

    // xf filter
    int nxf = LENGTH(xf_vals_);
    long xf_allow[16];
    for (int i = 0; i < nxf && i < 16; ++i)
        xf_allow[i] = (long)INTEGER(xf_vals_)[i];

    // Open BAM
    htsFile *in = sam_open(bam_path, "r");
    if (!in) error("Failed to open BAM: %s", bam_path);

    bam_hdr_t *hdr = sam_hdr_read(in);
    if (!hdr) { sam_close(in); error("Failed to read BAM header"); }

    hts_idx_t *idx = sam_index_load(in, bam_path);
    if (!idx) { bam_hdr_destroy(hdr); sam_close(in); error("No BAM index"); }

    // Initialize
    CBHash cb_hash;
    cb_hash_init(&cb_hash);

    TripletCOO triplets;
    triplet_init(&triplets, 500000);

    bam1_t *b = bam_init1();
    kstring_t ks_start = {0, 0, NULL};
    kstring_t ks_end = {0, 0, NULL};
    int coords[100];

    // Process all reads
    hts_itr_t *itr = sam_itr_queryi(idx, 0, 0, 0);

    while (sam_itr_next(in, itr, b) >= 0) {
        uint16_t flag = b->core.flag;
        if ((flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) != 0) continue;

        // xf filter for 10X
        if (strcmp(tech, "10X") == 0) {
            uint8_t *xfp = bam_aux_get(b, "xf");
            if (!xfp) continue;
            long xf = bam_aux2i(xfp);
            int ok = 0;
            for (int i = 0; i < nxf; ++i) if (xf == xf_allow[i]) { ok = 1; break; }
            if (!ok) continue;
        }

        // Get CB
        uint8_t *cbp = bam_aux_get(b, "CB");
        if (!cbp) continue;
        const char *CB = bam_aux2Z(cbp);
        int cell_id = cb_hash_get(&cb_hash, CB, strlen(CB));

        // Get blocks
        ks_start.l = ks_end.l = 0;
        int nblk = write_blocks_one(b, &ks_start, &ks_end, NULL);
        if (nblk <= 0) continue;

        // Parse coordinates
        int n_start = parse_coords(ks_start.s, coords, 100);
        int n_end = parse_coords(ks_end.s, coords + nblk, 100);

        // Map to clusters
        for (int bi = 0; bi < nblk && bi < n_start; bi++) {
            int blk_start = coords[bi];
            int blk_end = coords[nblk + bi];
            int start_idx = find_first_overlap(cluster_end, n_clusters, blk_start);
            for (int ci = start_idx; ci < n_clusters; ci++) {
                if (cluster_start[ci] > blk_end) break;
                if (cluster_end[ci] >= blk_start && cluster_start[ci] <= blk_end) {
                    triplet_add(&triplets, cell_id, ci);
                }
            }
        }
    }

    // Cleanup BAM
    ks_free(&ks_start);
    ks_free(&ks_end);
    bam_destroy1(b);
    hts_itr_destroy(itr);
    hts_idx_destroy(idx);
    bam_hdr_destroy(hdr);
    sam_close(in);

    int n_cells = cb_hash.count;
    int nnz = triplets.size;

    // ========================================================================
    // Write output files
    // ========================================================================

    char path[1024];

    // 1. barcodes.tsv.gz
    snprintf(path, 1024, "%s/barcodes.tsv.gz", out_dir);
    gzFile fp_barcodes = gzopen(path, "wb");
    if (!fp_barcodes) error("Cannot open %s", path);
    for (int i = 0; i < CB_HASH_SIZE; i++) {
        CBEntry *e = cb_hash.buckets[i];
        while (e) {
            if (e->cell_id >= 1 && e->cell_id <= n_cells) {
                gzprintf(fp_barcodes, "%s\n", e->key);
            }
            e = e->next;
        }
    }
    gzclose(fp_barcodes);

    // 2. features.tsv.gz (exon clusters)
    snprintf(path, 1024, "%s/features.tsv.gz", out_dir);
    gzFile fp_features = gzopen(path, "wb");
    if (!fp_features) error("Cannot open %s", path);
    gzprintf(fp_features, "feature_id\tfeature_name\tfeature_type\n");
    for (int i = 0; i < n_clusters; i++) {
        if (cluster_chr) {
            gzprintf(fp_features, "cluster_%d\t%s:%d-%d\texon_cluster\n",
                     i + 1, cluster_chr, cluster_start[i], cluster_end[i]);
        } else {
            gzprintf(fp_features, "cluster_%d\t%d-%d\texon_cluster\n",
                     i + 1, cluster_start[i], cluster_end[i]);
        }
    }
    gzclose(fp_features);

    // 3. matrix.mtx.gz (COO format)
    snprintf(path, 1024, "%s/matrix.mtx.gz", out_dir);
    gzFile fp_matrix = gzopen(path, "wb");
    if (!fp_matrix) error("Cannot open %s", path);

    // Matrix Market header
    gzprintf(fp_matrix, "%%%%MatrixMarket matrix coordinate integer general\n");
    gzprintf(fp_matrix, "%%%%metadata_json: {\"software_version\": \"exonBlocks\"}\n");
    gzprintf(fp_matrix, "%d %d %d\n", n_clusters, n_cells, nnz);

    // Write COO entries (feature_id cell_id count)
    for (int k = 0; k < nnz; k++) {
        gzprintf(fp_matrix, "%d %d %d\n",
                 triplets.cluster_idx[k] + 1,
                 triplets.cell_idx[k] + 1,
                 triplets.counts[k]);
    }
    gzclose(fp_matrix);

    // Cleanup
    cb_hash_free(&cb_hash);
    triplet_free(&triplets);

    // ========================================================================
    // Return summary as list
    // ========================================================================

    SEXP result = PROTECT(allocVector(VECSXP, 4));
    SEXP names = PROTECT(allocVector(STRSXP, 4));

    SET_VECTOR_ELT(result, 0, ScalarInteger(n_cells));
    SET_VECTOR_ELT(result, 1, ScalarInteger(n_clusters));
    SET_VECTOR_ELT(result, 2, ScalarInteger(nnz));

    snprintf(path, 1024, "%s/barcodes.tsv.gz", out_dir);
    SET_VECTOR_ELT(result, 3, mkString(path));

    SET_STRING_ELT(names, 0, mkChar("n_cells"));
    SET_STRING_ELT(names, 1, mkChar("n_features"));
    SET_STRING_ELT(names, 2, mkChar("nnz"));
    SET_STRING_ELT(names, 3, mkChar("barcodes_path"));

    setAttrib(result, R_NamesSymbol, names);
    UNPROTECT(2);
    return result;
}
