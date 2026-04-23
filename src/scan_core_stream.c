#include "exonblocks_shared.h"
#include <stdio.h>
#include <zlib.h>

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
            int start_idx = find_first_overlap_end(cluster_end, n_clusters, blk_start);
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
