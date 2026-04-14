#include "bam2db.h"
#include "mt19937ar.h"
#include "utils.h"

#include <R.h>
#include <Rinternals.h>

/* ── DNA encoding ────────────────────────────────────────────────────────── */

uint8_t encode_base(char base)
{
    switch (base) {
    case 'A': return 0b00;
    case 'C': return 0b01;
    case 'G': return 0b10;
    case 'T': return 0b11;
    default:  return BYTE_MASK;
    }
}

uint8_t *encode_DNA(const char *DNA_seq)
{
    uint16_t len_DNA = strlen(DNA_seq);
    uint8_t *encoded_DNA = (uint8_t *)calloc((len_DNA + BYTE_SIZE / BASE_BITS - 1) / BASE_BITS + 1, sizeof(uint8_t));
    if (!encoded_DNA) return NULL;

    uint8_t *p = encoded_DNA;
    int bit_index = BYTE_SIZE - BASE_BITS;
    for (size_t i = 0; i < len_DNA; i++) {
        uint8_t base = encode_base(DNA_seq[i]);
        if (base == BYTE_MASK) {
            free(encoded_DNA);
            return NULL;
        }
        *p |= (base & BASE_MASK) << bit_index;
        bit_index -= BASE_BITS;
        if (bit_index < 0) {
            bit_index = BYTE_SIZE - BASE_BITS;
            p++;
        }
    }
    return encoded_DNA;
}

char *decode_DNA(uint8_t *encoded_DNA, size_t n)
{
    char *decoded_DNA = (char *)calloc(n + 1, sizeof(char));
    if (!decoded_DNA) return NULL;
    decoded_DNA[n] = '\0';

    char *p = decoded_DNA;
    int bit_index = BYTE_SIZE - BASE_BITS;
    for (size_t i = 0; i < n; i++) {
        uint8_t base = (encoded_DNA[i / BASE_BITS] >> bit_index) & BASE_MASK;
        switch (base) {
        case 0b00: *p = 'A'; break;
        case 0b01: *p = 'C'; break;
        case 0b10: *p = 'G'; break;
        case 0b11: *p = 'T'; break;
        default:
            free(decoded_DNA);
            return NULL;
        }
        bit_index -= BASE_BITS;
        if (bit_index < 0)
            bit_index = BYTE_SIZE - BASE_BITS;
        p++;
    }
    return decoded_DNA;
}

/* ── Internal helpers ────────────────────────────────────────────────────── */

static uint64_t hash_djb2(const char *str, size_t len)
{
    uint64_t h = 5381;
    for (size_t i = 0; i < len; i++)
        h = ((h << 5) + h) + (uint64_t)(unsigned char)str[i];
    return h;
}

static int exec_sql(sqlite3 *db, const char *sql)
{
    char *err = NULL;
    int rc = sqlite3_exec(db, sql, NULL, 0, &err);
    if (rc != SQLITE_OK) {
        REprintf("SQL error: %s\n", err);
        sqlite3_free(err);
    }
    return rc;
}

/* ── Protocol filter ─────────────────────────────────────────────────────── */

static int protocol_filter(const bam1_t *rec, protocol_t proto)
{
    switch (proto) {
    case PROTO_10X: {
        uint8_t *xf = bam_aux_get(rec, "xf");
        if (!xf) return 0;
        int q = bam_aux2i(xf);
        return (q == 17 || q == 25);
    }
    case PROTO_ZUMIS:
        return 1;
    default:
        return 0;
    }
}

/* ── Identity extraction ─────────────────────────────────────────────────── */

static read_identity extract_identity(const bam1_t *rec, protocol_t proto)
{
    read_identity id = {NULL, NULL};
    switch (proto) {
    case PROTO_10X: {
        uint8_t *cb = bam_aux_get(rec, "CB");
        uint8_t *ub = bam_aux_get(rec, "UB");
        if (cb && ub) {
            id.cb = bam_aux2Z(cb);
            id.umi = bam_aux2Z(ub);
        }
        break;
    }
    case PROTO_ZUMIS: {
        uint8_t *bx = bam_aux_get(rec, "BX");
        uint8_t *ub = bam_aux_get(rec, "UB");
        if (bx && ub) {
            id.cb = bam_aux2Z(bx);
            id.umi = bam_aux2Z(ub);
        }
        break;
    }
    }
    return id;
}

/* ── CIGAR block extraction ──────────────────────────────────────────────── */

typedef struct {
    hts_pos_t start;
    hts_pos_t end;
} aligned_block;

static int extract_blocks(const bam1_t *rec, aligned_block *blocks, int max_blocks)
{
    int n = 0;
    uint32_t *cigar = bam_get_cigar(rec);
    int n_cigar = rec->core.n_cigar;
    hts_pos_t pos = rec->core.pos;

    for (int i = 0; i < n_cigar && n < max_blocks; i++) {
        int op = bam_cigar_op(cigar[i]);
        int len = bam_cigar_oplen(cigar[i]);

        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
            blocks[n].start = pos + 1;
            blocks[n].end = pos + len;
            n++;
            pos += len;
        } else if (op == BAM_CDEL || op == BAM_CREF_SKIP) {
            pos += len;
        } else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) {
            /* don't advance */
        } else if (op == BAM_CHARD_CLIP) {
            /* nothing */
        } else if (op == BAM_CPAD) {
            pos += len;
        }
    }
    return n;
}

/* ── Region overlap ──────────────────────────────────────────────────────── */

static int bsearch_chr(const region_registry *reg, const char *chr)
{
    int lo = 0, hi = (int)reg->n_chr - 1;
    while (lo <= hi) {
        int mid = lo + (hi - lo) / 2;
        int c = strcmp(chr, reg->chr_names[mid]);
        if (c == 0) return mid;
        if (c < 0) hi = mid - 1;
        else lo = mid + 1;
    }
    return -1;
}

static size_t find_region_overlaps(const bam1_t *rec,
                                    const bam_hdr_t *hdr,
                                    const region_registry *reg,
                                    size_t *out_indexes, size_t max_out)
{
    if (!reg || reg->n == 0) return 0;
    if (rec->core.tid < 0) return 0;

    const char *chr = hdr->target_name[rec->core.tid];
    int cid = bsearch_chr(reg, chr);
    if (cid < 0) return 0;

    /* Extract CIGAR blocks */
    aligned_block blocks[128];
    int n_blocks = extract_blocks(rec, blocks, 128);
    if (n_blocks == 0) return 0;

    size_t n_out = 0;
    size_t offset = reg->chr_first[cid];
    size_t count = reg->chr_count[cid];

    for (int b = 0; b < n_blocks; b++) {
        for (size_t i = 0; i < count; i++) {
            const region_t *r = &reg->items[offset + i];
            if (blocks[b].end < r->start) break;
            if (blocks[b].start > r->end) continue;

            /* Check for duplicate region_index */
            int dup = 0;
            for (size_t j = 0; j < n_out; j++) {
                if (out_indexes[j] == (offset + i + 1)) { dup = 1; break; }
            }
            if (dup) continue;

            if (n_out < max_out) {
                out_indexes[n_out] = offset + i + 1;
                n_out++;
            }
        }
    }
    return n_out;
}

/* ── Phase: load cell barcodes ───────────────────────────────────────────── */

static int load_cells(sqlite3 *db, const char *barcodes_file,
                      float rate_cell, unsigned int seed,
                      hash_table **out_ht, size_t *out_n_sampled)
{
    gzFile fp = gzopen(barcodes_file, "r");
    if (!fp) {
        REprintf("ERROR: Cannot open barcode file %s\n", barcodes_file);
        return 1;
    }

    size_t cap = 4096, n = 0;
    char **cells = (char **)malloc(cap * sizeof(char *));
    if (!cells) { gzclose(fp); return 1; }

    char buf[MAX_LINE_LENGTH];
    while (gzgets(fp, buf, sizeof(buf)) != NULL) {
        buf[strcspn(buf, "\n\r\t")] = '\0';
        if (strlen(buf) == 0) continue;
        if (n >= cap) {
            cap *= 2;
            cells = (char **)realloc(cells, cap * sizeof(char *));
            if (!cells) { gzclose(fp); return 1; }
        }
        cells[n] = strdup(buf);
        if (!cells[n]) { gzclose(fp); return 1; }
        n++;
    }
    gzclose(fp);
    REprintf("Total cell barcodes loaded: %zu\n", n);

    if (n == 0) { free(cells); return 1; }

    size_t n_sampled = (size_t)(n * rate_cell);
    if (n_sampled > n) n_sampled = n;
    if (n_sampled == 0) n_sampled = 1;

    size_t *seq = GetSeqInt(0, n - 1, 1);
    size_t *idx = SampleInt(seq, n, n_sampled, 0, seed);
    free(seq);
    qsort(idx, n_sampled, sizeof(size_t), vsI);

    hash_table *ht = hash_table_create((uint32_t)(n_sampled * 2), hash_djb2, free);
    if (!ht) { free(idx); free(cells); return 1; }

    if (exec_sql(db, "CREATE TABLE cell (cell_barcode TEXT);") != SQLITE_OK) {
        hash_table_destroy(ht); free(idx); free(cells); return 1;
    }

    sqlite3_stmt *stmt;
    if (sqlite3_prepare_v2(db, "INSERT INTO cell VALUES (?1);", -1, &stmt, NULL) != SQLITE_OK) {
        hash_table_destroy(ht); free(idx); free(cells); return 1;
    }

    exec_sql(db, "BEGIN TRANSACTION");
    size_t cell_index = 1;
    for (size_t i = 0; i < n_sampled; i++) {
        size_t ci = idx[i];
        size_t *idx_ptr = (size_t *)calloc(1, sizeof(size_t));
        if (!idx_ptr) continue;
        *idx_ptr = cell_index;

        if (hash_table_insert(ht, cells[ci], idx_ptr)) {
            sqlite3_bind_text(stmt, 1, cells[ci], -1, SQLITE_TRANSIENT);
            sqlite3_step(stmt);
            sqlite3_reset(stmt);
            cell_index++;
        } else {
            free(idx_ptr);
        }
    }
    exec_sql(db, "END TRANSACTION");
    sqlite3_finalize(stmt);

    *out_ht = ht;
    *out_n_sampled = cell_index - 1;

    for (size_t i = 0; i < n; i++) free(cells[i]);
    free(cells);
    free(idx);
    return 0;
}

/* ── Phase: load regions ─────────────────────────────────────────────────── */

static int load_regions(sqlite3 *db, const char *regions_file,
                        region_registry **out_reg)
{
    region_registry *reg = region_registry_build(regions_file);
    if (!reg) return 1;

    if (exec_sql(db, "CREATE TABLE region (region_index INTEGER PRIMARY KEY, id TEXT, symbol TEXT, chr TEXT, start INTEGER, end INTEGER, annotation1 TEXT, annotation2 TEXT);") != SQLITE_OK) {
        region_registry_destroy(reg);
        return 1;
    }

    sqlite3_stmt *stmt;
    if (sqlite3_prepare_v2(db, "INSERT INTO region VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8);", -1, &stmt, NULL) != SQLITE_OK) {
        region_registry_destroy(reg);
        return 1;
    }

    exec_sql(db, "BEGIN TRANSACTION");
    for (size_t i = 0; i < reg->n; i++) {
        sqlite3_bind_int(stmt, 1, (int)(i + 1));
        sqlite3_bind_text(stmt, 2, reg->items[i].id, -1, SQLITE_STATIC);
        sqlite3_bind_text(stmt, 3, reg->items[i].symbol, -1, SQLITE_STATIC);
        sqlite3_bind_text(stmt, 4, reg->items[i].chr, -1, SQLITE_STATIC);
        sqlite3_bind_int(stmt, 5, (int)reg->items[i].start);
        sqlite3_bind_int(stmt, 6, (int)reg->items[i].end);
        sqlite3_bind_text(stmt, 7, reg->items[i].annotation1 ? reg->items[i].annotation1 : "", -1, SQLITE_STATIC);
        sqlite3_bind_text(stmt, 8, reg->items[i].annotation2 ? reg->items[i].annotation2 : "", -1, SQLITE_STATIC);
        sqlite3_step(stmt);
        sqlite3_reset(stmt);
    }
    exec_sql(db, "END TRANSACTION");
    sqlite3_finalize(stmt);

    REprintf("Regions loaded: %zu\n", reg->n);
    *out_reg = reg;
    return 0;
}

/* ── Region query struct ─────────────────────────────────────────────────── */

typedef struct {
    int enabled;
    const char *chr;
    hts_pos_t start;
    hts_pos_t end;
} region_query;

/* ── Phase: process BAM reads ────────────────────────────────────────────── */

static int process_bam(sqlite3 *db, const char *bam_file,
                       hash_table *ht_cell,
                       const region_registry *reg,
                       float rate_depth, unsigned int seed,
                       protocol_t proto,
                       const region_query *roi,
                       size_t *out_total, size_t *out_sampled, size_t *out_valid)
{
    samFile *bam = sam_open(bam_file, "r");
    if (!bam) {
        REprintf("ERROR: Cannot open BAM file %s\n", bam_file);
        return 1;
    }

    bam_hdr_t *hdr = sam_hdr_read(bam);
    if (!hdr) {
        REprintf("ERROR: Cannot read BAM header from %s\n", bam_file);
        hts_close(bam);
        return 1;
    }

    bam1_t *rec = bam_init1();
    if (!rec) {
        bam_hdr_destroy(hdr);
        hts_close(bam);
        return 1;
    }

    init_genrand(seed);

    if (exec_sql(db, "CREATE TABLE umi (cell_index INTEGER, region_index INTEGER, encoded_umi BLOB, umi_len INTEGER);") != SQLITE_OK) {
        bam_destroy1(rec); bam_hdr_destroy(hdr); hts_close(bam); return 1;
    }

    sqlite3_stmt *stmt;
    if (sqlite3_prepare_v2(db, "INSERT INTO umi VALUES (?1, ?2, ?3, ?4);", -1, &stmt, NULL) != SQLITE_OK) {
        bam_destroy1(rec); bam_hdr_destroy(hdr); hts_close(bam); return 1;
    }

    exec_sql(db, "BEGIN TRANSACTION");

    size_t total = 0, sampled = 0, valid = 0;
    size_t region_indexes[128];

    /* Build chr name → tid map for ROI filtering */
    int roi_tid = -1;
    if (roi && roi->enabled) {
        roi_tid = bam_name2id(hdr, roi->chr);
        if (roi_tid < 0) {
            REprintf("ERROR: Unknown contig '%s' in BAM %s\n", roi->chr, bam_file);
            sqlite3_finalize(stmt);
            exec_sql(db, "END TRANSACTION");
            bam_destroy1(rec); bam_hdr_destroy(hdr); hts_close(bam);
            return 1;
        }
    }

    /* Decide iteration mode */
    int use_iterator = (roi && roi->enabled);
    hts_itr_t *itr = NULL;
    if (use_iterator) {
        hts_idx_t *idx = sam_index_load(bam, bam_file);
        if (!idx) {
            REprintf("ERROR: No index found for %s (region query requires .bai)\n", bam_file);
            sqlite3_finalize(stmt);
            exec_sql(db, "END TRANSACTION");
            bam_destroy1(rec); bam_hdr_destroy(hdr); hts_close(bam);
            return 1;
        }
        itr = sam_itr_queryi(idx, roi_tid, roi->start, roi->end);
        hts_idx_destroy(idx);
        if (!itr) {
            REprintf("ERROR: Failed to create iterator for %s:%ld-%ld\n",
                    roi->chr, (long)roi->start, (long)roi->end);
            sqlite3_finalize(stmt);
            exec_sql(db, "END TRANSACTION");
            bam_destroy1(rec); bam_hdr_destroy(hdr); hts_close(bam);
            return 1;
        }
        REprintf("Querying region %s:%ld-%ld\n", roi->chr,
                (long)roi->start, (long)roi->end);
    }

    int ret_code = 0;
    while (1) {
        int read_ok = use_iterator ? sam_itr_next(bam, itr, rec) : sam_read1(bam, hdr, rec);
        if (read_ok < 0) break;

        total++;

        /* ROI filter: skip reads outside target chr */
        if (roi && roi->enabled && rec->core.tid != roi_tid) continue;

        /* Protocol-specific quality filter */
        if (!protocol_filter(rec, proto)) continue;

        /* Extract cell barcode + UMI */
        read_identity id = extract_identity(rec, proto);
        if (!id.cb || !id.umi) continue;

        /* Cell barcode lookup */
        void *cell_lookup = hash_table_lookup(ht_cell, id.cb);
        if (!cell_lookup) continue;
        size_t cell_index = *(size_t *)cell_lookup;

        /* Depth subsampling */
        if (genrand_real1() >= rate_depth) continue;
        sampled++;

        /* Region overlap */
        size_t n_ov = find_region_overlaps(rec, hdr, reg, region_indexes, 128);
        if (n_ov == 0) continue;

        /* Encode UMI */
        uint8_t *encoded = encode_DNA(id.umi);
        if (!encoded) continue;
        size_t enc_size = (strlen(id.umi) + BYTE_SIZE / BASE_BITS - 1) / BASE_BITS;

        for (size_t k = 0; k < n_ov; k++) {
            sqlite3_bind_int(stmt, 1, (int)cell_index);
            sqlite3_bind_int(stmt, 2, (int)region_indexes[k]);
            sqlite3_bind_blob(stmt, 3, encoded, (int)enc_size, NULL);
            sqlite3_bind_int(stmt, 4, (int)strlen(id.umi));

            if (sqlite3_step(stmt) == SQLITE_DONE)
                valid++;

            sqlite3_reset(stmt);
        }

        free(encoded);
    }

    exec_sql(db, "END TRANSACTION");
    sqlite3_finalize(stmt);

    if (use_iterator && itr) hts_itr_destroy(itr);

    REprintf("BAM %s: total=%zu, sampled=%zu, valid=%zu\n",
            bam_file, total, sampled, valid);

    *out_total = total;
    *out_sampled = sampled;
    *out_valid = valid;

    bam_destroy1(rec);
    bam_hdr_destroy(hdr);
    hts_close(bam);
    return ret_code;
}

/* ── Phase: generate 10x output ──────────────────────────────────────────── */

static int generate_10x(sqlite3 *db, const char *path_out,
                        const char *bam_file,
                        float rate_cell, float rate_depth,
                        size_t total, size_t sampled, size_t valid,
                        int umi_copies_flag)
{
    char path[1024];

    snprintf(path, sizeof(path), "%s/barcodes.tsv.gz", path_out);
    gzFile f_bc = gzopen(path, "wb");
    if (!f_bc) { REprintf("ERROR: cannot open %s\n", path); return 1; }

    snprintf(path, sizeof(path), "%s/features.tsv.gz", path_out);
    gzFile f_ft = gzopen(path, "wb");
    if (!f_ft) { gzclose(f_bc); REprintf("ERROR: cannot open %s\n", path); return 1; }

    snprintf(path, sizeof(path), "%s/matrix.mtx.gz", path_out);
    gzFile f_mx = gzopen(path, "wb");
    if (!f_mx) { gzclose(f_bc); gzclose(f_ft); REprintf("ERROR: cannot open %s\n", path); return 1; }

    if (exec_sql(db,
        "CREATE TABLE mtx AS "
        "SELECT region_index, cell_index, COUNT(DISTINCT encoded_umi) AS expression_level "
        "FROM umi GROUP BY cell_index, region_index;") != SQLITE_OK) {
        gzclose(f_bc); gzclose(f_ft); gzclose(f_mx); return 1;
    }

    size_t n_mtx = nrow_sql_table(db, "mtx");
    size_t n_bc  = nrow_sql_table(db, "cell");
    size_t n_ft  = nrow_sql_table(db, "region");

    gzprintf(f_mx,
        "%%MatrixMarket matrix coordinate integer general\n"
        "%%%%\t\"parent_bam\": \"%s\",\n"
        "%%%%\t\"rate_cell\": %.3f,\n"
        "%%%%\t\"rate_depth\": %.3f,\n"
        "%%%%\t\"total_n_reads\": %zu,\n"
        "%%%%\t\"sampled_n_reads\": %zu,\n"
        "%%%%\t\"valid_n_reads\": %zu\n",
        bam_file, rate_cell, rate_depth, total, sampled, valid);
    gzprintf(f_mx, "%zu %zu %zu\n", n_ft, n_bc, n_mtx);

    table2gz(db, "mtx", f_mx, 0, " ");
    table2gz(db, "cell", f_bc, 0, "\n");
    table2gz(db, "region", f_ft, 0, "\t");

    gzclose(f_bc);
    gzclose(f_ft);
    gzclose(f_mx);

    if (umi_copies_flag) {
        if (exec_sql(db,
            "CREATE TABLE numi AS "
            "SELECT region_index, cell_index, encoded_umi, umi_len, COUNT(*) AS n_copy "
            "FROM umi GROUP BY cell_index, region_index, encoded_umi;") != SQLITE_OK)
            return 1;

        snprintf(path, sizeof(path), "%s/umi.tsv.gz", path_out);
        gzFile f_umi = gzopen(path, "wb");
        if (f_umi) {
            table2gz(db, "numi", f_umi, 0, "\t");
            gzclose(f_umi);
        }
    }

    return 0;
}

/* ── Public API ──────────────────────────────────────────────────────────── */

int bam2db(
    char *bam_file,
    char *db_file,
    char *path_out,
    char *barcodes_file,
    char *regions_file,
    float rate_cell,
    float rate_depth,
    unsigned int seed,
    int umi_copies_flag,
    protocol_t protocol,
    const char *region_chr,
    hts_pos_t region_start,
    hts_pos_t region_end)
{
    if (rate_cell <= 0.0f || rate_cell > 1.0f) {
        REprintf("ERROR: rate_cell must be in (0, 1], got %f\n", rate_cell);
        return 1;
    }
    if (rate_depth <= 0.0f || rate_depth > 1.0f) {
        REprintf("ERROR: rate_depth must be in (0, 1], got %f\n", rate_depth);
        return 1;
    }

    sqlite3 *db = NULL;
    hash_table *ht_cell = NULL;
    region_registry *reg = NULL;
    int ret = 1;

    region_query roi;
    roi.enabled = (region_chr != NULL && region_chr[0] != '\0');
    roi.chr = region_chr;
    roi.start = region_start;
    roi.end = region_end;

    if (sqlite3_open(db_file, &db) != SQLITE_OK) {
        REprintf("ERROR: Cannot open/create database: %s\n", sqlite3_errmsg(db));
        sqlite3_close(db);
        return 1;
    }

    size_t n_cells = 0;
    if (load_cells(db, barcodes_file, rate_cell, seed, &ht_cell, &n_cells) != 0)
        goto cleanup;

    if (load_regions(db, regions_file, &reg) != 0)
        goto cleanup;

    size_t total = 0, sampled = 0, valid = 0;
    if (process_bam(db, bam_file, ht_cell, reg, rate_depth, seed,
                    protocol, &roi, &total, &sampled, &valid) != 0)
        goto cleanup;

    if (generate_10x(db, path_out, bam_file, rate_cell, rate_depth,
                     total, sampled, valid, umi_copies_flag) != 0)
        goto cleanup;

    ret = 0;

cleanup:
    if (ht_cell)  hash_table_destroy(ht_cell);
    if (reg)      region_registry_destroy(reg);
    if (db)       sqlite3_close(db);
    return ret;
}

/* ── Utility: dump table to gzipped file ─────────────────────────────────── */

int table2gz(sqlite3 *db_handle, const char *table_name,
             gzFile gz_file_ptr, unsigned int header, const char *delim)
{
    char sql[1024];
    snprintf(sql, sizeof(sql), "SELECT * FROM %s;", table_name);

    sqlite3_stmt *stmt;
    if (sqlite3_prepare_v2(db_handle, sql, -1, &stmt, NULL) != SQLITE_OK) {
        REprintf("SQL error: %s\n", sqlite3_errmsg(db_handle));
        return 1;
    }

    int ncols = sqlite3_column_count(stmt);

    if (header) {
        for (int i = 0; i < ncols; i++) {
            gzprintf(gz_file_ptr, "%s", sqlite3_column_name(stmt, i));
            if (i < ncols - 1) gzprintf(gz_file_ptr, "%s", delim);
        }
        gzprintf(gz_file_ptr, "\n");
    }

    while (sqlite3_step(stmt) == SQLITE_ROW) {
        for (int i = 0; i < ncols; i++) {
            switch (sqlite3_column_type(stmt, i)) {
            case SQLITE_INTEGER:
                gzprintf(gz_file_ptr, "%d", sqlite3_column_int(stmt, i));
                break;
            case SQLITE_FLOAT:
                gzprintf(gz_file_ptr, "%.3f", sqlite3_column_double(stmt, i));
                break;
            case SQLITE_TEXT:
                gzprintf(gz_file_ptr, "%s", sqlite3_column_text(stmt, i));
                break;
            case SQLITE_BLOB: {
                int blob_len = sqlite3_column_bytes(stmt, i);
                const uint8_t *blob = (const uint8_t *)sqlite3_column_blob(stmt, i);
                if (strcmp(table_name, "numi") == 0 && i == 2) {
                    int umi_len = sqlite3_column_int(stmt, 3);
                    char *decoded = decode_DNA((uint8_t *)blob, (size_t)umi_len);
                    if (decoded) { gzprintf(gz_file_ptr, "%s", decoded); free(decoded); }
                    else gzprintf(gz_file_ptr, "NULL");
                } else if (strcmp(table_name, "umi") == 0 && i == 2) {
                    int umi_len = sqlite3_column_int(stmt, 3);
                    char *decoded = decode_DNA((uint8_t *)blob, (size_t)umi_len);
                    if (decoded) { gzprintf(gz_file_ptr, "%s", decoded); free(decoded); }
                    else gzprintf(gz_file_ptr, "NULL");
                } else {
                    gzprintf(gz_file_ptr, "<BLOB:%d>", blob_len);
                }
                break;
            }
            default:
                gzprintf(gz_file_ptr, "NULL");
                break;
            }
            if (i < ncols - 1) gzprintf(gz_file_ptr, "%s", delim);
        }
        gzprintf(gz_file_ptr, "\n");
    }

    sqlite3_finalize(stmt);
    return 0;
}

/* ── Utility: row count ──────────────────────────────────────────────────── */

size_t nrow_sql_table(sqlite3 *db_handle, const char *table_name)
{
    char sql[256];
    snprintf(sql, sizeof(sql), "SELECT COUNT(*) FROM %s;", table_name);

    sqlite3_stmt *stmt;
    if (sqlite3_prepare_v2(db_handle, sql, -1, &stmt, NULL) != SQLITE_OK) {
        REprintf("SQL error: %s\n", sqlite3_errmsg(db_handle));
        return 0;
    }

    size_t nrow = 0;
    if (sqlite3_step(stmt) == SQLITE_ROW)
        nrow = (size_t)sqlite3_column_int(stmt, 0);

    sqlite3_finalize(stmt);
    return nrow;
}
