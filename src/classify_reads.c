#include "exonblocks_shared.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <zlib.h>
#include <omp.h>

typedef enum {
    READ_SPLICED,
    READ_CHIMERIC,
    READ_UNSPLICED,
    READ_AMBIGUOUS,
    READ_UNASSIGNED
} ReadClass;

typedef struct {
    int32_t start;
    int32_t end;
    int32_t gene_id;
    int8_t is_rev;
} SortedGene;

typedef struct {
    SortedGene *fwd;
    SortedGene *rev;
    int *fwd_ends;
    int *rev_ends;
    int n_fwd;
    int n_rev;
} StrandGeneSpan;

typedef struct {
    ExonClusterIndex fwd_exon;
    ExonClusterIndex rev_exon;
    GeneIndex fwd_gene;
    GeneIndex rev_gene;
    StrandGeneSpan gene_span;
    int *cluster_start_sorted;
    int *cluster_end_sorted;
    int *cluster_orig_idx;
    int8_t *cluster_is_rev;
    int n_clusters;
    char *name;
} ChrAnnotNew;

static int cmp_sorted_gene(const void *a, const void *b) {
    const SortedGene *ga = (const SortedGene *)a;
    const SortedGene *gb = (const SortedGene *)b;
    if (ga->end != gb->end) return (ga->end > gb->end) - (ga->end < gb->end);
    if (ga->start != gb->start) return (ga->start > gb->start) - (ga->start < gb->start);
    return (ga->gene_id > gb->gene_id) - (ga->gene_id < gb->gene_id);
}

static void build_gene_span_c(int *starts, int *ends, int *gene_ids, char **strands, int n,
                                StrandGeneSpan *out) {
    int max_gid = 0;
    for (int i = 0; i < n; i++)
        if (gene_ids[i] > max_gid) max_gid = gene_ids[i];

    int32_t *g_min_start = (int32_t *)calloc(max_gid + 1, sizeof(int32_t));
    int32_t *g_max_end = (int32_t *)calloc(max_gid + 1, sizeof(int32_t));
    int8_t *g_strand = (int8_t *)calloc(max_gid + 1, sizeof(int8_t));
    int *g_set = (int *)calloc(max_gid + 1, sizeof(int));

    for (int i = 0; i <= max_gid; i++) g_min_start[i] = INT32_MAX;

    for (int i = 0; i < n; i++) {
        int gid = gene_ids[i];
        if (!g_set[gid]) {
            g_min_start[gid] = starts[i];
            g_max_end[gid] = ends[i];
            g_strand[gid] = (strands[i][0] == '-') ? 1 : 0;
            g_set[gid] = 1;
        } else {
            if (starts[i] < g_min_start[gid]) g_min_start[gid] = starts[i];
            if (ends[i] > g_max_end[gid]) g_max_end[gid] = ends[i];
        }
    }

    int n_fwd = 0, n_rev = 0;
    for (int i = 1; i <= max_gid; i++) {
        if (g_set[i]) {
            if (g_strand[i] == 0) n_fwd++;
            else n_rev++;
        }
    }

    out->fwd = NULL; out->rev = NULL;
    out->fwd_ends = NULL; out->rev_ends = NULL;
    out->n_fwd = 0; out->n_rev = 0;

    if (n_fwd > 0) {
        out->fwd = (SortedGene *)malloc(n_fwd * sizeof(SortedGene));
        out->fwd_ends = (int *)malloc(n_fwd * sizeof(int));
    }
    if (n_rev > 0) {
        out->rev = (SortedGene *)malloc(n_rev * sizeof(SortedGene));
        out->rev_ends = (int *)malloc(n_rev * sizeof(int));
    }

    int fi = 0, ri = 0;
    for (int i = 1; i <= max_gid; i++) {
        if (!g_set[i]) continue;
        SortedGene gs;
        gs.start = g_min_start[i];
        gs.end = g_max_end[i];
        gs.gene_id = i;
        gs.is_rev = g_strand[i];
        if (g_strand[i] == 0) {
            out->fwd[fi] = gs;
            fi++;
        } else {
            out->rev[ri] = gs;
            ri++;
        }
    }

    if (n_fwd > 0) {
        qsort(out->fwd, n_fwd, sizeof(SortedGene), cmp_sorted_gene);
        for (int i = 0; i < n_fwd; i++) out->fwd_ends[i] = out->fwd[i].end;
        out->n_fwd = n_fwd;
    }
    if (n_rev > 0) {
        qsort(out->rev, n_rev, sizeof(SortedGene), cmp_sorted_gene);
        for (int i = 0; i < n_rev; i++) out->rev_ends[i] = out->rev[i].end;
        out->n_rev = n_rev;
    }

    free(g_min_start); free(g_max_end); free(g_strand); free(g_set);
}

static void free_strand_gene_span(StrandGeneSpan *sgs) {
    if (sgs->fwd) { free(sgs->fwd); free(sgs->fwd_ends); }
    if (sgs->rev) { free(sgs->rev); free(sgs->rev_ends); }
}

static void build_exon_cluster_index(int *starts, int *ends, int *gene_ids,
                                      int *cluster_ids, int n,
                                      ExonClusterIndex *idx) {
    if (n == 0) { idx->n = 0; idx->start = NULL; idx->end = NULL; idx->gene_id = NULL; idx->exon_cluster_id = NULL; return; }
    idx->n = n;
    idx->start = (int32_t *)malloc(n * sizeof(int32_t));
    idx->end = (int32_t *)malloc(n * sizeof(int32_t));
    idx->gene_id = (int32_t *)malloc(n * sizeof(int32_t));
    idx->exon_cluster_id = (int32_t *)malloc(n * sizeof(int32_t));
    for (int i = 0; i < n; i++) {
        idx->start[i] = starts[i];
        idx->end[i] = ends[i];
        idx->gene_id[i] = gene_ids[i];
        idx->exon_cluster_id[i] = cluster_ids[i];
    }
}

static void free_exon_cluster_index(ExonClusterIndex *idx) {
    if (idx->n > 0) { free(idx->start); free(idx->end); free(idx->gene_id); free(idx->exon_cluster_id); }
}

static void build_gene_index(int *starts, int *ends, int *gene_ids, int n,
                               GeneIndex *idx) {
    if (n == 0) { idx->n = 0; idx->start = NULL; idx->end = NULL; idx->gene_id = NULL; return; }
    idx->n = n;
    idx->start = (int32_t *)malloc(n * sizeof(int32_t));
    idx->end = (int32_t *)malloc(n * sizeof(int32_t));
    idx->gene_id = (int32_t *)malloc(n * sizeof(int32_t));
    for (int i = 0; i < n; i++) {
        idx->start[i] = starts[i];
        idx->end[i] = ends[i];
        idx->gene_id[i] = gene_ids[i];
    }
}

static void free_gene_index(GeneIndex *idx) {
    if (idx->n > 0) { free(idx->start); free(idx->end); free(idx->gene_id); }
}

static void sort_exon_cluster_index_by_end(ExonClusterIndex *idx) {
    if (idx->n <= 1) return;
    for (size_t i = 1; i < idx->n; i++) {
        int32_t tmp_end = idx->end[i];
        int32_t tmp_start = idx->start[i];
        int32_t tmp_gid = idx->gene_id[i];
        int32_t tmp_cid = idx->exon_cluster_id[i];
        size_t j = i - 1;
        while (j >= 0 && idx->end[j] > tmp_end) {
            idx->end[j+1] = idx->end[j];
            idx->start[j+1] = idx->start[j];
            idx->gene_id[j+1] = idx->gene_id[j];
            idx->exon_cluster_id[j+1] = idx->exon_cluster_id[j];
            if (j == 0) break;
            j--;
        }
        if (idx->end[0] > tmp_end) {
            idx->end[0] = tmp_end;
            idx->start[0] = tmp_start;
            idx->gene_id[0] = tmp_gid;
            idx->exon_cluster_id[0] = tmp_cid;
        } else {
            idx->end[j+1] = tmp_end;
            idx->start[j+1] = tmp_start;
            idx->gene_id[j+1] = tmp_gid;
            idx->exon_cluster_id[j+1] = tmp_cid;
        }
    }
}

static void sort_gene_index_by_end(GeneIndex *idx) {
    if (idx->n <= 1) return;
    for (size_t i = 1; i < idx->n; i++) {
        int32_t tmp_end = idx->end[i];
        int32_t tmp_start = idx->start[i];
        int32_t tmp_gid = idx->gene_id[i];
        size_t j = i - 1;
        while (j >= 0 && idx->end[j] > tmp_end) {
            idx->end[j+1] = idx->end[j];
            idx->start[j+1] = idx->start[j];
            idx->gene_id[j+1] = idx->gene_id[j];
            if (j == 0) break;
            j--;
        }
        if (idx->end[0] > tmp_end) {
            idx->end[0] = tmp_end;
            idx->start[0] = tmp_start;
            idx->gene_id[0] = tmp_gid;
        } else {
            idx->end[j+1] = tmp_end;
            idx->start[j+1] = tmp_start;
            idx->gene_id[j+1] = tmp_gid;
        }
    }
}

static void split_exon_index_by_strand(ExonClusterIndex *src,
                                        char **strands, int n,
                                        ExonClusterIndex *fwd, ExonClusterIndex *rev) {
    int n_fwd = 0, n_rev = 0;
    for (int i = 0; i < n; i++) {
        if (strands[i][0] == '+') n_fwd++;
        else if (strands[i][0] == '-') n_rev++;
    }
    fwd->n = n_fwd;
    fwd->start = n_fwd ? (int32_t *)malloc(n_fwd * sizeof(int32_t)) : NULL;
    fwd->end = n_fwd ? (int32_t *)malloc(n_fwd * sizeof(int32_t)) : NULL;
    fwd->gene_id = n_fwd ? (int32_t *)malloc(n_fwd * sizeof(int32_t)) : NULL;
    fwd->exon_cluster_id = n_fwd ? (int32_t *)malloc(n_fwd * sizeof(int32_t)) : NULL;
    rev->n = n_rev;
    rev->start = n_rev ? (int32_t *)malloc(n_rev * sizeof(int32_t)) : NULL;
    rev->end = n_rev ? (int32_t *)malloc(n_rev * sizeof(int32_t)) : NULL;
    rev->gene_id = n_rev ? (int32_t *)malloc(n_rev * sizeof(int32_t)) : NULL;
    rev->exon_cluster_id = n_rev ? (int32_t *)malloc(n_rev * sizeof(int32_t)) : NULL;
    int fi = 0, ri = 0;
    for (int i = 0; i < n; i++) {
        if (strands[i][0] == '+') {
            fwd->start[fi] = src->start[i];
            fwd->end[fi] = src->end[i];
            fwd->gene_id[fi] = src->gene_id[i];
            fwd->exon_cluster_id[fi] = src->exon_cluster_id[i];
            fi++;
        } else if (strands[i][0] == '-') {
            rev->start[ri] = src->start[i];
            rev->end[ri] = src->end[i];
            rev->gene_id[ri] = src->gene_id[i];
            rev->exon_cluster_id[ri] = src->exon_cluster_id[i];
            ri++;
        }
    }
    sort_exon_cluster_index_by_end(fwd);
    sort_exon_cluster_index_by_end(rev);
}

static void split_gene_index_by_strand(GeneIndex *src,
                                        char **strands, int n,
                                        GeneIndex *fwd, GeneIndex *rev) {
    int n_fwd = 0, n_rev = 0;
    for (int i = 0; i < n; i++) {
        if (strands[i][0] == '+') n_fwd++;
        else if (strands[i][0] == '-') n_rev++;
    }
    fwd->n = n_fwd;
    fwd->start = n_fwd ? (int32_t *)malloc(n_fwd * sizeof(int32_t)) : NULL;
    fwd->end = n_fwd ? (int32_t *)malloc(n_fwd * sizeof(int32_t)) : NULL;
    fwd->gene_id = n_fwd ? (int32_t *)malloc(n_fwd * sizeof(int32_t)) : NULL;
    rev->n = n_rev;
    rev->start = n_rev ? (int32_t *)malloc(n_rev * sizeof(int32_t)) : NULL;
    rev->end = n_rev ? (int32_t *)malloc(n_rev * sizeof(int32_t)) : NULL;
    rev->gene_id = n_rev ? (int32_t *)malloc(n_rev * sizeof(int32_t)) : NULL;
    int fi = 0, ri = 0;
    for (int i = 0; i < n; i++) {
        if (strands[i][0] == '+') {
            fwd->start[fi] = src->start[i];
            fwd->end[fi] = src->end[i];
            fwd->gene_id[fi] = src->gene_id[i];
            fi++;
        } else if (strands[i][0] == '-') {
            rev->start[ri] = src->start[i];
            rev->end[ri] = src->end[i];
            rev->gene_id[ri] = src->gene_id[i];
            ri++;
        }
    }
    sort_gene_index_by_end(fwd);
    sort_gene_index_by_end(rev);
}

static int find_overlapping_gene_ids_exon(ExonClusterIndex *idx,
                                           int qstart, int qend,
                                           int *out_ids, int *out_clusters, int max_out) {
    int n_out = 0;
    if (idx->n == 0) return 0;
    int start_idx = find_first_overlap_end(idx->end, (int)idx->n, qstart);
    for (int i = start_idx; i < (int)idx->n; i++) {
        if (idx->start[i] > qend) break;
        if (idx->end[i] >= qstart && idx->start[i] <= qend) {
            int gid = idx->gene_id[i];
            int found = 0;
            for (int j = 0; j < n_out; j++) {
                if (out_ids[j] == gid) { found = 1; break; }
            }
            if (!found && n_out < max_out) {
                out_ids[n_out] = gid;
                out_clusters[n_out] = idx->exon_cluster_id[i];
                n_out++;
            }
        }
    }
    return n_out;
}

static int find_overlapping_gene_ids_gene(GeneIndex *idx,
                                           int qstart, int qend,
                                           int *out_ids, int max_out) {
    int n_out = 0;
    if (idx->n == 0) return 0;
    int start_idx = find_first_overlap_end(idx->end, (int)idx->n, qstart);
    for (int i = start_idx; i < (int)idx->n; i++) {
        if (idx->start[i] > qend) break;
        if (idx->end[i] >= qstart && idx->start[i] <= qend) {
            int gid = idx->gene_id[i];
            int found = 0;
            for (int j = 0; j < n_out; j++) {
                if (out_ids[j] == gid) { found = 1; break; }
            }
            if (!found && n_out < max_out) out_ids[n_out++] = gid;
            if (n_out >= max_out) break;
        }
    }
    return n_out;
}

static int intersect_sorted_ints(int *a, int na, int *b, int nb) {
    int i = 0, j = 0, k = 0;
    while (i < na && j < nb) {
        if (a[i] < b[j]) i++;
        else if (a[i] > b[j]) j++;
        else { a[k++] = a[i]; i++; j++; }
    }
    return k;
}

static int any_overlap_gene_index(GeneIndex *idx, int qstart, int qend) {
    if (idx->n == 0) return 0;
    int start_idx = find_first_overlap_end(idx->end, (int)idx->n, qstart);
    for (int i = start_idx; i < (int)idx->n; i++) {
        if (idx->start[i] > qend) break;
        if (idx->end[i] >= qstart && idx->start[i] <= qend) return 1;
    }
    return 0;
}

static int any_overlap_exon_index(ExonClusterIndex *idx, int qstart, int qend) {
    if (idx->n == 0) return 0;
    int start_idx = find_first_overlap_end(idx->end, (int)idx->n, qstart);
    for (int i = start_idx; i < (int)idx->n; i++) {
        if (idx->start[i] > qend) break;
        if (idx->end[i] >= qstart && idx->start[i] <= qend) return 1;
    }
    return 0;
}

static int collect_cluster_overlaps(ExonClusterIndex *idx, int qstart, int qend,
                                    int *out_clusters, int max_out) {
    int n_out = 0;
    if (idx->n == 0) return 0;
    int start_idx = find_first_overlap_end(idx->end, (int)idx->n, qstart);
    for (int i = start_idx; i < (int)idx->n && n_out < max_out; i++) {
        if (idx->start[i] > qend) break;
        if (idx->end[i] >= qstart && idx->start[i] <= qend) {
            int cid = idx->exon_cluster_id[i];
            int found = 0;
            for (int j = 0; j < n_out; j++) {
                if (out_clusters[j] == cid) { found = 1; break; }
            }
            if (!found) out_clusters[n_out++] = cid;
        }
    }
    return n_out;
}

static void write_category_output(const char *out_dir, const char *category,
                                  TripletCOO *triplets, CBHash *cb_hash,
                                  int n_clusters_global,
                                  int *cluster_start_global, int *cluster_end_global) {
    char path[2048];
    char subdir[1024];
    snprintf(subdir, 1024, "%s/%s", out_dir, category);
    mkdir_p(subdir);

    snprintf(path, 2048, "%s/barcodes.tsv.gz", subdir);
    gzFile fp_bc = gzopen(path, "wb");
    if (!fp_bc) error("Cannot open %s", path);
    for (int i = 0; i < CB_HASH_SIZE; i++) {
        CBEntry *e = cb_hash->buckets[i];
        while (e) {
            if (e->cell_id >= 1 && e->cell_id <= cb_hash->count)
                gzprintf(fp_bc, "%s\n", e->key);
            e = e->next;
        }
    }
    gzclose(fp_bc);

    snprintf(path, 2048, "%s/features.tsv.gz", subdir);
    gzFile fp_feat = gzopen(path, "wb");
    if (!fp_feat) error("Cannot open %s", path);
    gzprintf(fp_feat, "feature_id\tfeature_name\tfeature_type\n");
    for (int i = 0; i < n_clusters_global; i++) {
        gzprintf(fp_feat, "cluster_%d\t%d-%d\texon_cluster\n",
                 i + 1, cluster_start_global[i], cluster_end_global[i]);
    }
    gzclose(fp_feat);

    snprintf(path, 2048, "%s/matrix.mtx.gz", subdir);
    gzFile fp_mtx = gzopen(path, "wb");
    if (!fp_mtx) error("Cannot open %s", path);
    gzprintf(fp_mtx, "%%%%MatrixMarket matrix coordinate integer general\n");
    gzprintf(fp_mtx, "%%%%metadata_json: {\"software_version\": \"exonBlocks\", \"category\": \"%s\"}\n", category);
    gzprintf(fp_mtx, "%d %d %d\n", n_clusters_global, cb_hash->count, triplets->size);
    for (int k = 0; k < triplets->size; k++) {
        gzprintf(fp_mtx, "%d %d %d\n",
                 triplets->cluster_idx[k] + 1,
                 triplets->cell_idx[k] + 1,
                 triplets->counts[k]);
    }
    gzclose(fp_mtx);
}

SEXP _exonBlocks_classify_reads(
    SEXP bam_path_,
    SEXP exon_chr_,
    SEXP exon_starts_,
    SEXP exon_ends_,
    SEXP exon_gene_ids_,
    SEXP exon_strand_,
    SEXP intron_chr_,
    SEXP intron_starts_,
    SEXP intron_ends_,
    SEXP intron_gene_ids_,
    SEXP intron_strand_,
    SEXP cluster_chr_,
    SEXP cluster_starts_,
    SEXP cluster_ends_,
    SEXP cluster_strand_,
    SEXP cluster_ids_,
    SEXP tech_,
    SEXP xf_vals_,
    SEXP out_dir_,
    SEXP gene_span_pad_)
{
    const char *bam_path = CHAR(STRING_ELT(bam_path_, 0));
    const char *tech = CHAR(STRING_ELT(tech_, 0));
    const char *out_dir = CHAR(STRING_ELT(out_dir_, 0));
    int gene_span_pad = INTEGER(gene_span_pad_)[0];

    int n_exons_all = LENGTH(exon_starts_);
    int n_introns_all = LENGTH(intron_starts_);
    int n_clusters_all = LENGTH(cluster_starts_);

    int *exon_starts_all = INTEGER(exon_starts_);
    int *exon_ends_all = INTEGER(exon_ends_);
    int *exon_gene_ids_all = INTEGER(exon_gene_ids_);
    int *intron_starts_all = INTEGER(intron_starts_);
    int *intron_ends_all = INTEGER(intron_ends_);
    int *intron_gene_ids_all = INTEGER(intron_gene_ids_);
    int *cluster_starts_all = INTEGER(cluster_starts_);
    int *cluster_ends_all = INTEGER(cluster_ends_);
    int *cluster_ids_all = INTEGER(cluster_ids_);

    int nxf = LENGTH(xf_vals_);
    long xf_allow[16];
    for (int i = 0; i < nxf && i < 16; ++i)
        xf_allow[i] = (long)INTEGER(xf_vals_)[i];
    TechProtocol tech_id = parse_tech(tech);

    SEXP strand_r = PROTECT(allocVector(STRSXP, n_exons_all));
    for (int i = 0; i < n_exons_all; i++) SET_STRING_ELT(strand_r, i, STRING_ELT(exon_strand_, i));
    char **exon_strand_c = (char **)malloc(n_exons_all * sizeof(char *));
    for (int i = 0; i < n_exons_all; i++) exon_strand_c[i] = (char *)CHAR(STRING_ELT(strand_r, i));

    SEXP strand_ri = PROTECT(allocVector(STRSXP, n_introns_all));
    for (int i = 0; i < n_introns_all; i++) SET_STRING_ELT(strand_ri, i, STRING_ELT(intron_strand_, i));
    char **intron_strand_c = (char **)malloc(n_introns_all * sizeof(char *));
    for (int i = 0; i < n_introns_all; i++) intron_strand_c[i] = (char *)CHAR(STRING_ELT(strand_ri, i));

    SEXP strand_rc = PROTECT(allocVector(STRSXP, n_clusters_all));
    for (int i = 0; i < n_clusters_all; i++) SET_STRING_ELT(strand_rc, i, STRING_ELT(cluster_strand_, i));
    char **cluster_strand_c = (char **)malloc(n_clusters_all * sizeof(char *));
    for (int i = 0; i < n_clusters_all; i++) cluster_strand_c[i] = (char *)CHAR(STRING_ELT(strand_rc, i));

    char **exon_chr_c = (char **)malloc(n_exons_all * sizeof(char *));
    for (int i = 0; i < n_exons_all; i++) exon_chr_c[i] = (char *)CHAR(STRING_ELT(exon_chr_, i));
    char **intron_chr_c = (char **)malloc(n_introns_all * sizeof(char *));
    for (int i = 0; i < n_introns_all; i++) intron_chr_c[i] = (char *)CHAR(STRING_ELT(intron_chr_, i));
    char **cluster_chr_c = (char **)malloc(n_clusters_all * sizeof(char *));
    for (int i = 0; i < n_clusters_all; i++) cluster_chr_c[i] = (char *)CHAR(STRING_ELT(cluster_chr_, i));

    int n_unique_chrs = 0;
    char **unique_chrs = NULL;
    for (int i = 0; i < n_exons_all; i++) {
        int found = 0;
        for (int j = 0; j < n_unique_chrs; j++)
            if (strcmp(exon_chr_c[i], unique_chrs[j]) == 0) { found = 1; break; }
        if (!found) {
            unique_chrs = (char **)realloc(unique_chrs, (n_unique_chrs+1) * sizeof(char *));
            unique_chrs[n_unique_chrs] = exon_chr_c[i];
            n_unique_chrs++;
        }
    }

    ChrAnnotNew *chr_annots = (ChrAnnotNew *)calloc(n_unique_chrs, sizeof(ChrAnnotNew));
    for (int ci = 0; ci < n_unique_chrs; ci++) {
        chr_annots[ci].name = unique_chrs[ci];

        int ne = 0, ni = 0, nc = 0;
        for (int i = 0; i < n_exons_all; i++) if (strcmp(exon_chr_c[i], unique_chrs[ci]) == 0) ne++;
        for (int i = 0; i < n_introns_all; i++) if (strcmp(intron_chr_c[i], unique_chrs[ci]) == 0) ni++;
        for (int i = 0; i < n_clusters_all; i++) if (strcmp(cluster_chr_c[i], unique_chrs[ci]) == 0) nc++;

        int *es=NULL,*ee=NULL,*eg=NULL; char **est=NULL;
        int *is=NULL,*ie=NULL,*ig=NULL; char **ist=NULL;
        int *cs=NULL,*ce=NULL,*cid=NULL; char **cst=NULL;
        if (ne > 0) { es=(int*)malloc(ne*sizeof(int)); ee=(int*)malloc(ne*sizeof(int)); eg=(int*)malloc(ne*sizeof(int)); est=(char**)malloc(ne*sizeof(char*)); }
        if (ni > 0) { is=(int*)malloc(ni*sizeof(int)); ie=(int*)malloc(ni*sizeof(int)); ig=(int*)malloc(ni*sizeof(int)); ist=(char**)malloc(ni*sizeof(char*)); }
        if (nc > 0) { cs=(int*)malloc(nc*sizeof(int)); ce=(int*)malloc(nc*sizeof(int)); cid=(int*)malloc(nc*sizeof(int)); cst=(char**)malloc(nc*sizeof(char*)); }

        int je=0,ji=0,jc=0;
        for (int i = 0; i < n_exons_all; i++) if (strcmp(exon_chr_c[i],unique_chrs[ci])==0) { es[je]=exon_starts_all[i]; ee[je]=exon_ends_all[i]; eg[je]=exon_gene_ids_all[i]; est[je]=exon_strand_c[i]; je++; }
        for (int i = 0; i < n_introns_all; i++) if (strcmp(intron_chr_c[i],unique_chrs[ci])==0) { is[ji]=intron_starts_all[i]; ie[ji]=intron_ends_all[i]; ig[ji]=intron_gene_ids_all[i]; ist[ji]=intron_strand_c[i]; ji++; }
        for (int i = 0; i < n_clusters_all; i++) if (strcmp(cluster_chr_c[i],unique_chrs[ci])==0) { cs[jc]=cluster_starts_all[i]; ce[jc]=cluster_ends_all[i]; cid[jc]=cluster_ids_all[i]; cst[jc]=cluster_strand_c[i]; jc++; }

        ExonClusterIndex exon_idx;
        build_exon_cluster_index(es, ee, eg, cid, ne, &exon_idx);

        GeneIndex intron_idx;
        build_gene_index(is, ie, ig, ni, &intron_idx);
        sort_gene_index_by_end(&intron_idx);

        split_exon_index_by_strand(&exon_idx, est, ne, &chr_annots[ci].fwd_exon, &chr_annots[ci].rev_exon);
        split_gene_index_by_strand(&intron_idx, ist, ni, &chr_annots[ci].fwd_gene, &chr_annots[ci].rev_gene);

        build_gene_span_c(es, ee, eg, est, ne, &chr_annots[ci].gene_span);

        chr_annots[ci].n_clusters = nc;
        if (nc > 0) {
            chr_annots[ci].cluster_start_sorted = (int *)malloc(nc * sizeof(int));
            chr_annots[ci].cluster_end_sorted = (int *)malloc(nc * sizeof(int));
            chr_annots[ci].cluster_orig_idx = (int *)malloc(nc * sizeof(int));
            chr_annots[ci].cluster_is_rev = (int8_t *)malloc(nc * sizeof(int8_t));
            for (int i = 0; i < nc; i++) {
                chr_annots[ci].cluster_start_sorted[i] = cs[i];
                chr_annots[ci].cluster_end_sorted[i] = ce[i];
                chr_annots[ci].cluster_orig_idx[i] = cid[i];
                chr_annots[ci].cluster_is_rev[i] = (cst[i][0] == '-') ? 1 : 0;
            }
            for (int i = 1; i < nc; i++) {
                int tmp_end = chr_annots[ci].cluster_end_sorted[i];
                int tmp_start = chr_annots[ci].cluster_start_sorted[i];
                int tmp_idx = chr_annots[ci].cluster_orig_idx[i];
                int8_t tmp_rev = chr_annots[ci].cluster_is_rev[i];
                int j = i - 1;
                while (j >= 0 && chr_annots[ci].cluster_end_sorted[j] > tmp_end) {
                    chr_annots[ci].cluster_end_sorted[j+1] = chr_annots[ci].cluster_end_sorted[j];
                    chr_annots[ci].cluster_start_sorted[j+1] = chr_annots[ci].cluster_start_sorted[j];
                    chr_annots[ci].cluster_orig_idx[j+1] = chr_annots[ci].cluster_orig_idx[j];
                    chr_annots[ci].cluster_is_rev[j+1] = chr_annots[ci].cluster_is_rev[j];
                    j--;
                }
                chr_annots[ci].cluster_end_sorted[j+1] = tmp_end;
                chr_annots[ci].cluster_start_sorted[j+1] = tmp_start;
                chr_annots[ci].cluster_orig_idx[j+1] = tmp_idx;
                chr_annots[ci].cluster_is_rev[j+1] = tmp_rev;
            }
        } else {
            chr_annots[ci].cluster_start_sorted = NULL;
            chr_annots[ci].cluster_end_sorted = NULL;
            chr_annots[ci].cluster_orig_idx = NULL;
            chr_annots[ci].cluster_is_rev = NULL;
        }

        free_exon_cluster_index(&exon_idx);
        free_gene_index(&intron_idx);
        free(es); free(ee); free(eg); free(est);
        free(is); free(ie); free(ig); free(ist);
        free(cs); free(ce); free(cid); free(cst);
    }

    int n_threads = omp_get_max_threads();
    if (n_threads > n_unique_chrs) n_threads = n_unique_chrs;

    long h_patterns[16] = {0};

    CBHash *thread_cb = (CBHash *)malloc(n_unique_chrs * sizeof(CBHash));
    TripletCOO *thread_spliced = (TripletCOO *)malloc(n_unique_chrs * sizeof(TripletCOO));
    TripletCOO *thread_unspliced = (TripletCOO *)malloc(n_unique_chrs * sizeof(TripletCOO));
    TripletCOO *thread_chimeric = (TripletCOO *)malloc(n_unique_chrs * sizeof(TripletCOO));
    TripletCOO *thread_ambiguous = (TripletCOO *)malloc(n_unique_chrs * sizeof(TripletCOO));
    long *thread_n_spliced = (long *)calloc(n_unique_chrs, sizeof(long));
    long *thread_n_unspliced = (long *)calloc(n_unique_chrs, sizeof(long));
    long *thread_n_chimeric = (long *)calloc(n_unique_chrs, sizeof(long));
    long *thread_n_ambiguous = (long *)calloc(n_unique_chrs, sizeof(long));
    long *thread_n_unassigned = (long *)calloc(n_unique_chrs, sizeof(long));
    long *thread_flag_H = (long *)calloc(n_unique_chrs, sizeof(long));
    long *thread_flag_U = (long *)calloc(n_unique_chrs, sizeof(long));
    long *thread_flag_G = (long *)calloc(n_unique_chrs, sizeof(long));
    long *thread_flag_E = (long *)calloc(n_unique_chrs, sizeof(long));

    for (int ci = 0; ci < n_unique_chrs; ci++) {
        cb_hash_init(&thread_cb[ci]);
        triplet_init(&thread_spliced[ci], 100000);
        triplet_init(&thread_unspliced[ci], 100000);
        triplet_init(&thread_chimeric[ci], 100000);
        triplet_init(&thread_ambiguous[ci], 100000);
    }

    #pragma omp parallel for schedule(dynamic, 1) num_threads(n_threads)
    for (int ci = 0; ci < n_unique_chrs; ci++) {
        ChrAnnotNew *ca = &chr_annots[ci];
        const char *chr_name = ca->name;

        htsFile *in = sam_open(bam_path, "r");
        if (!in) continue;

        bam_hdr_t *hdr = sam_hdr_read(in);
        if (!hdr) { sam_close(in); continue; }

        hts_idx_t *idx = sam_index_load(in, bam_path);
        if (!idx) { bam_hdr_destroy(hdr); sam_close(in); continue; }

        int tid = bam_name2id(hdr, chr_name);
        if (tid < 0) { hts_idx_destroy(idx); bam_hdr_destroy(hdr); sam_close(in); continue; }

        hts_itr_t *itr = sam_itr_queryi(idx, tid, 0, INT_MAX);

        CBHash *cb = &thread_cb[ci];
        TripletCOO *sp = &thread_spliced[ci];
        TripletCOO *us = &thread_unspliced[ci];
        TripletCOO *ch = &thread_chimeric[ci];
        TripletCOO *am = &thread_ambiguous[ci];
        long *n_sp = &thread_n_spliced[ci];
        long *n_us = &thread_n_unspliced[ci];
        long *n_ch = &thread_n_chimeric[ci];
        long *n_am = &thread_n_ambiguous[ci];
        long *n_un = &thread_n_unassigned[ci];
        long *f_H = &thread_flag_H[ci];
        long *f_U = &thread_flag_U[ci];
        long *f_G = &thread_flag_G[ci];
        long *f_E = &thread_flag_E[ci];

        bam1_t *b = bam_init1();
        int blk_start[64], blk_end[64];
        int exon_gene_ids_buf[32], exon_cluster_ids_buf[32];
        int block_gene_ids_buf[32], block_cluster_ids_buf[32];
        int gene_intersect_buf[32];
        int intron_gene_ids_buf[32];
        int cluster_hits_buf[64];

        while (sam_itr_next(in, itr, b) >= 0) {
            if (!read_passes_filter(b, tech_id, xf_allow, nxf)) continue;

            const char *CB = read_cell_id(b, tech_id);
            if (!CB) continue;
            int cell_id = cb_hash_get(cb, CB, strlen(CB));

            int nblk = extract_blocks_int(b, blk_start, blk_end, 64);
            if (nblk <= 0) continue;

            int is_reverse = (b->core.flag & BAM_FREVERSE) != 0;
            ExonClusterIndex *exon_idx = is_reverse ? &ca->rev_exon : &ca->fwd_exon;
            GeneIndex *intron_idx = is_reverse ? &ca->rev_gene : &ca->fwd_gene;
            SortedGene *gspans = is_reverse ? ca->gene_span.rev : ca->gene_span.fwd;
            int *gs_ends = is_reverse ? ca->gene_span.rev_ends : ca->gene_span.fwd_ends;
            int n_gs = is_reverse ? ca->gene_span.n_rev : ca->gene_span.n_fwd;

            int has_N = has_N_cigar(b);
            int flag = 0;
            int all_exonic = 1;
            int all_in_gene = 1;

            int cluster_hit_count = 0;

            if (has_N) {
                flag |= HUGE_H;

                int n_gene_set = find_overlapping_gene_ids_exon(exon_idx,
                    blk_start[0], blk_end[0],
                    exon_gene_ids_buf, exon_cluster_ids_buf, 32);

                if (n_gene_set > 0) {
                    for (int bi = 1; bi < nblk && n_gene_set > 0; bi++) {
                        int n_bg = find_overlapping_gene_ids_exon(exon_idx,
                            blk_start[bi], blk_end[bi],
                            block_gene_ids_buf, block_cluster_ids_buf, 32);
                        if (n_bg == 0) {
                            n_gene_set = -1;
                        } else {
                            memcpy(gene_intersect_buf, exon_gene_ids_buf, n_gene_set * sizeof(int));
                            n_gene_set = intersect_sorted_ints(gene_intersect_buf, n_gene_set, block_gene_ids_buf, n_bg);
                            memcpy(exon_gene_ids_buf, gene_intersect_buf, n_gene_set * sizeof(int));
                        }
                    }
                }

                if (n_gene_set > 0) {
                    flag |= HUGE_U;
                }

                for (int bi = 0; bi < nblk; bi++) {
                    int exon_hit = any_overlap_exon_index(exon_idx, blk_start[bi], blk_end[bi]);
                    int gene_hit = find_overlapping_gene_ids_gene(intron_idx, blk_start[bi], blk_end[bi], intron_gene_ids_buf, 32);
                    if (!exon_hit) all_exonic = 0;
                    if (!exon_hit && !gene_hit) all_in_gene = 0;
                }

                for (int bi = 0; bi < nblk; bi++) {
                    int nh = collect_cluster_overlaps(exon_idx, blk_start[bi], blk_end[bi],
                                                      cluster_hits_buf + cluster_hit_count,
                                                      64 - cluster_hit_count);
                    cluster_hit_count += nh;
                }

                for (int a = 0; a < cluster_hit_count; a++) {
                    for (int b2 = a + 1; b2 < cluster_hit_count; b2++) {
                        if (cluster_hits_buf[a] == cluster_hits_buf[b2]) {
                            cluster_hits_buf[b2] = cluster_hits_buf[cluster_hit_count - 1];
                            cluster_hit_count--;
                            b2--;
                        }
                    }
                }

            } else {
                int read_start = blk_start[0];
                int read_end = blk_end[nblk - 1];

                int in_gene_span = 0;
                if (n_gs > 0) {
                    int start_idx = find_first_overlap_end(gs_ends, n_gs, read_start - gene_span_pad);
                    for (int i = start_idx; i < n_gs; i++) {
                        if (gspans[i].start > read_end + gene_span_pad) break;
                        if (gspans[i].end >= read_start - gene_span_pad && gspans[i].start <= read_end + gene_span_pad) {
                            in_gene_span = 1;
                            break;
                        }
                    }
                }

                for (int bi = 0; bi < nblk; bi++) {
                    int exon_hit = any_overlap_exon_index(exon_idx, blk_start[bi], blk_end[bi]);
                    int gene_hit = find_overlapping_gene_ids_gene(intron_idx, blk_start[bi], blk_end[bi], intron_gene_ids_buf, 32);
                    if (!exon_hit) all_exonic = 0;
                    if (!exon_hit && !gene_hit) all_in_gene = 0;
                }

                if (in_gene_span) {
                    flag |= HUGE_U;
                }

                for (int bi = 0; bi < nblk; bi++) {
                    int nh = collect_cluster_overlaps(exon_idx, blk_start[bi], blk_end[bi],
                                                      cluster_hits_buf + cluster_hit_count,
                                                      64 - cluster_hit_count);
                    cluster_hit_count += nh;
                }

                for (int a = 0; a < cluster_hit_count; a++) {
                    for (int b2 = a + 1; b2 < cluster_hit_count; b2++) {
                        if (cluster_hits_buf[a] == cluster_hits_buf[b2]) {
                            cluster_hits_buf[b2] = cluster_hits_buf[cluster_hit_count - 1];
                            cluster_hit_count--;
                            b2--;
                        }
                    }
                }
            }

            if (all_in_gene) flag |= HUGE_G;
            if (all_exonic) flag |= HUGE_E;

            if (flag & HUGE_H) (*f_H)++;
            if (flag & HUGE_U) (*f_U)++;
            if (flag & HUGE_G) (*f_G)++;
            if (flag & HUGE_E) (*f_E)++;

            ReadClass category;
            if ((flag & HUGE_H) && (flag & HUGE_U) && (flag & HUGE_G) && (flag & HUGE_E)) {
                category = READ_SPLICED;
            } else if ((flag & HUGE_H) && !(flag & HUGE_U) && (flag & HUGE_G) && (flag & HUGE_E)) {
                category = READ_CHIMERIC;
            } else if ((flag & HUGE_U) && (flag & HUGE_G) && !(flag & HUGE_E)) {
                category = READ_UNSPLICED;
            } else if (!(flag & HUGE_H) && (flag & HUGE_U) && (flag & HUGE_G) && (flag & HUGE_E)) {
                category = READ_AMBIGUOUS;
            } else {
                category = READ_UNASSIGNED;
            }

            if (category != READ_UNASSIGNED) {
                TripletCOO *coo;
                switch (category) {
                    case READ_SPLICED:   coo = sp; (*n_sp)++; break;
                    case READ_CHIMERIC:  coo = ch; (*n_ch)++; break;
                    case READ_UNSPLICED: coo = us; (*n_us)++; break;
                    case READ_AMBIGUOUS: coo = am; (*n_am)++; break;
                    default: coo = NULL; break;
                }
                if (coo) {
                    for (int hi = 0; hi < cluster_hit_count; hi++) {
                        triplet_add(coo, cell_id, cluster_hits_buf[hi]);
                    }
                }
            } else {
                (*n_un)++;
            }

                {
                    int pidx = ((flag & HUGE_H) ? 8 : 0) | ((flag & HUGE_U) ? 4 : 0) | ((flag & HUGE_G) ? 2 : 0) | ((flag & HUGE_E) ? 1 : 0);
                    #pragma omp atomic
                    h_patterns[pidx]++;
                }
        }

        bam_destroy1(b);
        hts_itr_destroy(itr);
        hts_idx_destroy(idx);
        bam_hdr_destroy(hdr);
        sam_close(in);
    }

    CBHash global_cb;
    cb_hash_init(&global_cb);

    for (int ci = 0; ci < n_unique_chrs; ci++) {
        int max_local_id = 0;
        for (int i = 0; i < CB_HASH_SIZE; i++) {
            CBEntry *e = thread_cb[ci].buckets[i];
            while (e) { if (e->cell_id > max_local_id) max_local_id = e->cell_id; e = e->next; }
        }
        int *remap = NULL;
        if (max_local_id > 0) {
            remap = (int *)malloc((max_local_id + 1) * sizeof(int));
            memset(remap, 0, (max_local_id + 1) * sizeof(int));
        }
        for (int i = 0; i < CB_HASH_SIZE; i++) {
            CBEntry *e = thread_cb[ci].buckets[i];
            while (e) {
                CBEntry *next = e->next;
                int global_id = cb_hash_get(&global_cb, e->key, e->key_len);
                remap[e->cell_id] = global_id;
                free(e->key); free(e);
                thread_cb[ci].buckets[i] = next;
                e = next;
            }
        }
        if (remap) {
            TripletCOO *t_arr[] = {&thread_spliced[ci], &thread_unspliced[ci], &thread_chimeric[ci], &thread_ambiguous[ci]};
            for (int cat = 0; cat < 3; cat++) {
                TripletCOO *t = t_arr[cat];
                for (int k = 0; k < t->size; k++) t->cell_idx[k] = remap[t->cell_idx[k]];
            }
            free(remap);
        }
    }

    TripletCOO global_spliced, global_unspliced, global_chimeric, global_ambiguous;
    triplet_init(&global_spliced, 1000000);
    triplet_init(&global_unspliced, 1000000);
    triplet_init(&global_chimeric, 1000000);
    triplet_init(&global_ambiguous, 1000000);

    for (int ci = 0; ci < n_unique_chrs; ci++) {
        TripletCOO *srcs[] = {&thread_spliced[ci], &thread_unspliced[ci], &thread_chimeric[ci], &thread_ambiguous[ci]};
        TripletCOO *dsts[] = {&global_spliced, &global_unspliced, &global_chimeric, &global_ambiguous};
        for (int cat = 0; cat < 3; cat++) {
            for (int k = 0; k < srcs[cat]->size; k++) {
                triplet_add_count(dsts[cat], srcs[cat]->cell_idx[k], srcs[cat]->cluster_idx[k], srcs[cat]->counts[k]);
            }
        }
    }

    long total_spliced = 0, total_unspliced = 0, total_chimeric = 0, total_ambiguous = 0, total_unassigned = 0;
    long total_flag_H = 0, total_flag_U = 0, total_flag_G = 0, total_flag_E = 0;
    for (int ci = 0; ci < n_unique_chrs; ci++) {
        total_spliced += thread_n_spliced[ci];
        total_unspliced += thread_n_unspliced[ci];
        total_chimeric += thread_n_chimeric[ci];
        total_ambiguous += thread_n_ambiguous[ci];
        total_unassigned += thread_n_unassigned[ci];
        total_flag_H += thread_flag_H[ci];
        total_flag_U += thread_flag_U[ci];
        total_flag_G += thread_flag_G[ci];
        total_flag_E += thread_flag_E[ci];
    }

    for (int ci = 0; ci < n_unique_chrs; ci++) {
        triplet_free(&thread_spliced[ci]);
        triplet_free(&thread_unspliced[ci]);
        triplet_free(&thread_chimeric[ci]);
        triplet_free(&thread_ambiguous[ci]);
    }
    free(thread_spliced); free(thread_unspliced); free(thread_chimeric); free(thread_ambiguous);
    free(thread_cb);

    SEXP chr_names = PROTECT(allocVector(STRSXP, n_unique_chrs));
    SEXP chr_spliced = PROTECT(allocVector(INTSXP, n_unique_chrs));
    SEXP chr_unspliced = PROTECT(allocVector(INTSXP, n_unique_chrs));
    SEXP chr_chimeric = PROTECT(allocVector(INTSXP, n_unique_chrs));
    SEXP chr_ambiguous = PROTECT(allocVector(INTSXP, n_unique_chrs));
    SEXP chr_unassigned = PROTECT(allocVector(INTSXP, n_unique_chrs));
    for (int ci = 0; ci < n_unique_chrs; ci++) {
        SET_STRING_ELT(chr_names, ci, mkChar(unique_chrs[ci]));
        INTEGER(chr_spliced)[ci] = (int)thread_n_spliced[ci];
        INTEGER(chr_unspliced)[ci] = (int)thread_n_unspliced[ci];
        INTEGER(chr_chimeric)[ci] = (int)thread_n_chimeric[ci];
        INTEGER(chr_ambiguous)[ci] = (int)thread_n_ambiguous[ci];
        INTEGER(chr_unassigned)[ci] = (int)thread_n_unassigned[ci];
    }

    for (int ci = 0; ci < n_unique_chrs; ci++) {
        free_exon_cluster_index(&chr_annots[ci].fwd_exon);
        free_exon_cluster_index(&chr_annots[ci].rev_exon);
        free_gene_index(&chr_annots[ci].fwd_gene);
        free_gene_index(&chr_annots[ci].rev_gene);
        free_strand_gene_span(&chr_annots[ci].gene_span);
        if (chr_annots[ci].n_clusters > 0) {
            free(chr_annots[ci].cluster_start_sorted);
            free(chr_annots[ci].cluster_end_sorted);
            free(chr_annots[ci].cluster_orig_idx);
            free(chr_annots[ci].cluster_is_rev);
        }
    }
    free(chr_annots);
    free(unique_chrs);
    free(exon_chr_c); free(intron_chr_c); free(cluster_chr_c);
    free(exon_strand_c); free(intron_strand_c); free(cluster_strand_c);
    free(thread_n_spliced); free(thread_n_unspliced);
    free(thread_n_chimeric); free(thread_n_ambiguous); free(thread_n_unassigned);

    SEXP chr_flag_H = PROTECT(allocVector(INTSXP, n_unique_chrs));
    SEXP chr_flag_U = PROTECT(allocVector(INTSXP, n_unique_chrs));
    SEXP chr_flag_G = PROTECT(allocVector(INTSXP, n_unique_chrs));
    SEXP chr_flag_E = PROTECT(allocVector(INTSXP, n_unique_chrs));
    for (int ci = 0; ci < n_unique_chrs; ci++) {
        INTEGER(chr_flag_H)[ci] = (int)thread_flag_H[ci];
        INTEGER(chr_flag_U)[ci] = (int)thread_flag_U[ci];
        INTEGER(chr_flag_G)[ci] = (int)thread_flag_G[ci];
        INTEGER(chr_flag_E)[ci] = (int)thread_flag_E[ci];
    }

    free(thread_flag_H); free(thread_flag_U); free(thread_flag_G); free(thread_flag_E);

    UNPROTECT(13);

    SEXP chr_df = PROTECT(allocVector(VECSXP, 10));
    SEXP chr_df_names = PROTECT(allocVector(STRSXP, 10));
    SET_VECTOR_ELT(chr_df, 0, chr_names);
    SET_VECTOR_ELT(chr_df, 1, chr_spliced);
    SET_VECTOR_ELT(chr_df, 2, chr_unspliced);
    SET_VECTOR_ELT(chr_df, 3, chr_chimeric);
    SET_VECTOR_ELT(chr_df, 4, chr_ambiguous);
    SET_VECTOR_ELT(chr_df, 5, chr_unassigned);
    SET_VECTOR_ELT(chr_df, 6, chr_flag_H);
    SET_VECTOR_ELT(chr_df, 7, chr_flag_U);
    SET_VECTOR_ELT(chr_df, 8, chr_flag_G);
    SET_VECTOR_ELT(chr_df, 9, chr_flag_E);
    SET_STRING_ELT(chr_df_names, 0, mkChar("chr"));
    SET_STRING_ELT(chr_df_names, 1, mkChar("spliced"));
    SET_STRING_ELT(chr_df_names, 2, mkChar("unspliced"));
    SET_STRING_ELT(chr_df_names, 3, mkChar("chimeric"));
    SET_STRING_ELT(chr_df_names, 4, mkChar("ambiguous"));
    SET_STRING_ELT(chr_df_names, 5, mkChar("unassigned"));
    SET_STRING_ELT(chr_df_names, 6, mkChar("flag_H"));
    SET_STRING_ELT(chr_df_names, 7, mkChar("flag_U"));
    SET_STRING_ELT(chr_df_names, 8, mkChar("flag_G"));
    SET_STRING_ELT(chr_df_names, 9, mkChar("flag_E"));
    setAttrib(chr_df, R_NamesSymbol, chr_df_names);
    SEXP chr_df_class = PROTECT(allocVector(STRSXP, 1));
    SET_STRING_ELT(chr_df_class, 0, mkChar("data.frame"));
    setAttrib(chr_df, R_ClassSymbol, chr_df_class);
    SEXP row_names = PROTECT(allocVector(INTSXP, n_unique_chrs));
    for (int ci = 0; ci < n_unique_chrs; ci++) INTEGER(row_names)[ci] = ci + 1;
    setAttrib(chr_df, R_RowNamesSymbol, row_names);

    write_category_output(out_dir, "spliced", &global_spliced, &global_cb,
                          n_clusters_all, cluster_starts_all, cluster_ends_all);
    write_category_output(out_dir, "unspliced", &global_unspliced, &global_cb,
                          n_clusters_all, cluster_starts_all, cluster_ends_all);
    write_category_output(out_dir, "chimeric", &global_chimeric, &global_cb,
                          n_clusters_all, cluster_starts_all, cluster_ends_all);
    write_category_output(out_dir, "ambiguous", &global_ambiguous, &global_cb,
                          n_clusters_all, cluster_starts_all, cluster_ends_all);

    triplet_free(&global_spliced);
    triplet_free(&global_unspliced);
    triplet_free(&global_chimeric);
    triplet_free(&global_ambiguous);
    cb_hash_free(&global_cb);

    SEXP result = PROTECT(allocVector(VECSXP, 14));
    SEXP names = PROTECT(allocVector(STRSXP, 14));
    SET_VECTOR_ELT(result, 0, ScalarInteger(global_cb.count));
    SET_VECTOR_ELT(result, 1, ScalarInteger(n_clusters_all));
    SET_VECTOR_ELT(result, 2, ScalarInteger((int)total_spliced));
    SET_VECTOR_ELT(result, 3, ScalarInteger((int)total_unspliced));
    SET_VECTOR_ELT(result, 4, ScalarInteger((int)total_chimeric));
    SET_VECTOR_ELT(result, 5, ScalarInteger((int)total_ambiguous));
    SET_VECTOR_ELT(result, 6, ScalarInteger((int)total_unassigned));
    SET_VECTOR_ELT(result, 7, mkString(out_dir));
    SET_VECTOR_ELT(result, 8, chr_df);
    SET_VECTOR_ELT(result, 9, ScalarInteger((int)total_flag_H));
    SET_VECTOR_ELT(result, 10, ScalarInteger((int)total_flag_U));
    SET_VECTOR_ELT(result, 11, ScalarInteger((int)total_flag_G));
    SET_VECTOR_ELT(result, 12, ScalarInteger((int)total_flag_E));

    SEXP h_patt = PROTECT(allocVector(INTSXP, 16));
    for (int i = 0; i < 16; i++) INTEGER(h_patt)[i] = (int)h_patterns[i];
    SET_VECTOR_ELT(result, 13, h_patt);

    SET_STRING_ELT(names, 0, mkChar("n_cells"));
    SET_STRING_ELT(names, 1, mkChar("n_features"));
    SET_STRING_ELT(names, 2, mkChar("n_spliced"));
    SET_STRING_ELT(names, 3, mkChar("n_unspliced"));
    SET_STRING_ELT(names, 4, mkChar("n_chimeric"));
    SET_STRING_ELT(names, 5, mkChar("n_ambiguous"));
    SET_STRING_ELT(names, 6, mkChar("n_unassigned"));
    SET_STRING_ELT(names, 7, mkChar("out_dir"));
    SET_STRING_ELT(names, 8, mkChar("by_chr"));
    SET_STRING_ELT(names, 9, mkChar("flag_H"));
    SET_STRING_ELT(names, 10, mkChar("flag_U"));
    SET_STRING_ELT(names, 11, mkChar("flag_G"));
    SET_STRING_ELT(names, 12, mkChar("flag_E"));
    SET_STRING_ELT(names, 13, mkChar("h_patterns"));

    setAttrib(result, R_NamesSymbol, names);
    UNPROTECT(7);
    return result;
}