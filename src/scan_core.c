#include <R.h>
#include <Rinternals.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/kstring.h>

// Convert BAM query sequence to char at qpos
static const char nt16_table[16] = "=-ACMGRSVTWYHKDBN";

/* Return the nucleotide character at query position qpos */
static inline char base_at(const uint8_t *s, int qpos)
{
    return nt16_table[bam_seqi(s, qpos)];
}

// Build semicolon-joined blocks from one alignment:
// - reference blocks for ops M/= /X
// - corresponding query slices (skip I/S in query cursor)
static int write_blocks_one(const bam1_t *b, kstring_t *ks_start,
                            kstring_t *ks_end, kstring_t *ks_seq)
{
    uint32_t *cig = bam_get_cigar(b);
    int n_cig = b->core.n_cigar;
    int32_t ref_pos = b->core.pos + 1; // 1-based
    int qpos = 0;                      // query cursor (0-based on query sequence)
    const uint8_t *qseq = bam_get_seq(b);

    int nblocks = 0;
    int first = 1;
    for (int i = 0; i < n_cig; ++i)
    {
        int op = bam_cigar_op(cig[i]);
        int len = bam_cigar_oplen(cig[i]);
        switch (op)
        {
        case BAM_CMATCH:
        case BAM_CEQUAL:
        case BAM_CDIFF:
        {
            // reference block
            int rs = ref_pos;
            int re = ref_pos + len - 1;
            if (!first)
            {
                kputc(';', ks_start);
                kputc(';', ks_end);
                kputc(';', ks_seq);
            }
            kputw(rs, ks_start);
            kputw(re, ks_end);

            // append query slice for this block
            for (int j = 0; j < len; ++j)
            {
                kputc(base_at(qseq, qpos + j), ks_seq);
            }

            ref_pos += len;
            qpos += len;
            nblocks++;
            first = 0;
            break;
        }
        case BAM_CINS:
            qpos += len;
            break; // consumes query
        case BAM_CSOFT_CLIP:
            qpos += len;
            break; // consumes query
        case BAM_CDEL:
            ref_pos += len;
            break; // consumes ref
        case BAM_CREF_SKIP:
            ref_pos += len;
            break;           // N (intron)
        case BAM_CHARD_CLIP: /*fallthrough*/
        case BAM_CPAD:       /* neither */
            break;
        default:
            break;
        }
    }
    return nblocks;
}

// Fetch Z (string) tag or return empty string
static const char *get_tagZ(const bam1_t *b, const char *tag)
{
    uint8_t *p = bam_aux_get(b, tag);
    if (!p)
        return NULL;
    return bam_aux2Z(p);
}

// Fetch integer tag (i/I types), return 0 and flag via ok*
static long get_tagi(const bam1_t *b, const char *tag, int *ok)
{
    *ok = 0;
    uint8_t *p = bam_aux_get(b, tag);
    if (!p)
        return 0;
    long v = bam_aux2i(p);
    *ok = 1;
    return v;
}

// [[.Call]] entry: scans BAM over [chr:start-end], filters, writes out.
// Returns integer(3): c(n_reads, n_umis, n_cells)
SEXP _exonBlocks_scan_bam_blocks_hts(
    SEXP bam_,
    SEXP chr_,
    SEXP start_,
    SEXP end_,
    SEXP out_bam_,
    SEXP tsv_,
    SEXP xf_vals_)
{
    const char *bam_path = CHAR(STRING_ELT(bam_, 0));
    const char *chr = CHAR(STRING_ELT(chr_, 0));
    int start = INTEGER(start_)[0];
    int end = INTEGER(end_)[0];
    const char *out_bam = CHAR(STRING_ELT(out_bam_, 0)); // "" -> skip
    const char *tsv = CHAR(STRING_ELT(tsv_, 0));

    // Build a quick lookup for xf filter
    int nxf = LENGTH(xf_vals_);
    long xf_allow[16]; // small set expected
    for (int i = 0; i < nxf && i < 16; ++i)
        xf_allow[i] = (long)INTEGER(xf_vals_)[i];

    // open BAM/CRAM
    htsFile *in = sam_open(bam_path, "r");
    if (!in)
        error("Failed to open input: %s", bam_path);
    bam_hdr_t *hdr = sam_hdr_read(in);
    if (!hdr)
    {
        sam_close(in);
        error("Failed to read header");
    }

    // create iterator over region (note: iterator uses 0-based, half-open)
    hts_idx_t *idx = sam_index_load(in, bam_path);
    if (!idx)
    {
        bam_hdr_destroy(hdr);
        sam_close(in);
        error("No index for input");
    }
    int tid = bam_name2id(hdr, chr);
    if (tid < 0)
    {
        hts_idx_destroy(idx);
        bam_hdr_destroy(hdr);
        sam_close(in);
        error("Unknown contig: %s", chr);
    }
    hts_itr_t *itr = sam_itr_queryi(idx, tid, start - 1, end); // 0-based, half-open

    // Optionally open BAM writer (compressed BAM)
    htsFile *out = NULL;
    if (out_bam && out_bam[0] != '\0')
    {
        out = sam_open(out_bam, "wb");
        if (!out)
        {
            hts_itr_destroy(itr);
            hts_idx_destroy(idx);
            bam_hdr_destroy(hdr);
            sam_close(in);
            error("Failed to open out_bam");
        }
        if (sam_hdr_write(out, hdr) < 0)
        {
            sam_close(out);
            out = NULL;
        }
    }

    // TSV
    FILE *fp = fopen(tsv, "wb");
    if (!fp)
    {
        if (out)
            sam_close(out);
        hts_itr_destroy(itr);
        hts_idx_destroy(idx);
        bam_hdr_destroy(hdr);
        sam_close(in);
        error("Failed to open TSV for write: %s", tsv);
    }
    fputs("CB\tUMI\tblock_start\tblock_end\tblock_seq\tnum_blocks\n", fp);

    // Simple hash sets for UMIs and CBs? To keep code short, we’ll count unique later in R,
    // but we can do a tiny on-the-fly linear cache (works well for 10x gene-window scans).
    // For robustness here, just count reads; return 0 for uniqs and let R recompute if needed.
    long n_reads = 0;

    bam1_t *b = bam_init1();
    kstring_t ks_start = {0, 0, NULL}, ks_end = {0, 0, NULL}, ks_seq = {0, 0, NULL};

    while (sam_itr_next(in, itr, b) >= 0)
    {
        // flag filter: unmapped/secondary/supplementary dropped
        uint16_t flag = b->core.flag;
        if ((flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) != 0)
            continue;

        // xf filter
        int okxf = 0, tagok = 0;
        long xf = get_tagi(b, "xf", &tagok);
        if (!tagok)
            continue;
        for (int i = 0; i < nxf; ++i)
            if (xf == xf_allow[i])
            {
                okxf = 1;
                break;
            }
        if (!okxf)
            continue;

        // CB / UB must exist
        const char *CB = get_tagZ(b, "CB");
        const char *UB = get_tagZ(b, "UB");
        if (!CB || !UB)
            continue;

        // Build blocks
        ks_start.l = ks_end.l = ks_seq.l = 0;
        int nblk = write_blocks_one(b, &ks_start, &ks_end, &ks_seq);
        if (nblk <= 0)
            continue;

        // Write TSV row
        // CB\tUMI\tstart\tend\tseq\tnumb_blocks\n
        fputs(CB, fp);
        fputc('\t', fp);
        fputs(UB, fp);
        fputc('\t', fp);
        fwrite(ks_start.s, 1, ks_start.l, fp);
        fputc('\t', fp);
        fwrite(ks_end.s, 1, ks_end.l, fp);
        fputc('\t', fp);
        fwrite(ks_seq.s, 1, ks_seq.l, fp);
        fputc('\t', fp);
        fprintf(fp, "%d\n", nblk);

        // Write to filtered BAM if requested
        if (out)
            if (sam_write1(out, hdr, b) < 0) {
                warning("Failed to write alignment to out_bam: %s", out_bam);
                /* Close and disable further writes to avoid repeated errors */
                sam_close(out);
                out = NULL;
            }

        ++n_reads;
    }

    // Cleanup
    if (ks_start.s)
        free(ks_start.s);
    if (ks_end.s)
        free(ks_end.s);
    if (ks_seq.s)
        free(ks_seq.s);

    bam_destroy1(b);
    fclose(fp);

    if (out)
    {
        sam_close(out);
        // Build index for out_bam (BAI). 0 -> default nthreads/format
        {
            /* Only try to build a BAI for BAM output files (skip for CRAM/other formats). */
            const char *dot = strrchr(out_bam, '.');
            if (dot && strcmp(dot, ".bam") == 0)
            {
            int rc = sam_index_build3(out_bam, NULL, 0, 0);
            if (rc != 0)
                warning("Failed to build BAI for %s (sam_index_build3 rc=%d)", out_bam, rc);
            }
        }
    }

    hts_itr_destroy(itr);
    hts_idx_destroy(idx);
    bam_hdr_destroy(hdr);
    sam_close(in);

    SEXP ans = PROTECT(allocVector(INTSXP, 1));
    INTEGER(ans)[0] = (int)n_reads;
    UNPROTECT(1);
    return ans;
}
