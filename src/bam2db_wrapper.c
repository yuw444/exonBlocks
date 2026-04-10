#include <R.h>
#include <Rinternals.h>
#include "bam2db.h"

SEXP _exonBlocks_bam2db(
    SEXP bam_file_,
    SEXP db_file_,
    SEXP path_out_,
    SEXP barcodes_file_,
    SEXP regions_file_,
    SEXP rate_cell_,
    SEXP rate_depth_,
    SEXP seed_,
    SEXP umi_copies_flag_,
    SEXP protocol_,
    SEXP region_chr_,
    SEXP region_start_,
    SEXP region_end_)
{
    const char *bam_file = CHAR(STRING_ELT(bam_file_, 0));
    const char *db_file = CHAR(STRING_ELT(db_file_, 0));
    const char *path_out = CHAR(STRING_ELT(path_out_, 0));
    const char *barcodes_file = CHAR(STRING_ELT(barcodes_file_, 0));
    const char *regions_file = CHAR(STRING_ELT(regions_file_, 0));
    float rate_cell = (float)REAL(rate_cell_)[0];
    float rate_depth = (float)REAL(rate_depth_)[0];
    unsigned int seed = (unsigned int)INTEGER(seed_)[0];
    int umi_copies_flag = INTEGER(umi_copies_flag_)[0];
    protocol_t protocol = (protocol_t)INTEGER(protocol_)[0];

    const char *region_chr = NULL;
    hts_pos_t region_start = 0, region_end = 0;
    if (LENGTH(region_chr_) > 0 && STRING_ELT(region_chr_, 0) != NA_STRING) {
        region_chr = CHAR(STRING_ELT(region_chr_, 0));
        region_start = (hts_pos_t)INTEGER(region_start_)[0];
        region_end = (hts_pos_t)INTEGER(region_end_)[0];
    }

    int ret = bam2db(
        (char *)bam_file, (char *)db_file, (char *)path_out,
        (char *)barcodes_file, (char *)regions_file,
        rate_cell, rate_depth, seed,
        umi_copies_flag, protocol,
        region_chr, region_start, region_end);

    if (ret != 0)
        error("bam2db failed with exit code %d", ret);

    return R_NilValue;
}
