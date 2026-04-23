#include <R.h>
#include <Rinternals.h>

SEXP _exonBlocks_scan_bam_blocks_hts(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP _exonBlocks_bam_to_cellranger(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP _exonBlocks_extract_unique_tags(SEXP, SEXP, SEXP);
SEXP _exonBlocks_bam2db(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP _exonBlocks_classify_reads(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_exonBlocks_scan_bam_blocks_hts", (DL_FUNC)&_exonBlocks_scan_bam_blocks_hts, 8},
    {"_exonBlocks_bam_to_cellranger", (DL_FUNC)&_exonBlocks_bam_to_cellranger, 8},
    {"_exonBlocks_extract_unique_tags", (DL_FUNC)&_exonBlocks_extract_unique_tags, 3},
    {"_exonBlocks_bam2db", (DL_FUNC)&_exonBlocks_bam2db, 13},
    {"_exonBlocks_classify_reads", (DL_FUNC)&_exonBlocks_classify_reads, 20},
    {NULL, NULL, 0}};

void R_init_exonBlocks(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
