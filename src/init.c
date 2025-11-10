#include <R.h>
#include <Rinternals.h>

SEXP _exonBlocks_scan_bam_blocks_hts(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_exonBlocks_scan_bam_blocks_hts", (DL_FUNC) &_exonBlocks_scan_bam_blocks_hts, 7},
  {NULL, NULL, 0}
};

void R_init_exonBlocks(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
