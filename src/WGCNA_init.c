#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void bicor1Fast(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void bicorFast(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void checkAvailableMemoryForR(void *);
extern void corFast(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void mean(void *, void *, void *, void *);
extern void minWhichMin(void *, void *, void *, void *, void *);
extern void minWhichMin_row(void *, void *, void *, void *, void *);
extern void quantileC(void *, void *, void *, void *, void *);
extern void rowQuantileC(void *, void *, void *, void *, void *);
extern void tomSimilarityFromAdj(void *, void *, void *, void *, void *, void *, void *);

/* .Call calls */
extern SEXP bucketOrder_R(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP cor1Fast_call(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP qorder(SEXP);
extern SEXP tomSimilarity_call(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
    {"bicor1Fast",               (DL_FUNC) &bicor1Fast,               14},
    {"bicorFast",                (DL_FUNC) &bicorFast,                20},
    {"checkAvailableMemoryForR", (DL_FUNC) &checkAvailableMemoryForR,  1},
    {"corFast",                  (DL_FUNC) &corFast,                  14},
    {"mean",                     (DL_FUNC) &mean,                      4},
    {"minWhichMin",              (DL_FUNC) &minWhichMin,               5},
    {"minWhichMin_row",          (DL_FUNC) &minWhichMin_row,           5},
    {"quantileC",                (DL_FUNC) &quantileC,                 5},
    {"rowQuantileC",             (DL_FUNC) &rowQuantileC,              5},
    {"tomSimilarityFromAdj",     (DL_FUNC) &tomSimilarityFromAdj,      7},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"bucketOrder_R",      (DL_FUNC) &bucketOrder_R,       5},
    {"cor1Fast_call",      (DL_FUNC) &cor1Fast_call,       9},
    {"qorder",             (DL_FUNC) &qorder,              1},
    {"tomSimilarity_call", (DL_FUNC) &tomSimilarity_call, 15},
    {NULL, NULL, 0}
};

// void R_init_WGCNA(DllInfo *dll)
// {
//     R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
//     R_useDynamicSymbols(dll, FALSE);
// }
