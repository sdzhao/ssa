#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP fsdr(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP fsdr_all(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP fsdr_p(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ldd(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"fsdr",     (DL_FUNC) &fsdr,     7},
    {"fsdr_all", (DL_FUNC) &fsdr_all, 5},
    {"fsdr_p",   (DL_FUNC) &fsdr_p,   6},
    {"ldd",      (DL_FUNC) &ldd,      4},
    {NULL, NULL, 0}
};

void R_init_ssa(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
