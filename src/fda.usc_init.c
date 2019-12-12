#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(adot)(void *, void *, void *);
extern void F77_NAME(pcvm_statistic)(void *, void *, void *, void *);
extern void F77_NAME(rp_stat)(void *, void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
    {"adot",           (DL_FUNC) &F77_NAME(adot),           3},
    {"pcvm_statistic", (DL_FUNC) &F77_NAME(pcvm_statistic), 4},
    {"rp_stat",        (DL_FUNC) &F77_NAME(rp_stat),        5},
    {NULL, NULL, 0}
};

void R_init_fda_usc(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
