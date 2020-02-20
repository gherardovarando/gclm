#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(dgelyp)(void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(grddsllc)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(pnllbc)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(prxgrdllb)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(prxgrdlsb)(void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
    {"dgelyp",    (DL_FUNC) &F77_NAME(dgelyp),     8},
    {"grddsllc",  (DL_FUNC) &F77_NAME(grddsllc),  12},
    {"pnllbc",    (DL_FUNC) &F77_NAME(pnllbc),    12},
    {"prxgrdllb", (DL_FUNC) &F77_NAME(prxgrdllb),  10},
    {"prxgrdlsb", (DL_FUNC) &F77_NAME(prxgrdlsb),  9},
    {NULL, NULL, 0}
};

void R_init_clggm(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}