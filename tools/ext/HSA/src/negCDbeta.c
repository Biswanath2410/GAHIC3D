#include <R.h>
#include <Rdefines.h>
#include <R_ext/Error.h>
#include "000.api.h"

SEXP negCDbeta(SEXP C, SEXP D, SEXP beta) {
  SEXP res;

  double *c = REAL(C);
  double *d = REAL(D);
  double b = REAL(beta)[0];

  SEXP dim = getAttrib(C, R_DimSymbol);
  int nrow = INTEGER(dim)[0];
  int ncol = INTEGER(dim)[1];

  dim = getAttrib(D, R_DimSymbol);
  if (INTEGER(dim)[0] != nrow)
    error("Argument \'C\' and \'D\' have different number of rows");
  if (INTEGER(dim)[1] != ncol)
    error("Argument \'C\' and \'D\' have different number of columns");

  PROTECT(res = allocMatrix(REALSXP, nrow, ncol));
  double *y = REAL(res);

  R_xlen_t ii = 0;
  for (int cc = 0; cc < ncol; ++cc) {
    for (int rr = 0; rr < nrow; ++rr) {
      y[ii] = - c[ii] * pow(d[ii], b);
      ii++;
    }
  }

  UNPROTECT(1);
  return res;
}

