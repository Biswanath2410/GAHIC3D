#include <R.h>
#include <Rdefines.h>
#include <R_ext/Error.h>
#include "000.api.h"

#define min(x, y) (((x) < (y)) ? (x) : (y))

SEXP upper_triangle(SEXP X, SEXP diag) {
  SEXP res;

  double *x = REAL(X);
  int d = LOGICAL(diag)[0];

  SEXP dim = getAttrib(X, R_DimSymbol);
  int nrow = INTEGER(dim)[0];
  int ncol = INTEGER(dim)[1];

  R_xlen_t ii = 0;
  if (d) {
    for (int cc = 0; cc < ncol; cc++) {
      for (int rr = 0; rr <= min(cc, nrow - 1); rr++) ++ii;
    }
  } else {
    for (int cc = 0; cc < ncol; cc++) {
      for (int rr = 0; rr < min(cc, nrow); rr++) ++ii;
    }
  }

  PROTECT(res = allocVector(REALSXP, ii));
  double *y = REAL(res);

  ii = 0;
  if (d) {
    for (int cc = 0; cc < ncol; cc++) {
      for (int rr = 0; rr <= min(cc, nrow - 1); rr++) {
	y[ii++] = x[cc * nrow + rr];
      }
    }
  } else {
    for (int cc = 0; cc < ncol; cc++) {
      for (int rr = 0; rr < min(cc, nrow); rr++) {
        y[ii++] = x[cc * nrow + rr];
      }
    }
  }

  UNPROTECT(1);
  return res;
}

