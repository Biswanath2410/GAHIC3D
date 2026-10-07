#include <R.h>
#include <Rdefines.h>
#include <R_ext/Error.h>
#include "000.api.h"

SEXP dist_matrix(SEXP X, SEXP size, SEXP square) {
  SEXP res;

  double *x = REAL(X);
  int sq = LOGICAL(square)[0];
  int sz = INTEGER(size)[0];
  
  PROTECT(res = allocMatrix(REALSXP, sz, sz));
  double *y = REAL(res);

  R_xlen_t ii = 0;
  R_xlen_t offset_cc = 0, offset_rr;

  if (sq) {
    for (int cc = 0; cc < sz; ++cc) {
      y[offset_cc + cc] = 0;
      offset_rr = cc * sz;
      for (int rr = cc + 1; rr < sz; ++rr) {
        offset_rr += sz;
	double value = x[ii] * x[ii];
        y[offset_cc + rr] = value;
        y[offset_rr + cc] = value;
        ii++;
      }
      offset_cc += sz;
    }
  } else {
    for (int cc = 0; cc < sz; ++cc) {
      y[offset_cc + cc] = 0;
      offset_rr = cc * sz;
      for (int rr = cc + 1; rr < sz; ++rr) {
        offset_rr += sz;
        y[offset_cc + rr] = x[ii];
        y[offset_rr + cc] = x[ii];
        ii++;
      }
      offset_cc += sz;
    }
  }

  UNPROTECT(1);
  return res;
}
