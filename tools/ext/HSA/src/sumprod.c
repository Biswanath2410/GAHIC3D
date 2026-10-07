#include <R.h>
#include <Rdefines.h>
#include <R_ext/Error.h>
#include "000.api.h"

SEXP sumprod(SEXP x, SEXP y) {
  SEXP res;
  R_xlen_t n = xlength(x);
  double *xx = REAL(x);
  double *yy = REAL(y);
  double s = 0;
  if (xlength(y) != n) error("Argument \'x\' and \'y\' are of different lengths");
  for (R_xlen_t ii = 0; ii < n; ++ii) s += xx[ii] * yy[ii];
  return ScalarReal(s);
}

SEXP sumsq(SEXP x) {
  SEXP res;
  R_xlen_t n = xlength(x);
  double *xx = REAL(x);
  double s = 0;
  for (R_xlen_t ii = 0; ii < n; ++ii) s += xx[ii] * xx[ii];
  return ScalarReal(s);
}

