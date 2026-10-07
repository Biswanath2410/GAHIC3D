#include <R.h>
#include <Rdefines.h>

SEXP negCDbeta(SEXP C, SEXP D, SEXP beta);
SEXP sumprod(SEXP x, SEXP y);
SEXP sumsq(SEXP x);
SEXP upper_triangle(SEXP X, SEXP diag);
SEXP dist_matrix(SEXP X, SEXP size, SEXP square);
