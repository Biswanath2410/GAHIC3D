## sumsq(x) is a faster version of sum(x^2)
sumsq <- function(x) {
  .Call(C_sumsq, x)
}

## sumprod(x, y) is a faster version of sum(x * y)
sumprod <- function(x, y) {
  .Call(C_sumprod, x, y)
}

## negCDbeta(C, D, beta) is a faster version of - C * D ^ beta
negCDbeta <- function(C, D, beta) {
  .Call(C_negCDbeta, C, D, beta)
}

## upper_triangle(X) is a faster version of X[upper.tri(X)]
upper_triangle <- function(x, diag = FALSE) {
  .Call(C_upper_triangle, x, diag)
}

## PERFORMANCE: dist_matrix(x) is a faster version of as.matrix(dist(x)),
## because it avoids the overhead from S3 method dispatching, handling of
## non-needed attributes etc.  Moreover, dist_matrix(x, square = TRUE) is
## avoids internal duplication of the distance matrix.
#' @importFrom stats dist
dist_matrix <- function(x, square = FALSE) {
  x <- dist(x)
  size <- attr(x, "Size")
  .Call(C_dist_matrix, x, size, square)
}


## PERFORMANCE: log_det(x) is a faster version of log(det(x)), because
## it avoids overhead from S3 method dispatching and skips an internal
## log(exp(t)) step.
log_det <- function(x) {
  determinant.matrix(x, logarithm = TRUE)[[1L]]
}

## PERFORMANCE: Avoid overhead from S3 dispatch on solve() and as.matrix()
#' @importFrom utils globalVariables
solve <- local({
  bs <- list()
  
  globalVariables("solve_matrix")
  solve_matrix <- NULL
  
  function(a, b = NULL) {
    if (is.null(solve_matrix)) {
      fcn <- base::solve.default
      body <- body(fcn)
      body <- body[[length(body)]]
      ## Did we get what we expected, otherwise ignore it.
      if (identical(body[[1]], as.symbol(".Internal"))) body(fcn) <- body
      solve_matrix <<- fcn
    }
    
    if (is.null(b)) {
      n <- nrow(a)
      ## Memoization of 'b'
      if (n <= length(bs)) b <- bs[[n]]
      if (is.null(b)) {
        b <- diag(1, nrow = n)
        bs[[n]] <<- b
      }
    }

    solve_matrix(a, b)
  }
})

## PERFORMANCE: Avoid overhead from setting dimnames in cbind() and rbind()
cbind <- local({
  base_cbind <- base::cbind
  function(...) base_cbind(..., deparse.level = 0L)
})

rbind <- local({
  base_rbind <- base::rbind
  function(...) base_rbind(..., deparse.level = 0L)
})

colSums <- function(x) {
  dim <- dim(x)
  .colSums(x, m = dim[1], n = dim[2], na.rm = FALSE)
}

rowSums <- function(x) {
  dim <- dim(x)
  .rowSums(x, m = dim[1], n = dim[2], na.rm = FALSE)
}

colMeans <- function(x) {
  dim <- dim(x)
  .colMeans(x, m = dim[1], n = dim[2], na.rm = FALSE)
}

## PERFORMANCE: Remove all unnecessary overhead from sapply()
sapply2 <- function(X, FUN, ...) {
  names(X) <- NULL
  x <- lapply(X = X, FUN = FUN, ...)
  n <- length(x) 
  if (n == 0L) return(x)
  
  ns <- lengths(x, use.names = FALSE)
  common.len <- unique(ns)
  if (length(common.len) > 1L) return(x)

  if (common.len == 0L) return(x)
  
  r <- unlist(x, recursive = FALSE, use.names = FALSE)
  if (common.len == 1L) return(r)

  d <- c(common.len, n)
  if (prod(d) != length(r)) return(x)
  
  dim(r) <- d
  r
}


## BACKWARD COMPATIBILITY
t_tx_OP_y <- local({
  if (packageVersion("matrixStats") <= "0.53.1") {
    fcn0 <- matrixStats::t_tx_OP_y
    function(x, y, OP) {
      OP <- c("+", "-", "*", "/")[OP]
      fcn0(x, y, OP)
    }
  } else {
    matrixStats::t_tx_OP_y
  }
})


## AD HOC: Trick cstruct1.R code to write files with 12 digits
## (still plenty) instead of 15 digits for easier 'diff' comparisons
#' @importFrom utils write.table
write_tsv <- function(x, ..., row.names = FALSE, col.names = FALSE,
                      sep = "\t", eol = "\n", digits = 12L) {
  for (kk in seq_along(x)) {
    if (is.double(x[[kk]])) x[[kk]] <- round(x[[kk]], digits = digits)
  }
  write.table(x, row.names = row.names, col.names = col.names,
              sep = sep, eol = eol, ...)
}
