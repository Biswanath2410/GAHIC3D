#' @importFrom stats optim rnorm
finital <- function(pbin, A0, b0, invSigma0, beta1, covmat0, mat, floglike, fdloglike) {
  N <- length(pbin)
  if (N <= 500L) {
    P <- Sis(d = 5, pbin = pbin, A = A0, b = b0, invSigma = invSigma0,
             beta = beta1, cx0 = covmat0, mat0 = mat, q0 = c(0, 0, 0),
             fL = function(x, ...) -loglikelihood(x, ...))
  } else {
    if (N > 500L && N <= 2000L) {
      m <- 100L
    } else {
      m <- 200L
    }
    index <- cbind(c(1, seq(from = m, to = N - m, by = m)), c(seq(from = m, to = N - m, by = m), N))
    index[-1, 1] <- index[-1, 1] + 1L
    lP <- lapply(1:dim(index)[1], FUN = function(x) {
      idxs <- index[x, 1]:index[x, 2]
      subinitial(pbin = pbin[idxs], A0 = A0, b0 = b0, invSigma0 = invSigma0,
                 beta1 = beta1, covmat0 = covmat0[idxs, idxs, ],
                 mat = mat[idxs, c(1L, idxs + 1L), ],
                 floglike = floglike, fdloglike = fdloglike)
    })
    P <- matrix(0, nrow = N, ncol = 3L)
    P[index[1, 1]:index[1, 2], ] <- lP[[1]]
    for (i in 2:dim(index)[1]) {
      index_ri_c1 <- index[i, 1]
      index_ri_c2 <- index[i, 2]
      idxs <- 1:index_ri_c2
      P_i <- P[1:(index[i - 1, 2]), ]
      pbin_i <- pbin[idxs]
      lP_i <- lP[[i]]
      covmat0_i <- covmat0[idxs, idxs, ]
      mat_i <- mat[idxs, 1:(index_ri_c2 + 1L), ]
      theta <- optim(as.vector(rbind(diag(3), P[index[i - 1, 2], ] - P[index_ri_c1, ] + rnorm(3L, sd = 1/100))), fn = function(x) {
        piece(theta0 = x, P1 = P_i, P2 = lP_i, pbin = pbin_i, A0 = A0, b0 = b0,
              invSigma0 = invSigma0, beta1 = beta1, covmat0 = covmat0_i,
              mat = mat_i, floglike = floglike)
      })
      theta <- matrix(theta[["par"]], nrow = 4L, ncol = 3L)
      P[index_ri_c1:index_ri_c2, ] <- lP_i %*% theta[-4, ] + rep(1, times = index_ri_c2 - index_ri_c1 + 1L) %*% t(theta[4, ])
    }
  }
  avsmth(pbin, P = P)
}
