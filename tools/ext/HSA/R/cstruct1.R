# Publications that use results obtained from this software please include a citation of the paper:
# Zou, C., Zhang, Y., Ouyang, Z. (2016) HSA: integrating multi-track Hi-C data for genome-scale reconstruction of 3D chromatin structure. Submitted.

normP <- function(P) {
  P1 <- t(t(P) - colMeans(P))
  5 / sqrt(max(rowSums(P1^2))) * P1
}

kinetic0 <- function(p0) {
  sumsq(p0) / dim(p0)[1]
}

kinetic0_1 <- function(p0) {
  sumsq(p0)
}

kinetic <- function(p0, N, rho) {
  (sumsq(p0) - sumprod(p0[-1, ], p0[-N, ]) * rho * 2) / N
}

kinetic_1 <- function(p0, N, rho) {
  sumsq(p0) - sumprod(p0[-1, ], p0[-N, ]) * rho * 2
}

momentum <- local({
  zeros_3 <- rep(0, times = 3L)
  function(p0, N, rho) {
    p0 - rbind(p0[-1, ] * rho, zeros_3) - rbind(zeros_3, p0[-N, ] * rho)
  }
})

momentum0 <- function(p0) p0

#' @importFrom stats rnorm splinefun
#' @importFrom matrixStats colSds
finistructure <- function(S0, bin) {
  n <- dim(S0)[1]
  N <- nrow(bin)
  if (n == N) {
    if (dim(S0)[2] == 3) {
      S <- S0
    } else {
      S <- S0[, 3:5]
    }
  } else {
    pts <- c(S0[, 1], S0[n, 2])
    S0_c35 <- S0[, 3:5]
    Y <- rbind(S0_c35, rnorm(3, mean = S0_c35[n, ], sd = colSds(S0_c35[-1, ] - S0_c35[-n, ])))
    bin_c1 <- bin[, 1]
    S <- normP(sapply2(1:3, FUN = function(x) splinefun(pts, Y[, x])(bin_c1)))
    S <- S + rnorm(3 * N, mean = 0, sd = sqrt(5/N))
  }
  S
}

## PERFORMANCE: fmkorder_temp(), fmkorder(), and fmkorder2() now returns a list of
## parameter estimates that can be used "as is" in the likelihood calculations.
## Previously these functions would stack these parameters up as columns in a matrix
## which then had to be extracted again.  Returning the parameters as is also makes
## it much clear what is calculated.
fmkorder_temp <- function(m, A, b, sigma, S) {
  if (m < 2L) {
    mu <- A %*% S + b
    Sigma <- sigma
  } else {
    tmp <- fmkorder_temp(m - 1L, A, b, sigma, S)
    mu <- A %*% tmp[[1L]] + b                       ## tmp[["mu"]]
    temp <- tmp[[3L]] %*% A                         ## tmp[["A"]]
    Sigma <- tmp[[2L]] + t(temp) %*% sigma %*% temp ## tmp[["Sigma"]]
    A <- temp
  }
  list(mu = mu, Sigma = Sigma, A = A)
}

fmkorder <- function(m, A, b, sigma, S) {
  temp <- fmkorder_temp(m, A, b, sigma, S)
  Sigma <- temp[[2L]]
  Sigma <- solve(Sigma)
  temp[[1L]] <- c(temp[[1L]])  ## temp[["mu"]]
  temp[[2L]] <- Sigma
  temp
}

fmkorder2 <- local({
  zeros_3 <- rep(0, times = 3L)
  function(m, A, b, sigma) fmkorder(m, A, b, sigma, S = zeros_3)
})

fnormvec <- function(a, b) {
  c(a[2] * b[3] - a[3] * b[2],
    a[3] * b[1] - a[1] * b[3],
    a[1] * b[2] - a[2] * b[1])
}

fangle <- function(a, b) {
  acos(sumprod(a, b) / sqrt(sumsq(a) * sumsq(b)))
}

frotanyvec <- function(x, v, theta) {
  n <- v / sqrt(sumsq(v))
  t <- sumprod(x, n) * n
  cos(theta) * (x - t) + sin(theta) * fnormvec(n, x) + t
}

fbead <- function(S1, S2) {
  m <- dim(S1)[1]
  S <- S1
  S_r1 <- S[1, ]
  S_rm <- S[m, ]
  S2_r1 <- S2[1, ]
  S2_r2 <- S2[2, ]
  n <- fnormvec(S_rm - S_r1, S2_r2 - S2_r1)
  theta <- fangle(S_rm - S_r1, S2_r2 - S2_r1)
  S_t <- t(S) - S_r1
  S_t <- apply(S_t, MARGIN = 2L, FUN = frotanyvec, v = n, theta = theta)
  S <- t(S_t - S_t[, 1] + S1[1, ])
  S
}

fmirror <- local({
  diag_3 <- diag(3)
  function(v) {
    if (any(v)) {
      diag_3 - v %*% t(v) / sumsq(v)
    } else {
      diag_3
    }
  }
})

#' @importFrom stats optim
#' @importFrom MASS ginv
tranS <- local({
  zeros_9 <- rep(0, times = 9L)
  function(S1, S2, I_scale = TRUE) {
    tmp <- cbind(1, S1)
    tmp_t <- t(tmp)
    beta <- ginv(tmp_t %*% tmp) %*% (tmp_t %*% S2)
    s <- svd(beta[-1, ])
    tmp2 <- s[["u"]] %*% (diag(sign(s[["d"]]))) %*% t(s[["v"]])
    beta[-1, ] <- tmp2
    S <- tmp %*% beta
    if (I_scale) {
      beta[-1, ] <- mean(abs(s[["d"]])) * tmp2
      tmp <- optim(c(1, 0, 0, 0, 0, 0, 0), fn = function(x) {
        sum(sqrt(rowSums(
          (t(t(x[1] * S %*% angle2mat(x[2:4])) + x[5:7]) - S2)^2
        )))
      })
      tmp <- tmp[["par"]]
      S <- t(t(tmp[1] * S %*% angle2mat(tmp[2:4])) + tmp[5:7])
    } else {
      tmp <- optim(zeros_9, fn = function(x) {
        sum(sqrt(rowSums(
          (t(t(S %*% angle2mat(x[4:6]) %*% fmirror(x[1:3])) + x[7:9]) - S2)^2
        )))
      })
      tmp <- tmp[["par"]]
      S <- t(t(S %*% angle2mat(tmp[4:6]) %*% fmirror(tmp[1:3])) + tmp[7:9])
    }
    S
  }
})

#' @importFrom stats median splinefun
rmol <- function(loci, P) {
  n <- dim(P)[1]
  m <- dim(P)[2]
  P1 <- P
  v <- !is.finite(P1)
  if (any(v)) {
    outlier <- v[, 1] | v[, 2] | v[, 3]
    P2 <- P[!outlier, ]
    loci_outlier <- loci[outlier]
    loci_nonoutlier <- loci[!outlier]
    tmp <- sapply2(1:m, FUN = function(x) {
      splinefun(loci_nonoutlier, P2[, x])(loci_outlier)
    })
    P1[outlier, ] <- tmp
  }
  
  d1 <- rowSums((P1[-1, ] - P1[-n, ])^2)
  cutoff <- 100 * median(d1)
  outlier <- (d1 >= cutoff)
  if (any(outlier)) {
    v <- which(outlier)
    P1_t <- t(P1)
    for (i in 1:length(v)) {
      v_i <- v[i]
      idxs <- 1:v_i
      P1_t[, idxs] <- P1_t[, idxs] + (P[v_i + 1L, ] - P[v_i, ]) * 0.8
    }
    P1 <- t(P1_t)
    P1 <- tranS(P1, S2 = P)
  }

  P1
}

#' @importFrom stats filter
avsmth <- local({
  thirds_3 <- rep(1/3, times = 3L)
  function(bin, P) {
    N <- length(bin)
    P0 <- filter(P, thirds_3, sides = 2L)
    P0[1, ] <- (P[1, ] + P[2, ]) / 2
    P0[N, ] <- (P[N, ] + P[N - 1L, ]) / 2
    matrix(P0, nrow = N, ncol = 3L)
  }
})

loglikelihood0 <- function(P0, A, b, invSigma, beta, cx, mat, pos = NULL, v = NULL, mak = NULL) {
  L <- 0.0
  P <- P0[, -1]
  sigma <- solve(invSigma)
  N <- dim(P)[1]
  C <- dim(cx)[3]
  if (is.na(C)) {
    C <- 1L
    cx <- array(cx, dim = c(dim(cx), 1L))
    mat <- array(mat, dim = c(dim(mat), 1L))
  }
  if (is.null(pos)) {
    pos <- apply(mat, MARGIN = 3L, FUN = function(x) which(!is.na(x[, 1])))
    if (!is.list(pos)) {
      pos <- lapply(1:ncol(pos), FUN = function(i) pos[, i])
    }
  }
  if (is.null(v)) {
    v <- lapply(1:C, FUN = function(i) which(mat[pos[[i]], pos[[i]] + 1L, i] > 0))
  }
  distmat <- dist_matrix(P)
  for (i in 1:C) {
    pos_i <- pos[[i]]
    v_i <- v[[i]]
    beta_i <- beta[i]
    
    distmat_i <- distmat[pos_i, pos_i]
    cx_i <- cx[pos_i, pos_i, i]
    temp <- negCDbeta(C = cx_i, D = distmat_i, beta = beta_i)
    temp[v_i] <- temp[v_i] + mat[pos_i, pos_i + 1L, i][v_i] * (beta_i * log(distmat_i[v_i]) + log(cx_i[v_i]))
    L_i <- sum(upper_triangle(temp)) / (3 * N)
    
    L <- L + L_i
  }
  L
}

#' @importFrom matrixStats t_tx_OP_y
dloglikelihood0 <- function(P0, A, b, invSigma, beta, cx, mat, pos, v = NULL, mak = NULL) {
  P <- P0[, -1]
  sigma <- solve(invSigma)
  N <- dim(P)[1]
  C <- dim(cx)[3]
  if (is.null(v)) {
    v <- lapply(1:C, FUN = function(i) which(mat[pos[[i]], pos[[i]] + 1L, i] > 0))
  }

  dL <- matrix(0, nrow = N, ncol = 3L)
  distmat <- dist_matrix(P, square = TRUE)
  temp <- matrix(0, nrow = N, ncol = N)
  for (i in 1:C) {
    pos_i <- pos[[i]]
    v_i <- v[[i]]
    beta_i <- beta[i]
    distmat_i <- distmat[pos_i, pos_i]
    temp[pos_i, pos_i] <- temp[pos_i, pos_i] + beta_i * negCDbeta(C = cx[pos_i, pos_i, i], D = distmat_i, beta = beta_i / 2 - 1)
    temp[pos_i, pos_i][v_i] <- temp[pos_i, pos_i][v_i] + beta_i * mat[pos_i, pos_i + 1L, i][v_i] / distmat_i[v_i]
  }
  diag(temp) <- 0
  tmp <- temp * P[, 1]
  dL[, 1] <- rowSums(t_tx_OP_y(tmp, tmp, OP = 2L)) / (3 * N)
  tmp <- temp * P[, 2]
  dL[, 2] <- rowSums(t_tx_OP_y(tmp, tmp, OP = 2L)) / (3 * N)
  tmp <- temp * P[, 3]
  dL[, 3] <- rowSums(t_tx_OP_y(tmp, tmp, OP = 2L)) / (3 * N)
  dL
}

loglikelihood <- function(P0, A, b, invSigma, beta, cx, mat, pos = NULL, v = NULL, mak = NULL) {
  L <- 0.0
  P <- P0[, -1]
  sigma <- solve(invSigma)
  N <- dim(P)[1]
  C <- dim(cx)[3]
  if (is.na(C)) {
    C <- 1L
    cx <- array(cx, dim = c(dim(cx), 1L))
    mat <- array(mat, dim = c(dim(mat), 1L))
  }
  if (is.null(pos)) {
    pos <- apply(mat, MARGIN = 3L, FUN = function(x) which(!is.na(x[, 1])))
    if (!is.list(pos)) {
      pos <- lapply(1:ncol(pos), FUN = function(i) pos[, i])
    }
  }
  if (is.null(v)) {
    v <- lapply(1:C, FUN = function(i) which(mat[pos[[i]], pos[[i]] + 1L, i] > 0))
  }
  distmat <- dist_matrix(P)
  for (i in 1:C) {
    pos_i <- pos[[i]]
    v_i <- v[[i]]
    beta_i <- beta[i]
    
    distmat_i <- distmat[pos_i, pos_i]
    cx_i <- cx[pos_i, pos_i, i]
    temp <- negCDbeta(C = cx_i, D = distmat_i, beta = beta_i)
    temp[v_i] <- temp[v_i] + mat[pos_i, pos_i + 1L, i][v_i] * (beta_i * log(distmat_i[v_i]) + log(cx_i[v_i]))
    L_i <- sum(upper_triangle(temp)) / (3 * N)
    
    L <- L + L_i
  }
  if (is.null(mak)) {
    L <- L + sum(unlist(lapply(2:N, FUN = function(ii) {
      nmp <- fmkorder(P0[ii, 1] - P0[ii - 1L, 1], A = A, b = b, sigma = sigma, S = P[ii - 1L, ])
      R_ii <- P[ii, ] - nmp[[1L]] ## npm[["mu"]]
      Sigma <- nmp[[2L]]          ## npm[["Sigma"]]
      -R_ii %*% Sigma %*% R_ii / 2 + log_det(Sigma) / 2
    }), recursive = FALSE, use.names = FALSE)) / (3 * N)
  } else {
    L <- L + sum(unlist(lapply(2:N, FUN = function(ii) {
      mak_ii1 <- mak[[ii - 1L]]
      mu <- mak_ii1[[1L]]     ## mak_ii1[["mu"]]
      Sigma <- mak_ii1[[2L]]  ## mak_ii1[["Sigma"]]
      A <- mak_ii1[[3]]       ## mak_ii1[["A"]]
      mu <- mu + A %*% P[ii - 1L, ]
      R_ii <- P[ii, ] - mu
      -t(R_ii) %*% Sigma %*% R_ii / 2 + log_det(Sigma) / 2
    }), recursive = FALSE, use.names = FALSE)) / (3 * N)
  }
  L
}

#' @importFrom matrixStats t_tx_OP_y
dloglikelihood <- function(P0, A, b, invSigma, beta, cx, mat, pos, v = NULL, mak = NULL) {
  P <- P0[, -1]
  sigma <- solve(invSigma)

  N <- dim(P)[1]
  C <- dim(cx)[3]
  dL <- matrix(0, nrow = N, ncol = 3L)
  if (is.null(v)) {
    v <- lapply(1:C, FUN = function(i) which(mat[pos[[i]], pos[[i]] + 1L, i] > 0))
  }
  distmat <- dist_matrix(P, square = TRUE)

  temp <- matrix(0, nrow = N, ncol = N)
  for (i in 1:C) {
    pos_i <- pos[[i]]
    v_i <- v[[i]]
    beta_i <- beta[i]
    distmat_i <- distmat[pos_i, pos_i]
    temp[pos_i, pos_i] <- temp[pos_i, pos_i] + beta_i * negCDbeta(C = cx[pos_i, pos_i, i], D = distmat_i, beta = beta_i / 2 - 1)
    temp[pos_i, pos_i][v_i] <- temp[pos_i, pos_i][v_i] + beta_i * mat[pos_i, pos_i + 1L, i][v_i] / distmat_i[v_i]
  }
  diag(temp) <- 0
  tmp <- temp * P[, 1]
  dL[, 1] <- rowSums(t_tx_OP_y(tmp, tmp, OP = 2L)) / (3 * N)
  tmp <- temp * P[, 2]
  dL[, 2] <- rowSums(t_tx_OP_y(tmp, tmp, OP = 2L)) / (3 * N)
  tmp <- temp * P[, 3]
  dL[, 3] <- rowSums(t_tx_OP_y(tmp, tmp, OP = 2L)) / (3 * N)
  if (is.null(mak)) {
    nmp <- lapply(2:N, FUN = function(ii) {
      fmkorder(P0[ii, 1] - P0[ii - 1L, 1], A = A, b = b, sigma = sigma, S = P[ii - 1L, ])
    })
    
    temp1 <- sapply2(2:N, FUN = function(ii) {
      nmp_iim1 <- nmp[[ii - 1L]]
      mu <- nmp_iim1[[1L]]    ## nmp_ii1[["mu"]]
      Sigma <- nmp_iim1[[2L]] ## nmp_ii1[["Sigma"]]
      -t(P[ii, ] - mu) %*% Sigma
    }) / (3 * N)
    dL[2:N, ] <- dL[2:N, ] + t(temp1)
    
    temp2 <- sapply2(1:(N - 1L), FUN = function(ii) {
      nmp_ii <- nmp[[ii]]
      mu <- nmp_ii[[1L]]    ## nmp_ii1[["mu"]]
      Sigma <- nmp_ii[[2L]] ## nmp_ii1[["Sigma"]]
      A <- nmp_ii[[3L]]     ## nmp_ii1[["A"]]
      T <- A %*% P[ii, ]
      U <- t(A) %*% Sigma
      -U %*% T - U %*% (mu - T - P[ii + 1L, ])
    }) / (3 * N)
    dL[1:(N - 1L), ] <- dL[1:(N - 1L), ] + t(temp2)
  } else {
    temp1 <- sapply2(2:N, FUN = function(ii) {
      mak_ii1 <- mak[[ii - 1L]]
      mu <- mak_ii1[[1L]]    ## mak_ii1[["mu"]]
      Sigma <- mak_ii1[[2L]] ## mak_ii1[["Sigma"]]
      A <- mak_ii1[[3L]]     ## mak_ii1[["A"]]
      mu <- mu + A %*% P[ii - 1L, ]
      -t(P[ii, ] - mu) %*% Sigma
    }) / (3 * N)
    dL[2:N, ] <- dL[2:N, ] + t(temp1)
    temp2 <- sapply2(1:(N - 1L), FUN = function(ii) {
      mak_ii <- mak[[ii]]
      mu <- mak_ii[[1L]]    ## mak_ii[["mu"]]
      Sigma <- mak_ii[[2L]] ## mak_ii[["Sigma"]]
      A <- mak_ii[[3L]]     ## mak_ii[["A"]]
      U <- t(A) %*% Sigma
      -U %*% A %*% P[ii, ] - U %*% (mu - P[ii + 1L, ])
    }) / (3 * N)
    dL[1:(N - 1L), ] <- dL[1:(N - 1L), ] + t(temp2)
  }
  dL
}

mkcloglikelihood <- function(theta, P0) {
  P <- P0[, -1]
  N <- dim(P)[1]
  A <- matrix(theta[1:9], nrow = 3L, ncol = 3L)
  b <- theta[10:12]
  invSigma <- matrix(0, nrow = 3L, ncol = 3L)
  invSigma[upper.tri(invSigma, diag = TRUE)] <- theta[-(1:12)]
  invSigma <- invSigma + t(invSigma) - diag(diag(invSigma))
  sigma <- solve(invSigma)
  temp <- sapply2(2:N, FUN = function(ii) {
    nmp <- fmkorder(P0[ii, 1] - P0[ii - 1L, 1], A = A, b = b, sigma = sigma, S = P[ii - 1L, ])
    R_ii <- P[ii, ] - nmp[[1L]] ## npm[["mu"]]
    Sigma <- nmp[[2L]]          ## npm[["Sigma"]]
    -R_ii %*% Sigma %*% R_ii / 2 + log_det(Sigma) / 2
  })
  mean(temp) / 3
}

#' @importFrom matrixStats t_tx_OP_y
dhllk <- function(index, theta, P0, A, b, invSigma, beta, cx, mat, pos, v = NULL) {
  ## HB: This is very convoluted; candidate for speed up
  P <- rbind(P0[index[1, 1]:index[1, 2], -1], do.call(rbind, args = sapply2(2:dim(index)[1], FUN = function(x) {
    t(t(P0[index[x, 1]:index[x, 2], -1] %*% matrix(theta[x - 1L, 1:9], nrow = 3L, ncol = 3L)) + theta[x - 1L, 10:12])
  })))

  N <- dim(P)[1]
  C <- dim(cx)[3]
  dL <- matrix(0, nrow = dim(theta)[1], ncol = 12L)
  if (is.null(v)) {
    v <- lapply(1:C, FUN = function(i) which(mat[pos[[i]], pos[[i]] + 1L, i] > 0))
  }
  distmat <- dist_matrix(P, square = TRUE)
  temp <- matrix(0, nrow = N, ncol = N)
  dD <- array(0, dim = c(N, N, 3L))
  for (i in 1:C) {
    pos_i <- pos[[i]]
    v_i <- v[[i]]
    beta_i <- beta[i]
    distmat_i <- distmat[pos_i, pos_i]
    temp[pos_i, pos_i] <- temp[pos_i, pos_i] + beta_i * negCDbeta(C = cx[pos_i, pos_i, i], D = distmat_i, beta = beta_i / 2 - 1)
    temp[pos_i, pos_i][v_i] <- temp[pos_i, pos_i][v_i] + beta_i * mat[pos_i, pos_i + 1L, i][v_i] / distmat_i[v_i]
  }

  diag(temp) <- 0
  tmp <- temp * P[, 1]
  dD[, , 1] <- t(tmp) - tmp
  tmp <- temp * P[, 2]
  dD[, , 2] <- t(tmp) - tmp
  tmp <- temp * P[, 3]
  dD[, , 3] <- t(tmp) - tmp
  index_rn1 <- index[-1, ]
  tmp <- t(apply(index_rn1, MARGIN = 1L, FUN = function(x) {
    idxs <- x[1]:x[2]    
    P0t <- P0[idxs, -1]
    c(colSums(dD[-idxs, idxs, 1] %*% P0t),
      colSums(dD[-idxs, idxs, 2] %*% P0t),
      colSums(dD[-idxs, idxs, 3] %*% P0t))
  })) / (3 * N)
  dL[, 1:9] <- tmp

  dL[, 1:9] <- dL[, 1:9] + t(apply(index_rn1, MARGIN = 1L, FUN = function(x) {
    idxs <- x[1]:x[2]
    sapply2(2:4, FUN = function(k) {
      P0t <- P0[idxs, k]
      T_1 <- dD[idxs, idxs, 1] * P0t
      T_2 <- dD[idxs, idxs, 2] * P0t
      T_3 <- dD[idxs, idxs, 3] * P0t
      c(sum(upper.tri(t(T_1)) - upper.tri(T_1)),
        sum(upper.tri(t(T_2)) - upper.tri(T_2)),
        sum(upper.tri(t(T_3)) - upper.tri(T_3)))
    })
  })) / (3 * N)
  
  dL[, 10] <- apply(index_rn1, MARGIN = 1L, FUN = function(x) {
    idxs <- x[1]:x[2]
    sum(dD[-idxs, idxs, 1])
  }) / (3 * N)
  dL[, 11] <- apply(index_rn1, MARGIN = 1L, FUN = function(x) {
    idxs <- x[1]:x[2]
    sum(dD[-idxs, idxs, 2])
  }) / (3 * N)
  dL[, 12] <- apply(index_rn1, MARGIN = 1L, FUN = function(x) {
    idxs <- x[1]:x[2]
    sum(dD[-idxs, idxs, 3])
  }) / (3 * N)
  dL
}

rotamat <- function(theta, full = TRUE) {
  sin_theta <- sin(theta)
  cos_theta <- cos(theta)
  R <- matrix(c(
    cos_theta[1] * cos_theta[3] - cos_theta[2] * sin_theta[1] * sin_theta[3],
   -cos_theta[2] * cos_theta[3] * sin_theta[1] - cos_theta[1] * sin_theta[3],
    sin_theta[1] * sin_theta[2],
    cos_theta[3] * sin_theta[1] + cos_theta[1] * cos_theta[2] * sin_theta[3],
    cos_theta[1] * cos_theta[2] * cos_theta[3] - sin_theta[1] * sin_theta[3],
   -cos_theta[1] * sin_theta[2],
    sin_theta[2] * sin_theta[3],
    cos_theta[3] * sin_theta[2],
    cos_theta[2]
  ), nrow = 3L, ncol = 3L)
  
  if (full) R <- cbind(R, sin_theta, cos_theta)

  R
}

dDtotheta <- function(p, b_t, temp) {
  l <- matrix(0, nrow = dim(b_t)[2], ncol = 3L)

  temp_11 <- temp[1, 1]
  temp_12 <- temp[1, 2]
  temp_13 <- temp[1, 3]
  temp_14 <- temp[1, 4]
  temp_15 <- temp[1, 5]
  temp_21 <- temp[2, 1]
  temp_22 <- temp[2, 2]
  temp_23 <- temp[2, 3]
  temp_31 <- temp[3, 1]
  temp_32 <- temp[3, 2]
  temp_33 <- temp[3, 3]
  temp_35 <- temp[3, 5]

  b_c1 <- b_t[1, ]
  b_c2 <- b_t[2, ]
  b_c3 <- b_t[3, ]
  
  l[, 1] <- -2 * (p[2] * -temp_22 - p[1] * temp_12 - temp_32 * p[3]) * (b_c1 - p[1] * temp_11 + p[2] * -temp_21 - temp_31 * p[3]) - 2 * (p[1] * temp_11 - p[2] * -temp_21 + temp_31 * p[3]) * (b_c2 - p[1] * temp_12 + p[2] * -temp_22 - temp_32 * p[3])

  l[, 2] <- 2 * (temp_15 * temp_33 * p[3] + temp_15 * temp_23 * p[2] + temp_15 * temp_13 * p[1]) * (b_c2 - p[1] * temp_12 + p[2] * -temp_22 - temp_32 * p[3]) - 2 * (temp_33 * temp_14 * p[3] + temp_35 * temp_31 * p[2] + temp_14 * temp_13 * p[1]) * (b_c1 - p[1] * temp_11 + p[2] * -temp_21 - temp_31 * p[3]) + 2 * (temp_33 * temp_35 * p[2] - temp[2, 4] * p[3] + temp_33 * temp[3, 4] * p[1]) * (temp_33 * p[3] - b_c3 + temp_23 * p[2] + temp_13 * p[1])

  l[, 3] <- 2 * (temp_23 * p[1] - temp_13 * p[2]) * (temp_33 * p[3] - b_c3 + temp_23 * p[2] + temp_13 * p[1]) + 2 * (p[1] * -temp_22 + p[2] * temp_12) * (b_c2 - p[1] * temp_12 + p[2] * -temp_22 - temp_32 * p[3]) + 2 * (p[1] * -temp_21 + p[2] * temp_11) * (b_c1 - p[1] * temp_11 + p[2] * -temp_21 - temp_31 * p[3])

  l
}

#' @importFrom matrixStats t_tx_OP_y
dhllk1 <- function(index, theta, P0, A, b, invSigma, beta, cx, mat, pos, v = NULL) {
  matheta <- lapply(2:dim(index)[1], FUN = function(x) rotamat(theta[x - 1L, ]))
  ## HB: This is very convoluted; candidate for speed up
  P <- rbind(P0[index[1, 1]:index[1, 2], -1], do.call(rbind, args = sapply2(2:dim(index)[1], FUN = function(x) {
    t(t(P0[index[x, 1]:index[x, 2], -1] %*% matheta[[x - 1L]][, 1:3]) + theta[x - 1L, 4:6])
  })))

  N <- dim(P)[1]
  C <- dim(cx)[3]
  dL <- matrix(0, nrow = dim(theta)[1], ncol = 6L)
  if (is.null(v)) {
    v <- lapply(1:C, FUN = function(i) which(mat[pos[[i]], pos[[i]] + 1L, i] > 0))
  }
  distmat <- dist_matrix(P, square = TRUE)
  temp <- matrix(0, nrow = N, ncol = N)
  dD <- array(0, dim = c(N, N, 3L))
  for (i in 1:C) {
    pos_i <- pos[[i]]
    v_i <- v[[i]]
    beta_i <- beta[i]
    distmat_i <- distmat[pos_i, pos_i]
    temp[pos_i, pos_i] <- temp[pos_i, pos_i] + beta_i * negCDbeta(C = cx[pos_i, pos_i, i], D = distmat_i, beta = beta_i / 2 - 1)
    temp[pos_i, pos_i][v_i] <- temp[pos_i, pos_i][v_i] + beta_i * mat[pos_i, pos_i + 1L, i][v_i] / distmat_i[v_i]
  }
  diag(temp) <- 0
  tmp <- temp * P[, 1]
  dD[, , 1] <- t(tmp) - tmp
  tmp <- temp * P[, 2]
  dD[, , 2] <- t(tmp) - tmp
  tmp <- temp * P[, 3]
  dD[, , 3] <- t(tmp) - tmp
  dL[, 1:3] <- t(sapply2(2:dim(index)[1], FUN = function(x) {
    idxs <- index[x, 1]:index[x, 2]
    rowSums(sapply2(idxs, FUN = function(i) colSums(dDtotheta(P0[i, -1], b_t = t(P[-idxs, ]) - theta[x - 1, 4:6], temp = matheta[[x - 1]]) * temp[-idxs, i])))
  }))
  index_rn1 <- index[-1, ]
  dL[, 4] <- apply(index_rn1, MARGIN = 1L, FUN = function(x) {
    idxs <- x[1]:x[2]
    sum(dD[-idxs, idxs, 1])
  })
  dL[, 5] <- apply(index_rn1, MARGIN = 1L, FUN = function(x) {
    idxs <- x[1]:x[2]
    sum(dD[-idxs, idxs, 2])
  })
  dL[, 6] <- apply(index_rn1, MARGIN = 1L, FUN = function(x) {
    idxs <- x[1]:x[2]
    sum(dD[-idxs, idxs, 3])
  })
  dL / (3 * N)
}

angle2mat <- function(theta) {
  rotamat(theta, full = FALSE)
}

v2mat <- function(theta) {
  matrix(theta[1:9], nrow = 3L, ncol = 3L)
}

Leapfrog <- function(grad_U, L, epsilon, p0, q0, fM) {
  N <- dim(p0)[1]
  if (L < 2L) {
    q <- q0 + epsilon * fM(p0)
    p <- p0 - epsilon * grad_U(q) / 2
  } else {
    temp <- Leapfrog(grad_U, L = L - 1L, epsilon, p0, q0, fM)
    q <- temp[[2]] + epsilon * fM(temp[[1]])
    p <- temp[[1]] - epsilon * grad_U(temp[[2]]) / 2
  }
  list(p, q)
}

#' @importFrom stats rnorm runif
HMC <- function(U, grad_U, epsilon, L, current_q0, T0, fK, fM, I_trans = FALSE) {
  N <- dim(current_q0)[1]
  m <- dim(current_q0)[2]
  current_q <- current_q0
  q <- current_q
  if (I_trans) {
    q <- rmol(1:N, P = current_q)
  }

  p <- array(rnorm(N * m, mean = 0, sd = 1 / sqrt(N)), dim = dim(q)) # independent standard normal variates
  current_p <- p
  # Make a half step for momentum at the beginning
  p <- p - epsilon * grad_U(q) / 2
  # Alternate full steps for position and momentum
  temp <- Leapfrog(grad_U, L, epsilon, p, q, fM)
  # Make a half step for momentum at the end.
  q <- temp[[2]]
  p <- temp[[1]] - epsilon * grad_U(temp[[2]]) / 2
  # Negate momentum at end of trajectory to make the proposal symmetric
  p <- -p

  # Evaluate potential and kinetic energies at start and end of trajectory
  current_U <- U(current_q)
  current_K <- fK(current_p) / 2
  proposed_U <- U(q)
  proposed_K <- fK(p) / 2
  const <- (current_U - proposed_U + current_K - proposed_K)

  # Accept or reject the state at end of trajectory, returning either
  # the position at the end of the trajectory or the initial position
  if (any(is.infinite(const) + is.na(const))) {
    cat("NA or inf produced!!", "\n")
    return(rmol(1:N, P = current_q))
  } else {
    if (runif(1L) < exp((const) / T0)) {
      if (I_trans) {
        q <- rmol(1:N, P = q)
        q <- tranS(q, S2 = current_q)
      }
      return(q) # accept
    } else {
      return(current_q) # reject
    }
  }
}

#' @importFrom stats rnorm runif
## HMC1(U, grad_U, epsilon, L, current_q0, T0, fK, fM, I_trans = FALSE)
HMC1 <- function(U, grad_U, epsilon, L, current_q0, T0, fK, fM, I_trans = FALSE) {
  N <- dim(current_q0)[1]
  m <- dim(current_q0)[2]
  current_q <- current_q0
  q <- current_q

  p <- array(rnorm(N * m, mean = 0, sd = 1 / sqrt(N)), dim = dim(q)) # independent standard normal variates
  current_p <- p
  # Make a half step for momentum at the beginning
  p <- p - epsilon * grad_U(q) / 2
  # Alternate full steps for position and momentum
  temp <- Leapfrog(grad_U, L, epsilon, p, q, fM)
  # Make a half step for momentum at the end.
  q <- temp[[2]]
  p <- temp[[1]] - epsilon * grad_U(temp[[2]]) / 2
  # Negate momentum at end of trajectory to make the proposal symmetric
  p <- -p
  # Evaluate potential and kinetic energies at start and end of trajectory
  current_U <- U(current_q)
  current_K <- fK(current_p) / 2
  proposed_U <- U(q)
  proposed_K <- fK(p) / 2
  const <- (current_U - proposed_U) / log(runif(1L) * 0.5)
  const <- ifelse(const > 0, const, 1)

  # Accept or reject the state at end of trajectory, returning either
  # the position at the end of the trajectory or the initial position
  if (any(c(is.infinite(const), is.na(const)))) {
    cat("NA or inf produced!!", "\n")
    return(rmol(1:N, P = current_q))
  } else {
    ## HB: (const) / T0 / const ?!?
    if (runif(1L) < exp((const) / T0 / const)) {
      if (I_trans) {
        q <- rmol(1:N, P = q)
        q <- tranS(q, S2 = current_q)
      }
      return(q) # accept
    } else {
      return(current_q) # reject
    }
  }
}


#' @importFrom stats optim rnorm
Qsis <- function(N, pbin, A, b, invSigma, beta, cx, mat, q0, fL) {
  if (N == 2L) {
    pbin_t <- pbin[1:2]
    cx_t <- cx[1:2, 1:2, ]
    mat_t <- mat[1:2, 1:3, ]
    Padd <- optim(b + rnorm(3L, sd = 1/5) + q0, fn = function(q) {
      fL(cbind(pbin_t, rbind(q0, q)), A, b, invSigma, beta, cx_t, mat_t)
    })
    return(rbind(q0, t(Padd[["par"]])))
  } else {
    temp <- Qsis(N = N - 1L, pbin = pbin, A = A, b = b, invSigma = invSigma,
                 beta = beta, cx = cx, mat = mat, q0 = q0, fL = fL)
    pbin_t <- pbin[1:N]
    cx_t <- cx[1:N, 1:N, ]
    mat_t <- mat[1:N, 1:(N + 1L), ]
    Padd <- optim(rnorm(3L, sd = 1/5) + A %*% temp[dim(temp)[1], ] + b, fn = function(q) {
      fL(cbind(pbin_t, rbind(temp, q)), A, b, invSigma, beta, cx_t, mat_t)
    })
    return(rbind(temp, t(Padd[["par"]])))
  }
}

#' @importFrom stats optim rnorm
Sis <- function(d, pbin, A, b, invSigma, beta, cx0, mat0, q0, fL) {
  N <- dim(mat0)[1]
  if (is.na(dim(cx0)[3])) {
    cx <- array(cx0, dim = c(dim(cx0), 1L))
    mat <- array(mat0, dim = c(dim(mat0), 1L))
  } else {
    cx <- cx0
    mat <- mat0
  }
  if (N == d) {
    return(Qsis(N = d, pbin = pbin, A = A, b = b, invSigma = invSigma,
                beta = beta, cx = cx, mat = mat, q0 = q0, fL = fL))
  } else {
    temp <- Sis(d = d, pbin = pbin[-N], A = A, b = b, invSigma = invSigma,
                beta = beta, cx0 = cx[-N, -N, ], mat0 = mat[-N, -(N + 1L), ],
                q0 = q0, fL = fL)
    idxs <- (N - d):N
    pbin_t <- pbin[idxs]
    nrow <- nrow(temp)
    temp_t <- temp[(nrow - d + 1L):nrow, ]
    cx_t <- cx[idxs, idxs, ]
    mat_t <- mat[idxs, c(1, (N - d + 1L):(N + 1L)), ]
    addq <- optim(rnorm(3L, sd = 1/5) + A %*% temp[nrow, ] + b, fn = function(x) {
      fL(cbind(pbin_t, rbind(temp_t, x)), A, b, invSigma, beta, cx_t, mat_t)
    })
    addq <- addq[["par"]]
    return(rbind(temp, t(addq)))
  }
}

subinitial <- function(pbin, A0, b0, invSigma0, beta1, covmat0, mat, floglike, fdloglike) {
  N <- length(pbin)
  if (is.na(dim(covmat0)[3])) {
    covmat0 <- array(covmat0, dim = c(dim(covmat0), 1L))
    mat <- array(mat, dim = c(dim(mat), 1L))
  }
  pos <- apply(mat, MARGIN = 3L, FUN = function(x) which(!is.na(x[, 1])))
  if (!is.list(pos)) {
    pos <- lapply(1:ncol(pos), FUN = function(i) pos[, i])
  }
  P0 <- Sis(d = 4, pbin = pbin, A = A0, b = b0, invSigma = invSigma0,
            beta = beta1, cx0 = covmat0, mat0 = mat, q0 = c(0, 0, 0),
            fL = function(x, ...) -loglikelihood(x, ...))
  u <- 0L
  while (u < 100L) {
    ## HMC(U, grad_U, epsilon, L, current_q0, T0, fK, fM, I_trans = FALSE)
    P <- HMC(
      U = function(x) -floglike(cbind(pbin, x), A0, b0, invSigma0, beta1, covmat0, mat),
      grad_U = function(y) -fdloglike(cbind(pbin, y), A0, b0, invSigma0, beta1, covmat0, mat, pos),
      epsilon = 0.002,
      L = 20L,
      current_q0 = P0,
      T0 = 10 * exp(-u / 20),
      fK = function(p) kinetic(p, N = N, rho = 0.1),
      fM = function(p) momentum(p, N = N, rho = 0.1),
      I_trans = TRUE)
    P0 <- P
    u <- u + 1L
  }
  P
}

suboptimz <- function(pbin, P0, A0, b0, invSigma0, beta1, covmat0, mat, floglike, fdloglike, I_mle = TRUE) {
  epslon <- 0.0005
  stp <- 35L
  M <- 100L
  if (is.na(dim(covmat0)[3])) {
    covmat0 <- array(covmat0, dim = c(dim(covmat0), 1L))
    mat <- array(mat, dim = c(dim(mat), 1L))
  }
  pos <- apply(mat, MARGIN = 3L, FUN = function(x) which(!is.na(x[, 1])))
  if (!is.list(pos)) {
    pos <- lapply(1:ncol(pos), FUN = function(i) pos[, i])
  }
  u <- 0L
  Po <- P0
  N <- dim(P0)[1]
  Lo <- floglike(cbind(pbin, P0), A0, b0, invSigma0, beta1, covmat0, mat)
  ## HMC1(U, grad_U, epsilon, L, current_q0, T0, fK, fM, I_trans = FALSE)
  P <- HMC1(
    U = function(x) -floglike(cbind(pbin, x), A0, b0, invSigma0, beta1, covmat0, mat),
    grad_U = function(y) -fdloglike(cbind(pbin, y), A0, b0, invSigma0, beta1, covmat0, mat, pos),
    epsilon = epslon,
    L = stp,
    current_q0 = P0,
    T0 = exp(-u * 3.5 / M),
    fK = kinetic0_1,
    fM = momentum0)
  P0 <- P
  if (I_mle) {
    while (u < M) {
      ## HMC1(U, grad_U, epsilon, L, current_q0, T0, fK, fM, I_trans = FALSE)
      P <- HMC1(
        U = function(x) -floglike(cbind(pbin, x), A0, b0, invSigma0, beta1, covmat0, mat),
        grad_U = function(y) -fdloglike(cbind(pbin, y), A0, b0, invSigma0, beta1, covmat0, mat, pos),
        epsilon = epslon,
        L = stp,
        current_q0 = P0,
        T0 = exp(-u * 3.5 / M),
        fK = kinetic0_1,
        fM = momentum0)
      P0 <- P
      L <- floglike(cbind(pbin, P), A0, b0, invSigma0, beta1, covmat0, mat)
      if (L >= Lo) {
        Po <- P
        Lo <- L
      }
      u <- u + 1L
    }
  } else {
    while (u < M) {
      ## HMC1(U, grad_U, epsilon, L, current_q0, T0, fK, fM, I_trans = FALSE)
      P <- HMC1(
        U = function(x) -floglike(cbind(pbin, x), A0, b0, invSigma0, beta1, covmat0, mat),
        grad_U = function(y) -fdloglike(cbind(pbin, y), A0, b0, invSigma0, beta1, covmat0, mat, pos),
        epsilon = epslon,
        L = stp,
        current_q0 = P0,
        T0 = exp(-u * 3.5 / M),
        fK = kinetic0_1,
        fM = momentum0)
      P0 <- P
      u <- u + 1L
    }
    Po <- P0
  }
  Po
}

piece <- function(theta0, P1, P2, pbin, A0, b0, invSigma0, beta1, covmat0, mat, floglike) {
  theta <- matrix(theta0, nrow = 4L, ncol = 3L)
  -floglike(cbind(pbin, rbind(P1, t(t(P2 %*% theta[-4, ]) + theta[4, ]))), A0, b0, invSigma0, beta1, covmat0, mat)
}
