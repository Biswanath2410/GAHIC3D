library(hsa)
args <- commandArgs(trailingOnly=TRUE)
data <- as.matrix(read.table(args[1], header=FALSE)); storage.mode(data) <- "double"
set.seed(as.integer(args[4]))
out <- fmain(lsmap0=list(data), lscov0=0L, outfile=args[2], Maxiter=10L, submaxiter=300L, lambda=25,
             Leapfrog=3L, epslon=0.003, mkfix=0, rho=0, mk=as.integer(args[3]))
