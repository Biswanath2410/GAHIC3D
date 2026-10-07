<a href="https://travis-ci.org/HenrikBengtsson/hsa"><img src="https://travis-ci.org/HenrikBengtsson/hsa.svg" alt="Build status"></a> <a href="https://codecov.io/gh/HenrikBengtsson/hsa"><img src="https://codecov.io/gh/HenrikBengtsson/hsa/branch/master/graph/badge.svg" alt="Coverage Status"/></a>



# hsa: A Non-Official R Package Version of hsa/cstruct1.R

This is an R package of the GPL-2 licensed [R scripts provided by the authors](http://ouyanglab.jax.org/hsa/) of:

* Zou C, Zhang Y, Ouyang Z. HSA: integrating multi-track Hi-C data for genome-scale reconstruction of 3D chromatin structure. Genome Biology, 2016, 17:40. [PMC4774023](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4774023/)

The code has been optimized for performance and memory usage, but has otherwise not been modified.  Validation of numerical reproducibility has been done based on public data sets.  See the R package tests for details.


## Installation

The R package **hsa** is only available on [GitHub](https://github.com/HenrikBengtsson/hsa).  It can be installed in R using:
```r
remotes::install_github("HenrikBengtsson/hsa")
```
This will install the package from source.
