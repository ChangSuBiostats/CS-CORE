# CSCORE 1.0.1 (2025-06-26)
## Features
* Added support for covariate adjustment in moment-based regressions for estimating mean, variance, and covariance.

## Improvements
* Introduced Rcpp implementations of `WLS_mean()` and `WLS_cov()` to enhance memory efficiency while preserving computational performance.
* Implemented `CSCORE_IRLS_cpp()`, the cpp implementation of the full IRLS procedure to replace the original base R implementation.
* Added tests to ensure consistency with the original base R implementation.
* Updated documentation.

# CSCORE 1.0.0
* The first official release of the CS-CORE R pakcage (with NC revision).

# CSCORE 0.0.0.9000

* The first version of the CS-CORE R package.
* Added a `NEWS.md` file to track changes to the package.
