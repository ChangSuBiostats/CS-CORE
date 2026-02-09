# Changelog

## CSCORE 1.0.2 (2025-09-13)

### Improvements

- Fixed BLAS/LAPACK linking on Windows by adding `src/Makevars.win`.
- Standardized dimension and loop variables to use `arma::uword`,
  removing signed/unsigned comparison warnings.
- Ensured compatibility with R 4.5.1 and RcppArmadillo 15.x (validated
  via win-builder and GitHub Actions CI).
- Improved consistency of header includes (`RcppArmadillo.h` first) for
  more robust builds.

### Infrastructure

- Added GitHub Actions workflow for Windows CI with latest R and
  RcppArmadillo.
- Added Windows build status badge to README.

## CSCORE 1.0.1 (2025-06-26)

### Features

- Added support for covariate adjustment in moment-based regressions for
  estimating mean, variance, and covariance.

### Improvements

- Introduced Rcpp implementations of
  [`WLS_mean()`](https://changsubiostats.github.io/CS-CORE/reference/WLS_mean.md)
  and
  [`WLS_cov()`](https://changsubiostats.github.io/CS-CORE/reference/WLS_cov.md)
  to enhance memory efficiency while preserving computational
  performance.
- Implemented
  [`CSCORE_IRLS_cpp()`](https://changsubiostats.github.io/CS-CORE/reference/CSCORE_IRLS_cpp.md),
  the cpp implementation of the full IRLS procedure to replace the
  original base R implementation.
- Added tests to ensure consistency with the original base R
  implementation.
- Updated documentation.

## CSCORE 1.0.0

- The first official release of the CS-CORE R pakcage (with NC
  revision).

## CSCORE 0.0.0.9000

- The first version of the CS-CORE R package.
- Added a `NEWS.md` file to track changes to the package.
