#' Weighted Least Squares Covariance Estimation
#'
#' Computes WLS estimates of gene-gene covariance and test statistics that assess
#' the statistical significance of co-expression using moment-based regressions.
#'
#' @param D Design matrix for WLS (size n x k)
#' @param X Gene expression matrix (n x p)
#' @param W Weight matrix (n x p)
#'
#' @return A list with:
#' \describe{
#'   \item{cov_hat}{A k x p x p array of WLS estimates}
#'   \item{test_stat}{A p x p matrix of test statistics}
#' }
#'
#' @export
WLS_cov <- function(D, X, W) {
  .Call(`_CSCORE_WLS_cov`, D, X, W)
}

#' Weighted Least Squares Mean and Variance Estimation
#'
#' Compute WLS estimates of mean or variance for all p genes
#'
#' @param D Design matrix for WLS (size n x k), where the first column represents baseline mean/var,
#'   and others represents covariates.
#' @param X Response matrix for WLS (size n x p)
#' @param W Weight matrix (n x p)
#'
#' @return A k x p matrix of WLS estimates
#'
#' @export
WLS_mean <- function(D, X, W) {
  .Call(`_CSCORE_WLS_mean`, D, X, W)
}
