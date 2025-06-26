#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;


/*
 * Compute WLS estimates of mean or variance for all p genes
 *
 * D: Design matrix for WLS (size n x k), where the first column represents baseline mean/var, and others represents covariates
 * X: Response matrix for WLS (size n x p)
 * W: Weight matrix (size n x p)
 *
 * Returns: numeric matrix of size k x p, where (, j) stores the WLS estimates for gene j
 */
// [[Rcpp::export]]
arma::mat WLS_mean(arma::mat D, arma::mat X, arma::mat W) {
  int n = D.n_rows;   // number of rows of D (i.e., number of cells)
  int k = D.n_cols;   // number of columns of D (i.e., number of covariates)
  int p = X.n_cols;   // number of genes

  if (X.n_rows != n) {
    stop("Dimensions of design matrix and gene expression data do not match.");
  }

  arma::mat result(k, p, arma::fill::zeros);  // k x p matrix, WLS estimates

  arma::mat D_T = D.t();  // precompute transpose once

  for (int j = 0; j < p; ++j) {
    // Compute D^T * diag(w) * D
    arma::mat DTD_w = D_T * (D.each_col() % W.col(j));

    // Inverse of DTD_w
    arma::mat inv_DTD_w = arma::inv(DTD_w);

    // Compute D^T * (w % x_prod)
    arma::vec Dt_w_x = D_T * (X.col(j) % W.col(j));

    // WLS estimator
    arma::vec res_vec = inv_DTD_w * Dt_w_x;

    // Save WLS estimates
    result.col(j) = res_vec;
  }

  return result;
}
