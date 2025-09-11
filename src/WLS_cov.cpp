#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

/*
 * Compute WLS estimates of covariance for all p by p gene pairs
 *
 * D: Design matrix for WLS (size n x k)
 * X: Gene expression matrix (size n x p)
 * W: Weight matrix (size n x p)
 *
 * Returns: a list with two components
 * cov_hat: numeric array of size k x p x p, where (, j ,j') stores the WLS estimates for gene pair j, j' across k covariates
 * test_stat: numeric matrix of size p x p, where (j, j') stores the test statistics for gene j, j' (H_0: independence)
 */
// [[Rcpp::export]]
Rcpp::List WLS_cov(arma::mat D, arma::mat X, arma::mat W) {
  const arma::uword n = D.n_rows;   // number of rows of D (i.e., number of cells)
  const arma::uword k = D.n_cols;   // number of columns of D (i.e., number of covariates)
  const arma::uword p = X.n_cols;   // number of genes

  if (X.n_rows != n) {
    stop("Dimensions of design matrix and gene expression datado not match.");
  }

  arma::cube result(k, p, p, arma::fill::zeros);  // k x p x p cube, WLS estimates
  //arma::cube inv_DTD_w_cube(k, k, p * p, arma::fill::zeros); // k x k x p x p cube, variance of WLS estimators under the null
  arma::mat ts_mat(p, p, arma::fill::zeros); // p x p matrix, test statistics

  arma::mat D_T = D.t();  // precompute transpose once

  for (arma::uword j = 0; j < p; ++j) {
    for (arma::uword jprime = j; jprime < p; ++jprime) {

      // Compute x_prod = X[, j] * X[, j']
      arma::vec x_prod = X.col(j) % X.col(jprime);  // element-wise product

      // Compute weights
      arma::vec w = W.col(j) % W.col(jprime);

      // Compute D^T * diag(w) * D
      arma::mat D_weighted = D;
      D_weighted.each_col() %= w;

      // Compute inverse of DTD_w
      arma::mat DTD_w = D_T * D_weighted;
      //arma::mat inv_DTD_w = arma::inv_sympd(DTD_w);  // safer inverse
      arma::mat inv_DTD_w;
      //try {
      //  inv_DTD_w = arma::inv_sympd(DTD_w);
      //} catch (...) {
        //Rcpp::Rcout << "j" << j << std::endl;
        //Rcpp::Rcout << "jprime" << jprime << std::endl;
        //Rcpp::Rcout << "DTD_w matrix:\n" << DTD_w << std::endl;
        // Rcpp::Rcout << "[WARNING] inv_sympd failed, using general inverse" << std::endl;
      //  inv_DTD_w = arma::inv(DTD_w);  // more stable but slower
      //}
      inv_DTD_w = arma::inv(DTD_w);  // more stable but slower

      // Compute D^T * (w % x_prod)
      arma::vec Dt_w_x = D_T * (w % x_prod);

      // WLS estimator
      arma::vec res_vec = inv_DTD_w * Dt_w_x;

      // Save covariance estimates
      result.slice(jprime).col(j) = res_vec;

      // Save inverse matrix
      //int pair_idx = j * p + jprime;
      //inv_DTD_w_cube.slice(pair_idx) = inv_DTD_w;

      // Save test statistics
      ts_mat(j, jprime) = res_vec(0) / std::sqrt(inv_DTD_w(0, 0));
    }
  }

  for (arma::uword j = 0; j < p; ++j) {
    for (arma::uword jprime = j + 1; jprime < p; ++jprime) {
      result.slice(j).col(jprime) = result.slice(jprime).col(j);
      ts_mat(jprime, j) = ts_mat(j, jprime);
    }
  }

  // Return both result and inverse cube
  return List::create(Named("cov_hat") = result,
                      //Named("var_cov_hat") = inv_DTD_w_cube,
                      Named("test_stat") = ts_mat);
}
