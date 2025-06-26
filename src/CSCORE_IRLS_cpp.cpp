#include "WLS_mean.h"
#include "WLS_cov.h"
#include <RcppArmadillo.h>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <string>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

arma::mat post_process_est(arma::mat est) {
  int p = est.n_rows;

  // Identify diagonal elements that are NA or Inf
  std::vector<int> neg_gene_inds;
  for (int i = 0; i < p; ++i) {
    double val = est(i, i);
    if (std::isnan(val) || !std::isfinite(val)) {
      neg_gene_inds.push_back(i);
    }
  }

  if (!neg_gene_inds.empty()) {
    Rcpp::Rcout << "[INFO] " << neg_gene_inds.size() << " among " << p
                << " genes have invalid variance estimates. Their co-expressions with other genes were set to 0." << std::endl;
  }

  // Set entire rows and columns of invalid genes to 0
  for (int idx : neg_gene_inds) {
    est.row(idx).zeros();
    est.col(idx).zeros();
  }

  // Set diagonal to 1
  est.diag().ones();

  // Clip values outside [-1, 1] and count them
  int upper_count = 0, lower_count = 0, total = 0;
  for (int i = 0; i < p; ++i) {
    for (int j = i + 1; j < p; ++j) {
      total++;
      if (est(i, j) > 1.0) {
        upper_count++;
        est(i, j) = 1.0;
        est(j, i) = 1.0;
      } else if (est(i, j) < -1.0) {
        lower_count++;
        est(i, j) = -1.0;
        est(j, i) = -1.0;
      }
    }
  }

  Rcpp::Rcout << "[INFO] "
              << std::fixed << std::setprecision(4)
              << 100.0 * upper_count / total
              << "% co-expression estimates were greater than 1 and were set to 1." << std::endl;
  Rcpp::Rcout << "[INFO] "
              << 100.0 * lower_count / total
              << "% co-expression estimates were smaller than -1 and were set to -1." << std::endl;

  return est;
}


// ---- inline quantile approximation ----
inline double arma_quantile(const arma::vec& v, double prob) {
  arma::uword n = v.n_elem;
  if (n == 0) return NA_REAL;

  arma::vec sorted = arma::sort(v);
  double index = (n - 1) * prob;
  arma::uword lo = std::floor(index);
  arma::uword hi = std::ceil(index);

  //Rcpp::Rcout << "[DEBUG] Input vector size: " << v.n_elem << std::endl;
  //Rcpp::Rcout << "[DEBUG] Sorted vector head: " << sorted.head(5).t();
  //Rcpp::Rcout << "[DEBUG] index: " << index << ", lo: " << lo << ", hi: " << hi << std::endl;

  if (lo == hi) {
    return sorted[lo];
  } else {
    double weight = index - lo;
    return (1.0 - weight) * sorted[lo] + weight * sorted[hi];
  }
}


// [[Rcpp::export]]
Rcpp::List CSCORE_IRLS_cpp_impl(const arma::mat& X,
                                const arma::vec& seq_depth_sq,
                                const arma::mat& D_mu,
                                const arma::mat& D_sigma2,
                                const arma::mat& D_sigma,
                                const bool post_process = true,
                                const int n_iter = 10,
                                const double eps = 0.05,
                                const bool verbose = false,
                                const std::string& conv = "q95",
                                const bool return_all = false) {

  int n = X.n_rows;
  int p = X.n_cols;

  arma::mat W = arma::ones<arma::mat>(X.n_rows, X.n_cols);
  arma::mat mu_beta = WLS_mean(D_mu, X, W);
  arma::mat M = D_mu * mu_beta;
  arma::rowvec mu = mu_beta.row(0);

  arma::mat X_centered = X - M;
  arma::mat sigma2_beta = WLS_mean(D_sigma2, (arma::square(X_centered) - M), W);
  arma::rowvec sigma2 = sigma2_beta.row(0);
  arma::rowvec theta = arma::square(mu) / sigma2;

  int iter = 0;
  double delta = std::numeric_limits<double>::infinity();

  while (delta > eps && iter <= n_iter) {
    arma::rowvec theta_prev = theta;
    arma::vec theta_pos = theta.elem(arma::find(theta > 0));
    double theta_median = arma::median(theta_pos);

    arma::rowvec mu_sq_over_theta = arma::square(mu) / theta_median;
    arma::mat weight = M + arma::repmat(seq_depth_sq, 1, p) % arma::repmat(mu_sq_over_theta, n, 1);
    weight.elem(arma::find_nonfinite(weight)).fill(1);
    weight.elem(arma::find(weight <= 0)).fill(1);

    mu_beta = WLS_mean(D_mu, X, 1 / weight);
    mu = mu_beta.row(0);
    M = D_mu * mu_beta;
    X_centered = X - M;

    arma::mat h = arma::square(arma::square(M) / theta_median + M);
    h.elem(arma::find(h <= 0)).fill(1);
    sigma2_beta = WLS_mean(D_sigma2, (arma::square(X_centered) - M), 1 / h);
    sigma2 = sigma2_beta.row(0);
    theta = arma::square(mu) / sigma2;

    arma::uvec valid_idx = arma::find(theta > 0 && theta_prev > 0);
    arma::vec ratio = arma::log(theta.elem(valid_idx) / theta_prev.elem(valid_idx));
    if (conv == "max") {
      delta = arma::max(arma::abs(ratio));
    } else if (conv == "q95") {
      delta = arma_quantile(arma::abs(ratio), 0.95);;
    }
    iter++;
    if (verbose) {
      Rcpp::Rcout << "[INFO] " << iter << ": delta=" << delta << "\n";
    }
  }

  if (iter > n_iter && delta > eps) {
    Rcpp::Rcout << "[WARNING] IRLS did not converge after " << n_iter << " iterations. \n Please increase the number of iterations or check your data.\n";
  } else {
    Rcpp::Rcout << "[INFO] IRLS converged after " << iter << " iterations.\n";
  }


  arma::rowvec mu_sq_over_theta = arma::square(mu) / arma::median(theta.elem(arma::find(theta > 0)));
  arma::mat w_cov = M + arma::repmat(seq_depth_sq, 1, p) % arma::repmat(mu_sq_over_theta, n, 1);
  w_cov.elem(arma::find_nonfinite(w_cov)).fill(1);
  w_cov.elem(arma::find(w_cov <= 0)).fill(1);
  X_centered = X - M;

  auto start = std::chrono::system_clock::now();
  std::time_t start_time = std::chrono::system_clock::to_time_t(start);
  Rcpp::Rcout << "[INFO] Starting WLS for covariance at " << std::ctime(&start_time);  // includes newline
  auto tic = std::chrono::high_resolution_clock::now();

  // Covariance estimation
  Rcpp::List covar_list = WLS_cov(D_sigma, X_centered, 1 / w_cov);
  arma::cube covar_beta = Rcpp::as<arma::cube>(covar_list["cov_hat"]);
  arma::mat covar = covar_beta.row(0);  // Extract k = 0 slice => p × p matrix

  // Statistical testing
  arma::mat Sigma = M + arma::repmat(seq_depth_sq, 1, p) % arma::repmat(sigma2, n, 1);
  arma::mat ele_inv_Sigma = 1 / Sigma;

  Rcpp::List ts_list = WLS_cov(D_sigma, X_centered, ele_inv_Sigma);
  arma::vec sigma = arma::sqrt(arma::clamp(sigma2.t(), 0, arma::datum::inf));
  //arma::cube covar_beta = Rcpp::as<arma::cube>(ts_list["cov_hat"]);
  //arma::mat covar = covar_beta.row(0);  // Extract k = 0 slice => p × p matrix
  arma::mat est = covar / (sigma * sigma.t());
  arma::mat test_stat = ts_list["test_stat"];
  arma::mat p_value = 2 * arma::normcdf(-arma::abs(test_stat));

  if (post_process) {
    est = post_process_est(est);
  }

  auto toc = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = toc - tic;
  Rcpp::Rcout << "[INFO] Finished WLS. Elapsed time: " << elapsed.count() << " seconds.\n";

  Rcpp::List result = List::create(Named("est") = est,
                                   Named("p_value") = p_value,
                                   Named("test_stat") = test_stat);

  if (return_all) {
    result["mu_beta"] = mu_beta;
    result["sigma2_beta"] = sigma2_beta;
    result["cov_beta"] = covar_beta;
  }

  return result;
}


