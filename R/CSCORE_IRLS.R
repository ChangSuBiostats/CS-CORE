#' Iteratively reweighted least squares (IRLS) procedure in CS-CORE
#'
#' This function implements the IRLS procedure used in CS-CORE for estimating and testing
#' cell-type-specific co-expression from single-cell RNA sequencing data.
#'
#' Let \eqn{x_{ij}} denote the UMI count of gene \eqn{j} in cell \eqn{i};
#' \eqn{s_i} denote the sequencing depth;
#' \eqn{\mu_j,\sigma_{jj}, \sigma_{jj'}} denote the mean, variance and covariance;
#' \eqn{c_{ik}} denote additional covariate \eqn{k} for cell \eqn{i} (e.g. disease status or cellular state).
#' The procedure consists of two main steps:
#' \enumerate{
#'   \item \strong{Mean and variance estimation:} Estimate gene-specific mean and variance parameters
#'         using two moment-based regressions:
#'     \itemize{
#'       \item \strong{Mean model:} \eqn{x_{ij} = s_i (\mu_j + \sum_k c_{ik} \beta_k) + \epsilon_{ij}}
#'       \item \strong{Variance model:} \eqn{(x_{ij} - s_i \mu_{ij})^2 = s_i \mu_{ij} + s_i^2 (\sigma_{jj} + \sum_k c_{ik} \gamma_k) + \eta_{ij}},
#'       where \eqn{\mu_{ij} = \mu_j + \sum_k c_{ik} \beta_k}
#'     }
#'
#'   \item \strong{Covariance estimation and hypothesis testing:} Estimate gene-gene covariance
#'         and compute test statistics to assess the statistical significance of gene co-expression using a third moment-based regression:
#'     \itemize{
#'       \item \strong{Covariance model:} \eqn{(x_{ij} - s_i \mu_{ij})(x_{ij'} - s_i \mu_{ij'}) = s_i^2 (\sigma_{jj'} + \sum_k c_{ik} \theta_k) + \xi_{ijj'}}
#'     }
#' }
#' We note that
#' \enumerate{
#'   \item The formulation above assumes that the covariates alter the mean / variance / covariance of underlying gene expression,
#'   rather than observed counts. If you believe that the covariates directly affect the observed counts independent of underlying gene expression
#'   (e.g. \eqn{x_{ij}|z_{ij} \sim \text{Poisson}(s_i z_{ij}+\sum_k c_{ik} \beta_k)} or \eqn{x_{ij} = s_i \mu_j  + \sum_k c_{ik} \beta_k + \epsilon_{ij}}),
#'   please specify \code{covariate_level="x"}.
#'   \item The original CS-CORE published in <https://doi.org/10.1038/s41467-023-40503-7>
#'   did not consider adjusting for covariates \eqn{c_{ik}}'s. This is equivalent to setting \code{covariates} to \code{NULL}.
#' }
#' @source Su, C., Xu, Z., Shan, X., Cai, B., Zhao, H., & Zhang, J. (2023).
#' Cell-type-specific co-expression inference from single cell RNA-sequencing data.
#' \emph{Nature Communications}.
#' doi: <https://doi.org/10.1038/s41467-023-40503-7>
#'
#' @param X A numeric matrix of UMI counts (\code{n x p}),
#'   where \code{n} is the number of cells and \code{p} is the number of genes.
#' @param seq_depth A numeric vector of sequencing depths of length \code{n}.
#' @param covariates Optional. A numeric matrix of covariates (\code{n x K}) to be adjusted for in the moment-based regressions.
#'   Can be a length n vector if \eqn{K=1}.
#'   If \code{NULL}, no covariates are adjusted for in the regression.
#'   Defaults to \code{NULL}.
#' @param post_process Optional. Logical; whether to rescale the estimated co-expressions to lie between –1 and 1.
#'   Defaults to \code{TRUE}.
#' @param covariate_level Optional. A character string indicating whether covariates are assumed to affect
#'   the underlying gene expression levels (\code{"z"}) or the observed counts (\code{"x"}).
#'   See the *Details* section for further explanation.
#'   Defaults to \code{"z"}.
#' @param adjust_setting Optional. A named logical vector of length 3;
#'   whether to adjust for covariates at the mean, variance, and covariance level.
#'   Must be named \code{c("mean", "var", "covar")}.
#'   Defaults to \code{c(mean = T,var = T, covar = T)}.
#' @param return_all Logical; whether to return all estimates, including the effect sizes for covariates.
#'   Defaults to \code{FALSE}.
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{est}{A \eqn{p \times p} matrix of co-expression estimates.}
#'   \item{p_value}{A \eqn{p \times p} matrix of p-values.}
#'   \item{test_stat}{A \eqn{p \times p} matrix of test statistics evaluating the statistical significance of co-expression.}
#'   \item{mu_beta}{A \eqn{k \times p} matrix of regression coefficients from the mean model. Returned if \code{return_all = TRUE}.}
#'   \item{sigma2_beta}{A \eqn{k \times p} matrix of regression coefficients from the variance model. Returned if \code{return_all = TRUE}.}
#'   \item{cov_beta}{A \eqn{k \times p \times p} array of regression coefficients from the covariance model. Returned if \code{return_all = TRUE}.}
#' }
#' @keywords internal
#'

CSCORE_IRLS <- function(X, seq_depth,
                        covariates = NULL,
                        post_process = TRUE,
                        covariate_level = 'z',
                        adjust_setting = c('mean' = T, 'var' = T, 'covar' = T),
                        return_all = FALSE){
  if(is.null(seq_depth)){
    seq_depth = apply(X, 1, sum, na.rm = T)
  }
  if(nrow(X) != length(seq_depth)){
    stop('The length of the sequencing depth must match the number of cells.')
  }
  if(has_non_integer(X)){
    stop("CS-CORE takes the UMI count matrix as input. The input matrix X contains non-integer values.")
  }
  if(!is.null(covariates)){
    if(is.vector(covariates)){
      if(nrow(X) != length(covariates)) stop('The length of the covariate vector should match the number of cells.')
    }else{
      if(nrow(X) != nrow(covariates)) stop('The number of rows in the covariate matrix should match the number of cells.')
      message(sprintf('The design matrix for regression has %i columns', ncol(covariates) + 1))
    }
  }
  n_cell = nrow(X)
  n_gene = ncol(X)
  seq_depth_sq = seq_depth^2
  seq_2 = sum(seq_depth_sq)
  seq_4 = sum(seq_depth^4)
  # Construct design matrix
  if(is.null(covariates)){
    D <- matrix(rep(1, n_cell), ncol = 1)
  }else{
    D <- cbind(1, covariates)
  }
  D_mu <- set_D(seq_depth, D, adjust_setting['mean'], covariate_level)
  D_sigma2 <- set_D(seq_depth_sq, D, adjust_setting['var'], covariate_level)
  D_sigma <- set_D(seq_depth_sq, D, adjust_setting['covar'], covariate_level)
  # initialize with OLS
  mu_beta <- WLS_mean(D_mu, X, matrix(1, nrow = n_cell, ncol = n_gene))
  M <- D_mu %*% mu_beta
  mu <- mu_beta[1,]
  X_centered = X - M
  sigma2_beta <- WLS_mean(D_sigma2, X_centered^2 - M, matrix(1, nrow = n_cell, ncol = n_gene))
  sigma2 <- sigma2_beta[1,]
  theta = mu^2/sigma2
  j = 0
  delta = Inf

  # IRLS for estimating mean and variance parameters
  while( delta > 0.05 & j <= 10 ){
    theta_previous = theta
    theta_median = stats::quantile(theta[theta > 0], na.rm = T, probs = 0.5)
    theta[theta < 0] = Inf
    w = M + outer(seq_depth_sq, mu^2/theta_median)
    w[is.na(w)|w <= 0] = 1
    mu_beta <- WLS_mean(D_mu, X, 1/w)
    mu <- mu_beta[1,]
    M <- D_mu %*% mu_beta
    X_centered = X - M
    h = (M^2/theta_median + M)^2
    h[h <= 0] = 1
    sigma2_beta <- WLS_mean(D_sigma2, X_centered^2 - M, 1/h)
    sigma2 <- sigma2_beta[1,]
    theta = mu^2/sigma2
    j = j+1
    delta = max(abs(log((theta/theta_previous)[theta > 0 & theta_previous > 0])), na.rm = T)
  }
  if(j == 10 & delta > 0.05){
    print('IRLS failed to converge after 10 iterations. Please check your data.')
  }else{
    print(sprintf('IRLS converged after %i iterations.', j))
  }

  # Update the weights for estimating covariance
  theta_median = stats::quantile(theta[theta > 0], na.rm = T, probs = 0.5)
  theta[theta < 0] = Inf
  w = M + outer(seq_depth_sq, mu^2/theta_median)
  w[is.na(w)|w <= 0] = 1

  X_centered <- X - M
  covar <- WLS_cov(D_sigma, X_centered, 1/w)$cov_hat[1,,]
  # Evaluate test statistics and p values
  Sigma <- M + outer(seq_depth_sq, sigma2)
  ele_inv_Sigma <- 1/Sigma
  ts_res <- WLS_cov(D_sigma, X_centered, ele_inv_Sigma)
  est <- ts_res$cov_hat[1,,]
  p_value <- 2 * stats::pnorm(abs(ts_res$test_stat), lower.tail = F)

  # Evaluate co-expression estimates
  neg_gene_inds <- which(sigma2 < 0)
  sigma2[neg_gene_inds] <- 0
  sigma <- sqrt(sigma2)
  diag(covar) <- sigma
  est <- covar/outer(sigma, sigma)
  test_stat <- ts_res$test_stat

  # Post-process the co-expression estimates
  #diag(est) <- 1
  if(post_process) est <- post_process_est(est)
  rownames(est) <- colnames(est) <- rownames(p_value) <- colnames(p_value) <- rownames(test_stat) <- colnames(test_stat) <- colnames(X)
  result_list <- list(est = est, p_value = p_value, test_stat = test_stat)
  if(return_all){
    result_list$mu_beta <- mu_beta
    result_list$sigma2_beta <- sigma2_beta
    result_list$cov_beta <- ts_res$cov_hat
  }
  return(result_list)
}

