#' Iteratively reweighted least squares (IRLS) procedure in CS-CORE (Rcpp)
#'
#' Estimate and test cell-type-specific co-expression using the IRLS procedure with optional covariate adjustment.
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
#' Note: This is an R wrapper for the `CSCORE_IRLS_cpp_impl()` function implemented in Rcpp.
#'
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
#' @param IRLS_par Optional. A named list of length 4 specifying parameters for the IRLS algorithm:
#'   \describe{
#'     \item{\code{n_iter}}{Maximum number of iterations.}
#'     \item{\code{eps}}{Convergence threshold for log-ratio change \code{delta}, computed as \code{abs(log(theta / theta_prev))}.}
#'     \item{\code{verbose}}{Logical; whether to print the convergence metric (\code{delta}) at each iteration.}
#'     \item{\code{conv}}{Character string; determine convergence based on \code{q95} (0.95 quantile) or \code{max} of \code{delta}.}
#'   }
#'   Defaults to \code{list(n_iter = 10, eps = 0.05, verbose = FALSE, cov = "q95")}.
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
#' @export
#'
#' @examples
#' ## Toy example:
#' ## run CSCORE on a simulated independent gene pair
#' cscore_example <- CSCORE_IRLS_cpp(ind_gene_pair$counts, ind_gene_pair$seq_depths)
#'
#' ## Estimated co-expression between two genes
#' cscore_example$est[1,2]
#' # close to 0: 0.007820124
#'
#' ## p-values
#' cscore_example$p_value[1,2]
#' # not significant: 0.961981
#'

CSCORE_IRLS_cpp <- function(X,
                            seq_depth,
                            covariates = NULL,
                            post_process = TRUE,
                            covariate_level = 'z',
                            adjust_setting = c('mean' = T, 'var' = T, 'covar' = T),
                            IRLS_par = list('n_iter' = 10, 'eps' = 0.05, 'verbose' = FALSE, conv = 'q95'),
                            return_all = FALSE) {
  if(is.null(seq_depth)){
    seq_depth = apply(X, 1, sum, na.rm = T)
  }
  if(nrow(X) != length(seq_depth)){
    stop('[ERROR] The length of the sequencing depth must match the number of cells.')
  }
  if(! all.equal(names(IRLS_par), c('n_iter', 'eps', 'verbose', 'conv'))){
    stop('[ERROR] IRLS_par must be a named list with integer n_iter, double eps, and logical verbose.')
  }
  if(has_non_integer(X)){
    stop("[ERROR] CS-CORE takes the UMI count matrix as input. The input matrix X contains non-integer values.")
  }
  if(!is.null(covariates)){
    if(is.vector(covariates)){
      if(nrow(X) != length(covariates)) stop('[ERROR] The length of the covariate vector should match the number of cells.')
    }else{
      if(nrow(X) != nrow(covariates)) stop('[ERROR] The number of rows in the covariate matrix should match the number of cells.')
      #cat(sprintf('[INFO] The design matrix for regression has %i columns', ncol(covariates) + 1))
    }
  }
  if (length(adjust_setting) != 3) stop("[ERROR] adjust_setting must be of length 3.")

  # Construct design matrix
  n_cell <- nrow(X)
  if(is.null(covariates)){
    D <- matrix(rep(1, n_cell), ncol = 1)
  }else{
    D <- cbind(1, covariates)
  }
  seq_depth_sq <- seq_depth^2
  D_mu <- set_D(seq_depth, D, adjust_setting['mean'], covariate_level)
  D_sigma2 <- set_D(seq_depth_sq, D, adjust_setting['var'], covariate_level)
  D_sigma <- set_D(seq_depth_sq, D, adjust_setting['covar'], covariate_level)

  res <- .Call(`_CSCORE_CSCORE_IRLS_cpp_impl`,
               X, seq_depth_sq, D_mu, D_sigma2, D_sigma, # regression parameters
               post_process, # post processing
               IRLS_par[['n_iter']], IRLS_par[['eps']], IRLS_par[['verbose']], IRLS_par[['conv']], # IRLS parameters
               return_all)

  # Assign row/colnames to square matrices
  gene_names <- colnames(X)
  rownames(res$est) <- colnames(res$est) <- gene_names
  rownames(res$p_value) <- colnames(res$p_value) <- gene_names
  rownames(res$test_stat) <- colnames(res$test_stat) <- gene_names

  return(res)
}

#' Set the design matrix for moment-based regressions
#'
#'
#' @param s A numeric vector of sequencing depths (for mean regression) or squared sequencing depths (for variance and covariance)
#' @param D A numeric matrix of intercept and covariates (\code{n x K})
#' @param adjust_setting Logical; whether to adjust for covariates
#' @param covariate_level A character string indicating whether covariates are assumed to affect
#'   the underlying gene expression levels (\code{"z"}) or the observed counts (\code{"x"}).
#'
#' @return Design matrix (n by K) for moment-based regressions
#'
#' @keywords internal
#'
set_D <- function(s, D, adjust_setting, covariate_level){
  if(!adjust_setting){
    return(matrix(s, ncol = 1))
  }else{
    if(covariate_level == 'z'){
      return(s * D)
    }else if(covariate_level == 'x'){
      return(cbind(s, D[,-1]))
    }
  }
}


#' Check for non-integer values in a matrix
#'
#' This function checks whether a numeric matrix contains any non-integer values.
#'
#' @param mat A numeric matrix.
#'
#' @return A logical value: \code{TRUE} if any element is non-integer, otherwise \code{FALSE}.
#'
#' @keywords internal
#'
has_non_integer <- function(mat) {
  for (val in mat) {
    if (abs(val - round(val)) > .Machine$double.eps^0.5) return(TRUE)
  }
  return(FALSE)
}


#' Post-process IRLS estimates
#'
#' The IRLS procedure does not guarantee the variance estimates to be postive nor the co-expression parameters to be bounded.
#' To address this, this function evaluates the percentage of genes with negative variance estimates;
#' sets the their co-expressions to 0 as these genes do not have sufficient biological variations.
#' This function also evaluates the percentage of gene pairs with out-of-bound co-expression estimates;
#' sets the co-expressions greater than 1 to 1; set the co-expressions smaller than -1 to -1.
#'
#' @param est Estimated co-expression matrix from IRLS
#'
#' @return Post-processed correlation matrix
#'
#' @export
#'
post_process_est <- function(est){
  p <- nrow(est)
  # Post-process CS-CORE estimates
  neg_gene_inds <- which(sapply(diag(est), function(x) is.infinite(x) | is.na(x)))
  #which(is.infinite(diag(est)))
  if(length(neg_gene_inds) > 0){
    print(sprintf('%i among %i genes have negative variance estimates. Their co-expressions with other genes were set to 0.',
                  length(neg_gene_inds), p))
  }
  # Negative variances suggest insufficient biological variation,
  # and also lack of correlation
  est[neg_gene_inds, ] <- 0
  est[, neg_gene_inds] <- 0
  # Set all diagonal values to 1
  diag(est) <- 1
  # Gene pairs with out-of-bound estimates
  print(sprintf('%.4f%% co-expression estimates were greater than 1 and were set to 1.',
                mean(est[upper.tri(est)] > 1, na.rm = T) * 100))
  print(sprintf('%.4f%% co-expression estimates were smaller than -1 and were set to -1.',
                mean(est[upper.tri(est)] < -1, na.rm = T) * 100))
  est[est > 1] <- 1
  est[est < -1] <- -1
  return(est)
}

