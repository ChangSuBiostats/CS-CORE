#' Iteratively reweighted least squares (IRLS) procedure in CS-CORE (base R, archived)
#'
#' This function was originally implemented in 2023 and included in the first release of the CS-CORE R package.
#' It has since been replaced by \code{\link{CSCORE_IRLS}} for two main reasons:
#'
#' \enumerate{
#'   \item It relies on base R for regression, which is slower and more memory-intensive than the Rcpp-based implementation in \code{CSCORE_IRLS}.
#'   \item It does not support covariate adjustment, a feature supported in \code{CSCORE_IRLS}.
#' }
#'
#' @note This function is retained for reference and backward compatibility, but users are encouraged to use \code{\link{CSCORE_IRLS}} for new analyses.
#'
#' @keywords internal
#'
#' @param X A n by p matrix of UMI counts, where n denotes the number of cells and p denotes the number of genes
#' @param seq_depth A length n vector of sequencing depths
#' @param post_process Whether to process the estimated co-expressions such that the estimates are between -1 and 1. Default to TRUE.
#'
#' @return A list of three p by p matrices:
#' \describe{
#'   \item{est}{co-expression estimates}
#'   \item{p_value}{p values}
#'   \item{test_stat}{test statistics}
#' }
#' @export
#'
#' @examples
#' ## Toy example:
#' ## run CSCORE on a simulated independent gene pair
#' cscore_example <- CSCORE_IRLS_base(ind_gene_pair$counts, ind_gene_pair$seq_depths)
#'
#' ## Estimated co-expression between two genes
#' cscore_example$est[1,2]
#' # close to 0: 0.007820124
#'
#' ## p-values
#' cscore_example$p_value[1,2]
#' # not significant: 0.961981
#'
#' @source Su, C., Xu, Z., Shan, X., Cai, B., Zhao, H., & Zhang, J. (2023).
#' Cell-type-specific co-expression inference from single cell RNA-sequencing data.
#' \emph{Nature Communications}.
#' doi: <https://doi.org/10.1038/s41467-023-40503-7>
#'
CSCORE_IRLS_base <- function(X, seq_depth, post_process = TRUE){
  if (is.null(seq_depth)) {
    seq_depth = apply(X, 1, sum, na.rm = T)
  }
  if(nrow(X) != length(seq_depth)){
    stop('The length of the sequencing depth must match the number of cells.')
  }
  n_cell = nrow(X)
  n_gene = ncol(X)
  seq_depth_sq = seq_depth^2
  seq_2 = sum(seq_depth_sq)
  seq_4 = sum(seq_depth^4)
  mu = colSums(X * seq_depth)/seq_2
  M = outer(seq_depth, mu)
  X_centered = X - M
  sigma2 = colSums(((X_centered^2 - M) * seq_depth_sq))/seq_4
  theta = mu^2/sigma2
  j = 0
  delta = Inf

  while( delta > 0.05 & j <= 10 ){
    theta_previous = theta
    theta_median = stats::quantile(theta[theta > 0], na.rm = T, probs = 0.5)
    theta[theta < 0] = Inf
    w = M + outer(seq_depth_sq, mu^2/theta_median)
    w[is.na(w)|w <= 0] = 1
    mu = colSums((X/w) * seq_depth)/colSums(seq_depth_sq/w)
    M = outer(seq_depth, mu)
    X_centered = X - M
    h = (M^2/theta_median + M)^2
    h[h <= 0] = 1
    sigma2 = colSums(((X_centered^2 - M)/h * seq_depth_sq))/colSums(seq_depth_sq^2/h)
    theta = mu^2/sigma2
    j = j+1
    # print(paste0("Iteration: ", j, ", Median: ", theta_median))
    delta = max(abs(log((theta/theta_previous)[theta > 0 & theta_previous > 0])), na.rm = T)
    # print(delta)
  }
  if(j == 10 & delta > 0.05){
    print('IRLS failed to converge after 10 iterations. Please check your data.')
  }else{
    print(sprintf('IRLS converged after %i iterations.', j))
  }

  theta_median = stats::quantile(theta[theta > 0], na.rm = T, probs = 0.5)
  theta[theta < 0] = Inf
  w = M + outer(seq_depth_sq, mu^2/theta_median)
  w[is.na(w)|w <= 0] = 1

  covar = matrix(NA, nrow = n_gene, ncol = n_gene)
  covar <- (t(seq_depth_sq * X_centered/w) %*% (X_centered/w))/(t(seq_depth_sq/w) %*% (seq_depth_sq/w))

  # Evaluate test statistics and p values
  Sigma <- M + outer(seq_depth_sq, sigma2)
  ele_inv_Sigma <- 1/Sigma
  X_centered_scaled <- X_centered * ele_inv_Sigma
  num <- t(seq_depth_sq * X_centered_scaled) %*% X_centered_scaled
  deno <- sqrt(t(seq_depth^4 * ele_inv_Sigma) %*% ele_inv_Sigma)
  test_stat <- num/deno
  p_value <- 2 * stats::pnorm(abs(test_stat), lower.tail = F)

  # Evaluate co-expression estimates
  neg_gene_inds <- which(sigma2 < 0)
  sigma2[neg_gene_inds] <- 0
  sigma <- sqrt(sigma2)
  est <- covar/outer(sigma, sigma)

  # Post-process the co-expression estimates
  if(post_process) est <- post_process_est(est)
  return(list(est = est, p_value = p_value, test_stat = test_stat))
}
