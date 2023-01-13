#' Iteratively reweighted least squares (IRLS) procedure in CS-CORE
#'
#' The iteratively reweighted least squares procedure in CS-CORE for estimating and testing cell-type-specific co-expressions using single cell RNA sequencing data.
#' More details on this procedure can be found in the CS-CORE paper: <https://doi.org/10.1101/2022.12.13.520181>.
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
#' cscore_example <- CSCORE_IRLS(ind_gene_pair$counts, ind_gene_pair$seq_depths)
#'
#' ## Estimated co-expression between two genes
#' cscore_example$est[1,2]
#' # close to 0: 0.007820124
#'
#' ## p-values
#' cscore_example$p_value[1,2]
#' # not significant: 0.961981
#'
#' @source
#' Cell-type-specific co-expression inference from single cell RNA-sequencing data
#' Chang Su, Zichun Xu, Xinning Shan, Biao Cai, Hongyu Zhao, Jingfei Zhang;
#' bioRxiv 2022.12.13.520181; doi: <https://doi.org/10.1101/2022.12.13.520181>
CSCORE_IRLS <- function(X, seq_depth, post_process = TRUE){
  if (is.null(seq_depth)) {
    seq_depth = apply(X, 1, sum, na.rm = T)
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

  neg_gene_inds <- which(sigma2 < 0)
  sigma2[neg_gene_inds] <- 0

  # Evaluate test statistics and p values
  Sigma <- M + outer(seq_depth_sq, sigma2)
  ele_inv_Sigma <- 1/Sigma
  X_centered_scaled <- X_centered * ele_inv_Sigma
  num <- t(seq_depth_sq * X_centered_scaled) %*% X_centered_scaled
  deno <- sqrt(t(seq_depth^4 * ele_inv_Sigma) %*% ele_inv_Sigma)
  test_stat <- num/deno
  p_value <- 2 * stats::pnorm(abs(test_stat), lower.tail = F)

  # Evaluate co-expression estimates
  sigma <- sqrt(sigma2)
  est <- covar/outer(sigma, sigma)
  # diag(est)[!is.na(diag(est))] <- 1

  # Post-process the co-expression estimates
  if(post_process) est <- post_process_est(est)
  return(list(est = est, p_value = p_value, test_stat = test_stat))
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
  neg_gene_inds <- which(is.infinite(diag(est)))
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
