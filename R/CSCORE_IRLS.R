CSCORE_IRLS <- function(X, seq_depth = NULL,  probs = 0.5, covar_weight = "regularized"){
  if (is.null(seq_depth)) {
    seq_depth = apply(X, 1, sum, na.rm = T)
  }
  n_cell = nrow(X)
  n_gene = ncol(X)
  seq_depth_sq = seq_depth^2
  seq_2 = sum(seq_depth_sq)
  seq_4 = sum(seq_depth^4)
  mu = colSums(X*seq_depth)/seq_2
  M = outer(seq_depth, mu)
  X_centered = X - M
  sigma2 = colSums(((X_centered^2 - M) * seq_depth_sq))/seq_4
  theta = mu^2/sigma2
  j = 0
  error = Inf

  while( error > 0.05 & j <= 10 ){
    theta_previous = theta
    theta_median = quantile(theta[theta > 0], na.rm = T, probs = probs)
    theta[theta < 0] = Inf
    w1 = M + outer(seq_depth_sq, mu^2/theta_median)
    w1[is.na(w1)|w1 <= 0] = 1
    mu = colSums((X/w1)*seq_depth)/colSums(seq_depth_sq/w1)
    M = outer(seq_depth, mu)
    X_centered = X - M
    w2 = (M^2/theta_median + M)^2
    w2[w2 <= 0] = 1
    sigma2 = colSums(((X_centered^2 - M)/w2 * seq_depth_sq))/colSums(seq_depth_sq^2/w2)
    theta = mu^2/sigma2
    j = j+1
    print(paste0("Iteration: ", j, ", Median: ", theta_median))
    error = max(abs(log((theta/theta_previous)[theta > 0 & theta_previous > 0])), na.rm = T)
    print(error)
  }

  theta_median = quantile(theta[theta > 0], na.rm = T, probs = probs)
  theta[theta < 0] = Inf
  if( covar_weight == "regularized" ){
    w1 = M + outer(seq_depth_sq, mu^2/theta_median)
  } else{
    w1 = M + outer(seq_depth_sq, mu^2/theta)
  }
  w1[is.na(w1)|w1 <= 0] = 1

  covar = matrix(NA, nrow = n_gene, ncol = n_gene)
  covar <- (t(seq_depth_sq * X_centered/w1) %*% (X_centered/w1))/(t(seq_depth_sq/w1) %*% (seq_depth_sq/w1))

  sigma <- sqrt(sigma2)
  est <- covar/outer(sigma, sigma)
  diag(est)[!is.na(diag(est))] <- 1
  Sigma <- M + outer(seq_depth_sq, sigma2)
  ele_inv_Sigma <- 1/Sigma
  X_centered_scaled <- X_centered * ele_inv_Sigma
  num <- t(seq_depth_sq * X_centered_scaled) %*% X_centered_scaled
  deno <- sqrt(t(seq_depth^4 * ele_inv_Sigma) %*% ele_inv_Sigma)
  test_stat <- num/deno
  p_value <- 2 * pnorm(abs(test_stat), lower.tail = F)
  return(list(est = est, p_value = p_value, test_stat = test_stat))
}
