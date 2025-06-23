test_that("The covaiates' effect sizes are computed correctly", {
  set.seed(202212)
  n <- 6000
  sim_seq_depths <- exp(rnorm(n, 7.3, sd=0.6))
  sim_seq_depths <- sim_seq_depths[!sim_seq_depths < 400]
  sim_seq_depths <- round(sim_seq_depths)

  n <- length(sim_seq_depths)
  mu <- 4.83 * 7.56e-05 # extracted from a specific gene
  sigma2 <- 4.83 * 7.56e-05^2
  # Simulate a covariate g
  g <- rbinom(n, 2, 0.4)
  effect_sizes <- matrix(c(0.05, 0.1, 0.2, 0.04), nrow = 2)
  rownames(effect_sizes) <- c('mu', 'sigma2')
  effect_sizes
  covar <- sigma2 * 0.1
  covar_effect_size <- 1

  z_mat <- matrix(nrow = n, ncol = 2)
  for(g_val in unique(g)){
    g_val_inds <- which(g == g_val)
    covar_g <- covar + covar * covar_effect_size * g_val
    sigma2_val_vec <- numeric(2)
    for(j in 1:2){
      sigma2_val_vec[j] <- sigma2 + sigma2 * effect_sizes['sigma2', j] * g_val
    }
    rho_g <- covar_g / sqrt(prod(sigma2_val_vec))
    z_mat_gaussian <- MASS::mvrnorm(length(g_val_inds), rep(0,2),
                              matrix(c(1, rho_g, rho_g, 1), nrow = 2))

    for(j in 1:2){
      mu_val <- mu + mu * effect_sizes['mu', j] * g_val
      sigma2_val <- sigma2 + sigma2 * effect_sizes['sigma2', j] * g_val
      z_mat[g_val_inds,j] <- qgamma(pnorm(z_mat_gaussian[,j]), shape = mu_val^2 / sigma2_val, scale = sigma2_val / mu_val)
    }
  }
  x_mat <- matrix(rpois(n = 2*n, lambda=z_mat * sim_seq_depths), nrow = n, ncol = 2)

  cscore_example <- CSCORE_IRLS(x_mat, sim_seq_depths,
                                covariates = g, return_all = TRUE)
  expect_equal(cscore_example$est[1,2], 0.1068046, tolerance = 1e-6)
  expect_equal(cscore_example$p_value[1,2], 0.5056709, tolerance = 1e-6)
  expect_equal(cscore_example$mu_beta,
               matrix(c(3.654416e-04, 3.484695e-04, 1.090923e-05, 7.330102e-05), 2, 2, byrow=T), tolerance = 1e-6)
  expect_equal(cscore_example$sigma2_beta,
               matrix(c(2.483174e-08, 2.850686e-08, 1.477678e-10, -6.450020e-09), 2, 2, byrow=T), tolerance = 1e-6)
  expect_equal(cscore_example$cov_beta[,,1],
               matrix(c(1.861900e-07, 2.855003e-09, 2.254392e-09, 8.866591e-09), 2, 2, byrow=T), tolerance = 1e-6)
})
