## code to prepare `ind_gene_pair` dataset goes here

# Generate highly variable seq depths from a log normal distribution
set.seed(202212)
sim_seq_depths <- exp(rnorm(2000, 7.3, sd=0.6))
sim_seq_depths <- sim_seq_depths[!sim_seq_depths < 400]
sim_seq_depths <- round(sim_seq_depths)
summary(sim_seq_depths)

# Generate an independent gene pair
n <- length(sim_seq_depths)
z_mat <- matrix(rgamma(2*n, shape=4.83, scale=7.56e-05),
                nrow = n, ncol = 2)
x_mat <- matrix(rpois(n = 2*n, lambda=z_mat * sim_seq_depths), nrow = n, ncol = 2)
summary(x_mat[,1])
summary(x_mat[,2])

ind_gene_pair <- list(counts = x_mat, seq_depths = sim_seq_depths)

usethis::use_data(ind_gene_pair, overwrite = TRUE)
