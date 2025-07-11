test_that("CSCORE_IRLS Rcpp version matches base version on a simulated independent gene pair", {
  expect_equal(CSCORE_IRLS_cpp(ind_gene_pair$counts, ind_gene_pair$seq_depths,
                               IRLS_par = list(n_iter = 10, eps = 0.05, verbose = FALSE, conv = 'max')),
               CSCORE_IRLS_base(ind_gene_pair$counts, ind_gene_pair$seq_depths), tolerance = 1e-8)
})

test_that("CSCORE_IRLS Rcpp version matches base version on real data from 3k cells and on 500 genes", {
  test_dat = readRDS("fixtures/mic_highly_expressed_500_AD_control")
  count = test_dat[[1]]
  seq_depth = test_dat[[2]]
  expect_equal(CSCORE_IRLS_cpp(count, seq_depth,
                               IRLS_par = list(n_iter = 10, eps = 0.05, verbose = FALSE, conv = 'max')),
               CSCORE_IRLS_base(count, seq_depth), tolerance = 1e-8)
})

test_that("CSCORE_IRLS Rcpp version matches base version on real data from 3k cells and on 500 genes", {
  test_dat = readRDS("fixtures/mic_highly_expressed_500_AD_control")
  count = test_dat[[1]]
  seq_depth = test_dat[[2]]
  expect_equal(CSCORE_IRLS_cpp(count, seq_depth,
                               IRLS_par = list(n_iter = 10, eps = 0.05, verbose = FALSE, conv = 'max')),
               CSCORE_IRLS_base(count, seq_depth), tolerance = 1e-8)
})
