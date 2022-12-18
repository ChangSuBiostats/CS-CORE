test_that("correct dimension", {
  test_dat = readRDS("fixtures/mic_highly_expressed_500_AD_control")
  count = test_dat[[1]]
  seq_depth = test_dat[[2]]
  CSCORE_result = CSCORE_IRLS(X = count, seq_depth = seq_depth)
  rho = CSCORE_result$est
  pvalues = CSCORE_result$p_value
  test_stats = CSCORE_result$test_stat
  p = c(dim(count)[2], dim(count)[2])
  expect_true(all(dim(rho) == p,
                  dim(pvalues) == p,
                  dim(test_stats) == p))
})


test_that("pvalues are in the correct range", {
  test_dat = readRDS("fixtures/mic_highly_expressed_500_AD_control")
  count = test_dat[[1]]
  seq_depth = test_dat[[2]]
  CSCORE_result = CSCORE_IRLS(X = count, seq_depth = seq_depth)
  pvalues = CSCORE_result$p_value
  expect_true(mean(pvalues <= 1, na.rm = T) == 1)
})

test_that("sequencing depth must be equal to the number of cells", {
  test_dat = readRDS("fixtures/mic_highly_expressed_500_AD_control")
  count = test_dat[[1]]
  seq_depth = c(test_dat[[2]], 0)
  expect_error(CSCORE_result = CSCORE_IRLS(X = count, seq_depth = seq_depth))
})

test_that("The values estimates and p values are computed correctly", {
  cscore_example <- CSCORE_IRLS(ind_gene_pair$counts, ind_gene_pair$seq_depths)
  expect_equal(cscore_example$est[1,2], 0.0078201236)
  expect_equal(cscore_example$p_value[1,2], 0.96198097)
})
