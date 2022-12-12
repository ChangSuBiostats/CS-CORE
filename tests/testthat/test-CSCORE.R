test_that("input must be a Seurat object", {
  test_dat = readRDS("fixtures/mic_highly_expressed_500_AD_control")
  count = test_dat[[1]]
  genes = colnames(count)
  expect_error(CSCORE(count, genes))
})

