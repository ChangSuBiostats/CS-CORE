test_that("input must be a Seurat object", {
  library(Seurat)
  test_dat = readRDS("fixtures/mic_highly_expressed_500_AD_control")
  count = test_dat[[1]]
  genes = colnames(count)
  expect_error(CSCORE(count, genes))
})

test_that("sequencing depth is calculated correctly", {
  library(Seurat)
  test_dat = readRDS("fixtures/mic_highly_expressed_500_AD_control")
  count = test_dat[[1]]
  genes = colnames(count)
  obj = CreateSeuratObject(t(count))
  rho_1 = CSCORE(obj, genes)[[1]]
  rho_2 = CSCORE(obj, genes, rowSums(count))[[1]]
  expect_true(mean(rho_1 == rho_2, na.rm = T) == 1)
})

