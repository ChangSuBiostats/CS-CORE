# CS-CORE for cell-type-specific co-expression network inference

Run CS-CORE on a Seurat object to infer the cell-type-specific
co-expression network for a specified set of genes, with optional
adjustment for covariates. For more details on the covariate adjustment
and the moment-based regression, please refer to
[CSCORE_IRLS](https://changsubiostats.github.io/CS-CORE/reference/CSCORE_IRLS.md).

## Usage

``` r
CSCORE(
  object,
  genes,
  seq_depth = NULL,
  covariate_names = NULL,
  adjust_setting = c(mean = T, var = T, covar = T),
  IRLS_version = "Rcpp",
  IRLS_par = list(n_iter = 10, eps = 0.05, verbose = FALSE)
)
```

## Arguments

- object:

  A Seurat object containing single-cell RNA-seq data. The object should
  be subsetted to cells of a single cell type to ensure
  cell-type-specific inference. CS-CORE requires raw UMI counts as
  input, and assumes that the raw count matrix is stored in the
  `"counts"` slot of the `"RNA"` assay (i.e., `object[["RNA"]]@counts`).

- genes:

  A character vector of gene names (length \\p\\) for which the
  co-expression network will be estimated.

- seq_depth:

  A numeric vector of sequencing depths (length \\n\\). If `NULL`,
  sequencing depth will be computed as the total UMI count per cell.
  Defaults to `NULL`.

- covariate_names:

  Optional. A character vector specifying the names of cell-level
  covariates to adjust for in the regression models. These variables
  will be extracted from `object@meta.data[, covariate_names]`. Defaults
  to `NULL`.

- adjust_setting:

  Optional. A named logical vector of length 3 indicating whether to
  adjust for covariates in the estimation of mean, variance, and
  covariance. Must be named `c("mean", "var", "covar")`. Defaults to
  `c(mean = TRUE, var = TRUE, covar = TRUE)`.

- IRLS_version:

  Optional. A character string specifying the IRLS implementation to
  use: `"Rcpp"` or `"base_R"`. Only the `"Rcpp"` version supports
  covariate adjustment. The `"base_R"` version does not. When
  applicable, `"Rcpp"` offers improved memory efficiency (~10-100 times)
  but may be slower (~10 times), while `"base_R"` is faster but more
  memory intensive. Defaults to `"Rcpp"`.

- IRLS_par:

  Optional. A named list of length 3 specifying parameters for the IRLS
  algorithm:

  `n_iter`

  :   Maximum number of iterations.

  `eps`

  :   Convergence threshold for log-ratio change `delta`, computed as
      `abs(log(theta / theta_prev))`.

  `verbose`

  :   Logical; whether to print the convergence metric (`delta`) at each
      iteration.

  Defaults to `list(n_iter = 10, eps = 0.05, verbose = FALSE)`.

## Value

A list of three p by p matrices:

- est:

  Matrix of co-expression estimates.

- p_value:

  Matrix of p-values for testing co-expression.

- test_stat:

  Matrix of test statistics for evaluating the significance of
  co-expression.

## See also

[CSCORE_IRLS](https://changsubiostats.github.io/CS-CORE/reference/CSCORE_IRLS.md)

[CS-CORE online
tutorial](https://changsubiostats.github.io/CS-CORE/articles/CSCORE.html)

## Examples

``` r
# See a full example at:
# https://changsubiostats.github.io/CS-CORE/articles/CSCORE.html
```
