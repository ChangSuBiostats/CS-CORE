# Iteratively reweighted least squares (IRLS) procedure in CS-CORE

This function implements the IRLS procedure used in CS-CORE for
estimating and testing cell-type-specific co-expression from single-cell
RNA sequencing data.

## Usage

``` r
CSCORE_IRLS(
  X,
  seq_depth,
  covariates = NULL,
  post_process = TRUE,
  covariate_level = "z",
  adjust_setting = c(mean = T, var = T, covar = T),
  return_all = FALSE
)
```

## Source

Su, C., Xu, Z., Shan, X., Cai, B., Zhao, H., & Zhang, J. (2023).
Cell-type-specific co-expression inference from single cell
RNA-sequencing data. *Nature Communications*. doi:
<https://doi.org/10.1038/s41467-023-40503-7>

## Arguments

- X:

  A numeric matrix of UMI counts (`n x p`), where `n` is the number of
  cells and `p` is the number of genes.

- seq_depth:

  A numeric vector of sequencing depths of length `n`.

- covariates:

  Optional. A numeric matrix of covariates (`n x K`) to be adjusted for
  in the moment-based regressions. Can be a length n vector if \\K=1\\.
  If `NULL`, no covariates are adjusted for in the regression. Defaults
  to `NULL`.

- post_process:

  Optional. Logical; whether to rescale the estimated co-expressions to
  lie between –1 and 1. Defaults to `TRUE`.

- covariate_level:

  Optional. A character string indicating whether covariates are assumed
  to affect the underlying gene expression levels (`"z"`) or the
  observed counts (`"x"`). See the *Details* section for further
  explanation. Defaults to `"z"`.

- adjust_setting:

  Optional. A named logical vector of length 3; whether to adjust for
  covariates at the mean, variance, and covariance level. Must be named
  `c("mean", "var", "covar")`. Defaults to
  `c(mean = T,var = T, covar = T)`.

- return_all:

  Logical; whether to return all estimates, including the effect sizes
  for covariates. Defaults to `FALSE`.

## Value

A list containing the following components:

- est:

  A \\p \times p\\ matrix of co-expression estimates.

- p_value:

  A \\p \times p\\ matrix of p-values.

- test_stat:

  A \\p \times p\\ matrix of test statistics evaluating the statistical
  significance of co-expression.

- mu_beta:

  A \\k \times p\\ matrix of regression coefficients from the mean
  model. Returned if `return_all = TRUE`.

- sigma2_beta:

  A \\k \times p\\ matrix of regression coefficients from the variance
  model. Returned if `return_all = TRUE`.

- cov_beta:

  A \\k \times p \times p\\ array of regression coefficients from the
  covariance model. Returned if `return_all = TRUE`.

## Details

Let \\x\_{ij}\\ denote the UMI count of gene \\j\\ in cell \\i\\;
\\s_i\\ denote the sequencing depth; \\\mu_j,\sigma\_{jj},
\sigma\_{jj'}\\ denote the mean, variance and covariance; \\c\_{ik}\\
denote additional covariate \\k\\ for cell \\i\\ (e.g. disease status or
cellular state). The procedure consists of two main steps:

1.  **Mean and variance estimation:** Estimate gene-specific mean and
    variance parameters using two moment-based regressions:

    - **Mean model:** \\x\_{ij} = s_i (\mu_j + \sum_k c\_{ik} \beta_k) +
      \epsilon\_{ij}\\

    - **Variance model:** \\(x\_{ij} - s_i \mu\_{ij})^2 = s_i
      \mu\_{ij} + s_i^2 (\sigma\_{jj} + \sum_k c\_{ik} \gamma_k) +
      \eta\_{ij}\\, where \\\mu\_{ij} = \mu_j + \sum_k c\_{ik} \beta_k\\

2.  **Covariance estimation and hypothesis testing:** Estimate gene-gene
    covariance and compute test statistics to assess the statistical
    significance of gene co-expression using a third moment-based
    regression:

    - **Covariance model:** \\(x\_{ij} - s_i \mu\_{ij})(x\_{ij'} - s_i
      \mu\_{ij'}) = s_i^2 (\sigma\_{jj'} + \sum_k c\_{ik} \theta_k) +
      \xi\_{ijj'}\\

We note that

1.  The formulation above assumes that the covariates alter the mean /
    variance / covariance of underlying gene expression, rather than
    observed counts. If you believe that the covariates directly affect
    the observed counts independent of underlying gene expression (e.g.
    \\x\_{ij}\|z\_{ij} \sim \text{Poisson}(s_i z\_{ij}+\sum_k c\_{ik}
    \beta_k)\\ or \\x\_{ij} = s_i \mu_j + \sum_k c\_{ik} \beta_k +
    \epsilon\_{ij}\\), please specify `covariate_level="x"`.

2.  The original CS-CORE published in
    <https://doi.org/10.1038/s41467-023-40503-7> did not consider
    adjusting for covariates \\c\_{ik}\\'s. This is equivalent to
    setting `covariates` to `NULL`.
