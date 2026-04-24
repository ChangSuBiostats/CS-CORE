# Iteratively reweighted least squares (IRLS) procedure in CS-CORE (base R, archived)

This function was originally implemented in 2023 and included in the
first release of the CS-CORE R package. It has since been replaced by
[`CSCORE_IRLS`](https://changsubiostats.github.io/CS-CORE/reference/CSCORE_IRLS.md)
for two main reasons:

## Usage

``` r
CSCORE_IRLS_base(
  X,
  seq_depth,
  post_process = TRUE,
  return_cov = FALSE,
  return_mu = FALSE
)
```

## Source

Su, C., Xu, Z., Shan, X., Cai, B., Zhao, H., & Zhang, J. (2023).
Cell-type-specific co-expression inference from single cell
RNA-sequencing data. *Nature Communications*. doi:
<https://doi.org/10.1038/s41467-023-40503-7>

## Arguments

- X:

  A n by p matrix of UMI counts, where n denotes the number of cells and
  p denotes the number of genes

- seq_depth:

  A length n vector of sequencing depths

- post_process:

  Whether to process the estimated co-expressions such that the
  estimates are between -1 and 1. Default to TRUE.

- return_cov:

  Whether to return covariance matrix (instead of correlation). Default
  to FALSE.

- return_mu_sigma_sq:

  Whether to return mu and sigma_sq estimates. Default to FLASE.

## Value

A list containing:

- est:

  A p x p matrix of co-expression estimates (if `return_cov = FALSE`) or
  covariance estimates (if `return_cov = TRUE`).

- p_value:

  A p x p matrix of p-values.

- test_stat:

  A p x p matrix of test statistics.

- mu:

  A length-p vector of estimated means (expectation of underlying gene
  expression). Returned only if `return_mu_sigma_sq = TRUE`.

- sigma_sq:

  A length-p vector of estimated variances (variance of underlying gene
  expression). Returned only if `return_mu_sigma_sq = TRUE`.

## Details

1.  It relies on base R for regression, which is slower and more
    memory-intensive than the Rcpp-based implementation in
    `CSCORE_IRLS`.

2.  It does not support covariate adjustment, a feature supported in
    `CSCORE_IRLS`.

## Note

This function is retained for reference and backward compatibility, but
users are encouraged to use
[`CSCORE_IRLS`](https://changsubiostats.github.io/CS-CORE/reference/CSCORE_IRLS.md)
for new analyses.

## Examples

``` r
## Toy example:
## run CSCORE on a simulated independent gene pair
cscore_example <- CSCORE_IRLS_base(ind_gene_pair$counts, ind_gene_pair$seq_depths)
#> [1] "IRLS converged after 2 iterations."
#> [1] "0.0000% co-expression estimates were greater than 1 and were set to 1."
#> [1] "0.0000% co-expression estimates were smaller than -1 and were set to -1."

## Estimated co-expression between two genes
cscore_example$est[1,2]
#> [1] 0.007820124
# close to 0: 0.007820124

## p-values
cscore_example$p_value[1,2]
#> [1] 0.961981
# not significant: 0.961981
```
