# Weighted Least Squares Covariance Estimation (Rcpp)

Computes WLS estimates of gene-gene covariance and test statistics that
assess the statistical significance of co-expression using moment-based
regressions.

## Usage

``` r
WLS_cov(D, X, W)
```

## Arguments

- D:

  Design matrix for WLS (size n x k)

- X:

  Gene expression matrix (n x p)

- W:

  Weight matrix (n x p)

## Value

A list with:

- cov_hat:

  A k x p x p array of WLS estimates

- test_stat:

  A p x p matrix of test statistics
