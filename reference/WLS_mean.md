# Weighted Least Squares Mean and Variance Estimation (Rcpp)

Compute WLS estimates of mean or variance for all p genes

## Usage

``` r
WLS_mean(D, X, W)
```

## Arguments

- D:

  Design matrix for WLS (size n x k), where the first column represents
  baseline mean/var, and others represents covariates.

- X:

  Response matrix for WLS (size n x p)

- W:

  Weight matrix (n x p)

## Value

A k x p matrix of WLS estimates
