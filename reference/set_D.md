# Set the design matrix for moment-based regressions

Set the design matrix for moment-based regressions

## Usage

``` r
set_D(s, D, adjust_setting, covariate_level)
```

## Arguments

- s:

  A numeric vector of sequencing depths (for mean regression) or squared
  sequencing depths (for variance and covariance)

- D:

  A numeric matrix of intercept and covariates (`n x K`)

- adjust_setting:

  Logical; whether to adjust for covariates

- covariate_level:

  A character string indicating whether covariates are assumed to affect
  the underlying gene expression levels (`"z"`) or the observed counts
  (`"x"`).

## Value

Design matrix (n by K) for moment-based regressions
