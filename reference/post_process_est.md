# Post-process IRLS estimates

The IRLS procedure does not guarantee the variance estimates to be
postive nor the co-expression parameters to be bounded. To address this,
this function evaluates the percentage of genes with negative variance
estimates; sets the their co-expressions to 0 as these genes do not have
sufficient biological variations. This function also evaluates the
percentage of gene pairs with out-of-bound co-expression estimates; sets
the co-expressions greater than 1 to 1; set the co-expressions smaller
than -1 to -1.

## Usage

``` r
post_process_est(est)
```

## Arguments

- est:

  Estimated co-expression matrix from IRLS

## Value

Post-processed correlation matrix
