# Differential co-expression analysis with CS-CORE

Test for differential co-expression between two groups using a
permutation test. CS-CORE is run on each group separately to obtain
co-expression estimates. Group labels are then permuted `n_permu` times
at the biological sample level to build a null distribution of
differences in co-expression estimates. Two-sided permutation p-values
are returned for each gene pair.

## Usage

``` r
diff_CSCORE(
  object,
  group_label,
  sample_label,
  genes,
  group_levels = NULL,
  n_permu = 100,
  seed = NULL,
  n_cores = 1L,
  verbose = FALSE,
  ...
)
```

## Arguments

- object:

  A Seurat object containing single-cell RNA-seq data, subsetted to a
  single cell type.

- group_label:

  A character string giving the name of the `object@meta.data` column
  that contains the group labels (e.g. `"Status"`).

- sample_label:

  A character string giving the name of the `object@meta.data` column
  that identifies biological samples / donors (e.g. `"Donor"`).
  Permutations are performed at the sample level to preserve
  within-sample cell correlations.

- genes:

  A character vector of gene names (length \\p\\) for which
  co-expression will be estimated.

- group_levels:

  Optional. A length-2 character vector specifying the two group labels
  and their order, e.g. `c("Healthy", "COVID")`. The observed difference
  is computed as group 2 minus group 1. Defaults to
  `sort(unique(groups))`.

- n_permu:

  Integer; number of permutations. Defaults to `100`.

- seed:

  Optional integer seed passed to
  [`set.seed()`](https://rdrr.io/r/base/Random.html) before selecting
  permutation combinations. Defaults to `NULL`.

- n_cores:

  Integer; number of cores for parallelizing the permutation loop via
  [`parallel::mclapply`](https://rdrr.io/r/parallel/mclapply.html).
  Defaults to `1L` (sequential). Values greater than 1 are not supported
  on Windows.

- verbose:

  Logical; if `FALSE` (default), suppresses all output from the
  underlying `CSCORE` calls and prints only a progress bar for the
  permutation loop. If `TRUE`, prints detailed `[INFO]` messages from
  each step.

- ...:

  Additional arguments passed to
  [`CSCORE`](https://changsubiostats.github.io/CS-CORE/reference/CSCORE.md).

## Value

A list with four elements:

- est_group1:

  A \\p \times p\\ matrix of co-expression estimates for group 1.

- est_group2:

  A \\p \times p\\ matrix of co-expression estimates for group 2.

- obs_diff:

  A \\p \times p\\ matrix of observed differences in co-expression
  (`est_group2 - est_group1`).

- p_value:

  A \\p \times p\\ matrix of two-sided permutation p-values. Diagonal
  entries are `NA`.

## See also

[Differential co-expression
tutorial](https://changsubiostats.github.io/CS-CORE/articles/differential_coexpression.html)

## Examples

``` r
# See a full example at:
# https://changsubiostats.github.io/CS-CORE/articles/differential_coexpression.html
```
