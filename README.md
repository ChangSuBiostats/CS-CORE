# CS-CORE
[![R-CMD-check](https://github.com/ChangSuBiostats/CS-CORE/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ChangSuBiostats/CS-CORE/actions/workflows/R-CMD-check.yaml)
[![DOI](https://zenodo.org/badge/576326164.svg)](https://zenodo.org/badge/latestdoi/576326164)

`CS-CORE` is a R package for cell-type-specific co-expression inference from single cell RNA-sequencing data. It provides an implementation for the statistical method CS-CORE proposed in this [paper](https://doi.org/10.1038/s41467-023-40503-7).

## Installation

`CS-CORE` is available on GitHub. You can install it using

``` r
## Load devtools for installing R packages from Github
library(devtools)
## Install CS-CORE from Github
install_github("ChangSuBiostats/CS-CORE")
```

## Vignettes

The [vignette](https://changsubiostats.github.io/CS-CORE/articles/CSCORE.html) shows an example of using CS-CORE for cell-type-specific co-expression analysis with single cell RNA-sequencing data, which includes inferring co-expressions, extracting co-expressed gene modules and functional enrichment analysis as performed in our manuscript.

## Contact us 

[Chang Su](www.changsu.org), <chang.su@emory.edu>

## Reference and more

**Cell-type-specific co-expression inference from single cell RNA-sequencing data.**
Su, Chang, et al. "Cell-type-specific co-expression inference from single cell RNA-sequencing data." Nature Communications 14.1 (2023): 4846.

A Python implementation of CS-CORE is also provided [here](https://github.com/ChangSuBiostats/CS-CORE_python).
