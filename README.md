# CS-CORE
[![R-CMD-check](https://github.com/ChangSuBiostats/CS-CORE/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ChangSuBiostats/CS-CORE/actions/workflows/R-CMD-check.yaml)
![Windows CI](https://github.com/ChangSuBiostats/CS-CORE/actions/workflows/windows-check.yaml/badge.svg)
[![DOI](https://zenodo.org/badge/576326164.svg)](https://zenodo.org/badge/latestdoi/576326164)

`CS-CORE` is a R package for cell-type-specific co-expression inference from single cell RNA-sequencing data. 

## Reference
Su, Chang, et al. "Cell-type-specific co-expression inference from single cell RNA-sequencing data." *Nature Communications* 14.1 (2023): 4846. (https://doi.org/10.1038/s41467-023-40503-7)

## Installation

`CS-CORE` is available on GitHub. You can install it using

``` r
## Load devtools for installing R packages from Github
library(devtools)
## Install CS-CORE from Github
install_github("ChangSuBiostats/CS-CORE")
```

## Vignettes

1. [Get started](https://changsubiostats.github.io/CS-CORE/articles/CSCORE.html) shows an example of using CS-CORE for cell-type-specific co-expression analysis with single cell RNA-sequencing data. 
It includes inferring co-expressions, extracting co-expressed gene modules and functional enrichment analysis.

2. [Covariate adjustment](https://changsubiostats.github.io/CS-CORE/articles/covariate_adjustment.html) shows how to adjust for 
covariates in co-expression inference with CS-CORE.


## Contact us 

For issues or feature requests, please visit GitHub Issues. If an issue remains unanswered for a while, 
you are welcome to email the maintainer at <chang.su@emory.edu>

## A Python version

A Python implementation of CS-CORE is also provided [here](https://github.com/ChangSuBiostats/CS-CORE_python).
