# CS-CORE
[![R-CMD-check](https://github.com/ChangSuBiostats/CS-CORE/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ChangSuBiostats/CS-CORE/actions/workflows/R-CMD-check.yaml)

`CS-CORE` is a R package for cell-type-specific co-expression inference from single cell RNA-sequencing data. It provides an implementation for the statistical method CS-CORE proposed in this [manuscript](https://doi.org/10.1101/2022.12.13.520181).

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

Feel free to contact Chang Su at <c.su@yale.edu> for questions on applying CS-CORE to your own singe cell data. 

## Reference
**Cell-type-specific co-expression inference from single cell RNA-sequencing data.**
Chang Su, Zichun Xu, Xinning Shan, Biao Cai, Hongyu Zhao, Jingfei Zhang.
bioRxiv 2022.12.13.520181; doi: https://doi.org/10.1101/2022.12.13.520181
