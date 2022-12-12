[![R-CMD-check](https://github.com/ChangSuBiostats/CS-CORE/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ChangSuBiostats/CS-CORE/actions/workflows/R-CMD-check.yaml)

# CS-CORE
`CS-CORE` is a R package for cell-type-specific co-expression inference from single cell RNA-sequencing data.

## Installation

`CS-CORE` is available on github. You can install it using

``` r
## Load devtools for installing R packages from Github
library(devtools)
## Install CS-CORE from Github
install_github("ChangSuBiostats/CS-CORE")
```

## A toy example

We present here a toy example of inferring cell-type-specific co-expression networks with CS-CORE. 

We use the single cell RNA sequencing data on Peripheral blood mononuclear cells (PBMC) from COVID patients and healthy controls from [Wilk et al.](https://www.nature.com/articles/s41591-020-0944-y), which can be downloaded via

```
wget https://hosted-matrices-prod.s3-us-west-2.amazonaws.com/Single_cell_atlas_of_per
ipheral_immune_response_to_SARS_CoV_2_infection-25/blish_covid.seu.rds
```

``` r
## Load CS-CORE
library(CSCORE) 

## Load single cell data
library(Seurat)
pbmc_covid <- readRDS('data/blish_covid.seu.rds')

## Specify the cell type to study
# Here, we select B cells from healthy control subjects 
# to infer B cell-specific co-expression networks in healthy control subjects
cells_selected <- (pbmc_covid$cell.type.coarse %in% 'B') & (pbmc_covid$Status == "Healthy")

## Specify the genes for which co-expression networks are inferred
# Here, we select the top 2,000 highly expressed genes to infer the co-expression network
pbmc_covid_B <- pbmc_covid[,cells_selected] 
pbmc_covid_B <- NormalizeData(pbmc_covid_B)
mean_exp <- rowMeans(GetAssayData(pbmc_covid_B, slot = 'data'))
genes_selected <- names(sort.int(mean_exp, decreasing = T))[1:2000]

## Run CS-CORE
cscore_network <- CSCORE(pbmc_covid_B, genes = genes_selected)
# CS-CORE returns a list of three 2000*2000 matrices:
# co-expression estimates, test statistics and p values
str(cscore_network) 
dim(cscore_network[[1]])
```

## Vignette
Coming soon!

## Contact us 

Feel free to contact <c.su@yale.edu> for any questions on applying CS-CORE to your own singe cell data. 

## Reference
Coming soon!
