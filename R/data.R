#' A simulated independent gene pair
#'
#' This data set is created only for illustrative purposes and is used to test CSCORE_IRLS.R.
#' The source code for generating this data set is in data-raw/ on Github.
#'
#' @format ## `ind_gene_pair`
#' A list with 2 elements:
#' \describe{
#'   \item{counts}{A 1972 by 2 count matrix with 1972 cells and 2 genes}
#'   \item{seq_depths}{A length 1972 vector of sequencing depths}
#' }
"ind_gene_pair"

#' Top GO enrichment result in the CS-CORE vignette
#'
#' This data set stores the top GO enrichment result from the analysis in CS-CORE vignette.
#' It is saved for efficiently rendering the vignette without re-evaluating the codes.
#'
#' @format ## `top_enrich_go`
#' A list with five elements, where each element holds the top 3 GO terms enriched in a gene module with strong enrichment signals
#'
"top_enrich_go"

