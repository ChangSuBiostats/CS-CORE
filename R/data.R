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

#' Gene modules from the differential co-expression vignette
#'
#' WGCNA module assignments from the differential co-expression analysis of
#' CD14 Monocytes (COVID-19 vs. healthy controls) in the CS-CORE vignette.
#' Saved for rendering the vignette without re-running the full analysis.
#'
#' @format ## `diff_module_list`
#' A list of character vectors, one per module, each containing the gene names assigned to that module.
#'
"diff_module_list"

#' GO enrichment results from the differential co-expression vignette
#'
#' GO enrichment results (from \code{clusterProfiler::enrichGO}) for all WGCNA modules
#' identified in the differential co-expression analysis of CD14 Monocytes
#' (COVID-19 vs. healthy controls).
#' Saved for rendering the vignette without re-running the full analysis.
#'
#' @format ## `diff_ego_result`
#' A list of \code{enrichResult} objects, one per module.
#'
"diff_ego_result"

#' Co-expression matrices for the example module in the differential co-expression vignette
#'
#' Subsetted co-expression matrices (healthy and COVID-19) for the example gene module
#' (enriched for defense response to virus) from the CS-CORE differential co-expression vignette.
#' Saved for rendering the vignette without re-running the full analysis.
#'
#' @format ## `diff_est_example`
#' A named list with two elements:
#' \describe{
#'   \item{Healthy}{Co-expression matrix for the example module in healthy controls.}
#'   \item{COVID}{Co-expression matrix for the example module in COVID-19 patients.}
#' }
#'
"diff_est_example"

