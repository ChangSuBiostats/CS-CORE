#' CS-CORE for cell-type-specific co-expression network inference
#'
#' Run CS-CORE on a Seurat object to infer the cell-type-specific co-expression network for a specified set of genes,
#' with optional adjustment for covariates.
#' For more details on the covariate adjustment and the moment-based regression,
#' please refer to \link{CSCORE_IRLS}.
#' @seealso \link{CSCORE_IRLS}
#' @seealso \href{https://changsubiostats.github.io/CS-CORE/articles/CSCORE.html}{CS-CORE online tutorial}
#'
#' @param object A Seurat object containing single-cell RNA-seq data.
#'   The object should be subsetted to cells of a single cell type to ensure cell-type-specific inference.
#'   CS-CORE requires raw UMI counts as input, and assumes that the raw count matrix is stored in the \code{"counts"} slot of the \code{"RNA"} assay
#'   (i.e., \code{object[["RNA"]]@counts}).
#' @param genes A character vector of gene names (length \eqn{p}) for which the co-expression network will be estimated.
#' @param seq_depth A numeric vector of sequencing depths (length \eqn{n}).
#'   If \code{NULL}, sequencing depth will be computed as the total UMI count per cell.
#'   Defaults to \code{NULL}.
#' @param covariate_names Optional. A character vector specifying the names of cell-level covariates to adjust for in the regression models.
#'   These variables will be extracted from \code{object@meta.data[, covariate_names]}. Defaults to \code{NULL}.
#' @param adjust_setting Optional. A named logical vector of length 3 indicating whether to adjust for covariates in the estimation of mean, variance, and covariance.
#'   Must be named \code{c("mean", "var", "covar")}. Defaults to \code{c(mean = TRUE, var = TRUE, covar = TRUE)}.
#' @param IRLS_version Optional. A character string specifying the IRLS implementation to use: \code{"Rcpp"} or \code{"base_R"}.
#'   Only the \code{"Rcpp"} version supports covariate adjustment. The \code{"base_R"} version does not.
#'   When applicable, \code{"Rcpp"} offers improved memory efficiency (~10-100 times) but may be slower (~10 times),
#'   while \code{"base_R"} is faster but more memory intensive.
#'   Defaults to \code{"Rcpp"}.
#' @param IRLS_par Optional. A named list of length 3 specifying parameters for the IRLS algorithm:
#'   \describe{
#'     \item{\code{n_iter}}{Maximum number of iterations.}
#'     \item{\code{eps}}{Convergence threshold for log-ratio change \code{delta}, computed as \code{abs(log(theta / theta_prev))}.}
#'     \item{\code{verbose}}{Logical; whether to print the convergence metric (\code{delta}) at each iteration.}
#'   }
#'   Defaults to \code{list(n_iter = 10, eps = 0.05, verbose = FALSE)}.
#'
#' @return A list of three p by p matrices:
#' \describe{
#'   \item{est}{Matrix of co-expression estimates.}
#'   \item{p_value}{Matrix of p-values for testing co-expression.}
#'   \item{test_stat}{Matrix of test statistics for evaluating the significance of co-expression.}
#' }
#' @export
#'
#' @examples
#' # See a full example at:
#' # https://changsubiostats.github.io/CS-CORE/articles/CSCORE.html
#'
CSCORE <- function(object, genes, seq_depth = NULL,
                   covariate_names = NULL,
                   adjust_setting = c('mean' = T, 'var' = T, 'covar' = T),
                   IRLS_version = 'Rcpp',
                   IRLS_par = list(n_iter = 10, eps = 0.05, verbose = FALSE)){
  # Extract the UMI count matrix from the single cell object
  # with RNA as the default assay
  # such that slot `counts` corresponds to UMI counts
  count_matrix <- t(as.matrix(.get_assay_data(object = object, assay = "RNA", layer ='counts')))
  # Extract / calculate the sequencing depths
  if(is.null(seq_depth)){
    if('nCount_RNA' %in% colnames(object@meta.data)){
      seq_depth <- object$nCount_RNA
    }else{
      seq_depth <- rowSums(count_matrix)
    }
  }else{
    if(length(seq_depth) != nrow(count_matrix)) stop("The length of the sequencing depth must match the number of cells.")
  }
  # Run CS-CORE
  if(is.null(covariate_names)){
    if(IRLS_version == 'base_R'){
      CSCORE_result <- CSCORE_IRLS_base(count_matrix[,genes], seq_depth)
    }else{
      IRLS_par[['conv']] <- 'max'
      CSCORE_result <- CSCORE_IRLS_cpp(count_matrix[,genes], seq_depth,
                                       IRLS_par = IRLS_par)
    }
  }else{
    if(any(!covariate_names %in% colnames(object@meta.data))){
      stop("[ERROR] Not all covariates are included in the Seurat object. Please check `colnames(object@meta.data)`")
    }
    covariates <- object@meta.data[, covariate_names]
    cat("[INFO] Adjust for covariates:", paste(colnames(covariates), collapse = ", "), "\n")
    covariate_matrix <- stats::model.matrix(~ ., data = covariates)[,-1]
    covariate_matrix <- scale(covariate_matrix, center = T, scale = T)
    cat("[INFO] Variables in the design matrix:", paste(colnames(covariate_matrix), collapse = ", "), "\n")

    IRLS_par[['conv']] <- 'q95'
    CSCORE_result <- CSCORE_IRLS_cpp(count_matrix[,genes], seq_depth,
                                     covariates = covariate_matrix,
                                     adjust_setting = adjust_setting,
                                     IRLS_par = IRLS_par)
  }
  return(CSCORE_result)
}
