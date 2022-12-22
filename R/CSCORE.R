#' Run CS-CORE to infer the cell-type-specific co-expression network
#'
#' Run CS-CORE on a Seurat object to infer the cell-type-specific co-expression network for a specified set of genes.
#' Note that the Seurat object should have already been subsetted to cells from the same cell type,
#' in order to infer cell-type-specific co-expressions.
#'
#' @param object A Seurat single cell object
#' @param genes A vector of gene names, for which the co-expression network will be estimated
#' @param seq_depth A length n vector of sequencing depths. If "NULL", then it will be calculated as the total number of UMI counts in each cell
#'
#' @return A list of three p by p matrices:
#' \describe{
#'   \item{est}{co-expression estimates}
#'   \item{p_value}{p values}
#'   \item{test_stat}{test statistics}
#' }
#' @export
#'
#' @examples
#' # to be added
CSCORE <- function(object, genes, seq_depth = NULL){
  # Extract the UMI count matrix from the single cell object
  # with RNA as the default assay
  # such that slot `counts` corresponds to UMI counts
  count_matrix <- t(as.matrix(Seurat::GetAssayData(object, slot = 'counts', assay = 'RNA')))
  # Extract / calculate the sequencing depths
  if(is.null(seq_depth)){
    if(!is.null(object$nCount_RNA)){
      seq_depth <- object$nCount_RNA
    }else{
      seq_depth <- rowSums(count_matrix)
    }
  }else{
    if(length(seq_depth) != nrow(count_matrix)) stop("The length of the sequencing depth must match the number of cells.")
  }
  # Run CS-CORE
  CSCORE_result <- CSCORE_IRLS(count_matrix[,genes], seq_depth)
  return(CSCORE_result)
}
