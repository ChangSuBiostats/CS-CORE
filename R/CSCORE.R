#' Run CS-CORE on a Seurat object
#'
#' @param object A Seurat single cell object
#' @param genes A vector of gene names, for which the co-expression network will be estimated
#' @param seq_depth A length n vector of sequencing depths. If "NULL", then it will be calculated as the total number of UMI counts in each cell
#'
#' @return A list of three p by p matrices: co-expression estimates, p values and test statistics
#' @export
#'
#' @examples
#' # to be added
CSCORE <- function(object, genes, seq_depth = NULL){
  # Extract the UMI count matrix from the single cell object
  count_matrix <- t(as.matrix(Seurat::GetAssayData(object, slot = 'counts')))
  # Extract / calculate the sequencing depths
  if(is.null(seq_depth)){
    seq_depth <- object$nCount_RNA
  }
  # Run CS-CORE
  cscore_network <- CSCORE_IRLS(count_matrix[,genes], seq_depth)
  return(cscore_network)
}
