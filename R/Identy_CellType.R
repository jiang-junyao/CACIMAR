#' Title
#'
#' @param object seurat object
#' @param Marker1 marker gene list rownames is symbol first column is ENS ID second column is cell type
#'
#' @return a lsit with cell type
#' @export
#'
#' @examples
Identy_CellType <- function(object, Marker1) {
  MarkerRoc1 <- Identify_CellTypes1(object, Marker1)
  ClusterCellT1 <- Identify_CellTypes2(MarkerRoc1)
  return(ClusterCellT1)
}
