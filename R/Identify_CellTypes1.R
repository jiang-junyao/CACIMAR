#' caculate ROC
#'
#' @param pbmc seurat object
#' @param Marker1 marker gene list
#'
#' @return ???
#'
#'
#' @examples
Identify_CellTypes1 <- function(pbmc, Marker1) {
  NumCellType <- apply(Marker1, 1, function(X1) {
    length(strsplit(as.character(X1[2]), ",")[[1]])
  })
  Marker2 <- cbind(Marker1, NumCellType)

  MarkerRoc1 <- Cal_MarkersRoc(pbmc, rownames(Marker1))
  MarkerRoc2 <- cbind(Marker2[order(rownames(Marker2)), ], MarkerRoc1[order(rownames(MarkerRoc1)), ])
  MarkerRoc2 <- MarkerRoc2[order(MarkerRoc2$CellType), ]

  return(MarkerRoc2)
}
