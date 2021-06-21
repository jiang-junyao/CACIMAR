

#' Calculate the power of each marker gene
#'
#' @param pbmc01 seurat object
#' @param marker01 marker gene list
#'
#' @return ???
#'
#'
#' @examples
Cal_MarkersRoc <- function(pbmc01, marker01) {
  clus01 <- levels(pbmc01@active.ident)
  for (i in 1:length(clus01)) {
    print(clus01[i])
    markerRoc01 <- Markertest(pbmc01, cells.1 = Seurat::WhichCells(pbmc01, idents = clus01[i]), cells.2 = Seurat::WhichCells(pbmc01, idents = clus01[i], invert = TRUE), genes.use = marker01)
    if (i == 1) {
      markerRoc02 <- markerRoc01[order(rownames(markerRoc01)), ]
    } else {
      markerRoc02 <- cbind(markerRoc02, markerRoc01[order(rownames(markerRoc01)), ])
    }
    colnames(markerRoc02)[(ncol(markerRoc02) - 2):ncol(markerRoc02)] <- paste0(paste(paste0("C", clus01[i]), collapse = ""), "_", colnames(markerRoc02)[(dim(markerRoc02)[2] - 2):dim(markerRoc02)[2]])
  }

  return(markerRoc02)
}
