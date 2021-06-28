#' internal function from seurat2.3.4 to test marker gene
#'
#' @param object seurat object
#' @param cells.1 cell in cluster 1
#' @param cells.2 cell in cluster 2
#' @param genes.use gene to be tested
#' @param print.bar use which apply function
#' @param assay.type ???
#' @importFrom Seurat GetAssayData
#' @return Returns a 'predictive power' (abs(AUC-0.5) * 2) ranked matrix of
# putative differentially expressed genes.
#'
#'
#' @examples
Markertest <- function(object, cells.1, cells.2, genes.use = NULL, print.bar = TRUE,
                       assay.type = "RNA") {
  data.test <- Seurat::GetAssayData(
    object = object,
    slot = "data"
  )
  if (is.null(genes.use)) {
    genes.use <- rownames(x = data.test)
  }
  to.return <- AUCMarkerTest(
    data1 = data.test[, cells.1],
    data2 = data.test[, cells.2], mygenes = genes.use, print.bar = print.bar
  )
  to.return$power <- abs(x = to.return$myAUC - 0.5) * 2
  return(to.return)
}
