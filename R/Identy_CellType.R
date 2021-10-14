#' Identify cell type of each cluster
#' @description This function has three steps to identify cell type of each cluster.
#' (1) Calculate the power of each known marker based on AUC
#' (area under the receiver operating characteristic curve of gene expression)
#' which indicates the capability of marker i from cell type m to distinguish
#' cluster j and the other clusters. (2) Calculate the united power (UP)
#' for cell type m across each cluster j. (3) For each cluster j we determine
#' the cell type according to UP. Generally, the cluster beongs to the cell
#' type which have the highest united power or higher than the threshold of
#' the united power (for example > 0.9 power).
#' @param object seurat object
#' @param Marker1 marker gene list rownames is symbol first column is ENS ID
#' second column is cell type
#'
#' @return Cell type with the highest power in each cluster
#' @export
#'
#' @examples
Identy_CellType <- function(object, Marker1) {
  MarkerRoc1 <- Identify_CellTypes1(object, Marker1)
  ClusterCellT1 <- Identify_CellTypes2(MarkerRoc1)
  return(ClusterCellT1)
}

Identify_CellTypes1 <- function(pbmc, Marker1) {
  NumCellType <- apply(Marker1, 1, function(X1) {
    length(strsplit(as.character(X1[3]), ";")[[1]])
  })
  Marker2 <- cbind(Marker1, NumCellType)

  MarkerRoc1 <- Cal_MarkersRoc(pbmc, rownames(Marker1))
  MarkerRoc2 <- cbind(Marker2[order(rownames(Marker2)), ], MarkerRoc1[order(rownames(MarkerRoc1)), ])
  MarkerRoc2 <- MarkerRoc2[order(MarkerRoc2$CellType), ]

  return(MarkerRoc2)
}

Identify_CellTypes2 <- function(MarkerRoc1) {
  print("Identify list of cell types")
  CellType1 <- strsplit(levels(as.factor(MarkerRoc1$CellType)), ",")
  CellType2 <- c()
  for (i in 1:length(CellType1)) {
    CellType2 <- c(CellType2, CellType1[[i]])
  }
  uCellType2 <- unique(CellType2)

  ## Revise the power according to the number of cell types
  MarkerRoc2 <- cbind(MarkerRoc1[, 1:3], t(apply(MarkerRoc1, 1, function(x1) {
    x2 <- as.numeric(x1[4:length(x1)])
    x2[seq(4, length(x2), 3)] <- x2[seq(4, length(x2), 3)] / x2[1]
    return(x2)
  })))
  colnames(MarkerRoc2) <- colnames(MarkerRoc1)

  MarkerRoc5 <- c()
  for (i in 1:length(uCellType2)) {
    uCellType3 <- uCellType2[i]
    Ind1 <- c()
    for (j in 1:length(uCellType3)) {
      Ind11 <- grep(paste0("^", uCellType3[j], "$"), MarkerRoc2$CellType)
      Ind12 <- grep(paste0(",", uCellType3[j], "$"), MarkerRoc2$CellType)
      Ind13 <- grep(paste0("^", uCellType3[j], ","), MarkerRoc2$CellType)
      Ind1 <- c(Ind1, Ind11, Ind12, Ind13)
    }
    MarkerRoc3 <- MarkerRoc2[unique(Ind1), ]
    MarkerRoc4 <- Cal_JointPower2(MarkerRoc3[, 4:ncol(MarkerRoc3)])
    MarkerRoc4 <- t(as.matrix(MarkerRoc4))
    colnames(MarkerRoc4) <- gsub("_power", "", colnames(MarkerRoc4))
    MarkerRoc5 <- rbind(MarkerRoc5, MarkerRoc4)
  }
  rownames(MarkerRoc5) <- uCellType2

  MarkerRoc8 <- c()
  for (i in 1:dim(MarkerRoc5)[2]) {
    print(colnames(MarkerRoc5)[i])
    MarkerRoc6 <- t(as.matrix(MarkerRoc5[, i]))
    MarkerRoc7 <- Sort_MarkersPower(MarkerRoc6, 0)
    MarkerRoc8 <- rbind(MarkerRoc8, MarkerRoc7)
  }
  rownames(MarkerRoc8) <- colnames(MarkerRoc5)
  colnames(MarkerRoc8) <- c("Cell.Type", "Power", "Difference")

  return(MarkerRoc8)
}
