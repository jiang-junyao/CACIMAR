

#' Calculate the relative power of cell types for each cell cluster
#'
#' @param MarkerRoc1 MarkerRoc1
#'
#' @return ???
#'
#'
#' @examples
Identify_CellTypes2 <- function(MarkerRoc1) {
  print("Identify list of cell types")
  CellType1 <- strsplit(levels(as.factor(MarkerRoc1$CellType)), ",")
  CellType2 <- c()
  for (i in 1:length(CellType1)) {
    CellType2 <- c(CellType2, CellType1[[i]])
  }
  uCellType2 <- unique(CellType2)

  ## Revise the power according to the number of cell types
  MarkerRoc2 <- cbind(MarkerRoc1[, 1:2], t(apply(MarkerRoc1, 1, function(x1) {
    x2 <- as.numeric(x1[3:length(x1)])
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
