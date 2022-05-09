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
#' @param seurat_object seurat object
#' @param Marker_gene_table data.frame, indicating marker gene and its
#' corresponding cell type. Marker_gene_table should contain two columns: 'CellType'
#' represent correseponding cell types of each marker and 'Marker' represent Markers
#'
#' @return Cell type with the highest power in each cluster
#' @export
#'
#' @examples \dontrun{ annotation <- Identify_CellType(seurat_object,Marker_gene_table)
#' }
Identify_CellType <- function(seurat_object, Marker_gene_table) {
  validInput(seurat_object,'seurat_object','seuratobject')
  if (!'CellType' %in% colnames(Marker_gene_table)) {
    stop('Marker_gene_table should contain CellType column ')
  }
  if (!'Marker' %in% colnames(Marker_gene_table)) {
    stop('Marker_gene_table should contain Marker column ')
  }
  MarkerRoc1 <- Identify_CellTypes1(seurat_object, Marker_gene_table)
  ClusterCellT1 <- Identify_CellTypes2(MarkerRoc1)

  return(ClusterCellT1)
}


Identify_CellTypes1 <- function(object, Marker1) {
  ## calculate the power of each marker on distinguishing cell types
  CelltypeIdx <- grep('CellType',colnames(Marker1))
  NumCellType <- apply(Marker1, 1, function(X1) {
    length(strsplit(as.character(X1[CelltypeIdx]),",")[[1]])
  })
  Marker2 <- cbind(Marker1, NumCellType)

  MarkerRoc1 <- Cal_MarkersRoc(object, Marker1$Marker)
  MarkerRoc2 <- cbind(Marker2[order(Marker2$Marker), ], MarkerRoc1[order(rownames(MarkerRoc1)), ])
  MarkerRoc2 <- MarkerRoc2[order(MarkerRoc2$CellType), ]

  return(MarkerRoc2)
}


Identify_CellTypes2 <- function(MarkerRoc1) {
  ## Get the list of cell types
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
  #write.table(MarkerRoc2, 'MarkerRoc2.txt')

  ## Calculate the joint power for each cluster
  MarkerRoc5 <- c()
  cell_type_out <- c()
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
    startidx <- grep('NumCellType',colnames(MarkerRoc3))
    ### add the power of gene in the same cell types, if gene difference <0 ,
    ### set power of that gene negative
    if (nrow(MarkerRoc3)>0) {
      MarkerRoc4 <- Cal_JointPower2(MarkerRoc3[, (startidx+1):ncol(MarkerRoc3)])
      MarkerRoc4 <- t(as.matrix(MarkerRoc4))
      colnames(MarkerRoc4) <- gsub("_power", "", colnames(MarkerRoc4))
      MarkerRoc5 <- rbind(MarkerRoc5, MarkerRoc4)
      cell_type_out <- c(cell_type_out,uCellType3)
    }

  }
  rownames(MarkerRoc5) <- cell_type_out
  exclusive_cell_type <- paste(setdiff(uCellType2,cell_type_out),collapse = ' ,')
  if (!is.null(exclusive_cell_type)) {
    warning(paste0('No genes express in',exclusive_cell_type))
  }


  ## Sort assigned cell types according to the joint power for each cluster
  MarkerRoc8 <- c()
  for (i in 1:dim(MarkerRoc5)[2]) {
    print(colnames(MarkerRoc5)[i])
    MarkerRoc6 <- t(as.matrix(MarkerRoc5[, i]))
    MarkerRoc7 <- Sort_MarkersPower(MarkerRoc6, 0)
    MarkerRoc8 <- rbind(MarkerRoc8, MarkerRoc7)
  }
  rownames(MarkerRoc8) <- colnames(MarkerRoc5)
  colnames(MarkerRoc8) <- c("Cell.types", "Corresponding.powers", 'Predicted.cell.type')

  return(MarkerRoc8)
}

