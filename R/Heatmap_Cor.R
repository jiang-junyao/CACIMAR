#' plot the pheatmap of marker genes across different species
#'
#'
#' @param RNA1 expression
#' @param CellType cell type
#' @param RowType1 ??
#' @param ColType1 ??
#' @param cluster_cols boolean values determining if columns should be clustered or hclust object
#' @param cluster_rows boolean values determining if rows should be clustered or hclust object
#' @param Color1 vector of colors used in heatmap
#' @export
#' @importFrom pheatmap pheatmap
#' @return return a list used in Plot_tree
#'
#' @examples
Heatmap_Cor <- function(RNA1, CellType, RowType1='', ColType1='', cluster_cols=F, cluster_rows=F, Color1=NULL){
  RNA1 <- Handle_CellType(RNA1, CellType)
  Ind21 <- c(); Ind22 <- c();
  if(RowType1==''){ Ind21 <- 1:dim(RNA1)[1];
  }else{
    for(i in 1:length(RowType1)){
      Ind1 <- grep(RowType1[i], rownames(RNA2))
      Ind21 <- c(Ind21, Ind1)
    }
  }

  if(ColType1==''){ Ind22 <- 1:dim(RNA1)[2]
  }else{
    for(i in 1:length(ColType1)){
      Ind1 <- grep(ColType1[i], colnames(RNA2))
      Ind22 <- c(Ind22, Ind1)
    }
  }
  RNA2 <- RNA1[Ind21, Ind22]
  #RNA2[RNA2==1] <- NA; RNA2[is.na(RNA2)] <- max(RNA2[!is.na(RNA2)])

  white1 <- rgb(230/255,230/255,230/255); purple1 <- rgb(192/255,103/255,169/255); purple2 <- rgb(148/255,43/255,112/255);
  blue1 <- rgb(72/255,85/255,167/255); red1 <- rgb(239/255,58/255,37/255)
  black1 <- rgb(71/255,71/255,71/255); yellow1 <- rgb(250/255,240/255,21/255);
  if(is.null(Color1)){ Color1 <- c(blue1, 'white', red1) }

  Hier1 <- pheatmap(as.matrix(RNA2), cluster_cols =cluster_cols, cluster_rows =cluster_rows, color = colorRampPalette(Color1)(50), border_color=rgb(200/255,200/255,200/255))

  return(Hier1)
}

Handle_CellType<-function(RNA1,CellT1){
  CellT2 <- cbind(apply(CellT1, 1, function(x1){ x2 <- paste0(x1[1],x1[2]) }), CellT1)
  colnames(RNA1) <- apply(CellT2[match(colnames(RNA1), CellT2[, 1]), ], 1, function(x1){ x2 <- paste(c(x1[2],'.',x1[4]),collapse='') })
  rownames(RNA1) <- apply(CellT2[match(rownames(RNA1), CellT2[, 1]), ], 1, function(x1){ x2 <- paste(c(x1[2],'.',x1[4]),collapse='') })
  RNA1[RNA1==1] <- NA
  return(RNA1)
}
