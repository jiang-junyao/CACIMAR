#' plot the heatmap of marker genes across different species
#' @param RNA1 correlation of expression in each cell type
#' @param RowType1 character, indicating the cell types that you want to show
#' on the row in heatmap. RowType1='' means show all cell types
#' @param ColType1 character, indicating the cell types that you want to show
#' on the column in heatmap. RowType1='' means show all cell types
#' @param cluster_cols boolean values determining if columns should be clustered
#' or hclust object
#' @param cluster_rows boolean values determining if rows should be clustered or
#'
#' hclust object
#' @param Color1 vector of colors used in heatmap
#' @export
#' @importFrom pheatmap pheatmap
#' @importFrom grDevices rgb
#' @importFrom grDevices colorRampPalette
#' @importFrom methods is
#' @return pheatmap object
#'
#' @examples load(system.file("extdata", "network_example.rda", package = "CACIMAR"))
#' n1 <- Identify_ConservedNetworks(OrthG_Mm_Zf,mmNetwork,zfNetwork,'mm','zf')
#' Heatmap_Cor(n1[[2]],cluster_cols=TRUE, cluster_rows=FALSE)
Heatmap_Cor <- function(RNA1, RowType1='', ColType1='', cluster_cols=T
                        , cluster_rows=F, Color1=NULL){
  validInput(RowType1,'RowType1','character')
  validInput(ColType1,'ColType1','character')
  validInput(cluster_cols,'cluster_cols','logical')
  validInput(cluster_rows,'cluster_rows','logical')
  Ind21 <- c(); Ind22 <- c();
  if(RowType1==''){ Ind21 <- 1:dim(RNA1)[1];
  }else{
    for(i in 1:length(RowType1)){
      Ind1 <- grep(RowType1[i], rownames(RNA1))
      Ind21 <- c(Ind21, Ind1)
    }
  }

  if(ColType1==''){ Ind22 <- 1:dim(RNA1)[2]
  }else{
    for(i in 1:length(ColType1)){
      Ind1 <- grep(ColType1[i], colnames(RNA1))
      Ind22 <- c(Ind22, Ind1)
    }
  }
  RNA2 <- RNA1[Ind21, Ind22]
  #RNA2[RNA2==1] <- NA; RNA2[is.na(RNA2)] <- max(RNA2[!is.na(RNA2)])

  white1 <- rgb(230/255,230/255,230/255); purple1 <- rgb(192/255,103/255,169/255)
  purple2 <- rgb(148/255,43/255,112/255);
  blue1 <- rgb(72/255,85/255,167/255); red1 <- rgb(239/255,58/255,37/255)
  black1 <- rgb(71/255,71/255,71/255); yellow1 <- rgb(250/255,240/255,21/255);
  if(is.null(Color1)){ Color1 <- c(blue1, 'white', red1) }

  Hier1 <- pheatmap::pheatmap(as.matrix(RNA2), cluster_cols =cluster_cols, cluster_rows =
                      cluster_rows, color = colorRampPalette(Color1)(50),
                    border_color=rgb(200/255,200/255,200/255))

  return(Hier1)
}



Seurat_SubsetData <- function(pbmc1, SubG1, SubS1=NULL, ExSubS1=NULL){

  if(!is.null(SubS1)){ CellN1 <- c()
  for(SubS2 in SubS1){
    CellN1 <- c(CellN1, rownames(pbmc1@meta.data[pbmc1@meta.data[, SubG1]==SubS2, ]))
  }
  pbmc1 <- subset(pbmc1, cells=unique(CellN1) )
  }

  if(!is.null(ExSubS1)){ CellN2 <- c()
  for(ExSubS2 in ExSubS1){
    CellN2 <- c(CellN2, rownames(pbmc1@meta.data[pbmc1@meta.data[, SubG1]==ExSubS1, ]))
  }
  CellN1 <- setdiff(rownames(pbmc1@meta.data), CellN2)
  pbmc1 <- subset(pbmc1, cells=unique(CellN1) )
  }

  return(pbmc1)
}

#' CACIMAR colors palette
#'
#' @param color_number numeric, indicating used colors number
#'
#' @return vector of colors
#' @export
#'
#' @examples CACIMAR_cols(10)
#' CACIMAR_cols(20)
CACIMAR_cols <- function(color_number){
  cols <- c("OrangeRed","SlateBlue3","DarkOrange","GreenYellow","Purple",
            "DarkSlateGray","Gold","DarkGreen","DeepPink2","Red4","#4682B4",
            "#FFDAB9","#708090","#836FFF","#CDC673","#CD9B1D","#FF6EB4","#CDB5CD"
            ,"#008B8B","#43CD80","#483D8B","#66CD00","#CDC673","#CDAD00","#CD9B9B"
            ,"#FF8247","#8B7355","#8B3A62","#68228B","#CDB7B5","#CD853F","#6B8E23"
            ,"#696969","#7B68EE","#9F79EE","#B0C4DE","#7A378B","#66CDAA","#EEE8AA"
            ,"#00FF00","#EEA2AD","#A0522D","#000080","#E9967A","#00CDCD","#8B4500"
            ,"#DDA0DD","#EE9572","#EEE9E9","#8B1A1A","#8B8378","#EE9A49","#EECFA1"
            ,"#8B4726","#8B8878","#EEB4B4","#C1CDCD","#8B7500","#0000FF","#EEEED1"
            ,"#4F94CD","#6E8B3D","#B0E2FF","#76EE00","#A2B5CD","#548B54","#BBFFFF"
            ,"#B4EEB4","#00C5CD","#008B8B","#7FFFD4","#8EE5EE","#43CD80","#68838B"
            ,"#00FF00","#B9D3EE","#9ACD32","#00688B","#FFEC8B","#1C86EE","#CDCD00"
            ,"#473C8B","#FFB90F","#EED5D2","#CD5555","#CDC9A5","#FFE7BA","#FFDAB9"
            ,"#CD661D","#CDC5BF","#FF8C69","#8A2BE2","#CD8500","#B03060","#FF6347"
            ,"#FF7F50","#CD0000","#F4A460","#FFB5C5","#DAA52")
  cols_return <- cols[1:color_number]
  return(cols_return)
}


#' Plot Markers in each cell type
#' @description This function integrate R package pheatmap to plot markers in each
#' cell type
#' @param ConservedMarker
#' @param start_col numeric, indicating the start column of marker power in each
#' cell type
#' @param module_colors vector, indicating colors of modules (annotation_colors)
#' @param heatmap_colors vector, indicating colors used in heatmap
#' @param cluster_rows boolean values determining if rows should be clustered or
#' hclust object
#' @param cluster_cols boolean values determining if columns should be clustered
#' or hclust object
#' @param show_rownames boolean specifying if column names are be shown
#' @param show_colnames boolean specifying if column names are be shown
#' @param cellwidth individual cell width in points. If left as NA, then the
#' values depend on the size of plotting window
#' @param cellheight individual cell height in points. If left as NA, then the
#' values depend on the size of plotting window
#' @param legend logical to determine if legend should be drawn or not
#' @param annotation_legend boolean value showing if the legend for annotation
#' tracks should be drawn
#' @param annotation_names_row boolean value showing if the names for row
#' annotation tracks should be drawn
#' @param ... parameter in pheatmap
#' @importFrom pheatmap pheatmap
#' @importFrom viridisLite viridis
#' @return pheatmap object
#' @export
#'
#' @examples data("pbmc_small")
#' all.markers <- Identify_Markers(pbmc_small)
#' all.markers <- Format_Markers_Frac(all.markers)
#' Plot_MarkersHeatmap(all.markers)
Plot_MarkersHeatmap <- function(ConservedMarker,start_col = 7,module_colors = NA,
                                heatmap_colors = NA, cluster_rows = F,cluster_cols = F,
                                show_rownames = F, show_colnames = F,cellwidth = NA,
                                cellheight = NA,legend = F,annotation_legend=F,
                                annotation_names_row = F, ...){
  if (is.na(module_colors)) {
    all_cell_type <- c(ConservedMarker$cluster)
    all_cell_type <- all_cell_type[!duplicated(all_cell_type)]
    module_colors <- CACIMAR::CACIMAR_cols(length(all_cell_type))
    names(module_colors) <- all_cell_type
    module_colors <- list(module_colors)
    names(module_colors) <- 'Celltype'
  }
  if (is.na(heatmap_colors)) {
    heatmap_colors <- viridisLite::viridis(100,option = 'D')
  }
  annotation_row <- as.data.frame(ConservedMarker$cluster)
  rownames(annotation_row) <- rownames(ConservedMarker)
  colnames(annotation_row) <- 'Celltype'
  gap_row <- c()
  for (i in levels(as.factor(ConservedMarker$cluster))) {
    row1 <- grep(i,ConservedMarker$cluster)
    gap_row <- c(gap_row,row1[length(row1)])
  }
  p1=pheatmap::pheatmap(ConservedMarker[,start_col:(ncol(ConservedMarker)-1)],annotation_row =
                annotation_row,cluster_rows = cluster_rows,
           show_rownames = show_rownames,cluster_cols = cluster_cols,
           gaps_row = gap_row,show_colnames = show_colnames, color = heatmap_colors,
           cellwidth = cellwidth,cellheight = cellheight,legend = legend,
           annotation_colors = module_colors,annotation_legend = annotation_legend,
           annotation_names_row = annotation_names_row,...)
  return(p1)
}



