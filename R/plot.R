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
#' @param ... parameter in pheatmap
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
                        , cluster_rows=F, Color1=NULL, ...){
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
                    border_color=rgb(200/255,200/255,200/255),...)

  return(Hier1)
}


#' Format Conserved markers to plot heatmap
#'
#' @param ConservedMarker Result from 'Identify_ConservedMarkers'
#'
#' @return list contains two data.frame
#' @export
#'
#' @examples load(system.file("extdata", "zf_mm_markers.rda", package = "CACIMAR"))
#' ConservedMarker <- Identify_ConservedMarkers(OrthG_Mm_Zf,Mm_marker,Zf_marker,
#' Species_name1 = 'mm',Species_name2 = 'zf')
#' MarkersPlot <- FormatConservedMarkers(ConservedMarker)
FormatConservedMarkers <- function(ConservedMarker){
  validInput(ConservedMarker,'ConservedMarker','df')
  Species_name <- colnames(ConservedMarker)[grep('Used',colnames(ConservedMarker))]
  Species_name1 <- strsplit(Species_name[1],'_')[[1]][2]
  Species_name2 <- strsplit(Species_name[2],'_')[[1]][2]
  Orth_cluster <- data.frame(ConservedMarker[,grep(paste0( Species_name1,'Allcluster'),colnames(ConservedMarker))],
                             ConservedMarker[,grep(paste0( Species_name2,'Allcluster'),colnames(ConservedMarker))])
  final_cluster1 <- c()
  for (i in 1:nrow(Orth_cluster)) {
    cluster1 <- unlist(strsplit(Orth_cluster[i,1],',')[[1]])
    cluster2 <- unlist(strsplit(Orth_cluster[i,2],',')[[1]])
    final_cluster2 <- intersect(cluster1,cluster2)[1]
    final_cluster1 <- c(final_cluster1,final_cluster2)
  }
  Orth_cluster_gene <- data.frame(final_cluster1,ConservedMarker$mmgene,ConservedMarker$zfgene)

  Sp1 <- ConservedMarker[,grep(Species_name1,colnames(ConservedMarker))]
  rownames(Sp1) <- Sp1[,grep(paste0(Species_name1,'gene'),colnames(Sp1))]
  Sp2 <- ConservedMarker[,grep(Species_name2,colnames(ConservedMarker))]
  rownames(Sp2) <- Sp2[,grep(paste0(Species_name2,'gene'),colnames(Sp2))]

  Species1_Poweridx <- (grep(paste0(Species_name1,'Pvalue'),colnames(Sp1))+1):ncol(Sp1)
  Species2_Poweridx <- (grep(paste0(Species_name2,'Pvalue'),colnames(Sp2))+1):ncol(Sp2)
  Species1_clusteridx <- grep(paste0(Species_name1,'cluster'),colnames(Sp1))
  Species2_clusteridx <- grep(paste0(Species_name2,'cluster'),colnames(Sp2))

  Species1ClusterPower <- Sp1[,c(Species1_clusteridx,Species1_Poweridx)]
  Species2ClusterPower <- Sp2[,c(Species2_clusteridx,Species2_Poweridx)]

  Species1ClusterPower[,1] <- Orth_cluster_gene[,1]
  Species2ClusterPower[,1] <- Orth_cluster_gene[,1]
  Species1ClusterPower[,1] <- paste0(Species_name1,Species1ClusterPower[,1])
  Species2ClusterPower[,1] <- paste0(Species_name2,Species2ClusterPower[,1])
  Species1ClusterPower <- orderCellType(Species1ClusterPower)
  Species2ClusterPower <- orderCellType(Species2ClusterPower)
  Species1ClusterPower[,1] <- gsub(Species_name1,'',Species1ClusterPower[,1])
  Species2ClusterPower[,1] <- gsub(Species_name2,'',Species2ClusterPower[,1])
  OrthGList <- list(Species1ClusterPower,Species2ClusterPower)
  names(OrthGList) <- c(paste0(Species_name1,'ClusterPower'),paste0(Species_name2,'ClusterPower'))
  return(OrthGList)
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
  validInput(color_number,'color_number','numeric')
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
#' @param ConservedMarker Markers table
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
#' @param legend logicalal to determine if legend should be drawn or not
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
#' Plot_MarkersHeatmap(all.markers[,c(2,6,7,8)])
Plot_MarkersHeatmap <- function(ConservedMarker,start_col = 2,module_colors = NA,
                                heatmap_colors = NA, cluster_rows = F,cluster_cols = F,
                                show_rownames = F, show_colnames = F,cellwidth = 6,
                                cellheight = 2,legend = F,annotation_legend=F,
                                annotation_names_row = F, ...){
  validInput(ConservedMarker,'ConservedMarker','df')
  validInput(start_col,'start_col','numeric')
  validInput(cellwidth,'cellwidth','numeric')
  validInput(cellheight,'cellheight','numeric')
  validInput(cluster_rows,'cluster_rows','logical')
  validInput(cluster_cols,'cluster_cols','logical')
  validInput(show_rownames,'show_rownames','logical')
  validInput(show_colnames,'show_colnames','logical')
  validInput(legend,'legend','logical')
  validInput(annotation_legend,'annotation_legend','logical')
  validInput(annotation_names_row,'annotation_names_row','logical')


  if (is.na(module_colors[1])) {
    all_cell_type <- c(ConservedMarker[,1])
    all_cell_type <- all_cell_type[!duplicated(all_cell_type)]
    module_colors <- CACIMAR::CACIMAR_cols(length(all_cell_type))
    names(module_colors) <- all_cell_type
    module_colors <- list(module_colors)
    names(module_colors) <- 'Celltype'
  }else{
    all_cell_type <- c(ConservedMarker[,1])
    all_cell_type <- all_cell_type[!duplicated(all_cell_type)]
    module_colors <- module_colors[1:length(all_cell_type)]
    names(module_colors) <- all_cell_type
    module_colors <- list(module_colors)
    names(module_colors) <- 'Celltype'
  }
  if (is.na(heatmap_colors)) {
    heatmap_colors <- viridisLite::viridis(100,option = 'D')
  }
  annotation_row <- as.data.frame(ConservedMarker[,start_col-1])
  rownames(annotation_row) <- rownames(ConservedMarker)
  colnames(annotation_row) <- 'Celltype'
  gap_row <- c()
  for (i in levels(as.factor(ConservedMarker[,1]))) {
    row1 <- which(ConservedMarker[,1]==i)
    gap_row <- c(gap_row,row1[length(row1)])
  }
  p1=pheatmap::pheatmap(ConservedMarker[,start_col:ncol(ConservedMarker)],annotation_row =
                annotation_row,cluster_rows = cluster_rows,
           show_rownames = show_rownames,cluster_cols = cluster_cols,
           gaps_row = gap_row,show_colnames = show_colnames, color = heatmap_colors,
           cellwidth = cellwidth,cellheight = cellheight,legend = legend,
           annotation_colors = module_colors,annotation_legend = annotation_legend,
           annotation_names_row = annotation_names_row,...)
  return(p1)
}


#' Plot communication between cell types
#'
#' @description  This function takes a dataframe containing information about communication between cell types,
#' and plots a graph to visualize the communication pattern. The graph represents sender and
#' receiver cell types as nodes, and the strength of communication as the width of edges connecting them.
#'
#' @param graph_df The input dataframe containing communication data, see the example in tutorial.
#' @param sendercolumn single column name, The name of the column containing sender cell types.
#' @param receivercolumn single column name, The name of the column containing receiver cell types.
#' @param widthcolumn single column name, The name of the column containing the width of edges (communication strength).
#' @param colors Optional vector of colors for nodes and edges.
#' @param useLabels Logical value indicating whether to display labels for nodes.
#' @param nodeSize Size of the nodes.
#' @param vertex.label.color Color of node labels.
#' @param vertex.label.cex Size of node labels.
#' @param ...
#'
#' @return The plotted graph as an invisible object
#' @import igraph
#' @importFrom igraph graph_from_data_frame layout.circle
#' @export
#'
#' @examples
#' plot_graph_network(long_df, sendercolumn = "celltype1", receivercolumn = "celltype2", widthcolumn = "score", useLabels = T, colors = colors, nodeSize = 10, vertex.label.color = "black", vertex.label.cex = 1.5)
#'
Plot_Celltype.Communication <- function(graph_df, sendercolumn, receivercolumn, widthcolumn, colors = NULL, useLabels = TRUE, nodeSize = 5, vertex.label.color, edge.label.color, vertex.label.cex, ...) {
  graph_df$sender <- as.character(graph_df[[sendercolumn]])
  graph_df$receiver <- as.character(graph_df[[receivercolumn]])
  graph_df$width <- graph_df[[widthcolumn]]

  vertices <- data.frame(id = unique(graph_df$sender), size = nodeSize,
                         stringsAsFactors = FALSE, label = unique(graph_df$sender))

  if (!is.null(colors)) {
    vertices$label.color <- colors[vertices$id]
    vertices$color <- colors[vertices$id]
    graph_df$color <- colors[graph_df$sender]
  }

  plot_graph <- igraph::graph_from_data_frame(graph_df, vertices = vertices)

  layout <- layout.circle(plot_graph)

  if (useLabels) {
    plot(plot_graph, layout = layout, vertex.frame.color = NA, edge.arrow.mode = "-", edge.curved = 0.3,
         edge.label = NULL,
         edge.alpha = 0.7,
         vertex.label = vertices$label,
         vertex.label.color = vertex.label.color,
         vertex.label.cex = vertex.label.cex,
         ...)
  } else {
    plot(plot_graph, layout = layout, vertex.frame.color = NA, edge.arrow.mode = "-", vertex.label = NA,
         edge.curved = 0.3,
         edge.label = NULL,
         edge.alpha = 0.7,
         vertex.label = vertices$label,
         vertex.label.color = vertex.label.color,
         vertex.label.cex = vertex.label.cex,
         ...)
  }
  invisible(plot_graph)
}



#' @title Plot Phylogenetic Tree of Different Species Cell Types
#' @description Plot a tree of different species cell types using the conserved score of cell types
#'
#' @param dist_matrix Distance matrix
#' @param hcluster.method chariter, Hierarchical clustering method, possible values are "average" (default), "single", "complete"
#' @param species.vector vector, Species vector
#' @param layout.tree Tree layout ("rectangular", "dendrogram", "fan", "circular"), default is "rectangular"
#' @param offset Tiplab offset, default is 0.01
#' @param tiplab.size Tiplab size, default is 5
#' @param tippoint.shape Tippoint shape,  default is 22
#' @param tippoint.shape.size Tippoint size, default is 4
#' @param geom_nodepoint Nodepoint size, default is 0
#' @param col.value Color vector for species, default is c(rgb(102/255,46/255,115/255), rgb(31/255,153/255,139/255))
#'
#' @return a list, which contain the Tree data and the Tree plot
#' @export
#'
#' @examples load(system.file("extdata", "zf_mm_markers.rda", package = "CACIMAR"))
#' expression <- Identify_ConservedCellTypes(OrthG_Mm_Zf,Zf_marker,Mm_marker,'zf','mm')
#' dist_matrix <- expression[[2]]
#' species.vector <- substr(rownames(dist_matrix), 1, 2)
#' Species_CellType_Tree <- Plot_Species_CellType_Tree(dist_matrix = dist_matrix, species.vector = species.vector, hcluster.method = "average", geom_nodepoint = 0, layout.tree = "rectangular", offset = 0.02, tiplab.size = 5, tippoint.shape = 22, tippoint.shape.size = 4)
#' Species_CellType_Tree$Plot
Plot_Species_CellType_Tree <- function(dist_matrix, hcluster.method = c("average", "single", "complete"),
                                       species.vector,
                                       layout.tree = c('rectangular', 'dendrogram', 'fan', 'circular'),
                                       offset = 0.01, tiplab.size = 5, tippoint.shape = c(22, 21),
                                       tippoint.shape.size = 4, geom_nodepoint = 0,
                                       col.value = c(rgb(102/255,46/255,115/255), rgb(31/255,153/255,139/255))) {
  stopifnot(is.matrix(dist_matrix))
  stopifnot(length(hcluster.method) == 1 && hcluster.method %in% c("average", "single", "complete"))
  stopifnot(is.vector(species.vector))
  stopifnot(length(layout.tree) == 1 && layout.tree %in% c('rectangular', 'dendrogram', 'fan', 'circular'))
  stopifnot(is.numeric(offset) && offset >= 0)
  stopifnot(is.numeric(tiplab.size) && tiplab.size > 0)
  stopifnot(is.numeric(tippoint.shape) && length(tippoint.shape) == 1)
  stopifnot(is.numeric(tippoint.shape.size) && tippoint.shape.size > 0)
  stopifnot(is.numeric(geom_nodepoint) && geom_nodepoint >= 0)
  stopifnot(length(col.value) >= 2 && is.character(col.value))


  # Perform hierarchical clustering
  hclust_tree <- hclust(as.dist(1 - dist_matrix), method = hcluster.method)

  # Convert the clustering result into a phylogenetic tree object
  tree <- ape::as.phylo(hclust_tree)
  tree$species <- species.vector

  # Generate the plot of the cell type tree
  tree_plot <- Plot_CellType_Tree(tree, layout.tree, offset, tiplab.size, tippoint.shape,
                                  tippoint.shape.size, geom_nodepoint = geom_nodepoint, col.value = col.value)

  # Return the tree object and the plot in a list
  return(list("Tree" = tree, "Plot" = tree_plot))
}

Plot_CellType_Tree <- function(tree, layout.tree, geom_nodepoint, col.value, offset, tiplab.size, tippoint.shape, tippoint.shape.size){
  tree.p <- ggtree::ggtree(tree, layout=layout.tree, size=0.7, col="black")
  data <- fortify(tree)
  tree.df <- data.frame("tip.label" = tree$tip.label, "species" = tree$species)
  tregraph <- tree.p %<+% tree.df +
    geom_tiplab(aes(col = species), offset = offset, size = tiplab.size) +
    geom_tippoint(aes(col = species, fill = species), shape = tippoint.shape, size = tippoint.shape.size) +
    geom_nodepoint(size = geom_nodepoint) +
    theme_tree2() +
    xlim(NA, max(data$x) * 1.3) +
    scale_fill_manual(values = col.value) +
    scale_color_manual(values = col.value) +
    guides(shape = ggplot2::guide_legend(override.aes = list(shape = c(tippoint.shape, tippoint.shape)))) +
    theme(legend.text = element_text(size = 16),  # 设置图例文本大小为12
          legend.title = element_text(size = 16))  # 设置图例标题大小为14

  return(tregraph)
}

