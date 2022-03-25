#' plot the heatmap of marker genes across different species
#' @param RNA1 correlation of expression in each cell type
#' @param CellType cell type. First column should be Species, second column should
#' be Cluster, third column should be CellType, fourth column should be Group
#' @param RowType1 character, indicating the cell types that you want to show
#' on the row in heatmap. RowType1='' means show all cell types
#' @param ColType1 character, indicating the cell types that you want to show
#' on the column in heatmap. RowType1='' means show all cell types
#' @param cluster_cols boolean values determining if columns should be clustered
#' or hclust object
#' @param cluster_rows boolean values determining if rows should be clustered or
#' hclust object
#' @param Color1 vector of colors used in heatmap
#' @export
#' @importFrom pheatmap pheatmap
#' @return return a list used in Plot_tree
#'
#' @examples
Heatmap_Cor <- function(RNA1, CellType, RowType1='', ColType1='', cluster_cols=T
                        , cluster_rows=F, Color1=NULL){
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

  white1 <- rgb(230/255,230/255,230/255); purple1 <- rgb(192/255,103/255,169/255)
  purple2 <- rgb(148/255,43/255,112/255);
  blue1 <- rgb(72/255,85/255,167/255); red1 <- rgb(239/255,58/255,37/255)
  black1 <- rgb(71/255,71/255,71/255); yellow1 <- rgb(250/255,240/255,21/255);
  if(is.null(Color1)){ Color1 <- c(blue1, 'white', red1) }

  Hier1 <- pheatmap(as.matrix(RNA2), cluster_cols =cluster_cols, cluster_rows =
                      cluster_rows, color = colorRampPalette(Color1)(50),
                    border_color=rgb(200/255,200/255,200/255))

  return(Hier1)
}

Handle_CellType<-function(RNA1,CellT1){
  CellT2 <- cbind(apply(CellT1, 1, function(x1){ x2 <- paste0(x1[1],x1[2]) })
                  , CellT1)
  colnames(RNA1) <- apply(CellT2[match(colnames(RNA1), CellT2[, 1]), ],
                          1, function(x1){ x2 <- paste(c(x1[2],'.',x1[4])
                                                       ,collapse='') })
  rownames(RNA1) <- apply(CellT2[match(rownames(RNA1), CellT2[, 1]), ],
                          1, function(x1){ x2 <- paste(c(x1[2],'.',x1[4])
                                                       ,collapse='') })
  RNA1[RNA1==1] <-NA
  return(RNA1)
}





#' Title
#'
#' @param Seurat_object Seurat object which should contain clustering inforamtion, reduction information(UMAP or TSNE)
#' @param Gene1 character, indicating the genes you want to plot
#' @param Marker_list Marker gene for each cell types
#' @param fileprefix the prefix of output figure
#' @param ident.include numeric, indicating which clusters include in the plot
#' @param width numeric, indicating width of the figure
#' @param height numeric, indicating height of the figure
#' @param PlotType1 types of plot, including: 'ggplot', 'Dotplot', 'VlnPlot' and 'FeaturePlot'
#' @importFrom ggplot2 ggplot
#' @return
#' @export
#'
#' @examples
Plot_MarkersExp<-function(Seurat_object,Gene1,Marker_list,File1,ident.include=NULL,width=12,height=10,PlotType1=c('ggPlot', 'DotPlot', 'VlnPlot', 'FeaturePlot')){
  print('Plot single-cell gene expression')
  if (is.null(ident.include)) {
    ident.include <- as.numeric(levels(Seurat_object@active.ident))
  }

  CellType1<-'VascularEndothelialCells';
  plot_MarkersExp(Seurat_object, File1=File1, Marker_list = Marker_list,Gene1=Gene1, PlotType1=PlotType1, ident.include=ident.include, width=width, height=height)
}







plot_MarkersExp <- function(pbmc1, File1, Marker_list, PlotType1=c('FeaturePlot', 'DotPlot', 'VlnPlot', 'ggPlot'), Gene1='', CellType1='', ident.include=NULL, group.by=NULL, width=12, height=10, nrow=3, RedMeth1='tSNE', no.legend=F, legend.text.size=12, legend.label=T, size.title.use=10, size.axis.text=10, point.size.use=2, Xtick1=T, dot.scale=10, IDtype1=c('EnsID', 'Symbol'), Col2=NULL){

  Col1 <- c(rgb(169/255,169/255,169/255), rgb(103/255,199/255,193/255), rgb(3/255,161/255,198/255), rgb(97/255,156/255,255/255), rgb(0/255,114/255,189/255), rgb(140/255,198/255,63/255), rgb(192/255,193/255,48/255), rgb(247/255,147/255,30/255), rgb(220/255,20/255,60/255), rgb(230/255,134/255,201/255), rgb(157/255,115/255,194/255), rgb(143/255,72/255,156/255), rgb(97/255,156/255,255/255), rgb(163/255,165/255,0/255), rgb(82/255,24/255,74/255), rgb(129/255,70/255,58/255))
  if(is.null(Col2)){ Col2 <- c(rgb(247/255,147/255,30/255), rgb(248/255,115/255,106/255), rgb(220/255,20/255,60/255), rgb(169/255,169/255,169/255), rgb(150/255,206/255,180/255), rgb(163/255,165/255,0/255), rgb(192/255,193/255,48/255), rgb(157/255,115/255,194/255), rgb(183/255,76/255,171/255), rgb(230/255,134/255,201/255), rgb(140/255,198/255,63/255), rgb(255/255,191/255,15/255), rgb(103/255,199/255,193/255), rgb(3/255,161/255,198/255), rgb(97/255,156/255,255/255), rgb(129/255,70/255,58/255), rgb(0/255,114/255,189/255), rgb(74/255,76/255,191/255), rgb(103/255,97/255,156/255), rgb(42/255,122/255,155/255)) }


  if(is.element(PlotType1[1], colnames(pbmc1@meta.data))){
    if(RedMeth1=='UMAP'){
      xlim1 <- c(min(pbmc1@reductions$umap@cell.embeddings[, 1]), max(pbmc1@reductions$umap@cell.embeddings[, 1]))
      ylim1 <- c(min(pbmc1@reductions$umap@cell.embeddings[, 2]), max(pbmc1@reductions$umap@cell.embeddings[, 2]))
    }else{
      xlim1 <- c(min(pbmc1@reductions$tsne@cell.embeddings[, 1]), max(pbmc1@reductions$tsne@cell.embeddings[, 1]))
      ylim1 <- c(min(pbmc1@reductions$tsne@cell.embeddings[, 2]), max(pbmc1@reductions$tsne@cell.embeddings[, 2]))
    }
    uGroup1 <- sort(unique(pbmc1@meta.data[, PlotType1[1]])); P1 <- list()
    for(i in 1:length(uGroup1)){ print(uGroup1[i])
      pbmc2 <- Seurat_SubsetData(pbmc1, PlotType1[1], SubS1=uGroup1[i])
      if(RedMeth1=='UMAP'){ P1[[i]] <- DimPlot(object = pbmc2, reduction.use = 'umap', no.legend=no.legend, group.by = PlotType1[1], cols.use =Col1[1:length(unique(pbmc1@meta.data$protocol))]) + labs(x='UMAP 1', y='UMAP 2') + theme(legend.text=element_text(size=legend.text.size), axis.line=element_line(size=0.5), axis.text = element_text(size=8), axis.title = element_text(size=8), panel.border = element_blank()) + coord_cartesian(xlim=xlim1, ylim=ylim1)
      }else{ P1[[i]] <- TSNEPlot(object = pbmc2, group.by = PlotType1[1], do.return = TRUE, bty = "n", pt.size = 0.4, label.size=8, no.legend=no.legend, colors.use =Col1[i%%length(Col1)+1]) + labs(x='tSNE 1', y='tSNE 2') + theme(legend.text=element_text(size=legend.text.size), axis.line=element_line(size=0.5), axis.text = element_text(size=8), axis.title = element_text(size=8), panel.border = element_blank()) + coord_cartesian(xlim=xlim1, ylim=ylim1) }
    }
    pdf(file=paste0(File1,'_',PlotType1[1],'_Each.pdf'), width=width, height=height)
    print(plot_grid(plotlist = P1, nrow = nrow, align='hv'))
    dev.off()

  }else{
    Marker1 <- Marker_list
    if(Gene1[1]=='' & CellType1==''){ Gene2 <- rownames(Marker1)
    }else if(Gene1[1]=='' & CellType1!=''){ Gene2 <- rownames(Marker1[grep(CellType1, Marker1$CellType), ])
    }else{ Gene2 <- Gene1 }

    if(PlotType1[1]=='DotPlot'){ GeneID01 <- apply(GeneID01, 2, rev) }
    GeneID1 <- Gene2; Symb1 <- Gene2

    pbmc2 <- pbmc1
    if(length(GeneID1)==1){
      if(IDtype1[1]=='EnsID'){
        pbmc2@assays$RNA@data <- t(as.matrix(pbmc2@assays$RNA@data[GeneID1, ])); rownames(pbmc2@assays$RNA@data) <- Symb1
      }else{ pbmc2@assays$RNA@data <- t(as.matrix(pbmc2@assays$RNA@data[Symb1, ])) }
    }else if(length(GeneID1)>1){
      if(IDtype1[1]=='EnsID'){
        pbmc2@assays$RNA@data <- pbmc2@assays$RNA@data[GeneID1, ]; rownames(pbmc2@assays$RNA@data) <- Symb1
      }else{ pbmc2@assays$RNA@data <- as.matrix(pbmc2@assays$RNA@data[Symb1, ]) }
    }else{ stop('Error: no such gene') }

    pbmcM1 <- rowMeans(as.matrix(pbmc2@assays$RNA@data))
    pbmc2@assays$RNA@data <- pbmc2@assays$RNA@data[pbmcM1!=0, ]
    print(dim(pbmc2@assays$RNA@data))
    if(is.null(dim(pbmc2@assays$RNA@data))){ pbmc2@assays$RNA@data <- t(as.matrix(pbmc2@assays$RNA@data)); rownames(pbmc2@assays$RNA@data) <- Symb1[pbmcM1!=0]
    }else if(nrow(pbmc2@assays$RNA@data)==0){ stop('Error: no such gene') }

    if(PlotType1[1]=='DotPlot'){

      pdf(file=paste0(File1,'_DotPlot.pdf'), width=width, height=height)
      if(is.null(group.by)){
        DotPlot(pbmc2, genes.plot = rownames(pbmc2@assays$RNA@data), cols.use = c("blue", "red"), plot.legend = TRUE, dot.scale=dot.scale)
      }else{ DotPlot(pbmc2, genes.plot = rownames(pbmc2@data), cols.use = c("blue", "red"), group.by=group.by, plot.legend = TRUE, dot.scale=dot.scale) }
      dev.off()

    }else if(PlotType1[1]=='FeaturePlot'){
      pdf(file=paste0(File1,'_FeaturePlot',CellType1,'.pdf'), width=width, height=height)
      FeaturePlot(pbmc2, features.plot = rownames(pbmc2@data), min.cutoff = "q9", cols.use = c("lightgrey", "red"), pt.size = 0.5)
      dev.off()

    }else if(PlotType1[1]=='VlnPlot'){
      print('Plot violin figure of gene expression')
      source('/users/jwang3/RetReg/Programs/VlnPlot_2.R')

      pdf(file=paste(File1,'_VlnPlot.pdf',sep=''), width=width, height=height)
      if(is.null(ident.include)){
        if(Xtick1==T){
          print(VlnPlot_2(object = pbmc2, features.plot = rownames(pbmc2@data), group.by=group.by, xlab='', size.title.use=size.title.use, size.axis.text=size.axis.text, do.return = F, nCol=1, point.size.use = point.size.use, xtick=levels(pbmc2@meta.data[, group.by]), cols.use = Col2[1:length(unique(pbmc2@meta.data[, group.by]))]))
        }else{ print(VlnPlot_2(object = pbmc2, features.plot = rownames(pbmc2@data), group.by=group.by, xlab='', size.title.use=size.title.use, size.axis.text=size.axis.text, do.return = F, nCol=1, point.size.use = point.size.use, xtick=c(), cols.use = Col2[1:length(unique(pbmc2@meta.data[, group.by]))])) }
      }else{
        if(Xtick1==T){ print(VlnPlot_2(object = pbmc2, features.plot = rownames(pbmc2@data), ident.include = ident.include, group.by=group.by, xlab='', size.title.use=size.title.use, size.axis.text=size.axis.text, do.return = F, nCol=1, point.size.use = point.size.use, xtick=levels(pbmc2@meta.data[, group.by]), cols.use = Col2[ident.include]))
        }else{ print(VlnPlot_2(object = pbmc2, features.plot = rownames(pbmc2@data), ident.include = ident.include, group.by=group.by, xlab='', size.title.use=size.title.use, size.axis.text=size.axis.text, do.return = F, nCol=1, point.size.use = point.size.use, xtick=c(), cols.use = Col2[ident.include])) }
      }
      dev.off()

    }else if(PlotType1[1]=='ggPlot'){
      print('Plot ggPlot figure of each gene expression')
      if(RedMeth1=='UMAP'){
        pbmc3 <- as.data.frame(cbind(pbmc2@reductions$umap@cell.embeddings[,1:2], t(as.matrix(pbmc2@assays$RNA@data))))
      }else{ pbmc3 <- as.data.frame(cbind(pbmc2@reductions$tsne@cell.embeddings, t(as.matrix(pbmc2@assays$RNA@data)))) }
      print(dim(pbmc3))
      for(i in 3:dim(pbmc3)[2]){ print(colnames(pbmc3)[i])
        pdf(file=paste0(File1,'_ggPlot_',colnames(pbmc3)[i],'.pdf'), width=width, height=height)
        if(legend.label==T){
          if(RedMeth1=='UMAP'){
            print(ggplot(pbmc3, aes(UMAP_1, UMAP_2)) + geom_point(aes(colour = pbmc3[,i]), size=0.6) + scale_colour_gradient(name=colnames(pbmc3)[i], low = rgb(169/255,169/255,169/255), high = rgb(220/255,20/255,60/255)) + labs(x=paste0(RedMeth1,' 1'), y=paste0(RedMeth1,' 2')) + theme(legend.text=element_text(size=legend.text.size), axis.line=element_line(size=0.75), axis.text = element_text(size=24), axis.title = element_text(size=24), panel.border = element_blank()) )
          }else{ print(ggplot(pbmc3, aes(tSNE_1, tSNE_2)) + geom_point(aes(colour = pbmc3[,i]), size=0.6) + scale_colour_gradient(name=colnames(pbmc3)[i], low = rgb(169/255,169/255,169/255), high = rgb(220/255,20/255,60/255)) + labs(x=paste0(RedMeth1,' 1'), y=paste0(RedMeth1,' 2')) + theme(legend.text=element_text(size=legend.text.size), axis.line=element_line(size=0.75), axis.text = element_text(size=24), axis.title = element_text(size=24), panel.border = element_blank()) ) }
        }else{
          if(RedMeth1=='UMAP'){
            print(ggplot(pbmc3, aes(UMAP_1, UMAP_2)) + geom_point(aes(colour = pbmc3[,i]), size=0.6) + scale_colour_gradient(name='', low = rgb(169/255,169/255,169/255), high = rgb(220/255,20/255,60/255)) + labs(x=paste0(RedMeth1,' 1'), y=paste0(RedMeth1,' 2')) + theme(legend.text=element_text(size=legend.text.size), axis.line=element_line(size=0.75), axis.text = element_text(size=24), axis.title = element_text(size=24), panel.border = element_blank()) )
          }else{ print(ggplot(pbmc3, aes(tSNE_1, tSNE_2)) + geom_point(aes(colour = pbmc3[,i]), size=0.6) + scale_colour_gradient(name='', low = rgb(169/255,169/255,169/255), high = rgb(220/255,20/255,60/255)) + labs(x=paste0(RedMeth1,' 1'), y=paste0(RedMeth1,' 2')) + theme(legend.text=element_text(size=legend.text.size), axis.line=element_line(size=0.75), axis.text = element_text(size=24), axis.title = element_text(size=24), panel.border = element_blank()) ) }
        }
        dev.off()
      }
    }
  }
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



#' Plot the Marker genes in each cluster
#'
#' @param scRNA1 dataframe, generated by
#' @param ModuleScale1 change the relative proportion of module bar and expression heatmap
#' @param Color1 vector of colors used in heatmap
#' @param ModuleColor1 colors for each moudle
#' @param Gene1 genes that you want to show in the pheatmap
#' @param NumRowBlank the blank between each group
#' @param show_colnames logic, indicating whether show column names
#' @param legend1 logic, indicating whether show the legend
#' @return
#' @export
#'
#' @examples
Plot_MarkersHeatmap<-function(scRNA1, ModuleScale1=15,
                              Color1=NULL, ModuleColor1=NULL, Gene1=NULL,
                              NumRowBlank=0,show_colnames=F, legend1=T,
                              Scale='none'){
  Scale1 <- Scale
  RevOrder1 <- -1
  scRNA12 <- scRNA1[, c(grep('Cluster\\d', colnames(scRNA1)))]
  RowGroup1 <- as.numeric(apply(scRNA1, 1, function(x1){ x2 <- gsub('Cluster','', x1[4]) }))
  FigWH1 <- c(ncol(scRNA12)/5, 10); NumRowBlank1=NumRowBlank
  white1 <- rgb(230/255,230/255,230/255); purple1 <- rgb(192/255,103/255,169/255)
  purple2 <- rgb(148/255,43/255,112/255);
  if (is.null(Color1)) {
    Color1 <- c(white1, purple1, purple2)
  }
  if (is.null(ModuleColor1)) {
    Col2 <- c(rgb(247/255,147/255,30/255), rgb(248/255,115/255,106/255),
              rgb(220/255,20/255,60/255), rgb(169/255,169/255,169/255),
              rgb(150/255,206/255,180/255), rgb(163/255,165/255,0/255),
              rgb(192/255,193/255,48/255), rgb(157/255,115/255,194/255),
              rgb(183/255,76/255,171/255), rgb(230/255,134/255,201/255),
              rgb(140/255,198/255,63/255), rgb(255/255,191/255,15/255),
              rgb(103/255,199/255,193/255), rgb(3/255,161/255,198/255),
              rgb(97/255,156/255,255/255), rgb(129/255,70/255,58/255),
              rgb(0/255,114/255,189/255), rgb(74/255,76/255,191/255),
              rgb(103/255,97/255,156/255), rgb(42/255,122/255,155/255))
  }else{Col2=ModuleColor1}
  scRNA2 <- scRNA12[, 1:ncol(scRNA12)]
  scRNA21<-apply(scRNA2, 2, as.numeric)
  rownames(scRNA21)<-rownames(scRNA2)
  scRNA3 <- Kmeans_heatmap(scRNA21, K1=1, RowGroup1=RowGroup1, Scale1=Scale1,
                           Gene1=Gene1, clustering_distance_rows='correlation',
                           NumRowBlank1=NumRowBlank1, show_colnames=show_colnames,
                           legend1=legend1, cluster_cols=F, Show.Module=T,
                           Color1=Color1, ModuleColor1=Col2, ModuleScale1=ModuleScale1,
                           RevOrder1=RevOrder1)
}



Kmeans_heatmap <- function(RNA1, K1=1, Gene1=NULL, RowGroup1=NULL, ColumnGroup1=NULL, Scale1='none', Range1=c(-Inf,Inf), Reorder1=T, RevOrder1=-1, clustering_distance_rows='euclidean', clustering_method='complete', Show.Module=F, Color1=NULL, show_colnames=T, NumRowBlank1=10, NumColumnBlank1=10, NAcolum1=NULL, ModuleColor1=NULL, ModuleScale1=10, cluster_cols=F, fontsize=10, legend1=T, border_color='white', ByPanel=F){
  library(pheatmap)

  if(Scale1=='row'){ RNA1 <- t(scale(t(RNA1))) }

  if(!is.null(ColumnGroup1)){
    print('sepate columns according to varible ColumnGroup1')
    ColumnBlank1 <- array(0, dim=c(nrow(RNA1), NumColumnBlank1))
    uColumnGroup1 <- unique(ColumnGroup1)
    if(length(uColumnGroup1)==1){
    }else{
      for(i in 1:length(uColumnGroup1)){
        RNA12 <- RNA1[, ColumnGroup1==uColumnGroup1[i]]
        if(i==1){ RNA21 <- cbind(RNA12, ColumnBlank1)
        }else if(i<length(uColumnGroup1)){
          RNA21 <- cbind(RNA21, RNA12, ColumnBlank1)
        }else{ RNA21 <- cbind(RNA21, RNA12) }
      }
    }
  }else{ RNA21 <- RNA1 }

  print('Perform k-means')
  if(is.null(RowGroup1)){
    if(!is.null(NAcolum1)){
      cRNA1 <- kmeans(RNA1[, -NAcolum1], K1)
      RNA02 <- cbind(cRNA1$cluster, RNA1[, -NAcolum1])
      Cluster1 <- cRNA1$cluster
    }else{
      if(K1==1){ Cluster1 <- rep(1, nrow(RNA1))
      }else{ cRNA1 <- kmeans(RNA1, K1); Cluster1 <- cRNA1$cluster }
      RNA02 <- cbind(Cluster1, RNA1)
    }
    RNA2 <- cbind(Cluster1, RNA21); NameInd1 <- 1:K1
    colnames(RNA2)[1] <- c('KmeansGroup'); print(table(RNA2[, 'KmeansGroup']))
  }else{
    if(is.factor(RowGroup1)){
      Name1 <- levels(RowGroup1); Name2 <- sort(levels(RowGroup1));
      NameInd1 <- match(Name1, Name2)
      RowGroup1 <- as.numeric(RowGroup1); uRowGroup1 <- NameInd1;
    }else{ uRowGroup1 <- sort(unique(RowGroup1)) }

    if(!is.null(NAcolum1)){ RNA02 <- cbind(RowGroup1, RNA1[, -NAcolum1])
    }else{ RNA02 <- cbind(RowGroup1, RNA1) }
    RNA2 <- cbind(RowGroup1, RNA21);
    K1 <- length(uRowGroup1)
    colnames(RNA2)[1] <- c('KmeansGroup')
    tRNA2 <- table(RNA2[, 'KmeansGroup']);
    if(is.factor(RowGroup1)){ names(tRNA2) <- Name2;
    tRNA2 <- tRNA2[NameInd1];
    }
    print(tRNA2)
  }

  print('Sort genes'); RNA22 <- c()
  RowBlank1 <- array(0, dim=c(NumRowBlank1, ncol(RNA2)))
  colnames(RowBlank1) <- colnames(RNA2)
  for(i in 1:K1){
    if(is.null(RowGroup1)){
      RNA20 <- RNA2[RNA2[, 'KmeansGroup']==i, ]
      RNA03 <- RNA02[RNA02[, 1]==i, ]
    }else{
      RNA20 <- RNA2[RNA2[, 'KmeansGroup']==uRowGroup1[i], ]
      RNA03 <- RNA02[RNA02[, 1]==uRowGroup1[i], ]
    }
    if(Reorder1==T){
      if(nrow(RNA03)>1){
        Hier1 <- hclust(as.dist((1 - cor(t(RNA03[, 2:ncol(RNA03)])))/2))
        if(RevOrder1[1]!=-1){ RevOrder2 <- F;
        for(j in 1:length(RevOrder1)){
          if(RevOrder1[j]==i){ RevOrder2 <- T; break }
        }
        if(RevOrder2==T){ Ind1 <- rev(Hier1$order)
        }else{ Ind1 <- Hier1$order }
        }else{ Ind1 <- Hier1$order }
      }else{ Ind1 <- 1:nrow(RNA20) }
    }else{ Ind1 <- 1:nrow(RNA20) }

    if(i!=K1){ RNA22 <- rbind(RNA22, RNA20[Ind1, ], RowBlank1)
    }else{ RNA22 <- rbind(RNA22, RNA20[Ind1, ]) }
  }

  print('Revise outlier')
  RNA23 <- RNA22[, 2:ncol(RNA22)];
  print(paste('Number of outlier:', c(length(RNA23[RNA23<Range1[1]]), length(RNA23[RNA23>Range1[2]]))))
  RNA23[RNA23<Range1[1]] <- Range1[1]; RNA23[RNA23>Range1[2]] <- Range1[2]
  RNA3 <- cbind(RNA22[,1], RNA23)
  colnames(RNA3)[1] <- 'Module'
  RNA4 <- RNA22[RNA22[,1]!=0, ]

  if(is.null(Gene1)){ show_rownames <- F; Gene2 <- ''
  }else if(Gene1[1]=='all'){ show_rownames <- T; Gene2 <- rownames(RNA3)
  }else{ show_rownames <- T; Gene2 <- rownames(RNA3)
  Gene2[!is.element(Gene2, Gene1)] <- '';
  Gene3 <- rownames(RNA3)[is.element(rownames(RNA3), Gene1)]; print(sort(Gene3))
  }

  if(is.null(Color1)){ Color1 <- c(rgb(72/255,85/255,167/255), rgb(255/255,255/255,255/255), rgb(239/255,58/255,37/255)) }

  if(Show.Module==T){
    library(gridExtra)

    if(is.null(ModuleColor1)){
      ModuleColor1 <- c(rgb(255/255,255/255,255/255), rgb(152/255,152/255,152/255), rgb(72/255,85/255,167/255))
    }else{
      if(!is.null(RowGroup1) & is.factor(RowGroup1)){ ModuleColor1 <- ModuleColor1[c(1,NameInd1+1)]; }
    }

    if(K1==1){
      ph1 <- pheatmap(RNA3[, 1], clustering_distance_rows=clustering_distance_rows, show_rownames=F, show_colnames=show_colnames, border_color='white', cluster_rows=F, cluster_cols=cluster_cols, color = ModuleColor1, legend=F, silent=T, breaks=seq(1, nrow(RNA3)))
    }else{
      ph1 <- pheatmap(RNA3[, 1], clustering_distance_rows=clustering_distance_rows, show_rownames=F, show_colnames=show_colnames, border_color='white', cluster_rows=F, cluster_cols=cluster_cols, color = ModuleColor1, legend=F, silent=T)
    }
    if(ByPanel==T){
      uColumnGroup1 <- unique(ColumnGroup1)
      RNA31 <- RNA3[, 2:(length(ColumnGroup1[ColumnGroup1==uColumnGroup1[1]])+1)]
      ph12 <- pheatmap(RNA3[, 2:(length(ColumnGroup1[ColumnGroup1==uColumnGroup1[1]])+1)], clustering_distance_rows=clustering_distance_rows, clustering_method=clustering_method, show_rownames=show_rownames, show_colnames=show_colnames, labels_row=Gene2, border_color=border_color, cluster_rows=F, cluster_cols=cluster_cols, color = colorRampPalette(Color1)(100), fontsize=fontsize, legend=legend1, silent=T)
      ph2 <- pheatmap(RNA31, clustering_distance_rows=clustering_distance_rows, clustering_method=clustering_method, show_rownames=show_rownames, show_colnames=show_colnames, labels_row=Gene2, border_color=border_color, cluster_rows=F, cluster_cols=cluster_cols, color = colorRampPalette(Color1)(100), fontsize=fontsize, legend=legend1, silent=T)
      plot_list <- list(); plot_list[[1]] <- ph1[[4]]; plot_list[[2]] <- ph2[[4]]
      ModuleScale2 <- rep(2, ceiling(ncol(RNA31)/ncol(RNA3)*ModuleScale1))
      if(length(uColumnGroup1)>1){
        for(i in 2:length(uColumnGroup1)){
          GroupLeng1 <- length(ColumnGroup1[ColumnGroup1<uColumnGroup1[i]])
          GroupLeng2 <- length(ColumnGroup1[ColumnGroup1==uColumnGroup1[i]])
          RNA32 <- RNA3[, (GroupLeng1+NumColumnBlank1*(i-1)+2):(GroupLeng1+NumColumnBlank1*(i-1)+GroupLeng2+1)]
          ph3 <- pheatmap(RNA32, clustering_distance_rows=clustering_distance_rows, clustering_method=clustering_method, show_rownames=show_rownames, show_colnames=show_colnames, labels_row=Gene2, border_color=border_color, cluster_rows=F, cluster_cols=cluster_cols, color = colorRampPalette(Color1)(100), fontsize=fontsize, legend=legend1, silent=T)
          plot_list[[i+1]] <- ph3[[4]]
          ModuleScale2 <- c(ModuleScale2, rep(i+1, ceiling(ncol(RNA32)/ncol(RNA3)*ModuleScale1)))
        } }
      layout_ncol <- length(uColumnGroup1)+1; layout_matrix <- matrix(c(1, ModuleScale2), nrow=1)
    }else{ ph2 <- pheatmap(RNA3[, 2:ncol(RNA3)], clustering_distance_rows=clustering_distance_rows, clustering_method=clustering_method, show_rownames=show_rownames, show_colnames=show_colnames, labels_row=Gene2, border_color=border_color, cluster_rows=F, cluster_cols=cluster_cols, color = colorRampPalette(Color1)(100), fontsize=fontsize, legend=legend1, silent=T)
    plot_list <- list(ph1[[4]], ph2[[4]]); layout_ncol <- 2
    layout_matrix <- matrix(c(1, rep(2,ModuleScale1)), nrow=1)
    }

    g <- gridExtra::grid.arrange(arrangeGrob(grobs= plot_list, ncol=layout_ncol, layout_matrix=layout_matrix))
  }else{
    pheatmap::pheatmap(RNA3[, 2:ncol(RNA3)], clustering_distance_rows=clustering_distance_rows, clustering_method=clustering_method, show_rownames=show_rownames, show_colnames=show_colnames, labels_row=Gene2, border_color=border_color, cluster_rows=F, cluster_cols=cluster_cols, color = colorRampPalette(Color1)(100), fontsize=fontsize, legend=legend1)
  }

  return(RNA4)
}


#' Title
#'
#' @param Tree1 list,
#' @param layout1
#' @param TipFont1
#' @param PdfFile1
#' @importFrom ggtree groupOTU
#' @importFrom ggtree ggtree
#' @importFrom ape as.phylo
#' @export
#'
#' @examples
Plot_tree <- function(Tree1, layout1='circular', TipFont1=5, PdfFile1=NULL){

  Tree1 <- ape::as.phylo(Tree1$tree_col)
  Col1 <- c(rgb(242/255,101/255,33/255), rgb(0,114/255,189/255), rgb(90/255,90/255,90/255))
  GroupInfo <- split(Tree1$tip.label, gsub("\\..*", "", Tree1$tip.label))
  Tree1 <- ggtree::groupOTU(Tree1, GroupInfo)
  print(Tree1$tip.label)
  Tree1$tip.label <- apply(as.matrix(Tree1$tip.label), 1, function(x1){
    x2 <- strsplit(x1, '\\.')[[1]][2]
  } )
  print(Tree1$tip.label)

  if(layout1=='circular'){
    Ggtree1 <- ggtree::ggtree(Tree1, ggtree::aes(color=group), layout='circular') + ggtree::scale_colour_manual(values=Col1) + ggtree::geom_tiplab(size=TipFont1, ggtree::aes(angle=angle))
  }else{ Ggtree1 <- ggtree::ggtree(Tree1) + ggtree::geom_tiplab(size=TipFont1, hjust=1) }

  if(!is.null(PdfFile1)){
    pdf(file=paste0(PdfFile1,'.pdf'), width=12, height=12)
    print(Ggtree1)
    dev.off()
  }else{
    Ggtree1
  }
  return(Ggtree1)
}
