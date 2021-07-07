


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






#' Inner function to plot MarkersExp

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







#' Inner function to subset seurat object

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
