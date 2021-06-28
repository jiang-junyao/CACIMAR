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
