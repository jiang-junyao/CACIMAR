Seurat_Markers <- function(pbmc1, Spec1=NULL, Method1='wilcox', PadjThr1=0.05, VarG1=0){
  if(Method1 == 'wilcox' | Method1 == 'roc'){
      ObjectM1 <- FindAllMarkers(object = pbmc1, test.use = Method1)
      if(Method1 == 'wilcox'){ ObjectM1 <- ObjectM1[ObjectM1$p_val_adj < PadjThr1, ] }
  }else if(Method1 == 'variable'){
      pbmc1 <- FindVariableGenes(object = pbmc1, do.plot = FALSE)
    if(VarG1==0){ ObjectM1 <- pbmc1@hvg.info[pbmc1@var.genes, ]
    }else{ ObjectM1 <- pbmc1@hvg.info[1:min(nrow(pbmc1@hvg.info), VarG1), ] }
      ObjectM1[, 'gene'] <- rownames(ObjectM1)
  }else{ stop('The method is not wilcox, roc, or variable') }
  if(Method1 == 'wilcox' | Method1 == 'roc'){
    ObjectM1 <- ObjectM1[order(abs(ObjectM1$avg_log2FC), decreasing=T), ]
    ObjectM1 <- ObjectM1[order(ObjectM1$p_val_adj), ]
  }
  if (Spec1=='Mm') {
    GeneSymb1 = MmscRNA_genes
  }else if (Spec1=='Hs') {
    GeneSymb1 = HsscRNA_genes
  }else if (Spec1=='Zf') {
    GeneSymb1 = ZfscRNA_genes
  }else if (Spec1=='Ch') {
    GeneSymb1 = ChscRNA_genes
  }else if (is.null(Spec1)){
    print('please input correct ')
  }
  if(is.null(Spec1)){
    ObjectM2 <- ObjectM1
  }else{
    if(grepl('ENS', rownames(ObjectM1)[1])){
      GeneID1 <- Converse_GeneIDSymbol(ObjectM1[, 'gene'], GeneInf1 = ZfscRNA_genes)
      ObjectM2 <- cbind(GeneID1, ObjectM1)
    }else{
      ObjectM2 <- ObjectM1
    }
  }

  return(ObjectM2)
}
