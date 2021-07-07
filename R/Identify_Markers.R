#' Identify markers of each cluster
#'
#' @param Seurat_object Seurat object, should contain cluster information
#' @param PowerThr1 numeric, indicating
#' @param Spec1 character, When GeneSymb1 == NULL, inner file will be used for iD changing, this parameter indicate the species of
#' @param GeneSymb1 Gene correspondence that is used for ID chaning. first column should be ENSEMBLE ID, second column should be Symbol, third column should be NCBIID and fourth column should be officalsymbol.
#' @param FracThr1 numeric, indicating
#'
#' @return
#' @export
#'
#' @examples
Identify_Markers<-function(Seurat_object, PowerThr1=0.4, Spec1=NULL, GeneSymb1=NULL, FracThr1=3 ){
  if (is.null(GeneSymb1)) {
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
  }
  MarkerRoc<-Identify_Markers1(Seurat_object,PowerThr1)
  MarkerRoc<-as.data.frame(MarkerRoc)
  Marker<-Identify_Markers2(Seurat_object,MarkerRoc,GeneSymb1=GeneSymb1,PowerThr1=PowerThr1)
  final_Markers<-Refine_Markers(Seurat_object,Spec1,GeneSymb1,Marker,FracThr1)
  return(final_Markers)
}
