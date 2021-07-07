#' inner function to refine markers
#'
#' @param Seurat_obejct seurat object
#' @param Spec1 character, indicating the species of data, including: 'Mm', 'Ch', 'Zf', 'Hs'
#' @param GeneSymb1 Gene correspondence that is used for ID chaning. first column should be ENSEMBLE ID, second column should be Symbol, third column should be NCBIID and fourth column should be officalsymbol.
#' @param Marker Marker gene list
#' @param FracThr1 ????
#'
#' @return
#'
#' @examples
Refine_Markers<-function(Seurat_obejct, Spec1, GeneSymb1, Marker, FracThr1=3){
  if(Spec1=='Zf'){ Cluster1 <- 'res.0.2' }else{ Cluster1 <- 'res.0.1' }
  MarkerRoc3 <- Identify_Markers3(Seurat_obejct, Marker, Cluster1, GeneSymb1, FracThr1=FracThr1)
  Ind1 <- grep('Symbol', colnames(MarkerRoc3))
  if(length(Ind1)>1){ MarkerRoc3 <- MarkerRoc3[, -Ind1[2]] }
  MarkerRoc4 <- MarkerRoc3[MarkerRoc3[,'Pvalue'] < 0.01, -grep('EnsemblID', colnames(MarkerRoc3))]
  print(c(nrow(MarkerRoc3), nrow(MarkerRoc4)))
  return(MarkerRoc4)
}


#' Inner function to refine markers
Identify_Markers3 <- function(pbmc, MarkerRoc2, Cluster1, Symb1, FracThr1=5){

  pbmc@assays$RNA@data <- GetAssayData(pbmc)[rownames(MarkerRoc2), ]
  Frac1 <- Get_scRNA_AggExp(pbmc, Symb1, Cluster1=Cluster1, ExpType1='fraction')
  colnames(Frac1) <- gsub('_0', '', colnames(Frac1)); colnames(Frac1) <- gsub('Frac', '', colnames(Frac1))
  MarkerRoc3 <- cbind(MarkerRoc2, Frac1)
  AllCluster1 <- as.matrix(table(pbmc@active.ident)); rownames(AllCluster1) <- paste0('C', rownames(AllCluster1))

  MarkerRoc31 <- apply(MarkerRoc3, 1, function(x1){
    x11 <- strsplit(x1[1], ',')[[1]]
    x12 <- x11[1:(length(x11)-1)]; x13 <- setdiff(rownames(AllCluster1), x12)
    TNum11 <- as.numeric(AllCluster1[match(x12, rownames(AllCluster1)),1])
    TNum12 <- as.numeric(AllCluster1[match(x13, rownames(AllCluster1)),1])
    Frac11 <- as.numeric(x1[match(x12, names(x1))]); mFrac11 <- which.min(Frac11)
    Frac12 <- as.numeric(x1[match(x13, names(x1))]); mFrac12 <- which.max(Frac12)
    TNum21 <- TNum11[mFrac11]; TNum22 <- TNum12[mFrac12]
    Frac21 <- Frac11[mFrac11]; Frac22 <- Frac12[mFrac12]
    TNum31 <- round(TNum21*Frac21); TNum32 <- round(TNum22*Frac22)

    Num4 <- matrix(c(TNum31, TNum21-TNum31, TNum32, TNum22-TNum32), nrow=2)
    Test1 <- fisher.test(Num4)
    if(Frac22!=0){
      if(as.numeric(Test1$estimate)/Frac22 < FracThr1){ return(Frac22)
      }else{ return(Test1$p.value) }
    }else{ return(Test1$p.value) }

  })
  MarkerRoc4 <- cbind(MarkerRoc3, MarkerRoc31); colnames(MarkerRoc4)[ncol(MarkerRoc4)] <- 'Pvalue'

  return(MarkerRoc4)
}


#' Inner function to refine markers

Get_scRNA_AggExp <- function(pbmc, Symb1, Cluster1=c('all', 'res.0.1'), Group1=c('all', 'protocol'), ExpType1=c('expression', 'fraction')){

  source('../Markers_Programs/Programs/Converse_GeneIDSymbol.R')

  if(Cluster1[1]=='all'){ pbmc@meta.data[,'all'] <- 0
  }else if(!is.element(Cluster1, colnames(pbmc@meta.data))){
    stop(paste('please input one of', paste(colnames(pbmc@meta.data)[grep('res',colnames(pbmc@meta.data))], collapse=',') ))
  }
  if(Group1[1]=='all'){ pbmc@meta.data[,'all'] <- 0
  }else if(!is.element(Group1, colnames(pbmc@meta.data))){
    stop(paste('please input one of', paste(colnames(pbmc@meta.data), collapse=',') ))
  }

  Meta01 <- pbmc@meta.data
  Cluster2 <- sort(unique(Meta01[, Cluster1[1]]))
  Group2 <- sort(unique(Meta01[, Group1[1]]))

  Name01 <- c(); flag01 <- 0;
  for(i0 in 1:length(Cluster2)){
    for(j0 in 1:length(Group2)){
      Cell01 <- Meta01[Meta01[, Cluster1[1]]==Cluster2[i0] & Meta01[, Group1[1]]==Group2[j0], ]
      Object02 <- as.matrix(GetAssayData(pbmc)[, rownames(Cell01)]); print(c(Cluster2[i0], as.character(Group2[j0]), ncol(Object02)))

      if(ncol(Object02)==0){
      }else{
        Frac1 <- t(apply(Object02, 1, function(x1){
          x12 <- length(x1[x1!=0])
          x2 <- x12/length(x1)
          return(c(x12, x2))
        }))
        if(is.element('expression', ExpType1) & is.element('fraction', ExpType1)){ mObject03 <- cbind(rowMeans(Object02), Frac1[,2])
        Name01 <- c(Name01, paste0(c('Exp','Frac'),paste0('C',Cluster2[i0],'_',Group2[j0])))
        }else if(is.element('expression', ExpType1)){ mObject03 <- rowMeans(Object02)
        Name01 <- c(Name01, paste0('ExpC',Cluster2[i0],'_',Group2[j0]))
        }else if(is.element('fraction', ExpType1)){ mObject03 <- Frac1[,2]
        Name01 <- c(Name01, paste0('FracC',Cluster2[i0],'_',Group2[j0]))
        }else if(is.element('number', ExpType1)){ mObject03 <- Frac1[,1]
        Name01 <- c(Name01, paste0('NumC',Cluster2[i0],'_',Group2[j0]))
        }else{ stop(paste0('Please choose ExpType1 for expression or fraction')) }

        if(flag01==0){ mObject04 <- mObject03; flag01 <- 1;
        }else{ mObject04 <- cbind(mObject04, mObject03) }
      } } }
  mObject04 <- as.matrix(mObject04); colnames(mObject04) <- Name01
  if(is.element('expression', ExpType1) & is.element('fraction', ExpType1)){
    mObject04 <- mObject04[, c(seq(1,ncol(mObject04),2), seq(2,ncol(mObject04),2))]
  }

  if(grepl('^ENS', rownames(GetAssayData(pbmc))[1])){
    GeneID1 <- Converse_GeneIDSymbol(rownames(mObject04), Symb1)
    mObject05 <- cbind(GeneID1, mObject04)
  }else{ mObject05 <- mObject04 }

  return(mObject05)
}
