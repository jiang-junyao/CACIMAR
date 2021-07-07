#' Inner function to make ID change
#'

Converse_GeneIDSymbol <- function(Gene1, GeneInf1=NULL){

  Gene1 <- as.character(Gene1)
  if(grepl('ENS', Gene1[1])){ GeneID01 <- Gene1;
  Ind1 <- match(GeneID01, GeneInf1[,1]);
  if(length(Ind1[!is.na(Ind1)])>0){ GeneID1 <- GeneID01[!is.na(Ind1)]
  Symb1 <- as.character(GeneInf1[Ind1[!is.na(Ind1)], 'OfficialSymbol'])
  NCBI1 <- GeneInf1[Ind1[!is.na(Ind1)], 'NCBIID']
  if(length(GeneID01[is.na(Ind1)])>0){
    message(paste('Warning: no', paste(GeneID01[is.na(Ind1)], collapse=' ')))
  }
  }else{ stop('Error: no such gene') }

  }else{ Symb01 <- Gene1;
  Ind1 <- apply(as.matrix(Symb01), 1, function(x1){
    if(is.element(x1, GeneInf1[, 1])){
      x2 <- match(x1, GeneInf1[, 1])
    }else if(is.element(x1, GeneInf1[, 4])){
      x2 <- match(x1, GeneInf1[, 4])
    }else{ x2 <- NA }
    return(x2)
  } )

  if(length(Ind1[!is.na(Ind1)])>0){ GeneID1 <- rownames(GeneInf1[Ind1[!is.na(Ind1)], ])
  Symb1 <- Symb01[!is.na(Ind1)]
  NCBI1 <- GeneInf1[Ind1[!is.na(Ind1)], 'NCBIID']
  if(length(Symb01[is.na(Ind1)])>0){
    message(paste('Warning: no', paste(Symb01[is.na(Ind1)], collapse=' ')))
  }
  }else{ stop('Error: no such gene') }
  }

  GeneID2 <- cbind(GeneID1, Symb1, NCBI1); colnames(GeneID2) <- c('EnsemblID', 'Symbol', 'NCBIID')

  return(GeneID2)
}
