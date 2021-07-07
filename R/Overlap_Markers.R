Overlap_Markers_Cond<-function(RNA2,CellT1,Spec2){
  Cluster1 <- apply(CellT1, 1, function(x1){ x2 <- paste(x1[1:2], collapse='.') })
  CellT12 <- cbind(CellT1, Cluster1)
  uCellT1 <- unique(CellT12$CellType)
  for(i in 1:length(uCellT1)){
    CellT2 <- CellT12[CellT12$CellType==uCellT1[i], ]
    Ind2 <- match(CellT2[, 'Cluster1'], colnames(RNA2))
    RNA3 <- RNA2[, Ind2]
    if(length(Ind2)==1){ RNA31 <- RNA3
    }else if(length(Ind2)>1){
      RNA31 <- apply(RNA3, 1, function(x1){
        x12 <- x1[!is.na(x1)]; x2 <- mean(x12)
        return(x2)
      })
    }else{ print(uCellT1[i]) }
    if(i==1){ RNA4 <- RNA31
    }else{ RNA4 <- cbind(RNA4, RNA31) }
  }
  colnames(RNA4) <- paste0(Spec2, '.',uCellT1)
  Ind3 <- c(match(c('Symbol','Type'), colnames(RNA2)), grep('Power', colnames(RNA2)))
  RNA41 <- RNA2[, Ind3]
  RNA42 <- RNA2[, grep('Cluster', colnames(RNA2))]; dRNA42 <- dim(RNA42); Name1 <- c()
  for(i in 1:dRNA42[2]){ Name1[i] <- strsplit(colnames(RNA42)[i], '\\.')[[1]][1] }
  RNA43 <- t(apply(RNA42, 1, function(x1){
    for(i in 1:dRNA42[2]){
      if(!is.na(x1[i])){
        x12 <- strsplit(x1[i], ',')[[1]]; Cluster2 <- c()
        for(j in 1:length(x12)){
          Cluster2[j] <- as.character(CellT12[CellT12[,'Cluster1']==paste0(Name1[i],'.',x12[j]), 'CellType'])
        }
        x1[i] <- paste(Cluster2, collapse=',');
      } }
    return(x1)
  }))
  RNA5 <- cbind(RNA41, RNA43, RNA4)
  return(RNA5)
}
