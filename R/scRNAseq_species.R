#' Inner function to get wilcox markers cond
#'

get_Wilcox_Markers_Cond<-function(Marker1, Cor1, CorThr1 = 0.2, CorType1 = 'Cor_mmNMDA_mmLD'){
  ShMarker1 <- c(); Marker3 <- list()
  Marker2 <- Marker1[Marker1[, grep('avg_logFC',colnames(Marker1))]>0 & Marker1[, grep('p_val_adj',colnames(Marker1))]<0.01, ]
  Marker2<-na.omit(Marker2)
  Marker3[[1]] <- Marker2
  ShMarker1 <- union(ShMarker1, rownames(Marker2))
  Cor2 <- Cor1[match(ShMarker1, Cor1[,1]),]; Marker5 <- Cor2
  Marker4 <- Marker3[[1]][match(ShMarker1, rownames(Marker3[[1]])), ]
  Marker5 <- cbind(Marker5, Marker4[,3:ncol(Marker4)])
  Marker6 <- Marker5[Marker5[,CorType1]>CorThr1, ]; print(nrow(Marker6))
  Ind1 <- grep('p_val_adj', colnames(Marker6)); Ind2 <- c()
  for(i in 1:length(Ind1)){ Ind2 <- c(Ind2, Ind1[i]-3, Ind1[i]-4, Ind1[i]) }
  Marker7 <- t(apply(Marker6[, Ind2], 1, function(x1){
    x11 <- x1[seq(3,length(x1),by=3)]; x11[is.na(x11)] <- 1; flag1=0;
    if(length(x11)>1){
      for(i in 2:length(x11)){
        if(x11[i]!=x11[1]){ flag1=1; break }
      } }

    if(flag1){ x12 <- which.min(x11);
    }else{ x13 <- x1[seq(1,length(x1),by=3)]; x13[is.na(x13)] <- 0
    x12 <- which.max(abs(x13))
    }
    x2 <- x1[(x12*3-2):(x12*3)]
    return(x2)
  }))
  colnames(Marker7) <- c('Topavg_logFC','Topp_val','Topp_val_adj')
  Marker8 <- cbind(Marker6, Marker7)
  return(Marker8)
}



