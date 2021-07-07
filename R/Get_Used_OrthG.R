Get_Used_OrthG<-function(OrthG1,mmExp01,zfExp01,Species){
  colnames(mmExp01)[1:2] <- paste0(Species[1],colnames(mmExp01)[1:2])
  colnames(zfExp01)[1:2] <- paste0(Species[2],colnames(zfExp01)[1:2])
  for(i in 1:2){
    if(i==1){ Exp01 <- mmExp01
    }else{ Exp01 <- zfExp01 }
    Ind1 <- grep('Power', colnames(Exp01))
    Power1 <- t(apply(Exp01, 1, function(x1){
      x11 <- x1[Ind1]; x12 <- x11[!is.na(x11)]; x13 <- c()
      for(i in 1:length(x12)){
        x13 <- c(x13, strsplit(x12[i],',')[[1]][1])
      }
      x2 <- mean(as.numeric(x13))
      return(c(x2,x2))
    }))
    if(i==1){ mmExp1 <- as.data.frame(Power1);
    }else{ zfExp1 <- as.data.frame(Power1) }
  }
  Exp2 <- Get_OrthG(OrthG1, mmExp1, zfExp1, Species)
  Type1 <- paste0('Used_',Species[1],'_ID'); Type2 <- paste0('Used_',Species[2],'_ID');
  mmExp2 <- mmExp1[match(Exp2[, Type1], rownames(mmExp1)), ]
  zfExp2 <- zfExp1[match(Exp2[, Type2], rownames(zfExp1)), ]
  Exp3 <- cbind(Exp2, mmExp2, zfExp2)
  Exp4 <- Exp3[!is.na(Exp3[, dim(Exp2)[2]+1]) & !is.na(Exp3[, dim(Exp2)[2]+dim(mmExp2)[2]+1]), ]
  Exp5 <- cbind(Exp4[,1:7], mmExp01[match(Exp4[,Type1], rownames(mmExp01)), ], zfExp01[match(Exp4[,Type2], rownames(zfExp01)), ])
  Exp5 <- Refine_Used_OrthG(Exp5,Species)
  return(Exp5)
}

#' Inner function to get orthology

Get_OrthG <- function(OrthG1, MmRNA1, ZfRNA1, Spec1, MmPattern1='', ZfPattern1=''){
  tOrthG1 <- table(OrthG1$Type); print(tOrthG1)
  Ind1 <- c(grep(paste0(Spec1[1],'_ID'), colnames(OrthG1)), grep(paste0(Spec1[2],'_ID'), colnames(OrthG1)))

  OrthG21 <- list()
  for(i in 1:length(tOrthG1)){ print(names(tOrthG1)[i])
    if(names(tOrthG1)[i]==paste(c(Spec1,'0T1'),collapse='_')){ OrthG2 <- OrthG1[OrthG1$Type==names(tOrthG1)[i], ]
    OrthG3 <- cbind(OrthG2, OrthG2[,Ind1])
    }else if(names(tOrthG1)[i]==paste(c(Spec1,'1T0'),collapse='_')){ OrthG2 <- OrthG1[OrthG1$Type==names(tOrthG1)[i], ]
    OrthG3 <- cbind(OrthG2, OrthG2[,Ind1])
    }else if(names(tOrthG1)[i]==paste(c(Spec1,'1T1'),collapse='_')){ OrthG2 <- OrthG1[OrthG1$Type==names(tOrthG1)[i], ]
    OrthG3 <- cbind(OrthG2, OrthG2[,Ind1])
    }else if(grepl(paste(c(Spec1,'.*N'),collapse='_'),names(tOrthG1)[i])){ OrthG2 <- OrthG1[OrthG1$Type==names(tOrthG1)[i], ]
    OrthG21 <- apply(OrthG2, 1 ,function(x1){
      for(j in 1:length(Ind1)) { x2 <- strsplit(as.character(x1[Ind1[j]]),'[;,]')[[1]]
      if(j==1){ RNA21 <- MmRNA1[match(x2, rownames(MmRNA1)), ]
      if(MmPattern1[1]!=''){ Pattern21 <- MmPattern1[match(x2, rownames(MmRNA1)), ] }
      }else{ RNA21 <- ZfRNA1[match(x2, rownames(ZfRNA1)), ]
      if(MmPattern1[1]!=''){ Pattern21 <- ZfPattern1[match(x2, rownames(ZfRNA1)), ] }
      }
      RNA2 <- RNA21[!is.na(RNA21[,1]), ];
      if(MmPattern1[1]!=''){ Pattern2 <- Pattern21[!is.na(RNA21[,1])] }

      if(is.null(dim(RNA2)) | dim(RNA2)[1]==1){
        RNA3 <- RNA2; Gene1 <- rownames(RNA2)
      }else if(dim(RNA2)[1]==0){
        RNA3 <- RNA2[1, ]; Gene1 <- rownames(RNA2)[1]
      }else{
        if(MmPattern1[1]!=''){ RNA31 <- RNA2[grepl('[UD]', Pattern2), ]
        if(dim(RNA31)[1]!=0){ RNA32 <- RNA31
        }else{ RNA32 <- RNA2 }
        }else{ RNA32 <- RNA2 }
        mRNA32 <- which.max(rowMeans(RNA32))
        RNA3 <- RNA32[mRNA32, ]; Gene1 <- rownames(RNA32)[mRNA32]
      }

      if(j==1){ RNA4 <- as.numeric(RNA3); Gene2 <- Gene1
      }else{ RNA4 <- c(RNA4, as.numeric(RNA3)); Gene2 <- c(Gene2, Gene1) }
      }
      return(Gene2)
    } )
    rownames(OrthG21) <- colnames(OrthG1)[Ind1]
    OrthG3 <- cbind(OrthG2, t(OrthG21))
    }else{ print(paste0('Not process ',tOrthG1[i])) }

    if(i==1){ OrthG4 <- as.matrix(OrthG3);
    }else{ OrthG4 <- rbind(OrthG4, as.matrix(OrthG3)) }
  }
  colnames(OrthG4)[(ncol(OrthG1)+1):(ncol(OrthG1)+2)] <- paste0('Used_',colnames(OrthG4)[(ncol(OrthG1)+1):(ncol(OrthG1)+2)])

  return(OrthG4)
}

#' Inner function to refine orthology
Refine_Used_OrthG<-function(ShMarker1,Species){
  Type1 <- paste0(Species[1],'.*\\.AllCluster'); Type2 <- paste0(Species[2], '.*\\.AllCluster')
  Spec1Type1 <- grep(Type1, colnames(ShMarker1))
  Spec1Type2 <- grep(Type2, colnames(ShMarker1))
  ShMarker2 <- apply(ShMarker1, 1, function(x1){
    x11 <- x1[Spec1Type1]; x12 <- x1[Spec1Type2]; x2 <- F
    for(i in 1:length(x11)){
      if(!is.na(x11[i])){
        x112 <- strsplit(x11[i], ',')[[1]]
        for(i1 in 1:length(x112)){
          for(j in 1:length(x12)){
            if(!is.na(x12[j])){
              x122 <- strsplit(x12[j], ',')[[1]]
              for(j1 in 1:length(x122)){
                if(x112[i1]==x122[j1] | grepl('BC',x112[i1]) & grepl('BC',x122[j1])){
                  x2 <- T
                } } } } } } }
    return(x2)
  })
  ShMarker3 <- ShMarker1[ShMarker2, ]
  ShMarker4 <- ShMarker3[order(ShMarker3[, grep('Power', colnames(ShMarker3))[1]], decreasing=T), ]
  ShMarker4 <- ShMarker4[order(ShMarker4[, Spec1Type1[1]]), ]
  return(ShMarker4)
}
