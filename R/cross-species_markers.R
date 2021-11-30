#' Overlap markers from different conditions of the same species
#'
#' @param markers_expression marker genes information with expression
#' @param celltype cell type data.frame, first column should be samples, second
#' column should be cluster, thrid column should be cell type
#' @param Spec1 character, indicating the species of your data
#'
#' @return marker genes information in different samples
#' @export
#'
#' @examples
Overlap_Markers_Cond<-function(markers_expression, celltype, Spec1){
  CellT1 = celltype  ; Spec2 = Spec1 ; RNA2 <- markers_expression
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

#' Identify orthologs genes based on orthologs database
#'
#' @param OrthG1
#' @param Species1_expression
#' @param Species2_expression
#' @param Species
#'
#' @return
#' @export
#'
#' @examples
Get_Used_OrthG<-function(OrthG1,Species1_expression,Species2_expression,Species){
  colnames(Species1_expression)[1:2] <- paste0(Species[1],
                                               colnames(Species1_expression)[1:2])
  colnames(Species2_expression)[1:2] <- paste0(Species[2],
                                               colnames(Species2_expression)[1:2])
  Species1_expression2 <- Species1_expression
  Species2_expression2 <- Species2_expression
  for(i in 1:2){##average the marker gene power in each condition
    if(i==1){ Exp01 <- Species1_expression
    }else{ Exp01 <- Species2_expression }
    Ind1 <- grep('Power', colnames(Exp01))
    Power1 <- t(apply(Exp01, 1, function(x1){
      x11 <- x1[Ind1]; x12 <- x11[!is.na(x11)]; x13 <- c()
      for(i in 1:length(x12)){
        x13 <- c(x13, strsplit(x12[i],',')[[1]][1])
      }
      x2 <- mean(as.numeric(x13))
      return(c(x2,x2))
    }))
    if(i==1){ Species1_expression <- as.data.frame(Power1);
    }else{ Species2_expression <- as.data.frame(Power1) }
  }
  ### only marker genes in orthG database will be retained
  Exp2 <- Get_OrthG(OrthG1, Species1_expression, Species2_expression, Species)
  Type1 <- paste0('Used_',Species[1],'_ID')
  Type2 <- paste0('Used_',Species[2],'_ID')
  Species1_expression <- Species1_expression[match(Exp2[, Type1],
                                                   rownames(Species1_expression)), ]
  Species2_expression <- Species2_expression[match(Exp2[, Type2],
                                                   rownames(Species2_expression)), ]
  Exp3 <- cbind(Exp2, Species1_expression, Species2_expression)
  Exp4 <- Exp3[!is.na(Exp3[, dim(Exp2)[2]+1]) & !is.na(Exp3[, dim(Exp2)[2]+dim(Species1_expression)[2]+1]), ]
  Exp5 <- cbind(Exp4[,1:7], Species1_expression2[match(Exp4[,Type1],
                                                      rownames(Species1_expression2)), ],Species2_expression2[match(Exp4[,Type2], rownames(Species2_expression2)), ])
 ###For each marker gene, if appear in the same cluster of two species, it will be retained
   Exp6 <- Refine_Used_OrthG(Exp5,Species)
  return(Exp6)
}


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


#' Refine cluster markers between two species
#'
#' @param ShMarker11
#' @param CellT11
#' @param CellT12
#' @param Species
#'
#' @return
#' @export
#'
#' @examples
Refine_TwoSpecies<-function(ShMarker11,CellT11,CellT12,Species){
  uCellT11 <- intersect(CellT11$CellType, CellT12$CellType)
  uCellT12 <- 'BC'
  uCellT2 <- c(uCellT11, uCellT12)
  SpecInd11 <- grep(paste0(Species[1],'.*\\.AllCluster'), colnames(ShMarker11))
  SpecInd12 <- grep(paste0(Species[2],'.*\\.AllCluster'), colnames(ShMarker11))
  UsedID1 <- paste0('Used_',Species,'_ID'); ShMarker4 <- c()
  for(i in 1:length(uCellT2)){ ShMarker22 <- list()
  for(i1 in 1:1){
    if(i1==1){ ShMarker1 <- ShMarker11
    SpecInd31 <- SpecInd11; SpecInd32 <- SpecInd12;
    }else{ ShMarker1 <- ShMarker12
    SpecInd31 <- SpecInd21; SpecInd32 <- SpecInd22;
    }

    ShMarker21 <- apply(ShMarker1, 1, function(x1){
      x11 <- x1[SpecInd31]; x12 <- x1[SpecInd32]; x2 <- F
      for(j in 1:length(x11)){
        if(!is.na(x11[j])){
          x112 <- strsplit(x11[j], ',')[[1]]
          for(j1 in 1:length(x112)){
            for(k in 1:length(x12)){
              if(!is.na(x12[k])){
                x122 <- strsplit(x12[k], ',')[[1]]
                for(k1 in 1:length(x122)){
                  if(x112[j1]==x122[k1] & x112[j1]==uCellT2[i] |
                     i==length(uCellT2) & grepl(uCellT2[i], x112[j1]) & grepl(uCellT2[i], x122[k1])){
                    x2 <- T
                  } } } } } } }
      return(x2)
    })
    ShMarker22[[i1]] <- ShMarker1[ShMarker21, ]
  }


  ShMarker2 <- intersect(ShMarker22[[1]][, UsedID1[1]],ShMarker22[[1]][, UsedID1[1]])
  if(length(ShMarker2)>0){
    ShMarker31 <- ShMarker22[[1]][match(ShMarker2, ShMarker22[[1]][, UsedID1[1]]), ]

    ShMarker4 <- rbind(ShMarker4, ShMarker31)
  }
  }
  ShMarker5<-Refine_Markers_Species(ShMarker4,Species)
  return(ShMarker5)
}

#' Inner function to refine markers

Refine_Markers_Species<-function(Marker1,Species){
  PowerTh1 <- 0.4; PowerTh2 <- gsub('\\.','',PowerTh1)
  SpecInd1 <- length(Species);
  Ind1 <- list(); Ind2 <- list();
  for(i in 1:length(Species)){
    Ind1[[i]] <- grep(paste0(Species[i],'.*\\.AllCluster'), colnames(Marker1))
    Ind2[[i]] <- grep(paste0(Species[i],'.*\\.Power'), colnames(Marker1))
  }

  Marker2 <- apply(Marker1, 1, function(x1){
    x2 <- F; x14 <- list();
    for(i in 1:length(Ind1)){
      x12 <- x1[Ind1[[i]]]
      for(i1 in 1:length(x12)){
        if(!is.na(x12[i1])){
          x13 <- strsplit(x12[i1], ',')[[1]]
          for(i2 in 1:length(x13)){
            x14[[x13[i2]]] <- x13[i2]
          } }
      } }
    x41 <- list()
    for(i in 1:length(x14)){
      x15 <- rep(0, length(Ind1))
      for(i1 in 1:length(Ind1)){
        x12 <- x1[Ind1[[i1]]]; x22 <- x1[Ind2[[i1]]]
        for(i2 in 1:length(x12)){
          if(!is.na(x12[i2])){
            x13 <- strsplit(x12[i2], ',')[[1]]
            x23 <- strsplit(x22[i2], ',')[[1]]
            for(i3 in 1:length(x13)){
              if(x14[[i]]==x13[i3] | grepl('MG',x14[[i]]) & grepl('MG',x13[i3]) | grepl('BC',x14[[i]]) & grepl('BC',x13[i3])){
                x15[i1] <- max(x15[i1], x23[i3])
              } } } } }

      x31 <- T; x32 <- F
      for(j in 1:length(Ind1)){
        if(as.numeric(x15[j])==0){ x31 <- F }
        if(as.numeric(x15[j])>PowerTh1){ x32 <- T
        } }

      if(x31==T & x32==T){ x41[[x14[[i]]]] <- x15
      } }

    if(length(x41)>0){ x2 <- T }

    return(x2)
  })

  Marker3 <- Marker1[Marker2, ]
  print(c(nrow(Marker1), nrow(Marker3)))
  return(Marker3)
}
