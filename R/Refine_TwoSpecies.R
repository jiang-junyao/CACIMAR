#' Title
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
