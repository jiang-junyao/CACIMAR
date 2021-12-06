#' Calculate fraction of shared markers between species
#'
#' @param Marker_list each element should be marker genes table
#' @param Species_names character, indicating the species of each marker genes
#' table in Marker_list
#' @param PowerTh1 numeric, indicating the threshold to filter marker genes
#'
#' @return
#' @export
#'
#' @examples
Identify_SharedMarkers <- function(Marker_list, Species_names, PowerTh1=0.1){
  Exp1 <- list(); ExpInd1 <- list(); AllCluster1 <- c()
  PowerTh12 <- gsub('\\.','',PowerTh1)
  for(i in 1:length(Marker_list)){
    Exp01 <- Marker_list[[i]]
    print(nrow(Exp01))
    Exp02 <- Exp01[Exp01$avg_diff>0 & Exp01$power>PowerTh1, ]
    Exp02[, 'cluster'] <- paste0(Species_names[i], 'C', Exp02[, 'cluster'])
    print(table(Exp02[, 'cluster'])); AllCluster1 <- c(AllCluster1, unique(Exp02[, 'cluster']))
    Exp1[[i]] <- Exp02

    ShMarker1 <- Cal_SharedMarkers(Exp1[[i]], Species_names[i])
    if(i==1){ ShMarker3 <- ShMarker1[[1]]
    Frac1 <- ShMarker1[[3]]
    }else{ ShMarker3 <- rbind(ShMarker3, ShMarker1[[1]])
    Frac1 <- rbind(Frac1, ShMarker1[[3]])
    } }

  print('Calculate fraction of shared markers between species')
  cFile1 <- combn(length(Marker_list), 2)
  for(i in 1:ncol(cFile1)){ print(cFile1[,i])
    Spec1 <- Species_names[cFile1[1,i]]
    Spec2 <- Species_names[cFile1[2,i]]
    if (Spec1 == 'ch' & Spec2 =='mm' | Spec2 == 'ch' & Spec1 =='mm') {
      OrthG01 <- OrthG_Mm_Ch
    }else if (Spec1 == 'zf' & Spec2 =='mm' | Spec2 == 'zf' & Spec1 =='mm') {
      OrthG01 <- OrthG_Mm_Zf
    }else if (Spec1 == 'ch' & Spec2 =='zf' | Spec2 == 'ch' & Spec1 =='zf') {
      OrthG01 <- OrthG_Zf_Ch
    }
    OrthG1 <- OrthG01[-grep('0', OrthG01$Type), ]

    for(j in 1:2){
      uExp1 <- unique(Exp1[[cFile1[j,i]]][, 'gene'])
      ExpInd01 <- t(apply(as.matrix(uExp1), 1, function(x1){
        x2 <- grep(x1, OrthG1[,j*2])
        if(length(x2)==1){ x3 <- c(x1, x2) }else{ x3 <- c(x1, NA) }
      }) )
      ExpInd1[[j]] <- ExpInd01[!is.na(ExpInd01[, 2]), ]
      colnames(ExpInd1[[j]]) <- c('ID', 'RowNum')
    }

    mmExp1 <- Exp1[[cFile1[1,i]]]
    zfExp1 <- Exp1[[cFile1[2,i]]]
    mmExpInd1 <- ExpInd1[[1]]
    zfExpInd1 <- ExpInd1[[2]]
    ShMarker2 <- Cal_SharedMarkers_Species(mmExp1, zfExp1, mmExpInd1, zfExpInd1, Species_names[cFile1[,i]])
    ShMarker3 <- rbind(ShMarker3, ShMarker2[[1]])
    Frac1 <- rbind(Frac1, ShMarker2[[3]])
  }

  AllCluster12 <- sort(AllCluster1); Frac3 <- c()
  for(i in 1:length(AllCluster1)){ Frac21 <- c()
  for(j in 1:length(AllCluster1)){
    Frac11 <- Frac1[Frac1[,1]==AllCluster1[i] & Frac1[,2]==AllCluster1[j], 3]
    Frac12 <- Frac1[Frac1[,1]==AllCluster1[j] & Frac1[,2]==AllCluster1[i], 3]
    if(length(Frac11)==1){ Frac2 <- as.numeric(Frac11)
    }else if(length(Frac12)==1){ Frac2 <- as.numeric(Frac12)
    }else{ print(paste('No',AllCluster1[i],AllCluster1[j])) }
    Frac21 <- c(Frac21, Frac2)
  }
  Frac3 <- rbind(Frac3, Frac21)
  }
  rownames(Frac3) <- AllCluster1; colnames(Frac3) <- AllCluster1

  ShMarker4 <- list(); ShMarker4[[1]] <- ShMarker3; ShMarker4[[2]] <- Frac3; ShMarker4[[3]] <- Frac1

  return(ShMarker4)
}

Cal_SharedMarkers <- function(mmExp1, Spec1='mm'){
  print('Calculate fraction of shared markers across clusters within species')

  mmCluster1 <- sort(unique(mmExp1$cluster))
  ShMarker1 <- c(); Fraction3 <- c(); Fraction4 <- c()
  for(i in 1:length(mmCluster1)){
    zfExp2 <- mmExp1[mmExp1$cluster==mmCluster1[i], ]
    zfPower1 <- sum(zfExp2$power)
    Fraction2 <- c()
    for(j in 1:length(mmCluster1)){
      mmExp2 <- mmExp1[mmExp1$cluster==mmCluster1[j], ]
      mmPower1 <- sum(mmExp2$power)
      Exp3 <- intersect(zfExp2[, 'gene'], mmExp2[, 'gene'])

      zfShMarker1 <- 'none'; zfShPower1 <-0;
      mmShMarker1 <- 'none'; mmShPower1 <-0;
      zfPower2 <- 0; mmPower2 <- 0
      if(length(Exp3)>0){
        zfExp3 <- zfExp2[match(Exp3, zfExp2[,'gene']), ]
        mmExp3 <- mmExp2[match(Exp3, mmExp2[,'gene']), ]
        zfPower2 <- sum(zfExp3$power); mmPower2 <- sum(mmExp3$power)
      }
      Fraction1 <- (zfPower2 + mmPower2)/(zfPower1 + mmPower1)
      Fraction2 <- c(Fraction2, Fraction1); Fraction4 <- rbind(Fraction4, c(mmCluster1[i], mmCluster1[j], Fraction1))
      ShMarker1 <- rbind(ShMarker1, c(mmCluster1[i], mmCluster1[j], dim(zfExp2)[1], dim(mmExp2)[1], dim(zfExp2)[1], dim(mmExp2)[1], length(Exp3), Fraction1, zfShMarker1, zfShPower1, mmShMarker1, mmShPower1))
    }
    Fraction3 <- rbind(Fraction3, Fraction2)
  }
  rownames(Fraction3) <- paste0('C',mmCluster1); colnames(Fraction3) <- paste0('C',mmCluster1);
  colnames(Fraction4) <- c('Cluster1', 'Cluster2', 'Fraction')
  colnames(ShMarker1) <- paste0(Spec1, c('Cluster', 'Cluster', 'MarkersNum', 'MarkersNum', 'OrthMarkersNum', 'OrthMarkersNum', 'ShMarkersNum', 'SharedFraction', 'SharedMarkers', 'ShMarkersPower', 'SharedMarkers', 'ShMarkersPower'))
  ShMarker2 <- list(); ShMarker2[[1]] <- ShMarker1; ShMarker2[[2]] <- Fraction3; ShMarker2[[3]] <- Fraction4;

  return(ShMarker2)
}

Cal_SharedMarkers_Species <- function(mmExp1, zfExp1, mmExpInd1, zfExpInd1, Spec1=c('mm','zf')){
  zfCluster1 <- sort(unique(zfExp1$cluster))
  mmCluster1 <- sort(unique(mmExp1$cluster))
  ShMarker1 <- c(); Fraction3 <- c(); Fraction4 <- c()
  for(i in 1:length(zfCluster1)){
    zfExp2 <- zfExp1[zfExp1$cluster==zfCluster1[i], ]
    zfPower1 <- sum(zfExp2$power)
    zfExp22 <- zfExpInd1[match(zfExp2[,'gene'], zfExpInd1[,'ID']), ]

    Fraction2 <- c()
    for(j in 1:length(mmCluster1)){
      mmExp2 <- mmExp1[mmExp1$cluster==mmCluster1[j], ]
      mmPower1 <- sum(mmExp2$power)
      mmExp22 <- mmExpInd1[match(mmExp2[, 'gene'], mmExpInd1[,'ID']), ]
      Exp3 <- intersect(zfExp22[, 2], mmExp22[, 2])
      Exp32 <- Exp3[!is.na(Exp3)]

      zfShMarker1 <- 'none'; zfShPower1 <-0;
      mmShMarker1 <- 'none'; mmShPower1 <-0;
      zfPower2 <- 0; mmPower2 <- 0
      if(length(Exp32)>0){
        Exp33 <- apply(as.matrix(Exp32), 1, function(x1){ x2 <- c()
        for(k in 1:2){
          if(k==1){ Exp23 <- zfExp22; Exp2 <- zfExp2
          }else{ Exp23 <- mmExp22; Exp2 <- mmExp2 }
          Exp41 <- Exp23[is.element(Exp23[,2], x1), ]
          if(is.null(dim(Exp41))){ Exp42 <- Exp2[is.element(Exp2[,'gene'], Exp41[1]), ]
          }else{ Exp42 <- Exp2[is.element(Exp2[,'gene'], Exp41[,1]), ] }
          ShMarker1 <- paste(Exp42[, 'gene'], collapse=',')
          ShPower1 <- paste(Exp42[, 'power'], collapse=',')
          Power2 <- sum(Exp42[, 'power'])
          x2 <- c(x2, ShMarker1, ShPower1, Power2)
        }
        return(x2)
        })
        zfShMarker1 <- paste(Exp33[1, ], collapse=';')
        zfShPower1 <- paste(Exp33[2, ], collapse=';')
        mmShMarker1 <- paste(Exp33[4, ], collapse=';')
        mmShPower1 <- paste(Exp33[5, ], collapse=';')
        zfPower2 <- sum(as.numeric(Exp33[3, ])); mmPower2 <- sum(as.numeric(Exp33[6, ]))
      }
      Fraction1 <- (zfPower2 + mmPower2)/(zfPower1 + mmPower1)
      Fraction2 <- c(Fraction2, Fraction1); Fraction4 <- rbind(Fraction4, c(zfCluster1[i], mmCluster1[j], Fraction1))
      ShMarker1 <- rbind(ShMarker1, c(zfCluster1[i], mmCluster1[j], dim(zfExp22)[1], dim(mmExp22)[1], dim(zfExp22)[1], dim(mmExp22)[1], length(Exp32), Fraction1, zfShMarker1, zfShPower1, mmShMarker1, mmShPower1))
    }
    Fraction3 <- rbind(Fraction3, Fraction2)
  }
  rownames(Fraction3) <- paste0('C',zfCluster1); colnames(Fraction3) <- paste0('C',mmCluster1);
  colnames(Fraction4) <- c('Cluster1', 'Cluster2', 'Fraction')
  colnames(ShMarker1) <- c(paste0(Spec1[1],'Cluster'), paste0(Spec1[2],'Cluster'), paste0(Spec1[1],'MarkersNum'), paste0(Spec1[2],'MarkersNum'), paste0(Spec1[1],'OrthMarkersNum'), paste0(Spec1[2],'OrthMarkersNum'), 'ShMarkersNum', 'SharedFraction', paste0(Spec1[1],'SharedMarkers'), paste0(Spec1[1],'ShMarkersPower'), paste0(Spec1[2],'mmSharedMarkers'), paste0(Spec1[2],'ShMarkersPower'))
  ShMarker2 <- list(); ShMarker2[[1]] <- ShMarker1; ShMarker2[[2]] <- Fraction3; ShMarker2[[3]] <- Fraction4;

  return(ShMarker2)
}
