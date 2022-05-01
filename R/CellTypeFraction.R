#' Identify conserved cell types based on power of genes and orthologs database
#'
#' @param OrthG ortholog genes database
#' @param Species1_Marker_table data.frame of species 1, should contain three column:
#' 'gene', 'cluster' and 'power'
#' @param Species2_Marker_table data.frame of species 2, should contain three column:
#' 'gene', 'cluster' and 'power'
#' @param Species_name1 character, indicating the species names of Species1_Marker_table
#' @param Species_name2 character, indicating the species names of Species2_Marker_table
#'
#' @return list contains two elements: first one is details of conserved cell types,
#' second one is matrix of cell types conserved score
#' @export
#'
#' @examples load(system.file("extdata", "CellTypeAllMarkers.rda", package = "CACIMAR"))
#' expression <- Identify_ConservedCellTypes(OrthG_Mm_Zf,mm_Marker,zf_Marker,'mm','zf')
Identify_ConservedCellTypes <- function(OrthG,Species1_Marker_table,Species2_Marker_table,
                                        Species_name1,Species_name2){
  ### inputs validation
  if (!'power' %in% colnames(Species1_Marker_table)) {
    stop('Species1_Marker_table should contain power column!')
  }
  if (!'power' %in% colnames(Species2_Marker_table)) {
    stop('Species2_Marker_table should contain power column!')
  }
  if (!'cluster' %in% colnames(Species1_Marker_table)) {
    stop('Species1_Marker_table should contain cluster column!')
  }
  if (!'cluster' %in% colnames(Species2_Marker_table)) {
    stop('Species2_Marker_table should contain cluster column!')
  }
  if (!'gene' %in% colnames(Species1_Marker_table)) {
    stop('Species1_Marker_table should contain gene column!')
  }
  if (!'gene' %in% colnames(Species2_Marker_table)) {
    stop('Species2_Marker_table should contain gene column!')
  }
  validInput(OrthG,'OrthG','df')
  validInput(Species_name1,'Species_name1','character')
  validInput(Species_name2,'Species_name2','character')
  Species_name1 <- tolower(Species_name1)
  Species_name2 <- tolower(Species_name2)
  Spec1 <- colnames(OrthG)[2]
  Spec2 <- colnames(OrthG)[4]
  Spec1 <- gsub('_ID','',Spec1)
  Spec2 <- gsub('_ID','',Spec2)
  if (Spec1 == Species_name1 & Spec2 == Species_name2) {
    Species_name <- c(Spec1,Spec2)
    Species1_Marker <- Species1_Marker_table
    Species2_Marker <- Species2_Marker_table
  }else if(Spec2 == Species_name1 & Spec1 == Species_name2){
    Species_name <- c(Spec2,Spec1)
    Species2_Marker <- Species1_Marker_table
    Species1_Marker <- Species2_Marker_table
  }else{stop('please input correct Species name')}
  Species1_Marker_table$power <- as.numeric(Species1_Marker_table$power)
  Species2_Marker_table$power <- as.numeric(Species2_Marker_table$power)

  ### fraction in same species
  Species1_Marker_table$cluster <- paste0(Species_name1, Species1_Marker_table[, 'cluster'])
  Species2_Marker_table$cluster <- paste0(Species_name2, Species2_Marker_table[, 'cluster'])
  ShMarker1 <- Cal_SharedMarkers(Species1_Marker_table, Species_name1)
  ShMarker2 <- Cal_SharedMarkers(Species2_Marker_table, Species_name2)
  Frac1 <- rbind(ShMarker1[[3]], ShMarker2[[3]])
  ShMarker3 <- rbind(ShMarker1[[1]],ShMarker2[[1]])
  ### fraction between different species
  Sp1Gene <- Species1_Marker_table$gene
  ExpInd01 <- t(apply(as.matrix(Sp1Gene), 1, function(x1){
    if(length(grep(x1, OrthG[,2]))==1){ x3 <- c(x1, grep(x1, OrthG[,2])) }
    else if(length(grep(x1, OrthG[,4]))==1){ x3 <- c(x1, grep(x1, OrthG[,4])) }else{ x3 <- c(x1, NA) }
  }) )
  Sp1Ind2 <- ExpInd01[!is.na(ExpInd01[, 2]), ]
  colnames(Sp1Ind2) <- c('ID', 'RowNum')
  Sp2Gene <- Species2_Marker_table$gene
  ExpInd02 <- t(apply(as.matrix(Sp2Gene), 1, function(x1){
    if(length(grep(x1, OrthG[,2]))==1){ x3 <- c(x1, grep(x1, OrthG[,2])) }
    else if(length(grep(x1, OrthG[,4]))==1){ x3 <- c(x1, grep(x1, OrthG[,4])) }else{ x3 <- c(x1, NA) }
  }) )
  Sp2Ind2 <- ExpInd02[!is.na(ExpInd02[, 2]), ]
  colnames(Sp2Ind2) <- c('ID', 'RowNum')
  ShMarker2 <- Cal_SharedMarkers_Species(Species1_Marker_table, Species2_Marker_table,
                                         Sp1Ind2, Sp2Ind2,
                                         c(Species_name1,Species_name2))
  ShMarker3 <- rbind(ShMarker3, ShMarker2[[1]])
  Frac1 <- rbind(Frac1, ShMarker2[[3]])

  AllCluster1 <- unique(c(Species1_Marker_table$cluster,Species2_Marker_table$cluster))
  AllCluster12 <- sort(AllCluster1); Frac3 <- c()
  for(i in 1:length(AllCluster1)){ Frac21 <- c()
  for(j in 1:length(AllCluster1)){
    Frac11 <- Frac1[Frac1[,1]==AllCluster1[i] & Frac1[,2]==AllCluster1[j], 3]
    Frac12 <- Frac1[Frac1[,1]==AllCluster1[j] & Frac1[,2]==AllCluster1[i], 3]
    if(length(Frac11)==1){ Frac2 <- as.numeric(Frac11)
    }else if(length(Frac12)==1){ Frac2 <- as.numeric(Frac12)
    }else{
      print(paste('No',AllCluster1[i],AllCluster1[j])) }
    Frac21 <- c(Frac21, Frac2)
  }
  Frac3 <- rbind(Frac3, Frac21)
  }
  rownames(Frac3) <- AllCluster1; colnames(Frac3) <- AllCluster1

  ShMarker4 <- list(); ShMarker4[[1]] <- ShMarker3
  Frac3[Frac3==1] <- NA
  ShMarker4[[2]] <- Frac3
  ShMarker4[[3]] <- Frac1
  return(ShMarker4)
}


Cal_SharedMarkers <- function(mmExp1, Spec1='mm'){

  mmCluster1 <- sort(unique(mmExp1$cluster))
  ShMarker1 <- c(); Fraction3 <- c(); Fraction4 <- c()
  for(i in 1:length(mmCluster1)){
    zfExp2 <- mmExp1[mmExp1$cluster==mmCluster1[i], ]
    zfPower1 <- sum(as.numeric(zfExp2$power))
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
    zfPower1 <- sum(as.numeric(zfExp2$power))
    zfExp22 <- zfExpInd1[match(zfExp2[,'gene'], zfExpInd1[,'ID']), ]

    Fraction2 <- c()
    for(j in 1:length(mmCluster1)){
      mmExp2 <- mmExp1[mmExp1$cluster==mmCluster1[j], ]
      mmPower1 <- sum(mmExp2$power)
      mmExp22 <- mmExpInd1[match(mmExp2[, 'gene'], mmExpInd1[,'ID']), ]
      if (class(mmExp22) == 'character') {
        mmExp22 <- matrix(mmExp22,ncol=2)
      }
      if (class(zfExp22) == 'character') {
        zfExp22 <- matrix(zfExp22,ncol=2)
      }
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
