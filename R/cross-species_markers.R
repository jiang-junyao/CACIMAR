#' Identify orthologs marker genes for two species
#' @description Identify orthologs marker genes for two species based on orthologs database
#' @param OrthG ortholog genes database
#' @param Species1_Marker_table data.frame of species 1, should contain 'gene' and
#' 'Allcluster' columns.
#' @param Species2_Marker_table data.frame of species 2, should contain 'gene' and
#' 'Allcluster' columns.
#' @param Species_name1 character, indicating the species names of Species1_Marker_table.
#' @param Species_name2 character, indicating the species names of Species2_Marker_table
#' @param match_cell_name characters contained in both cell names
#'  to match similar cell types
#' @param filter_marker logical, indicating whether filter markers
#'
#' @return Data frame of conserved markers
#' @export
#'
#' @examples load(system.file("extdata", "zf_mm_markers.rda", package = "CACIMAR"))
#' ConservedMarker <- Identify_ConservedMarkers(OrthG_Mm_Zf,Mm_marker,Zf_marker,
#' Species_name1 = 'mm',Species_name2 = 'zf')
Identify_ConservedMarkers <- function(OrthG,Species1_Marker_table,Species2_Marker_table,
                           Species_name1,Species_name2,
                           match_cell_name=NULL,filter_marker =TRUE){
  validInput(OrthG,'OrthG','df')
  validInput(Species_name1,'Species_name1','character')
  validInput(Species_name2,'Species_name2','character')
  if (!'gene' %in% colnames(Species1_Marker_table)) {
    stop('Species1_Marker_table should contain gene column')
  }
  if (!'gene' %in% colnames(Species2_Marker_table)) {
    stop('Species2_Marker_table should contain gene column')
  }
  if (!'Allcluster' %in% colnames(Species1_Marker_table)) {
    stop('Species1_Marker_table should contain Allcluster column')
  }
  if (!'Allcluster' %in% colnames(Species2_Marker_table)) {
    stop('Species2_Marker_table should contain Allcluster column')
  }
  geneidx1 <- grep('gene',colnames(Species1_Marker_table))
  geneidx2 <- grep('gene',colnames(Species2_Marker_table))
  Species1_Marker_table <- Species1_Marker_table[!duplicated(Species1_Marker_table[,geneidx1]),]
  Species2_Marker_table <- Species2_Marker_table[!duplicated(Species2_Marker_table[,geneidx2]),]
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
  colnames(Species1_Marker) <- paste0(Species_name[1],colnames(Species1_Marker))
  colnames(Species2_Marker) <- paste0(Species_name[2],colnames(Species2_Marker))
  Species12 <- Species1_Marker
  Species22 <- Species2_Marker
  Spec1_gene <- data.frame(rep(0,nrow(Species12)),
                           rep(1,nrow(Species12)))
  rownames(Spec1_gene) <- Species12[,geneidx1]
  Spec2_gene <- data.frame(rep(0,nrow(Species22)),
                           rep(1,nrow(Species22)))
  rownames(Spec2_gene) <- Species22[,geneidx2]
  Exp2 <- Get_OrthG(OrthG, Spec1_gene, Spec2_gene, Species_name)
  if (grepl('ENS',rownames(Spec1_gene)[1])) {
    Type1 <- paste0('Used_',Species_name[1],'_ID')
    Type2 <- paste0('Used_',Species_name[2],'_ID')
  }else{
    Type1 <- paste0('Used_',Species_name[1],'_Symbol')
    Type2 <- paste0('Used_',Species_name[2],'_Symbol')
  }
  Species1 <- Species1_Marker[match(Exp2[, Type1],Species1_Marker[,paste0(Species_name[1],'gene')]), ]
  Species2 <- Species2_Marker[match(Exp2[, Type2],Species2_Marker[,paste0(Species_name[2],'gene')]), ]
  Exp3 <- cbind(Exp2, Species1, Species2)
  Exp4 <- Exp3[!is.na(Exp3[, dim(Exp2)[2]+1]) & !is.na(Exp3[,dim(Exp2)[2]+dim(Species1)[2]+1]), ]
  if (nrow(Exp4)==0) {
    stop('No homologous genes appear!')
  }
  if (filter_marker) {
      Exp5 <- cbind(Exp4[,1:7], Species12[match(Exp4[,Type1],
                                            Species12[,paste0(Species_name[1],'gene')]), ],
                Species22[match(Exp4[,Type2], Species22[,paste0(Species_name[2],'gene')]), ])
      Exp6 <- Refine_Used_OrthG(Exp5,Species_name,match_cell_name)
      return(Exp6)
  }else{
    return(Exp4)
  }

}

#' Identify conserved markers in conserved celltype
#'
#' Identify orthologs marker genes for two species
#' @description Identify orthologs marker genes for two species based on orthologs database
#' @param OrthG ortholog genes database
#' @param Species1_Marker_table data.frame of species 1, should contain 'gene' and
#' 'Allcluster' columns.
#' @param Species2_Marker_table data.frame of species 2, should contain 'gene' and
#' 'Allcluster' columns.
#' @param Species_name1 character, indicating the species names of Species1_Marker_table.
#' @param Species_name2 character, indicating the species names of Species2_Marker_table
#' @param conserved_celltype_pair character, indicating the conserved celltypes
#'
#' @return
#' @export
#'
#' @examples
identify_conserved_marker <- function(OrthG,Species1_Marker_table,
                                      Species2_Marker_table,
                                      Species_name1,Species_name2,
                                      conserved_celltype_pair){
  ### check species name
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

  conserved_list = list()
  for (i in conserved_celltype_pair) {
    spc1_cluster = unlist(strsplit(i,'-'))[1]
    spc1_cluster = gsub(Species_name1,'',spc1_cluster)
    spc2_cluster = unlist(strsplit(i,'-'))[2]
    spc2_cluster = gsub(Species_name2,'',spc2_cluster)
    spc1_marker = Species1_Marker_table[grep(spc1_cluster,Species1_Marker_table$Allcluster),]
    spc2_marker = Species2_Marker_table[grep(spc2_cluster,Species2_Marker_table$Allcluster),]
    conserved_gene = identify_conserved_gene(OrthG,spc1_marker$gene,
                                             spc2_marker$gene,
                                             Species_name1,Species_name2)
    if (is.matrix(conserved_gene)) {
      conserved_gene = as.data.frame(conserved_gene[,6:7])
    }else{
      conserved_gene = as.data.frame(t(as.data.frame(conserved_gene[6:7])))
    }

    conserved_gene$spc1_cluster = spc1_cluster
    conserved_gene$spc2_cluster = spc2_cluster
    conserved_list[[i]] = conserved_gene
  }
  conserved_df = do.call(bind_rows,conserved_list)
  colnames(conserved_df) = c(paste0(Species_name1,'_conserved_marker'),
                             paste0(Species_name2,'_conserved_marker'),
                             paste0(Species_name1,'_cluster'),
                             paste0(Species_name2,'_cluster'))
  return(conserved_df)
}

#' Identify conserved gene
#'
#' @param OrthG ortholog genes database
#' @param spc1_marker vector, indicating the gene of species 1
#' @param spc2_marker vector, indicating the gene of species 2
#' @param Species_name1 character, indicating the species names of Species1_Marker_table.
#' @param Species_name2 character, indicating the species names of Species2_Marker_table
#'
#' @return
#' @export
#'
#' @examples
identify_conserved_gene <- function(OrthG,spc1_marker,spc2_marker,Species_name1,
                                    Species_name2){
  ### check species name
  Species_name1 <- tolower(Species_name1)
  Species_name2 <- tolower(Species_name2)
  Spec1 <- colnames(OrthG)[2]
  Spec2 <- colnames(OrthG)[4]
  Spec1 <- gsub('_ID','',Spec1)
  Spec2 <- gsub('_ID','',Spec2)
  if (Spec1 == Species_name1 & Spec2 == Species_name2) {
    Species_name <- c(Species_name1,Species_name2)
    Species1_Marker <- spc1_marker
    Species2_Marker <- spc2_marker
  }else if(Spec2 == Species_name1 & Spec1 == Species_name2){
    Species_name <- c(Species_name2,Species_name1)
    Species2_Marker <- spc1_marker
    Species1_Marker <- spc2_marker
  }else{stop('please input correct Species name')}
  spc1_marker = Species1_Marker
  spc2_marker = Species2_Marker
  Species_name1 = Species_name[1]
  Species_name2 = Species_name[2]

  Species_name1 <- tolower(Species_name1)
  Species_name2 <- tolower(Species_name2)
  Spec1 <- colnames(OrthG)[2]
  Spec2 <- colnames(OrthG)[4]
  Spec1 <- gsub('_ID','',Spec1)
  Spec2 <- gsub('_ID','',Spec2)
  Spec1_gene <- data.frame(rep(0,length(spc1_marker)),
                           rep(1,length(spc1_marker)))
  rownames(Spec1_gene) <- spc1_marker
  Spec2_gene <- data.frame(rep(0,length(spc2_marker)),
                           rep(1,length(spc2_marker)))
  rownames(Spec2_gene) <- spc2_marker
  Species_name = c(Species_name1,Species_name2)
  Exp2 <- Get_OrthG(OrthG, Spec1_gene, Spec2_gene, Species_name)
  if (grepl('ENS',rownames(Spec1_gene)[1])) {
    Type1 <- paste0('Used_',Species_name[1],'_ID')
    Type2 <- paste0('Used_',Species_name[2],'_ID')
  }else{
    Type1 <- paste0('Used_',Species_name[1],'_Symbol')
    Type2 <- paste0('Used_',Species_name[2],'_Symbol')
  }
  Exp2 = Exp2[Exp2[,6] %in% spc1_marker,]
  Exp2 = Exp2[Exp2[,7] %in% spc2_marker,]
  if (!is_empty(Exp2)) {
      return(Exp2)
  }else{
    return('no marker')
  }

}



Get_OrthG <- function(OrthG1, MmRNA1, ZfRNA1, Spec1, MmPattern1='', ZfPattern1=''){
  tOrthG1 <- table(OrthG1$Type)
  if (grepl('ENS',rownames(MmRNA1)[1])) {
    Ind1 <- c(grep(paste0(Spec1[1],'_ID'), colnames(OrthG1)), grep(paste0(Spec1[2],'_ID'), colnames(OrthG1)))
  }else{
    Ind1 <- c(grep(paste0(Spec1[1],'_Symbol'), colnames(OrthG1)), grep(paste0(Spec1[2],'_Symbol'), colnames(OrthG1)))
  }


  OrthG21 <- list()
  for(i in 1:length(tOrthG1)){
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
    }else{ print(paste0('Not process ',tOrthG1[i]))
     stop('Can not match any  orthologs marker genes, please whether input
           correct species names ') }
    if(i==1){ OrthG4 <- as.matrix(OrthG3);
    }else{ OrthG4 <- rbind(OrthG4, as.matrix(OrthG3)) }
  }
  colnames(OrthG4)[(ncol(OrthG1)+1):(ncol(OrthG1)+2)] <- paste0('Used_',colnames(OrthG4)[(ncol(OrthG1)+1):(ncol(OrthG1)+2)])

  return(OrthG4)
}


Refine_Used_OrthG<-function(ShMarker1,Species,smiliar_cell_name){
  Type1 <- paste0(Species[1],'.*\\Allcluster'); Type2 <- paste0(Species[2], '.*\\Allcluster')
  Spec1Type1 <- grep(Type1, colnames(ShMarker1))
  Spec1Type2 <- grep(Type2, colnames(ShMarker1))
  if (is.null(smiliar_cell_name)) {
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
                  if(x112[i1]==x122[j1]){
                    x2 <- T
                  } } } } } } }
      return(x2)
    })
  }else{
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
                if(x112[i1]==x122[j1] | grepl(smiliar_cell_name,x112[i1]) & grepl(smiliar_cell_name,x122[j1])){
                  x2 <- T
                } } } } } } }
    return(x2)
  })}
  ShMarker3 <- ShMarker1[ShMarker2, ]
  ShMarker4 <- ShMarker3[order(ShMarker3[, Spec1Type1[1]]), ]
  return(ShMarker4)
}



Refine_Markers_Species<-function(Marker1,Species){
  PowerTh1 <- 0.4; PowerTh2 <- gsub('\\.','',PowerTh1)
  SpecInd1 <- length(Species);
  Ind1 <- list(); Ind2 <- list();
  for(i in 1:length(Species)){
    Ind1[[i]] <- grep(paste0(Species[i],'.*\\.luster'), colnames(Marker1))
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

orderCellType <- function(df){
  ct <- colnames(df)[-1]
  idx <- c()
  for (i in ct) {
    idx1 <- grep(i,df[,1])
    idx <- c(idx,idx1)
  }
  df1 <- df[idx,]
  return(df1)
}
