
#' Identify conserved regulatory networks
#' @description Use Score of Conserved network to identify conserved regulatory
#' network modules based on homologous genes databased and topology of networks
#' @param OrthG ortholog genes database
#' @param Species1_GRN gene regulatory network of species 1
#' @param Species2_GRN gene regulatory network of species 2
#' @param Species_name1 character, indicating the species names of Species1_GRN
#' @param Species_name2 character, indicating the species names of Species2_GRN
#' @importFrom reshape2 dcast
#' @return list contains two df. First df contains details of conserved regulatory
#' network, second df contains NCS between module pairs
#' @export
#'
#' @examples load(system.file("extdata", "gene_network.rda", package = "CACIMAR"))
#' n1 <- Identify_ConservedNetworks(OrthG_Mm_Zf,mm_gene_network,zf_gene_network,'mm','zf')
Identify_ConservedNetworks <- function(OrthG,Species1_GRN,Species2_GRN,
                                       Species_name1,Species_name2){
  ### input check
  validInput(Species_name1,'Species_name1','character')
  validInput(Species_name2,'Species_name2','character')

  if (!'Source' %in% colnames(Species1_GRN)) {
    stop('Species1_GRN should contain Source column')
  }
  if (!'Target' %in% colnames(Species1_GRN)) {
    stop('Species1_GRN should contain Target column')
  }
  if (!'SourceGroup' %in% colnames(Species1_GRN)) {
    stop('Species1_GRN should contain SourceGroup column')
  }
  if (!'TargetGroup' %in% colnames(Species1_GRN)) {
    stop('Species1_GRN should contain TargetGroup column')
  }
  if (!'Source' %in% colnames(Species2_GRN)) {
    stop('Species2_GRN should contain Source column')
  }
  if (!'Target' %in% colnames(Species2_GRN)) {
    stop('Species2_GRN should contain Target column')
  }
  if (!'SourceGroup' %in% colnames(Species2_GRN)) {
    stop('Species2_GRN should contain SourceGroup column')
  }
  if (!'TargetGroup' %in% colnames(Species2_GRN)) {
    stop('Species2_GRN should contain TargetGroup column')
  }
  Species_name1 <- tolower(Species_name1)
  Species_name2 <- tolower(Species_name2)
  ### species check
  Spec1 <- colnames(OrthG)[2]
  Spec2 <- colnames(OrthG)[4]
  Spec1 <- gsub('_ID','',Spec1)
  Spec2 <- gsub('_ID','',Spec2)
  if (Spec1 == tolower(Species_name1) & Spec2 == tolower(Species_name2)) {
    Species_name <- c(Spec1,Spec2)
    RnList <- list(Species1_GRN,Species2_GRN)
  }else if(Spec2 == tolower(Species_name1) & Spec1 == tolower(Species_name2)){
    Species_name <- c(Spec1,Spec2)
    RnList <- list(Species2_GRN,Species1_GRN)
  }else{stop('please input correct Species name')}


  ### Extract genes in each group
  RnGeneList <- list()
  for (i in 1:length(RnList)) {
    Species <- Species_name[i]
    Rn <- RnList[[i]]
    RnTFGroup <- dplyr::select(Rn,c('Source','SourceGroup'))
    colnames(RnTFGroup) <- paste0(Species,c('Gene','Group'))
    RnTargetGroup <- dplyr::select(Rn,c('Target','TargetGroup'))
    colnames(RnTargetGroup) <- paste0(Species,c('Gene','Group'))
    GeneGroup <- rbind(RnTFGroup,RnTargetGroup)
    GeneGroup <- GeneGroup[!duplicated(GeneGroup),]
    RnGeneList[[Species]] <- GeneGroup
  }

  ### get orthology genes
  Sp1Gene <- RnGeneList[[Species_name[1]]]
  Sp2Gene <- RnGeneList[[Species_name[2]]]
  Sp1Gene <- Sp1Gene[!duplicated(Sp1Gene[,1]),]
  Sp2Gene <- Sp2Gene[!duplicated(Sp2Gene[,1]),]

  Spec1_gene <- data.frame(rep(0,nrow(Sp1Gene)),
                           rep(1,nrow(Sp1Gene)))
  rownames(Spec1_gene) <- Sp1Gene[,1]
  Spec2_gene <- data.frame(rep(0,nrow(Sp2Gene)),
                           rep(1,nrow(Sp2Gene)))
  rownames(Spec2_gene) <- Sp2Gene[,1]
  Exp2 <- Get_OrthG(OrthG, Spec1_gene, Spec2_gene, Species_name)
  if (str_sub(Sp1Gene$mmGene[1],1,3) == 'ENS') {
    Type1 <- paste0('Used_',Species_name[1],'_ID')
    Type2 <- paste0('Used_',Species_name[2],'_ID')
  }else{
    Type1 <- paste0('Used_',Species_name[1],'_Symbol')
    Type2 <- paste0('Used_',Species_name[2],'_Symbol')
  }

  Species1 <- Sp1Gene[match(Exp2[, Type1],Sp1Gene$mmGene), ]
  Species2 <- Sp2Gene[match(Exp2[, Type2],Sp2Gene$zfGene), ]
  Exp3 <- cbind(Exp2, Species1, Species2)
  Exp4 <- Exp3[!is.na(Exp3[, dim(Exp2)[2]+1]) &
                 !is.na(Exp3[,dim(Exp2)[2]+dim(Species1)[2]+1]), ]
  if (nrow(Exp4)==0) {
    NCS_df <- matrix(0,nrow = length(unique(Sp1Gene$SourceGroup)),
           ncol = length(unique(Sp2Gene$SourceGroup)))
    rownames(NCS_df) <- unique(Sp1Gene$SourceGroup)
    colnames(NCS_df) <- unique(Sp2Gene$SourceGroup)
    OrthG_list <- list(NA,NCS_df)
    return(OrthG_list)
  }

  ### calculate orthology fraction
  Sp1Freq <- as.data.frame(table(Sp1Gene$mmGroup))
  Sp2Freq <- as.data.frame(table(Sp2Gene$zfGroup))
  OrthG_df <- as.data.frame(table(Exp4$mmGroup,Exp4$zfGroup))
  OrthG_df <- OrthG_df[OrthG_df$Freq>0,]
  Sp1GeneNum <- Sp1Freq[match(OrthG_df$Var1,Sp1Freq$Var1),2]
  Sp2GeneNum <- Sp2Freq[match(OrthG_df$Var2,Sp2Freq$Var1),2]
  OrthG_df$Sp1GeneNum <- Sp1GeneNum
  OrthG_df$Sp2GeneNum <- Sp2GeneNum
  OrthGeneIdx <- match(paste0(Exp4[,9],Exp4[,11]),paste0(OrthG_df[,1],OrthG_df[,2]))
  Sp1OrthGenesList <- list()
  Sp2OrthGenesList <- list()
  for (i in 1:length(OrthGeneIdx)) {
    x1 <- OrthGeneIdx[i]
    Sp1OrthGenesList[[x1]] <- paste0(Sp1OrthGenesList[x1],';',Exp4[i,8])
    Sp2OrthGenesList[[x1]] <- paste0(Sp2OrthGenesList[x1],';',Exp4[i,10])
  }
  OrthG_df$Sp1OrthGenes <- gsub('NULL;','',unlist(Sp1OrthGenesList))
  OrthG_df$Sp2OrthGenes <- gsub('NULL;','',unlist(Sp2OrthGenesList))
  OrthG_df$SharedGenesFraction <- (OrthG_df$Freq*2)/(OrthG_df$Sp1GeneNum+OrthG_df$Sp2GeneNum)
  colnames(OrthG_df)[1:5] <- c('Group1','Group2','OrthGeneNum','Group1AllGenes','Group2AllGenes')

  ### calculate orthology edges fraction and topological orthology edges fraction
  OrthGEdgeSum <- c()
  AllEdgeNumSum <- c()
  OrthGEdgeFraction <- c()
  TopoOrthGEdgeSum <- c()
  OrthGTopoEdgeFraction <- c()
  for (i in 1:nrow(OrthG_df)) {
    Sp1GeneRow <- unlist(strsplit(OrthG_df$Sp1OrthGenes[i],';'))
    Sp2GeneRow <- unlist(strsplit(OrthG_df$Sp2OrthGenes[i],';'))

    NetworkGroup1 <- RnList[[1]][RnList[[1]]$SourceGroup==OrthG_df[i,1] &
                                   RnList[[1]]$TargetGroup==OrthG_df[i,1],]

    NetworkGroup2 <- RnList[[2]][RnList[[2]]$SourceGroup==OrthG_df[i,2] &
                                   RnList[[2]]$TargetGroup==OrthG_df[i,2],]

    Sp1OrthGEdge <- NetworkGroup1[NetworkGroup1$Source%in%Sp1GeneRow &
                                    NetworkGroup1$Target%in%Sp1GeneRow,]

    Sp2OrthGEdge <- NetworkGroup2[NetworkGroup2$Source%in%Sp2GeneRow &
                                    NetworkGroup2$Target%in%Sp2GeneRow,]

    OrthEdgeGNum <- nrow(Sp1OrthGEdge)+nrow(Sp2OrthGEdge)
    AllEdgeNum <- nrow(NetworkGroup1)+nrow(NetworkGroup2)

    OrthGEdgeSum <- c(OrthGEdgeSum,OrthEdgeGNum)
    AllEdgeNumSum <- c(AllEdgeNumSum,AllEdgeNum)
    OrthGEdgeFraction <- c(OrthGEdgeFraction,(OrthEdgeGNum/AllEdgeNum))
    ### topological fraction
    if (OrthG_df[i,3]>1) {

      GeneRowmatch <- data.frame(Sp1GeneRow,Sp2GeneRow)

      Sp1Regulation <- NetworkGroup1[NetworkGroup1$Source%in%Sp1GeneRow &
                                       NetworkGroup1$Target%in%Sp1GeneRow,]

      Sp2Regulation <- NetworkGroup2[NetworkGroup2$Source%in%Sp2GeneRow &
                                       NetworkGroup2$Target%in%Sp2GeneRow,]

      if (nrow(Sp1Regulation) > 0 & nrow(Sp2Regulation) > 0) {
        OrthGRelationshipsNum <- 0
        for (j in 1:nrow(Sp1Regulation)) {
          Sp2RelatedTF <- GeneRowmatch[GeneRowmatch$Sp1GeneRow == Sp1Regulation$Source[j],2]
          Sp2RelatedTarget <- GeneRowmatch[GeneRowmatch$Sp1GeneRow == Sp1Regulation$Target[j],2]
          if (paste(Sp2RelatedTF,Sp2RelatedTarget)%in%paste(Sp2Regulation$Source,Sp2Regulation$Target)) {
            OrthGRelationshipsNum <- OrthGRelationshipsNum+1
          }
        }
        Fraction <- (OrthGRelationshipsNum*2)/(nrow(Sp1Regulation)+nrow(Sp2Regulation))
        TopoOrthGEdgeSum <- c(TopoOrthGEdgeSum,OrthGRelationshipsNum*2)
        OrthGTopoEdgeFraction <- c(OrthGTopoEdgeFraction,Fraction)
      }else{
        TopoOrthGEdgeSum <- c(TopoOrthGEdgeSum,0)
        OrthGTopoEdgeFraction <- c(OrthGTopoEdgeFraction,0)
        }
    }else{
      TopoOrthGEdgeSum <- c(TopoOrthGEdgeSum,0)
      OrthGTopoEdgeFraction <- c(OrthGTopoEdgeFraction,0)
      }
  }
  OrthG_df$OrthGEdgeNum <- OrthGEdgeSum
  OrthG_df$EdgeNum <- AllEdgeNumSum
  OrthG_df$OrthGEdgeFraction <- OrthGEdgeFraction
  OrthG_df$TopologicalOrthGEdgeNum<-TopoOrthGEdgeSum
  OrthG_df$TopologicalOrthGEdgeFraction <- OrthGTopoEdgeFraction
  OrthG_df <- OrthG_df[,c(1:5,8:13,6,7)]

  ### calculate Score of Conserved Networks (SCN)
  Candidate_moudle <- OrthG_df[OrthG_df$OrthGEdgeNum > 0,]
  HNi <- mean(Candidate_moudle$OrthGeneNum)
  HNj <- mean(Candidate_moudle$OrthGeneNum)
  Ni <- mean(Candidate_moudle$Group1AllGenes)
  Nj <- mean(Candidate_moudle$Group2AllGenes)
  MIU <- ((Ni+Nj-1)/(HNi+HNj-1))
  #OrthG_df$NCS <- OrthG_df[,6] + (MIU*OrthG_df[,9]) + (MIU*2*OrthG_df[,11])
  OrthG_df$NCS <- OrthG_df[,6]+ (MIU*OrthG_df[,9])
  OrthG_df[,1] <- paste0(Species_name1,OrthG_df[,1])
  OrthG_df[,2] <- paste0(Species_name2,OrthG_df[,2])
  NCS_df <- reshape2::dcast(OrthG_df[,c(1,2,14)],Group1~Group2,fill = 0)
  rownames(NCS_df) <- NCS_df[,1]
  NCS_df <- NCS_df[,-1]
  OrthG_list <- list(OrthG_df,NCS_df)
  return(OrthG_list)
}


#' Identify conserved cell type specific regulatory networks
#'
#' @param OrthG ortholog genes database
#' @param Species1_GRN gene regulatory network of species 1
#' @param Species2_GRN gene regulatory network of species 2
#' @param Species_name1 character, indicating the species names of Species1_GRN
#' @param Species_name2 character, indicating the species names of Species2_GRN
#'
#' @return
#' @export
#'
#' @examples
identify_ct_ConservedNetworks <- function(OrthG,Species1_GRN,Species2_GRN,
                                          Species_name1,Species_name2){
  df_final <- list()
  for (i in unique(Species1_GRN$SourceGroup)) {
    for (j in unique(Species2_GRN$SourceGroup)) {

      spe1_network_use <- Species1_GRN[Species1_GRN$SourceGroup==i,]
      spe2_network_use <- Species2_GRN[Species2_GRN$SourceGroup==j,]
      rownames(spe1_network_use) <- 1:nrow(spe1_network_use)
      rownames(spe2_network_use) <- 1:nrow(spe2_network_use)
      n1 <- Identify_ConservedNetworks(OrthG_Mm_Zf,
                                       spe1_network_use,
                                       spe2_network_use,
                                       'mm','zf')
      if (!is.na(n1[[1]])) {
        homo_summary <- n1[[1]]
        n1[[1]]$pvalue <- chisq.test(data.frame(c(homo_summary[,3],homo_summary[,5]-homo_summary[,3]),
                              c(homo_summary[,3],homo_summary[,4]-homo_summary[,3])))$p.value
      }

      if (n1 == 'No conserved relationships') {

      }
      df_final[[paste0(i,'-',j)]] <- n1[[1]]
    }
  }

  df_final <- df_final[!is.na(df_final)]
  df_final <- do.call(bind_rows,df_final)
  ###NCS
  NCS_df <- reshape2::dcast(df_final[,c(1,2,14)],Group1~Group2,fill = 0)
  rownames(NCS_df) <- NCS_df[,1]
  NCS_df <- NCS_df[,-1]
  ###pvalue
  pv_df <- reshape2::dcast(df_final[,c(1,2,15)],Group1~Group2,fill = 0)
  rownames(pv_df) <- pv_df[,1]
  pv_df <- pv_df[,-1]
  OrthG_list <- list(df_final,NCS_df,pv_df)
  return(OrthG_list)
}




#' Identify conserved node and edge in the networks
#' Based on the Conservation analysis result from Identify_ConservedNetworks,
#' this function further identify the conserved node and edge in the network(graph)
#' @param Species1_GRN gene regulatory network of species 1
#' @param Species2_GRN gene regulatory network of species 2
#' @param ConservedNetworkTable result from Identify_ConservedNetworks
#' @param Species1_group character,indicating interested GRN group in species 1
#' @param Species2_group character,indicating interested GRN group in species 2
#'
#' @return
#' @export
#'
#' @examples
identify_network_relationships <- function(Species1_GRN,Species2_GRN,
                                           ConservedNetworkTable,
                                           Species1_group = '',
                                           Species2_group = ''){
  mm_con_gene = ConservedNetworkTable[ConservedNetworkTable$Group1==Species1_group &
                                        ConservedNetworkTable$Group2==Species2_group,12]
  zf_con_gene = ConservedNetworkTable[ConservedNetworkTable$Group1==Species1_group &
                                        ConservedNetworkTable$Group2==Species2_group,13]
  mm_con_gene = unlist(strsplit(mm_con_gene,';'))
  zf_con_gene = unlist(strsplit(zf_con_gene,';'))
  mm_source = c()
  mm_target = c()
  for (i in 1:nrow(Species1_GRN)) {
    if (Species1_GRN$Source[i] %in% mm_con_gene) {
      mm_source = c(mm_source,'conserved')
    }else{
      mm_source = c(mm_source,'unconserved')
    }
    if (Species1_GRN$Target[i] %in% mm_con_gene) {
      mm_target = c(mm_target,'conserved')
    }else{
      mm_target = c(mm_target,'unconserved')
    }
  }

  zf_source = c()
  zf_target = c()
  for (i in 1:nrow(Species2_GRN)) {
    if (Species2_GRN$Source[i] %in% zf_con_gene) {
      zf_source = c(zf_source,'conserved')
    }else{
      zf_source = c(zf_source,'unconserved')
    }
    if (Species2_GRN$Target[i] %in% zf_con_gene) {
      zf_target = c(zf_target,'conserved')
    }else{
      zf_target = c(zf_target,'unconserved')
    }
  }
  Species1_GRN$spec1_source_type = mm_source
  Species1_GRN$spec1_target_type = mm_target
  Species2_GRN$spec2_source_type = zf_source
  Species2_GRN$spec2_target_type = zf_target

  Species1_GRN$index = paste0(Species1_GRN$spec1_source_type,
                              Species1_GRN$spec1_target_type)
  Species1_GRN$edge_type = ifelse(Species1_GRN$index=='conservedconserved'
                                  ,'conserved','unconserved')
  Species2_GRN$index = paste0(Species2_GRN$spec2_source_type,
                              Species2_GRN$spec2_target_type)
  Species2_GRN$edge_type = ifelse(Species2_GRN$index=='conservedconserved'
                                  ,'conserved','unconserved')
  return(list(Species1_GRN,Species2_GRN))
}
