
#' Identify conserved regulatory networks
#' @description Use Score of Conserved network to identify conserved regulatory
#' network modules based on homologous genes databased and topology of networks
#' @param OrthG ortholog genes database
#' @param Species1_GRN gene regulatory network of species 1
#' @param Species2_GRN gene regulatory network of species 2
#' @param Species_name1 character, indicating the species names of Species1_GRN
#' @param Species_name2 character, indicating the species names of Species2_GRN
#'
#' @return list contains two df. First df contains details of conserved regulatory
#' network, second df contains NCS between module pairs
#' @export
#'
#' @examples
Identify_ConservedNetworks <- function(OrthG,Species1_GRN,Species2_GRN,Species_name1,Species_name2){
  ### input check
  validInput(Species_name1,'Species_name1','character')
  validInput(Species_name2,'Species_name2','character')
  validInput(k1,'k1','numeric')
  validInput(k1,'k2','numeric')
  validInput(k1,'k3','numeric')

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
    Species_name <- c(Species2_GRN,Species1_GRN)
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
  Sp1Gene <- RnGeneList[[1]]
  Sp2Gene <- RnGeneList[[2]]
  Sp1Group <- sort(Sp1Gene$Group[!duplicated(Sp1Gene$Group)])
  Sp2Group <- sort(Sp2Gene$Group[!duplicated(Sp2Gene$Group)])
  Spec1_gene <- data.frame(rep(0,nrow(Sp1Gene)),
                           rep(1,nrow(Sp1Gene)))
  rownames(Spec1_gene) <- Sp1Gene$mmGene
  Spec2_gene <- data.frame(rep(0,nrow(Sp2Gene)),
                           rep(1,nrow(Sp2Gene)))
  rownames(Spec2_gene) <- Sp2Gene$zfGene
  Exp2 <- Get_OrthG(OrthG, Spec1_gene, Spec2_gene, Species_name)
  Type1 <- paste0('Used_',Species_name[1],'_ID')
  Type2 <- paste0('Used_',Species_name[2],'_ID')
  Species1 <- Sp1Gene[match(Exp2[, Type1],Sp1Gene$mmGene), ]
  Species2 <- Sp2Gene[match(Exp2[, Type2],Sp2Gene$zfGene), ]
  Exp3 <- cbind(Exp2, Species1, Species2)
  Exp4 <- Exp3[!is.na(Exp3[, dim(Exp2)[2]+1]) &
                 !is.na(Exp3[,dim(Exp2)[2]+dim(Species1)[2]+1]), ]

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

    Sp2OrthGEdge <- NetworkGroup2[NetworkGroup2$Source%in%Sp1GeneRow &
                                    NetworkGroup2$Target%in%Sp1GeneRow,]

    OrthEdgeGNum <- nrow(Sp1OrthGEdge)+nrow(Sp2OrthGEdge)
    AllEdgeNum <- nrow(NetworkGroup1)+nrow(NetworkGroup2)

    OrthGEdgeSum <- c(OrthGEdgeSum,OrthEdgeGNum)
    AllEdgeNumSum <- c(AllEdgeNumSum,AllEdgeNum)
    OrthGEdgeFraction <- c(OrthGEdgeFraction,(OrthEdgeGNum/AllEdgeNum))
    ### topological fraction
    if (OrthG_df[,3]>1) {

      GeneRowmatch <- data.frame(Sp1GeneRow,Sp2GeneRow)

      Sp1Regulation <- NetworkGroup1[NetworkGroup1$Source%in%Sp1GeneRow &
                                       NetworkGroup1$Target%in%Sp1GeneRow,]

      Sp2Regulation <- NetworkGroup2[NetworkGroup2$Source%in%Sp1GeneRow &
                                       NetworkGroup2$Target%in%Sp1GeneRow,]

      if (nrow(Sp1Regulation) > 0 & nrow(Sp2Regulation) > 0) {
        OrthGRelationshipsNum <- 0
        for (j in 1:nrow(Sp1Regulation)) {
          Sp2RelatedTF <- GeneRowmatch[GeneRowmatch$Sp1GeneRow == Sp1Regulation$Source,2]
          Sp2RelatedTarget <- GeneRowmatch[GeneRowmatch$Sp1GeneRow == Sp1Regulation$Target,2]
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
  OrthG_df$NCS <- OrthG_df[,6] + (MIU*OrthG_df[,9]) + (MIU*2*OrthG_df[,11])
  OrthG_df[,1] <- paste0(Species_name1,OrthG_df[,1])
  OrthG_df[,2] <- paste0(Species_name2,OrthG_df[,2])
  NCS_df <- reshape2::dcast(OrthG_df[,c(1,2,14)],Group1~Group2,fill = 0)
  rownames(NCS_df) <- NCS_df[,1]
  NCS_df <- NCS_df[,-1]
  OrthG_list <- list(OrthG_df,NCS_df)
  return(OrthG_list)
}






