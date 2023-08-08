
#' Identify conserved cell-cell interactions in conserved cell types
#' @description Identify conserved cell-cell interactions in conserved cell types
#' @param OrthG ortholog genes database
#' @param Species1_CCI cell-cell interactions in species1. First column should
#' be ligand, second column should be target, third column should be corresponding
#' cell type of ligand, fourth column should be corresponding cel type of target
#' @param Species2_CCI cell-cell interactions in species2. First column should
#' be ligand, second column should be target, third column should be corresponding
#' cell type of ligand, fourth column should be corresponding cel type of target
#' @param ConservedCellType
#' @param Species_name1 character, indicating the species names of Species1_GRN
#' @param Species_name2 character, indicating the species names of Species2_GRN
#' @importFrom dplyr bind_rows
#' @importFrom  stringr str_sub
#' @return
#' @export
#'
#' @examples load(system.file("extdata", "cci_test.rda", package = "CACIMAR"))
#' Identify_ConservedCCI(OrthG_Hs_Mm,hs_cci_test,mm_cci_test,celltype,'hs','mm')
Identify_ConservedCCI <- function(OrthG,
                                  Species1_CCI,
                                  Species2_CCI,
                                  ConservedCellType=NULL,
                                  Species_name1,
                                  Species_name2){
  colnames(Species1_CCI) = paste0('sp1',colnames(Species1_CCI))
  colnames(Species2_CCI) = paste0('sp2',colnames(Species2_CCI))
  ### input check
  validInput(Species_name1,'Species_name1','character')
  validInput(Species_name2,'Species_name2','character')
  ### species name check
  Spec1 <- colnames(OrthG)[2]
  Spec2 <- colnames(OrthG)[4]
  Spec1 <- gsub('_ID','',Spec1)
  Spec2 <- gsub('_ID','',Spec2)
  if (Spec1 == tolower(Species_name1) & Spec2 == tolower(Species_name2)) {
    Species1_CCI = Species1_CCI
    Species2_CCI = Species2_CCI
    Species_name <- c(Species_name1,Species_name2)
  }else if(Spec2 == tolower(Species_name1) & Spec1 == tolower(Species_name2)){
    sp1 = Species1_CCI
    sp2 = Species2_CCI
    Species1_CCI = sp2
    Species2_CCI = sp1
    Species_name <- c(Species_name2,Species_name1)

  }else{stop('please input correct Species name')}
  ### check gene name format
  if (stringr::str_sub(Species1_CCI[1,1],1,3)=='ENS') {
    gene_name_type = 'ENS'
  }else{gene_name_type = 'Symbol'}
  ### set orthg idx
  if (gene_name_type=='ENS') {
    db_idx = c(2,4)
  }else if(gene_name_type=='Symbol'){
    db_idx = c(3,5)
  }
  ### retina same clusters
  sp1_all_clusters = c(Species1_CCI[,3],Species1_CCI[,4])
  sp1_all_clusters = sp1_all_clusters[!duplicated(sp1_all_clusters)]
  sp2_all_clusters = c(Species2_CCI[,3],Species2_CCI[,4])
  sp2_all_clusters = sp2_all_clusters[!duplicated(sp2_all_clusters)]
  clusters_use = intersect(sp1_all_clusters,sp2_all_clusters)
  Species1_CCI = Species1_CCI[Species1_CCI[,3]%in%clusters_use &
                                Species1_CCI[,4]%in%clusters_use,]
  Species2_CCI = Species2_CCI[Species2_CCI[,3]%in%clusters_use &
                                Species2_CCI[,4]%in%clusters_use,]

  ### retina same cluster pairs
  sp1_cluster_pair = paste0(Species1_CCI[,3],Species1_CCI[,4])
  sp1_cluster_pair_un = sp1_cluster_pair[!duplicated(sp1_cluster_pair )]
  sp2_cluster_pair = paste0(Species2_CCI[,3],Species2_CCI[,4])
  sp2_cluster_pair_un = sp2_cluster_pair[!duplicated(sp2_cluster_pair )]
  overlap_pairs = intersect(sp1_cluster_pair_un,sp2_cluster_pair_un)
  Species1_CCI = Species1_CCI[sp1_cluster_pair %in% overlap_pairs,]
  Species2_CCI = Species2_CCI[sp2_cluster_pair %in% overlap_pairs,]

  ### select cell types may contain interactions
  if (!is.null(ConservedCellType)) {
    sp1_ct_use = ConservedCellType[,1]
    sp2_ct_use = ConservedCellType[,2]
    ### filter CCI celltype
    Species1_CCI = Species1_CCI[Species1_CCI[,3] %in% sp1_ct_use,]
    Species1_CCI = Species1_CCI[Species1_CCI[,4] %in% sp1_ct_use,]
    Species2_CCI = Species2_CCI[Species2_CCI[,3] %in% sp2_ct_use,]
    Species2_CCI = Species2_CCI[Species2_CCI[,4] %in% sp2_ct_use,]
    ### initial CCI filtering
    Species1_CCI = Species1_CCI[Species1_CCI[,1] %in%
                                  unlist(strsplit(OrthG[,db_idx[1]],';')),]
    Species1_CCI = Species1_CCI[Species1_CCI[,2] %in%
                                  unlist(strsplit(OrthG[,db_idx[1]],';')),]
    Species2_CCI = Species2_CCI[Species2_CCI[,1] %in%
                                  unlist(strsplit(OrthG[,db_idx[2]],';')),]
    Species2_CCI = Species2_CCI[Species2_CCI[,2] %in%
                                  unlist(strsplit(OrthG[,db_idx[2]],';')),]
  }
  ### filter unconserved ligand
  Sp1Gene <- Species1_CCI[,1]
  Sp2Gene <- Species2_CCI[,1]
  Sp1Gene <- Sp1Gene[!duplicated(Sp1Gene)]
  Sp2Gene <- Sp2Gene[!duplicated(Sp2Gene)]
  Spec1_gene <- data.frame(rep(0,length(Sp1Gene)),
                           rep(1,length(Sp1Gene)))
  rownames(Spec1_gene) <- Sp1Gene
  Spec2_gene <- data.frame(rep(0,length(Sp2Gene)),
                           rep(1,length(Sp2Gene)))
  rownames(Spec2_gene) <- Sp2Gene
  Exp2 <- Get_OrthG(OrthG, Spec1_gene, Spec2_gene, Species_name)
  sp1idx = grep(paste0('Used_',Species_name[1]),colnames(Exp2))
  sp2idx = grep(paste0('Used_',Species_name[2]),colnames(Exp2))
  Exp3 = Exp2[!is.na(Exp2[,sp1idx])&
                !is.na(Exp2[,sp2idx]),]
  sp1_lr_used = Sp1Gene[Sp1Gene %in% Exp3[,sp1idx]]
  sp2_lr_used = Sp2Gene[Sp2Gene %in% Exp3[,sp2idx]]
  Species1_CCI = Species1_CCI[Species1_CCI[,1]%in% sp1_lr_used,]
  Species2_CCI = Species2_CCI[Species2_CCI[,1]%in% sp2_lr_used,]

  ### filter unconserved receptor
  Sp1Gene <- Species1_CCI[,2]
  Sp2Gene <- Species2_CCI[,2]
  Sp1Gene <- Sp1Gene[!duplicated(Sp1Gene)]
  Sp2Gene <- Sp2Gene[!duplicated(Sp2Gene)]
  Spec1_gene <- data.frame(rep(0,length(Sp1Gene)),
                           rep(1,length(Sp1Gene)))
  rownames(Spec1_gene) <- Sp1Gene
  Spec2_gene <- data.frame(rep(0,length(Sp2Gene)),
                           rep(1,length(Sp2Gene)))
  rownames(Spec2_gene) <- Sp2Gene
  Exp2 <- Get_OrthG(OrthG, Spec1_gene, Spec2_gene, Species_name)
  sp1idx = grep(paste0('Used_',Species_name[1]),colnames(Exp2))
  sp2idx = grep(paste0('Used_',Species_name[2]),colnames(Exp2))
  Exp3 = Exp2[!is.na(Exp2[,sp1idx])&
                !is.na(Exp2[,sp2idx]),]
  sp1_lr_used = Sp1Gene[Sp1Gene %in% Exp3[,sp1idx]]
  sp2_lr_used = Sp2Gene[Sp2Gene %in% Exp3[,sp2idx]]
  Species1_CCI = Species1_CCI[Species1_CCI[,2]%in% sp1_lr_used,]
  Species2_CCI = Species2_CCI[Species2_CCI[,2]%in% sp2_lr_used,]


  ### identify conserved ligand receptor pairs
  sp1_lr_pair = paste0(Species1_CCI[,1],'-',Species1_CCI[,2])
  sp1_lr_pair_un = sp1_lr_pair[!duplicated(sp1_lr_pair)]
  sp2_lr_pair = paste0(Species2_CCI[,1],'-',Species2_CCI[,2])
  sp2_lr_pair_un = sp2_lr_pair[!duplicated(sp2_lr_pair)]
  for (i in sp1_lr_pair_un) {
    print(Species1_CCI[sp1_lr_pair==i,c(3,4)])

  }
  sp1_orthg_cci = c()
  sp2_orthg_cci = c()
  for (i in sp1_lr_pair_un) {
    cci1 = unlist(strsplit(i,'-'))
    for (j in sp2_lr_pair_un) {
      cci2 = unlist(strsplit(j,'-'))
      orthg_check = filter_orthg(cci1,cci2,OrthG,db_idx)
      if (orthg_check) {
        sp1_orthg_cci = c(sp1_orthg_cci,i)
        sp2_orthg_cci = c(sp2_orthg_cci,j)
      }
    }
  }
  Species1_CCI = Species1_CCI[sp1_lr_pair %in% sp1_orthg_cci,]
  Species2_CCI = Species2_CCI[sp2_lr_pair %in% sp2_orthg_cci,]
  sp1_lr_pair = paste0(Species1_CCI[,1],'-',Species1_CCI[,2])
  sp2_lr_pair = paste0(Species2_CCI[,1],'-',Species2_CCI[,2])
  sp1_final_list = list()
  sp2_final_list = list()
  for (i in 1:length(sp1_orthg_cci)) {
    sp1_cci_use = Species1_CCI[sp1_lr_pair%in%sp1_orthg_cci[i],]
    sp2_cci_use = Species2_CCI[sp2_lr_pair%in%sp2_orthg_cci[i],]
    sp1_cci_use_cluster = paste(sp1_cci_use[,3],sp1_cci_use[,4])
    sp2_cci_use_cluster = paste(sp2_cci_use[,3],sp2_cci_use[,4])
    intersect_cluster = intersect(sp1_cci_use_cluster,sp2_cci_use_cluster)
    print(intersect_cluster)
    sp1_final_list[[i]]=sp1_cci_use[sp1_cci_use_cluster %in% intersect_cluster,]
    sp2_final_list[[i]]=sp2_cci_use[sp2_cci_use_cluster %in% intersect_cluster,]
  }
  sp1_final_df = do.call(dplyr::bind_rows,sp1_final_list)
  sp2_final_df = do.call(dplyr::bind_rows,sp2_final_list)
  final_list = list(sp1_final_df,sp2_final_df)
  names(final_list) = c('sp1_orthg_ccc_df','sp2_orthg_ccc_df')
  return(final_list)
}


filter_orthg <- function(cci1,cci2,OrthG,db_idx){
  source = FALSE
  target = FALSE
  ### identify 1T1
  idx = grep('*1T1',OrthG$Type)
  db = OrthG[idx,]
  orthg_genes1 = db[db[,db_idx[1]] %in% cci1[1],db_idx[2]]
  orthg_genes2 = db[db[,db_idx[1]] %in% cci1[2],db_idx[2]]
  if (cci2[1] %in% orthg_genes1) {
    source = TRUE
  }
  if (cci2[2] %in% orthg_genes2) {
    target = TRUE
  }
  ### identify 1TN
  idx = grep('*1TN',OrthG$Type)
  db = OrthG[idx,]
  orthg_genes1 = db[db[,db_idx[1]] %in% cci1[1],db_idx[2]]
  orthg_genes1 = unlist(strsplit(orthg_genes1,','))
  orthg_genes2 = db[db[,db_idx[1]] %in% cci1[2],db_idx[2]]
  orthg_genes2 = unlist(strsplit(orthg_genes2,','))
  if (cci2[1] %in% orthg_genes1) {
    source = TRUE
  }
  if (cci2[2] %in% orthg_genes2) {
    target = TRUE
  }
  ### identify NT1
  idx = grep('*NT1',OrthG$Type)
  db = OrthG[idx,]
  check_source = unlist(apply(db, 1, check_orthg_NT1,
                              gene1=cci1[1],gene2=cci2[1],db_idx=db_idx))
  check_target = unlist(apply(db, 1, check_orthg_NT1,
                              gene1=cci1[2],gene2=cci2[2],db_idx=db_idx))
  if (1 %in% check_source) {
    source = TRUE
  }
  if (1 %in% check_target) {
    target = TRUE
  }
  ### identify NTN
  idx = grep('*NTN',OrthG$Type)
  db = OrthG[idx,]
  check_source = unlist(apply(db, 1, check_orthg_NTN,
                              gene1=cci1[1],gene2=cci2[1],db_idx=db_idx))
  check_target = unlist(apply(db, 1, check_orthg_NTN,
                              gene1=cci1[2],gene2=cci2[2],db_idx=db_idx))
  if (1 %in% check_source) {
    source = TRUE
  }
  if (1 %in% check_target) {
    target = TRUE
  }
  if (source & target) {
    return(TRUE)
  }else{return(FALSE)}
}


check_orthg_NT1 <- function(db,gene1,gene2,db_idx){
  querry = unlist(strsplit(db[db_idx[1]],';'))
  if (gene1 %in% querry) {
    if (gene2 == db[db_idx[2]]) {
      return(1)
    }
  }
}

check_orthg_NTN <- function(db,gene1,gene2,db_idx){
  querry = unlist(strsplit(db[db_idx[1]],';'))
  target = unlist(strsplit(db[db_idx[2]],';'))
  if (gene1 %in% querry) {
    if (gene2 %in% target) {
      return(1)
    }
  }
}


#' filter cellchat CCC networks
#' @description Use conserved CCC to filter cellchat CCC networks
#' @param conserved_ccc_df
#' @param cellchat_df
#'
#' @return
#' @export
#'
#' @examples
add_cellchat_prob = function(conserved_ccc_df,cellchat_df){
  ta_idx <- paste0(cellchat_df[,3],cellchat_df[,4],
                   cellchat_df[,1],cellchat_df[,2])
  qu_idx <- paste0(conserved_ccc_df[,1],conserved_ccc_df[,2],
                   conserved_ccc_df[,3],conserved_ccc_df[,4])
  filtered_ta <- cellchat_df[ta_idx %in% qu_idx,]
  conserved_ccc_df$prob = filtered_ta$prob
  return(conserved_ccc_df)
}

