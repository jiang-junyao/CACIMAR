
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
#' @examples
Identify_ConservedCCI <- function(OrthG,
                                  Species1_CCI,
                                  Species2_CCI,
                                  ConservedCellType,
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
  }else if(Spec2 == tolower(Species_name1) & Spec1 == tolower(Species_name2)){
    sp1 = Species1_CCI
    sp2 = Species2_CCI
    Species1_CCI = sp2
    Species2_CCI = sp1

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
  ### select cell types may contain interactions
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
  ### identify conserved ligand
  orthg_cci = c()
  for (i in 1:nrow(Species1_CCI)) {
    cci1 = c(Species1_CCI[i,1],Species1_CCI[i,2])
    for (j in 1:nrow(Species2_CCI)) {
      cci2 = c(Species2_CCI[j,1],Species2_CCI[j,2])
      orthg_check = filter_orthg(cci1,cci2,OrthG,db_idx)
      if (orthg_check) {
        orthg_cci = c(orthg_cci,paste0(i,j))
      }
    }
  }
  conLRlist = list()
  for (i in orthg_cci) {
    LR_idx = unlist(strsplit(i,''))
    conLR = cbind(Species1_CCI[LR_idx[1],],Species1_CCI[LR_idx[2],])
    conLRlist[[i]] = conLR
  }
  return(do.call(dplyr::bind_rows,conLRlist))
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







