
#' Build homologous gene database
#' @description Build homologous gene database according to Vertebrate Homology
#' data in MGI database. This function currently supports ten species: cattle,
#' chicken, chimpanzee, dog, frog, human, macaque, mouse, rat, and zebrafish. After
#' building the database, this function also integrates biomaRt to add ENSEMBEL ID
#' for each gene in the database.
#' @param MGI MGI database, download from http://www.informatics.jax.org/
#' @param Species_name1 The name of the first species. input 'mm' for mouse, 'hs'
#' for human, 'zf' for zebrafish, 'ch' for chicken, 'cf' for dog, 'pt' for chimpanzee,
#' 'xt' for frog, 'rn' for rat, 'bt' for cattle, and 'rh' for macaque.
#' @param Species_name2 The name of the second species.  input 'mm' for mouse, 'hs'
#' for human, 'zf' for zebrafish, 'ch' for chicken, 'cf' for dog, 'pt' for chimpanzee,
#' 'xt' for frog, 'rn' for rat, 'bt' for cattle, and 'rh' for macaque.
#' @param
#' @importFrom  stats na.omit
#' @importFrom dplyr group_by
#' @importFrom dplyr group_map
#' @importFrom biomaRt useMart
#' @importFrom biomaRt getBM
#' @importFrom purrr map_chr
#' @return homologous gene database of two species
#' @export
#'
#' @example load(system.file("extdata", "MGIrda", package = "CACIMAR"))
#' buildHomDatabase(MGI,'mm','pt')
buildHomDatabase <- function(MGI,Species_name1,Species_name2){
  ### set information for species 1
  if (Species_name1=='mm') {
    Sp1name <- 'mouse, laboratory'
    use_name1 <- 'mm'
    Mart_ID1 <- "mmusculus_gene_ensembl"
  }
  if (Species_name1=='hs') {
    Sp1name <- 'human'
    use_name1 <- 'hs'
    Mart_ID1 <- "hsapiens_gene_ensembl"
  }
  if (Species_name1=='zf') {
    Sp1name <- 'zebrafish'
    use_name1 <- 'zf'
    Mart_ID1 <- "drerio_gene_ensembl"
  }
  if (Species_name1=='ch') {
    Sp1name <- 'chicken'
    use_name1 <- 'ch'
    Mart_ID1 <- "ggallus_gene_ensembl"
  }
  if (Species_name1=='cf') {
    Sp1name <- 'dog, domestic'
    use_name1 <- 'cf'
    Mart_ID1 <- "clfamiliaris_gene_ensembl"
  }
  if (Species_name1=='pt') {
    Sp1name <- 'chimpanzee'
    use_name1 <- 'pt'
    Mart_ID1 <- "ptroglodytes_gene_ensembl"
  }
  if (Species_name1=='xt') {
    Sp1name <- 'frog, western clawed'
    use_name1 <- 'xt'
    Mart_ID1 <- "xtropicalis_gene_ensembl"
  }
  if (Species_name1=='rn') {
    Sp1name <- 'rat'
    use_name1 <- 'rn'
    Mart_ID1 <- "rnorvegicus_gene_ensembl"
  }
  if (Species_name1=='bt') {
    Sp1name <- 'cattle'
    use_name1 <- 'bt'
    Mart_ID1 <- "btaurus_gene_ensembl"
  }
  if (Species_name1=='rh') {
    Sp1name <- 'macaque, rhesus'
    use_name1 <- 'rh'
    Mart_ID1 <- "mmulatta_gene_ensembl"
  }


  ### set information for species 2
  if (Species_name2=='mm') {
    Sp2name <- 'mouse, laboratory'
    use_name2 <- 'mm'
    Mart_ID2 <- "mmusculus_gene_ensembl"
  }
  if (Species_name2=='hs') {
    Sp2name <- 'human'
    use_name2 <- 'hs'
    Mart_ID2 <- "hsapiens_gene_ensembl"
  }
  if (Species_name2=='zf') {
    Sp2name <- 'zebrafish'
    use_name2 <- 'zf'
    Mart_ID2 <- "drerio_gene_ensembl"
  }
  if (Species_name2=='ch') {
    Sp2name <- 'chicken'
    use_name2 <- 'ch'
    Mart_ID2 <- "ggallus_gene_ensembl"
  }
  if (Species_name2=='cf') {
    Sp2name <- 'dog, domestic'
    use_name2 <- 'cf'
    Mart_ID2 <- "clfamiliaris_gene_ensembl"
  }
  if (Species_name2=='pt') {
    Sp2name <- 'chimpanzee'
    use_name2 <- 'pt'
    Mart_ID2 <- "ptroglodytes_gene_ensembl"
  }
  if (Species_name2=='xt') {
    Sp2name <- 'frog, western clawed'
    use_name2 <- 'xt'
    Mart_ID2 <- "xtropicalis_gene_ensembl"
  }
  if (Species_name2=='rn') {
    Sp2name <- 'rat'
    use_name2 <- 'rn'
    Mart_ID2 <- "rnorvegicus_gene_ensembl"
  }
  if (Species_name2=='bt') {
    Sp2name <- 'cattle'
    use_name2 <- 'bt'
    Mart_ID2 <- "btaurus_gene_ensembl"
  }
  if (Species_name2=='rh') {
    Sp2name <- 'macaque, rhesus'
    use_name2 <- 'rh'
    Mart_ID2 <- "mmulatta_gene_ensembl"
  }


  ### build database
  HOM_AllOrganism <- MGI
  colnames(HOM_AllOrganism)[1] = 'HomoloGene.ID'
  typeName <-paste0(use_name1,'_',use_name2,'_')
  Sp1 <- HOM_AllOrganism[HOM_AllOrganism$Common.Organism.Name==Sp1name,]
  Sp1 <- Sp1[order(Sp1$HomoloGene.ID),]
  Sp2 <- HOM_AllOrganism[HOM_AllOrganism$Common.Organism.Name==Sp2name,]
  Sp2 <- Sp2[order(Sp2$HomoloGene.ID),]
  orthg_group1 <- dplyr::group_by(Sp1,HomoloGene.ID)
  orthg_gene1 <- t(as.data.frame(dplyr::group_map(orthg_group1,~get_genes(.x))))
  OrthID1 <- Sp1[!duplicated(Sp1$HomoloGene.ID),1]
  Sp1_orthg <- as.data.frame(cbind(OrthID1,orthg_gene1))

  orthg_group2 <- dplyr::group_by(Sp2,HomoloGene.ID)
  orthg_gene2 <- t(as.data.frame(dplyr::group_map(orthg_group2,~get_genes(.x))))
  OrthID2 <- Sp2[!duplicated(Sp2$HomoloGene.ID),1]
  Sp2_orthg <- as.data.frame(cbind(OrthID2,orthg_gene2))

  Sp1_orthg$Species2Symbol <- Sp2_orthg[match(Sp1_orthg$OrthID1,Sp2_orthg$OrthID2),2]
  Sp1_orthg$Species2Type <- Sp2_orthg[match(Sp1_orthg$OrthID1,Sp2_orthg$OrthID2),3]
  Sp1_orthg <- na.omit(Sp1_orthg)
  OrthType <- paste0(typeName,Sp1_orthg[,4],'T',Sp1_orthg[,4])
  Sp1_orthg$Type <- OrthType
  final_orthg <- Sp1_orthg[,c(7,3,2,6,5)]
  colnames(final_orthg)[2:5] = c(paste0(Species_name1,'_ID'),
                            paste0(Species_name1,'_Symbol'),
                            paste0(Species_name2,'_ID'),
                            paste0(Species_name2,'_Symbol'))
  return(final_orthg)
  # ### add ens
  # if (AddENS == TRUE) {
  #   mart1 <- biomaRt::useMart("ensembl",Mart_ID1,host = host)
  #   gene_id1 <- biomaRt::getBM(attributes=c("external_gene_name","ensembl_gene_id"),
  #                            filters = "external_gene_name",values = unlist(strsplit(final_orthg$V2,';')),
  #                            mart = mart1)
  #   ENS_ID1 <- purrr::map_chr(final_orthg$V2,~match_gene1(.x,gene_id1))
  #   final_orthg$Species1ENS <- ENS_ID1
  #
  #
  #   mart2 <- biomaRt::useMart("ensembl",Mart_ID2,host = host)
  #   gene_id2<-biomaRt::getBM(attributes=c("external_gene_name","ensembl_gene_id"),
  #                            filters = "external_gene_name",values = unlist(strsplit(final_orthg$Species2Symbol,';')),
  #                            mart = mart2)
  #   ENS_ID2 <- purrr::map_chr(final_orthg$Species2Symbol,~match_gene2(.x,gene_id2))
  #   final_orthg$Species2ENS <- ENS_ID2
  #   final_orthg <- final_orthg[,c(1,4,2,5,3)]
  #   colnames(final_orthg)[2:5] <- c(paste0(use_name1,'_ID'),paste0(use_name1,'_Symbol')
  #                                   ,paste0(use_name2,'_ID'),paste0(use_name2,'_Symbol'))
  #   return(final_orthg)
  # }else{
  #   colnames(final_orthg) <- c('Type',paste0(use_name1,'_Symbol'),
  #                              paste0(use_name2,'_Symbol'))
  #   return(final_orthg)
  #   }
}





get_genes <- function(data1){
  gene<-paste(as.character(data1$Symbol),collapse = ',')
  ent<-paste(as.character(data1$EntrezGene.ID),collapse = ',')
  if (length(data1$Symbol)>1) {
    type <- 'N'
  }else{type<-'1'}
  df=c(gene,ent,type)
  return(df)
}

match_gene1 <- function(x,gene_id1){
  x1 <- strsplit(x,';')[[1]]
  x2 <- c()
  for (i in x1) {
    if (i %in% gene_id1$external_gene_name) {
      x2 <- c(x2,gene_id1[gene_id1$external_gene_name==i,2])
    }else{x2 <- c(x2,NA)}
  }
  x3 <- paste(x2,collapse = ';')
  return(x3)
}
match_gene2 <- function(x,gene_id2){
  x1 <- strsplit(x,';')[[1]]
  x2 <- c()
  for (i in x1) {
    if (i %in% gene_id2$external_gene_name) {
      x2 <- c(x2,gene_id2[gene_id2$external_gene_name==i,2])
    }else{x2 <- c(x2,NA)}
  }
  x3 <- paste(x2,collapse = ';')
  return(x3)
}
