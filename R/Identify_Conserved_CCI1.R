#' Identify conserved cell-cell interaction from two species
#' @description get conserved cell-cell interaction between two species
#'
#' @param species1_cci result of SingCellSignalR of species1
#' @param species2_cci result of SingCellSignalR of species2
#' @param species_names_ref dataframe, if you construct new homolog database and use new species homolog information in this analysis, you should provide the names corresponding relationship. it should be a dataframe contain at lease two column full_name and sy_name.
#' @param conserved_cell_types_df dataframe, like conserved_cell_types_mm_zf <- data.frame("mm" = c("mmRods", "mmRod BC", "mmCones", "mmPericytes", "mmV/E cells", "mmRGC", "mmGABAergic AC", "mmRPE", "mmResting MG", "mmActivated MG", "mmMicroglia"), "zf" = c("zfRods", "zfCone BC", "zfCones", "zfPericytes", "zfV/E cells", "zfRGC", "zfGABAergic AC", "zfRPE", "zfResting MG", "zfActivated MG", "zfMicroglia"))
#' @param species_name1 two character to represent species1, like "mm". You should set this value from dataframe species_names_ref.rda
#' @param species_name2 two character to represent species2, like "zf". You should set this value from dataframe species_names_ref.rda
#' @param HOM_matrix dataframe, default is NULL, if you used new created homolog database, you should put it here
#'
#' @return list for conserved result for species1 and species2
#' @export
#'
#' @examples
#' conserved_cell_types_mm_zf <- data.frame("mm" = c("mmRods", "mmRod BC", "mmCones", "mmPericytes", "mmV/E cells", "mmRGC", "mmGABAergic AC", "mmRPE", "mmResting MG", "mmActivated MG", "mmMicroglia"), "zf" = c("zfRods", "zfCone BC", "zfCones", "zfPericytes", "zfV/E cells", "zfRGC", "zfGABAergic AC", "zfRPE", "zfResting MG", "zfActivated MG", "zfMicroglia"))
#' conserved_result_mm_zf <- Identify_Conserved_CCI1(species1_cci=SingleCellSignalR_mouse_result,
#'                                                   species2_cci=SingleCellSignalR_zebrafish_result,
#'                                                   conserved_cell_types_df=conserved_cell_types_mm_zf,
#'                                                   species_name1 = "mm",
#'                                                   species_name2 = "zf")
Identify_Conserved_CCI1 <- function(species1_cci,
                                    species2_cci,
                                    species_names_ref = NULL,
                                    conserved_cell_types_df,
                                    species_name1 = "mm",
                                    species_name2 = "zf",
                                    HOM_matrix = NULL) {
  # load("/data2/lijinlian/CACIMAR2/ChordDiagram_three_species/species_names_ref.rda")
  if (is.null(species_names_ref)) {
    load(system.file("extdata", "species_names_ref.rda", package = "CACIMAR"))
  }
  if (is.null(HOM_matrix)) {
    load(system.file("extdata", "HOM_matrix.rda", package = "CACIMAR"))
  }
  db_species1 = species_names_ref$full_name[which(species_names_ref$sy_name == species_name1)]
  db_species2 = species_names_ref$full_name[which(species_names_ref$sy_name == species_name2)]
  species1_cci$source <- paste0(species_name1, species1_cci$source)
  species1_cci$target <- paste0(species_name1, species1_cci$target)
  sp1_cci <- subset(species1_cci, source %in% conserved_cell_types_df[[species_name1]])
  sp1_cci <- subset(sp1_cci, target %in% conserved_cell_types_df[[species_name1]])
  species2_cci$source <- paste0(species_name2, species2_cci$source)
  species2_cci$target <- paste0(species_name2, species2_cci$target)
  sp2_cci <- subset(species2_cci, source %in% conserved_cell_types_df[[species_name2]])
  sp2_cci <- subset(sp2_cci, target %in% conserved_cell_types_df[[species_name2]])
  sp1_cci$c_pair <- paste(sp1_cci$source, sp1_cci$target, sep = "_")
  sp2_cci$c_pair <- paste(sp2_cci$source, sp2_cci$target, sep = "_")
  sp1_cci$lr_p <- paste(sp1_cci$ligand, sp1_cci$receptor, sep = "_")
  sp2_cci$lr_p <- paste(sp2_cci$ligand, sp2_cci$receptor, sep = "_")
  sp1_cci_conserved_df_list <- list()
  sp2_cci_conserved_df_list <- list()
  for (i in seq_len(length(unique(sp1_cci$c_pair)))) {
    sp1_cci_sub <- subset(sp1_cci, c_pair %in% unique(sp1_cci$c_pair)[i])
    sp2_cci_sub <- subset(sp2_cci, source == conserved_cell_types_df[[species_name2]][grep(paste0("^", unique(sp1_cci_sub$source), "$"), conserved_cell_types_df[[species_name1]])] & target == conserved_cell_types_df[[species_name2]][grep(paste0("^", unique(sp1_cci_sub$target), "$"), conserved_cell_types_df[[species_name1]])])
    sp1_cci_sub_conserved_df <- data.frame()
    sp2_cci_sub_conserved_df <- data.frame()
    for (j in seq_len(nrow(sp1_cci_sub))) {
      sp2_l <- unlist(strsplit(HOM_matrix[[db_species2]][grep(paste0("^", sp1_cci_sub$ligand[j], "$"), HOM_matrix[[db_species1]])], ","))
      sp2_r <- unlist(strsplit(HOM_matrix[[db_species2]][grep(paste0("^", sp1_cci_sub$receptor[j], "$"), HOM_matrix[[db_species1]])], ","))
      sp2_l <- as.vector(na.omit(sp2_l))
      sp2_r <- as.vector(na.omit(sp2_r))
      if (length(sp2_l) > 0 & length(sp2_r) > 0) {
        combinations <- do.call(rbind, lapply(1:length(sp2_l), function(ll_rr) cbind(sp2_l[ll_rr], sp2_r)))
        sp2_exp_lr_p <- paste(combinations[,1], combinations[,2], sep = "_")
        sp2_conserved_id <- unlist(sapply(sp2_exp_lr_p, function(x) grep(x, sp2_cci_sub$lr_p)))
        names(sp2_conserved_id) <- NULL
        if (length(sp2_conserved_id) > 0) {
          if (!is.na(HOM_matrix[[db_species2]][grep(paste0("^", sp1_cci_sub$ligand[j], "$"), HOM_matrix[[db_species1]])]) & !is.na(HOM_matrix[[db_species2]][grep(paste0("^", sp1_cci_sub$receptor[j], "$"), HOM_matrix[[db_species1]])])) {
            sp2_cci_sub_conserved_df <- rbind(sp2_cci_sub_conserved_df, sp2_cci_sub[sp2_conserved_id, ])
            sp1_cci_sub_conserved_df <- rbind(sp1_cci_sub_conserved_df, sp1_cci_sub[j, ])
          }
        }
      }
    }
    sp1_cci_conserved_df_list[[unique(sp1_cci_sub$c_pair)]] <- sp1_cci_sub_conserved_df
    sp2_cci_conserved_df_list[[unique(sp2_cci_sub$c_pair)]] <- sp2_cci_sub_conserved_df
  }
  sp1_cci_conserved_df <- do.call(rbind, sp1_cci_conserved_df_list)
  sp2_cci_conserved_df <- do.call(rbind, sp2_cci_conserved_df_list)
  result_list <- list("sp1_cci_conserved_df"=sp1_cci_conserved_df, "sp2_cci_conserved_df"=sp2_cci_conserved_df, "sp1_cci_conserved_df_list"=sp1_cci_conserved_df_list, "sp2_cci_conserved_df_list"=sp2_cci_conserved_df_list)
  saveRDS(result_list, paste(species_name1, species_name2, "conserved_result.Rds", sep = "_"))
  return(result_list)
}
#' Caculate conserved score of cell-cell interaction
#' @description caculate the conserved score of cell-cell interaction with summed weights of ligand-receptor interactions for two species
#'
#' @param conserved_result_df dataframe, result from function Identify_Conserved_CCI1
#' @param species1_cci dataframe, result from SingleCellSignalR of species1
#' @param species2_cci dataframe, result from SingleCellSignalR of species2
#' @param conserved_cell_types_df dataframe, contain the conserved cell type for each species, like conserved_cell_types_mm_zf <- data.frame("mm" = c("mmRods", "mmRod BC", "mmCones", "mmPericytes", "mmV/E cells", "mmRGC", "mmGABAergic AC", "mmRPE", "mmResting MG", "mmActivated MG", "mmMicroglia"), "zf" = c("zfRods", "zfCone BC", "zfCones", "zfPericytes", "zfV/E cells", "zfRGC", "zfGABAergic AC", "zfRPE", "zfResting MG", "zfActivated MG", "zfMicroglia"))
#' @param species_name1 two character to represent species1, like "mm". You should set this value from dataframe species_names_ref.rda
#' @param species_name2 two character to represent species2, like "zf". You should set this value from dataframe species_names_ref.rda
#'
#' @return dataframe, conserved Weights table
#' @export
#'
#' @examples
#' conserved_cell_types_mm_zf <- data.frame("mm" = c("mmRods", "mmRod BC", "mmCones", "mmPericytes", "mmV/E cells", "mmRGC", "mmGABAergic AC", "mmRPE", "mmResting MG", "mmActivated MG", "mmMicroglia"), "zf" = c("zfRods", "zfCone BC", "zfCones", "zfPericytes", "zfV/E cells", "zfRGC", "zfGABAergic AC", "zfRPE", "zfResting MG", "zfActivated MG", "zfMicroglia"))
#' cci_conserved_Weights_table_mm_zf <- Caculate_cell_pair_cci_score(conserved_result_df=conserved_result_mm_zf,
#' species1_cci = SingleCellSignalR_mouse_result,
#' species2_cci = SingleCellSignalR_zebrafish_result,
#' conserved_cell_types_df = conserved_cell_types_mm_zf,
#' species_name1 = "mm",
#' species_name2 = "zf")
Caculate_cell_pair_cci_score <- function(conserved_result_df, species1_cci, species2_cci, conserved_cell_types_df, species_name1, species_name2) {
  sp1_cci_conserved_Weights <- calculate_cell_pair_Weights1(conserved_result_df[["sp1_cci_conserved_df"]], used_col = c("source", "target", "LRscore_scale"))
  sp2_cci_conserved_Weights <- calculate_cell_pair_Weights1(conserved_result_df[["sp2_cci_conserved_df"]], used_col = c("source", "target", "LRscore_scale"))
  sp1_cci_cell_pair_Weights <- calculate_cell_pair_Weights1(species1_cci, used_col = c("source", "target", "LRscore_scale"), value_name = "T_Weight")
  sp2_cci_cell_pair_Weights <- calculate_cell_pair_Weights1(species2_cci, used_col = c("source", "target", "LRscore_scale"), value_name = "T_Weight")
  sp1_cci_cell_pair_Weights$source <- paste0(species_name1, sp1_cci_cell_pair_Weights$source)
  sp1_cci_cell_pair_Weights$target <- paste0(species_name1, sp1_cci_cell_pair_Weights$target)
  sp2_cci_cell_pair_Weights$source <- paste0(species_name2, sp2_cci_cell_pair_Weights$source)
  sp2_cci_cell_pair_Weights$target <- paste0(species_name2, sp2_cci_cell_pair_Weights$target)
  cci_conserved_Weights_table <- sp1_cci_conserved_Weights
  for (k in seq_len(nrow(cci_conserved_Weights_table))) {
    cci_conserved_Weights_table$sp1_T_Weight[k] <- sp1_cci_cell_pair_Weights$T_Weight[which(sp1_cci_cell_pair_Weights$source == cci_conserved_Weights_table$source[k] & sp1_cci_cell_pair_Weights$target == cci_conserved_Weights_table$target[k])]
    cci_conserved_Weights_table$sp2_source[k] <- conserved_cell_types_df[[species_name2]][which(conserved_cell_types_df[[species_name1]] == cci_conserved_Weights_table$source[k])]
    cci_conserved_Weights_table$sp2_target[k] <- conserved_cell_types_df[[species_name2]][which(conserved_cell_types_df[[species_name1]] == cci_conserved_Weights_table$target[k])]
  }
  for (m in seq_len(nrow(cci_conserved_Weights_table))) {
    cci_conserved_Weights_table$sp2_conserved_weight[m] <- sp2_cci_conserved_Weights$conserved_weight[which(sp2_cci_conserved_Weights$source == cci_conserved_Weights_table$sp2_source[m] & sp2_cci_conserved_Weights$target == cci_conserved_Weights_table$sp2_target[m])]
    cci_conserved_Weights_table$sp2_T_Weight[m] <- sp2_cci_cell_pair_Weights$T_Weight[which(sp2_cci_cell_pair_Weights$source == cci_conserved_Weights_table$sp2_source[m] & sp2_cci_cell_pair_Weights$target == cci_conserved_Weights_table$sp2_target[m])]
  }
  for (n in seq_len(nrow(cci_conserved_Weights_table))) {
    cci_conserved_Weights_table$conserved_score[n] <- sum(cci_conserved_Weights_table[["conserved_weight"]][n], cci_conserved_Weights_table[["sp2_conserved_weight"]][n])/sum(cci_conserved_Weights_table[["sp1_T_Weight"]][n], cci_conserved_Weights_table[["sp2_T_Weight"]][n])
  }
  return(cci_conserved_Weights_table)
}
calculate_cell_pair_Weights1 <- function(data, used_col = c("source", "target", "LRscore_scale"), value_name = "conserved_weight") {
  ccc_cweight <- reshape2::dcast(data[, used_col],
                                 source ~ target, sum,
                                 value.var = 'LRscore_scale')
  rownames(ccc_cweight) <- ccc_cweight[, 1]
  all_cweight_df_long <- reshape2::melt(ccc_cweight,
                                        id.vars = "source",
                                        variable.name = "target",
                                        value.name = value_name)
  return(all_cweight_df_long)
}
#' Make data fit for ChordDiagram
#'
#' @param cci_conserved_Weights_table dataframe, result from function Caculate_cell_pair_cci_score
#' @param species_name1 two character to represent species1, like "mm"
#' @param species_name2 two character to represent species1, like "zf"
#'
#' @return dataframe, data for ChordDiagram
#' @export
#'
#' @examples
#' cci_data_mm_zf <- Make_ChordDiagram_data(cci_conserved_Weights_table = cci_conserved_Weights_table_mm_zf, species_name1 = "mm", species_name2 = "zf")
Make_ChordDiagram_data <- function(cci_conserved_Weights_table, species_name1, species_name2) {
  cci_data_sp1_sp2 <- cci_conserved_Weights_table[, c("source", "target", "sp2_source", "sp2_target", "conserved_score")]
  cci_data_sp1_sp2$source <- sub(species_name1, "", cci_data_sp1_sp2$source)
  cci_data_sp1_sp2$target <- sub(species_name1, "", cci_data_sp1_sp2$target)
  cci_data_sp1_sp2$sp2_source <- sub(species_name2, "", cci_data_sp1_sp2$sp2_source)
  cci_data_sp1_sp2$sp2_target <- sub(species_name2, "", cci_data_sp1_sp2$sp2_target)
  cci_data_sp1_sp2$rep_source <- ifelse(cci_data_sp1_sp2$source == cci_data_sp1_sp2$sp2_source, cci_data_sp1_sp2$source, paste(cci_data_sp1_sp2$source, cci_data_sp1_sp2$sp2_source, sep = "/"))
  cci_data_sp1_sp2$rep_target <- ifelse(cci_data_sp1_sp2$target == cci_data_sp1_sp2$sp2_target, cci_data_sp1_sp2$target, paste(cci_data_sp1_sp2$target, cci_data_sp1_sp2$sp2_target, sep = "/"))
  write.table(cci_data_sp1_sp2, paste(species_name1, species_name2, "cci_data_for_ChordDiagram.txt", sep = "_"))
  cci_data <- cci_data_sp1_sp2[, c("rep_source", "rep_target", "conserved_score")]
  return(cci_data)
}
