#' Identify conserved cell-cell interaction from three species
#' @description
#' find the conserved cell-cell interaction from three species and caculate the conserved score of cell-cell interaction
#'
#' @param species1_cci result of SingCellSignalR of species1
#' @param species2_cci result of SingCellSignalR of species2
#' @param species3_cci result of SingCellSignalR of species3
#' @param species_names_ref dataframe, if you construct new homolog database and use new species homolog information in this analysis, you should provide the names corresponding relationship. it should be a dataframe contain at lease two column full_name and sy_name.
#' @param conserved_cell_types_df dataframe, like conserved_cell_types_mm_zf_ch <- data.frame('mm' = c("mmResting MG","mmGABAergic AC", "mmRGC", "mmCones"), "zf" = c("zfResting MG", "zfGABAergic AC", "zfRGC", "zfCones"), 'ch' = c("chResting MG","chGABAergic AC", "chRGC", "chCones"))
#' @param species_name1 two character to represent species1, like "mm". You should set this value from dataframe species_names_ref.rda
#' @param species_name2 two character to represent species2, like "zf". You should set this value from dataframe species_names_ref.rda
#' @param species_name3 two character to represent species3, like "ch". You should set this value from dataframe species_names_ref.rda
#' @param HOM_matrix dataframe, default is NULL, if you used new created homolog database, you should put it here
#'
#' @return list for conserved result for each species
#' @export
#'
#' @examples
Identify_Conserved_CCI2 <- function(species1_cci,
                                    species2_cci,
                                    species3_cci,
                                    species_names_ref = NULL,
                                    conserved_cell_types_df,
                                    species_name1 = "mm",
                                    species_name2 = "zf",
                                    species_name3 = "ch",
                                    HOM_matrix = NULL) {
  if (is.null(species_names_ref)) {
    load(system.file("extdata", "species_names_ref.rda", package = "CACIMAR"))
  }
  if (is.null(HOM_matrix)) {
    load(system.file("extdata", "HOM_matrix.rda", package = "CACIMAR"))
  }
  db_species1 = species_names_ref$full_name[which(species_names_ref$sy_name == species_name1)]
  db_species2 = species_names_ref$full_name[which(species_names_ref$sy_name == species_name2)]
  db_species3 = species_names_ref$full_name[which(species_names_ref$sy_name == species_name3)]
  species1_cci$source <- paste0(species_name1, species1_cci$source)
  species1_cci$target <- paste0(species_name1, species1_cci$target)
  sp1_cci <- subset(species1_cci, source %in% conserved_cell_types_df[[species_name1]])
  sp1_cci <- subset(sp1_cci, target %in% conserved_cell_types_df[[species_name1]])
  species2_cci$source <- paste0(species_name2, species2_cci$source)
  species2_cci$target <- paste0(species_name2, species2_cci$target)
  sp2_cci <- subset(species2_cci, source %in% conserved_cell_types_df[[species_name2]])
  sp2_cci <- subset(sp2_cci, target %in% conserved_cell_types_df[[species_name2]])
  species3_cci$source <- paste0(species_name3, species3_cci$source)
  species3_cci$target <- paste0(species_name3, species3_cci$target)
  sp3_cci <- subset(species3_cci, source %in% conserved_cell_types_df[[species_name3]])
  sp3_cci <- subset(sp3_cci, target %in% conserved_cell_types_df[[species_name3]])
  sp1_cci$c_pair <- paste(sp1_cci$source, sp1_cci$target, sep = "_")
  sp2_cci$c_pair <- paste(sp2_cci$source, sp2_cci$target, sep = "_")
  sp3_cci$c_pair <- paste(sp3_cci$source, sp3_cci$target, sep = "_")
  sp1_cci$lr_p <- paste(sp1_cci$ligand, sp1_cci$receptor, sep = "_")
  sp2_cci$lr_p <- paste(sp2_cci$ligand, sp2_cci$receptor, sep = "_")
  sp3_cci$lr_p <- paste(sp3_cci$ligand, sp3_cci$receptor, sep = "_")
  sp1_cci_conserved_df_list <- list()
  sp2_cci_conserved_df_list <- list()
  sp3_cci_conserved_df_list <- list()
  for (i in seq_len(length(unique(sp1_cci$c_pair)))) {
    sp1_cci_sub <- subset(sp1_cci, c_pair %in% unique(sp1_cci$c_pair)[i])
    sp2_cci_sub <- subset(sp2_cci, source == conserved_cell_types_df[[species_name2]][grep(paste0("^", unique(sp1_cci_sub$source), "$"), conserved_cell_types_df[[species_name1]])] & target == conserved_cell_types_df[[species_name2]][grep(paste0("^", unique(sp1_cci_sub$target), "$"), conserved_cell_types_df[[species_name1]])])
    sp3_cci_sub <- subset(sp3_cci, source == conserved_cell_types_df[[species_name3]][grep(paste0("^", unique(sp1_cci_sub$source), "$"), conserved_cell_types_df[[species_name1]])] & target == conserved_cell_types_df[[species_name3]][grep(paste0("^", unique(sp1_cci_sub$target), "$"), conserved_cell_types_df[[species_name1]])])
    sp1_cci_sub_conserved_df <- data.frame()
    sp2_cci_sub_conserved_df <- data.frame()
    sp3_cci_sub_conserved_df <- data.frame()
    for (j in seq_len(nrow(sp1_cci_sub))) {
      sp2_l <- unlist(strsplit(HOM_matrix[[db_species2]][grep(paste0("^", sp1_cci_sub$ligand[j], "$"), HOM_matrix[[db_species1]])], ","))
      sp2_r <- unlist(strsplit(HOM_matrix[[db_species2]][grep(paste0("^", sp1_cci_sub$receptor[j], "$"), HOM_matrix[[db_species1]])], ","))
      sp2_l <- as.vector(na.omit(sp2_l))
      sp2_r <- as.vector(na.omit(sp2_r))
      if (length(sp2_l) > 0 & length(sp2_r) > 0) {
        combinations_sp2 <- do.call(rbind, lapply(1:length(sp2_l), function(ll_rr) cbind(sp2_l[ll_rr], sp2_r)))
        sp2_exp_lr_p <- paste(combinations_sp2[,1], combinations_sp2[,2], sep = "_")
        sp2_conserved_id <- unlist(sapply(sp2_exp_lr_p, function(x) grep(x, sp2_cci_sub$lr_p)))
        names(sp2_conserved_id) <- NULL
        sp2_conserved_id <- as.vector(na.omit(sp2_conserved_id))
      }
      sp3_l <- unlist(strsplit(HOM_matrix[[db_species3]][grep(paste0("^", sp1_cci_sub$ligand[j], "$"), HOM_matrix[[db_species1]])], ","))
      sp3_r <- unlist(strsplit(HOM_matrix[[db_species3]][grep(paste0("^", sp1_cci_sub$receptor[j], "$"), HOM_matrix[[db_species1]])], ","))
      sp3_l <- as.vector(na.omit(sp3_l))
      sp3_r <- as.vector(na.omit(sp3_r))
      if (length(sp3_l) > 0 & length(sp3_r) > 0) {
        combinations_sp3 <- do.call(rbind, lapply(1:length(sp3_l), function(lll_rrr) cbind(sp3_l[lll_rrr], sp3_r)))
        sp3_exp_lr_p <- paste(combinations_sp3[,1], combinations_sp3[,2], sep = "_")
        sp3_conserved_id <- unlist(sapply(sp3_exp_lr_p, function(x) grep(x, sp3_cci_sub$lr_p)))
        names(sp3_conserved_id) <- NULL
        sp3_conserved_id <- as.vector(na.omit(sp3_conserved_id))
      }
      if (length(sp2_conserved_id) > 0 & length(sp3_conserved_id) > 0) {
        if (!is.na(HOM_matrix[[db_species2]][grep(paste0("^", sp1_cci_sub$ligand[j], "$"), HOM_matrix[[db_species1]])]) & !is.na(HOM_matrix[[db_species2]][grep(paste0("^", sp1_cci_sub$receptor[j], "$"), HOM_matrix[[db_species1]])]) & !is.na(HOM_matrix[[db_species3]][grep(paste0("^", sp1_cci_sub$ligand[j], "$"), HOM_matrix[[db_species1]])]) & !is.na(HOM_matrix[[db_species3]][grep(paste0("^", sp1_cci_sub$receptor[j], "$"), HOM_matrix[[db_species1]])])) {
          sp3_cci_sub_conserved_df <- rbind(sp3_cci_sub_conserved_df, sp3_cci_sub[sp3_conserved_id, ])
          sp2_cci_sub_conserved_df <- rbind(sp2_cci_sub_conserved_df, sp2_cci_sub[sp2_conserved_id, ])
          sp1_cci_sub_conserved_df <- rbind(sp1_cci_sub_conserved_df, sp1_cci_sub[j, ])
        }
      }
    }
    sp1_cci_conserved_df_list[[unique(sp1_cci_sub$c_pair)]] <- sp1_cci_sub_conserved_df
    sp2_cci_conserved_df_list[[unique(sp2_cci_sub$c_pair)]] <- sp2_cci_sub_conserved_df
    sp3_cci_conserved_df_list[[unique(sp3_cci_sub$c_pair)]] <- sp3_cci_sub_conserved_df
  }
  sp1_cci_conserved_df <- do.call(rbind, sp1_cci_conserved_df_list)
  sp1_cci_conserved_df <- sp1_cci_conserved_df[complete.cases(sp1_cci_conserved_df),]
  sp2_cci_conserved_df <- do.call(rbind, sp2_cci_conserved_df_list)
  sp2_cci_conserved_df <- sp2_cci_conserved_df[complete.cases(sp2_cci_conserved_df),]
  sp3_cci_conserved_df <- do.call(rbind, sp3_cci_conserved_df_list)
  sp3_cci_conserved_df <- sp3_cci_conserved_df[complete.cases(sp3_cci_conserved_df),]
  result_list <- list("sp1_cci_conserved_df"=sp1_cci_conserved_df, "sp2_cci_conserved_df"=sp2_cci_conserved_df, "sp3_cci_conserved_df"=sp3_cci_conserved_df, "sp1_cci_conserved_df_list"=sp1_cci_conserved_df_list, "sp2_cci_conserved_df_list"=sp2_cci_conserved_df_list, "sp3_cci_conserved_df_list"=sp3_cci_conserved_df_list)
  saveRDS(result_list, paste(species_name1, species_name2, species_name3, "conserved_result.Rds", sep = "_"))
  return(result_list)
}
#' Caculate conserved score of cell-cell interaction for three species
#' @description caculate the conserved score of cell-cell interaction with summed weights of ligand-receptor interactions for three species
#' @param conserved_result_df dataframe, result from function Identify_Conserved_CCI2
#' @param species1_cci dataframe, result from SingleCellSignalR of species1
#' @param species2_cci dataframe, result from SingleCellSignalR of species2
#' @param species3_cci dataframe, result from SingleCellSignalR of species3
#' @param conserved_cell_types_df dataframe, contain the conserved cell type for each species, like conserved_cell_types_mm_zf_ch <- data.frame('mm' = c("mmResting MG","mmGABAergic AC", "mmRGC", "mmCones"), "zf" = c("zfResting MG", "zfGABAergic AC", "zfRGC", "zfCones"), 'ch' = c("chResting MG","chGABAergic AC", "chRGC", "chCones"))
#' @param species_name1 two character to represent species1, like "mm". You should set this value from dataframe species_names_ref.rda
#' @param species_name2 two character to represent species2, like "zf". You should set this value from dataframe species_names_ref.rda
#' @param species_name3 two character to represent species3, like "ch". You should set this value from dataframe species_names_ref.rda
#'
#' @return dataframe, conserved Weights table
#' @export
#'
#' @examples
Caculate_cell_pair_cci_score2 <- function(conserved_result_df, species1_cci, species2_cci, species3_cci, conserved_cell_types_df, species_name1, species_name2, species_name3) {
  sp1_cci_conserved_Weights <- calculate_cell_pair_Weights1(conserved_result_df[["sp1_cci_conserved_df"]], used_col = c("source", "target", "LRscore_scale"))
  sp2_cci_conserved_Weights <- calculate_cell_pair_Weights1(conserved_result_df[["sp2_cci_conserved_df"]], used_col = c("source", "target", "LRscore_scale"))
  sp3_cci_conserved_Weights <- calculate_cell_pair_Weights1(conserved_result_df[["sp3_cci_conserved_df"]], used_col = c("source", "target", "LRscore_scale"))
  sp1_cci_cell_pair_Weights <- calculate_cell_pair_Weights1(species1_cci, used_col = c("source", "target", "LRscore_scale"), value_name = "T_Weight")
  sp2_cci_cell_pair_Weights <- calculate_cell_pair_Weights1(species2_cci, used_col = c("source", "target", "LRscore_scale"), value_name = "T_Weight")
  sp3_cci_cell_pair_Weights <- calculate_cell_pair_Weights1(species3_cci, used_col = c("source", "target", "LRscore_scale"), value_name = "T_Weight")
  sp1_cci_cell_pair_Weights$source <- paste0(species_name1, sp1_cci_cell_pair_Weights$source)
  sp1_cci_cell_pair_Weights$target <- paste0(species_name1, sp1_cci_cell_pair_Weights$target)
  sp2_cci_cell_pair_Weights$source <- paste0(species_name2, sp2_cci_cell_pair_Weights$source)
  sp2_cci_cell_pair_Weights$target <- paste0(species_name2, sp2_cci_cell_pair_Weights$target)
  sp3_cci_cell_pair_Weights$source <- paste0(species_name3, sp3_cci_cell_pair_Weights$source)
  sp3_cci_cell_pair_Weights$target <- paste0(species_name3, sp3_cci_cell_pair_Weights$target)
  cci_conserved_Weights_table <- sp1_cci_conserved_Weights
  for (k in seq_len(nrow(cci_conserved_Weights_table))) {
    cci_conserved_Weights_table$sp1_T_Weight[k] <- sp1_cci_cell_pair_Weights$T_Weight[which(sp1_cci_cell_pair_Weights$source == cci_conserved_Weights_table$source[k] & sp1_cci_cell_pair_Weights$target == cci_conserved_Weights_table$target[k])]
    cci_conserved_Weights_table$sp2_source[k] <- conserved_cell_types_df[[species_name2]][which(conserved_cell_types_df[[species_name1]] == cci_conserved_Weights_table$source[k])]
    cci_conserved_Weights_table$sp2_target[k] <- conserved_cell_types_df[[species_name2]][which(conserved_cell_types_df[[species_name1]] == cci_conserved_Weights_table$target[k])]
    cci_conserved_Weights_table$sp3_source[k] <- conserved_cell_types_df[[species_name3]][which(conserved_cell_types_df[[species_name1]] == cci_conserved_Weights_table$source[k])]
    cci_conserved_Weights_table$sp3_target[k] <- conserved_cell_types_df[[species_name3]][which(conserved_cell_types_df[[species_name1]] == cci_conserved_Weights_table$target[k])]
  }
  for (m in seq_len(nrow(cci_conserved_Weights_table))) {
    cci_conserved_Weights_table$sp2_conserved_weight[m] <- sp2_cci_conserved_Weights$conserved_weight[which(sp2_cci_conserved_Weights$source == cci_conserved_Weights_table$sp2_source[m] & sp2_cci_conserved_Weights$target == cci_conserved_Weights_table$sp2_target[m])]
    cci_conserved_Weights_table$sp2_T_Weight[m] <- sp2_cci_cell_pair_Weights$T_Weight[which(sp2_cci_cell_pair_Weights$source == cci_conserved_Weights_table$sp2_source[m] & sp2_cci_cell_pair_Weights$target == cci_conserved_Weights_table$sp2_target[m])]
  }
  for (m in seq_len(nrow(cci_conserved_Weights_table))) {
    cci_conserved_Weights_table$sp3_conserved_weight[m] <- sp3_cci_conserved_Weights$conserved_weight[which(sp3_cci_conserved_Weights$source == cci_conserved_Weights_table$sp3_source[m] & sp3_cci_conserved_Weights$target == cci_conserved_Weights_table$sp3_target[m])]
    cci_conserved_Weights_table$sp3_T_Weight[m] <- sp3_cci_cell_pair_Weights$T_Weight[which(sp3_cci_cell_pair_Weights$source == cci_conserved_Weights_table$sp3_source[m] & sp3_cci_cell_pair_Weights$target == cci_conserved_Weights_table$sp3_target[m])]
  }
  for (n in seq_len(nrow(cci_conserved_Weights_table))) {
    cci_conserved_Weights_table$conserved_score[n] <- sum(cci_conserved_Weights_table[["conserved_weight"]][n], cci_conserved_Weights_table[["sp2_conserved_weight"]][n], cci_conserved_Weights_table[["sp3_conserved_weight"]][n])/sum(cci_conserved_Weights_table[["sp1_T_Weight"]][n], cci_conserved_Weights_table[["sp2_T_Weight"]][n], cci_conserved_Weights_table[["sp3_T_Weight"]][n])
  }
  return(cci_conserved_Weights_table)
}
#' Make data fit for ChordDiagram for conserved cell types from three species
#'
#' @param cci_conserved_Weights_table  dataframe, result from function Caculate_cell_pair_cci_score2
#' @param species_name1 two character to represent species1, like "mm"
#' @param species_name2 two character to represent species2, like "zf"
#' @param species_name3 two character to represent species3, like "ch"
#'
#' @return dataframe, data for ChordDiagram
#' @export
#'
#' @examples
#' # cci_conserved_Weights_table_sp1_sp2_sp3 is the result of function Caculate_cell_pair_cci_score2
#'  cci_data_sp1_sp2_sp3 <- Make_ChordDiagram_data2(cci_conserved_Weights_table = cci_conserved_Weights_table_sp1_sp2_sp3,
#'  species_name1 = "mm",
#'  species_name2 = "zf",
#'  species_name3 = "ch")
Make_ChordDiagram_data2 <- function(cci_conserved_Weights_table, species_name1, species_name2, species_name3) {
  cci_data_sp1_sp2 <- cci_conserved_Weights_table[, c("source", "target", "sp2_source", "sp2_target", "sp3_source", "sp3_target", "conserved_score")]
  cci_data_sp1_sp2$source <- sub(species_name1, "", cci_data_sp1_sp2$source)
  cci_data_sp1_sp2$target <- sub(species_name1, "", cci_data_sp1_sp2$target)
  cci_data_sp1_sp2$sp2_source <- sub(species_name2, "", cci_data_sp1_sp2$sp2_source)
  cci_data_sp1_sp2$sp2_target <- sub(species_name2, "", cci_data_sp1_sp2$sp2_target)
  cci_data_sp1_sp2$sp3_source <- sub(species_name3, "", cci_data_sp1_sp2$sp3_source)
  cci_data_sp1_sp2$sp3_target <- sub(species_name3, "", cci_data_sp1_sp2$sp3_target)
  cci_data_sp1_sp2$rep_source <- ifelse(cci_data_sp1_sp2$source == cci_data_sp1_sp2$sp2_source & cci_data_sp1_sp2$source == cci_data_sp1_sp2$sp3_source & cci_data_sp1_sp2$sp3_source == cci_data_sp1_sp2$sp2_source, cci_data_sp1_sp2$source, paste(cci_data_sp1_sp2$source, cci_data_sp1_sp2$sp2_source, cci_data_sp1_sp2$sp3_source, sep = "/"))
  cci_data_sp1_sp2$rep_target <- ifelse(cci_data_sp1_sp2$target == cci_data_sp1_sp2$sp2_target & cci_data_sp1_sp2$target == cci_data_sp1_sp2$sp3_target & cci_data_sp1_sp2$sp3_target == cci_data_sp1_sp2$sp2_target, cci_data_sp1_sp2$target, paste(cci_data_sp1_sp2$target, cci_data_sp1_sp2$sp2_target, cci_data_sp1_sp2$sp3_target, sep = "/"))
  write.table(cci_data_sp1_sp2, paste(species_name1, species_name2, species_name3, "cci_data_for_ChordDiagram.txt", sep = "_"))
  cci_data <- cci_data_sp1_sp2[, c("rep_source", "rep_target", "conserved_score")]
  return(cci_data)
}
#' Caculate conserved interaction score and Make data for ChordDiagram
#'
#' @param conserved_result_species dataframe, result from function Identify_Conserved_CCI2
#' @param SingleCellSignalR_sp1_result dataframe, result from SingleCellSignalR of species1
#' @param SingleCellSignalR_sp2_result dataframe, result from SingleCellSignalR of species2
#' @param SingleCellSignalR_sp3_result dataframe, result from SingleCellSignalR of species3
#' @param conserved_cell_types_df dataframe, contain the conserved cell type for each species, like conserved_cell_types_mm_zf_ch <- data.frame('mm' = c("mmResting MG","mmGABAergic AC", "mmRGC", "mmCones"), "zf" = c("zfResting MG", "zfGABAergic AC", "zfRGC", "zfCones"), 'ch' = c("chResting MG","chGABAergic AC", "chRGC", "chCones"))
#' @param species_name1 two character to represent species1, like "mm". You should set this value from dataframe species_names_ref.rda
#' @param species_name2 two character to represent species2, like "zf". You should set this value from dataframe species_names_ref.rda
#' @param species_name3 two character to represent species3, like "ch". You should set this value from dataframe species_names_ref.rda
#'
#' @return list of conserved Weights table and data for chordDiagram
#' @export
#'
#' @examples
#' #not run
#' #conserved_result_mm_zf_ch is the result from function Identify_Conserved_CCI2
#' conserved_cci_result <- conserved_interaction_score(conserved_result_species = conserved_result_mm_zf_ch,
#' SingleCellSignalR_sp1_result = SingleCellSignalR_mouse_result,
#' SingleCellSignalR_sp2_result = SingleCellSignalR_zebrafish_result,
#' SingleCellSignalR_sp3_result = SingleCellSignalR_chick_result,
#' conserved_cell_types_df = conserved_cell_types_mm_zf_ch,
#' species_name1 = "mm",
#' species_name2 = "zf",
#' species_name3 = "ch")
conserved_interaction_score <- function(conserved_result_species,
                                        SingleCellSignalR_sp1_result,
                                        SingleCellSignalR_sp2_result,
                                        SingleCellSignalR_sp3_result,
                                        conserved_cell_types_df,
                                        species_name1 = "mm",
                                        species_name2 = "zf",
                                        species_name3 = "ch"
){
  cci_conserved_Weights_table_sp1_sp2_sp3 <- Caculate_cell_pair_cci_score2(conserved_result_df=conserved_result_species,
                                                                           species1_cci = SingleCellSignalR_sp1_result,
                                                                           species2_cci = SingleCellSignalR_sp2_result,
                                                                           species3_cci = SingleCellSignalR_sp3_result,
                                                                           conserved_cell_types_df = conserved_cell_types_df,
                                                                           species_name1 = species_name1,
                                                                           species_name2 = species_name2,
                                                                           species_name3 = species_name3)
  cci_data_sp1_sp2_sp3 <- Make_ChordDiagram_data2(cci_conserved_Weights_table = cci_conserved_Weights_table_sp1_sp2_sp3,
                                                  species_name1 = species_name1,
                                                  species_name2 = species_name2,
                                                  species_name3 = species_name3)
  return(list(cci_conserved_Weights = cci_conserved_Weights_table_sp1_sp2_sp3, chordDiagram_data = cci_data_sp1_sp2_sp3))
}
