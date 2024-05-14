#' Make ligand-receptor data for pheatmap
#' @description prepared average expression data for pheatmap
#'
#' @param cci_conserved_results dataframe, the result from function Identify_Conserved_CCI2
#' @param seurat_object_sp1 seurat object of species1
#' @param seurat_object_sp2 seurat object of species2
#' @param seurat_object_sp3 seurat object of species3
#' @param species_name1 character, two character the representation the species1, like "mm"
#' @param species_name2 character, two character the representation the species2, like "zf"
#' @param species_name3 character, two character the representation the species2, like "ch"
#' @param avg_group character, which group used to calculate the average expression, like "celltype" in metadata
#' @param cci_conserved_Weights dataframe, the result from function Caculate_cell_pair_cci_score
#' @param subset_quantile numeric, the value to filter the conserved score of the cell-cell interaction, default is 0.75 of the consered score of cell-cell interaction
#' @param assay_set which assay used for average expression, like "RNA", or "SCT"
#'
#' @return list, contains the dataframe for heatmap, and the average expression for each species
#' @export
#'
#' @examples
#' merge_avg_width_df <- make_pheatmap_LR_data(cci_conserved_results = conserved_result_mm_zf_ch,
#' seurat_object_sp1 = Mm_seurat,
#' seurat_object_sp2 = Zf_seurat,
#' seurat_object_sp3 = ch_seurat,
#' avg_group = "cellname",
#' species_name1 = "mm",
#' species_name2 = "zf",
#' species_name3 = "ch",
#' cci_conserved_Weights = conserved_cci_result$cci_conserved_Weights)
make_pheatmap_LR_data <- function(cci_conserved_results,
                                  seurat_object_sp1,
                                  seurat_object_sp2,
                                  seurat_object_sp3,
                                  species_name1 = "mm",
                                  species_name2 = "zf",
                                  species_name3 = "ch",
                                  avg_group = "cellname",
                                  cci_conserved_Weights,
                                  subset_quantile = 0.75,
                                  assay_set = "RNA"
){
  specie_conserved_CCC_sp1 <- get_average_expression(cci_conserved_results$sp1_cci_conserved_df, seurat_object = seurat_object_sp1, species_set = species_name1, avg_group = avg_group, assay_set = assay_set)
  specie_conserved_CCC_sp2 <- get_average_expression(cci_conserved_results$sp2_cci_conserved_df, seurat_object = seurat_object_sp2, species_set = species_name2, avg_group = avg_group, assay_set = assay_set)
  specie_conserved_CCC_sp3 <- get_average_expression(cci_conserved_results$sp3_cci_conserved_df, seurat_object = seurat_object_sp3, species_set = species_name3, avg_group = avg_group, assay_set = assay_set)

  cci_conserved_Weights_q0.75 <- cci_conserved_Weights[!cci_conserved_Weights$conserved_score < quantile(cci_conserved_Weights$conserved_score, subset_quantile),]
  cci_conserved_Weights_q0.75_c_pair_sp1 <- paste(cci_conserved_Weights_q0.75$source, cci_conserved_Weights_q0.75$target, sep = "_")
  cci_conserved_Weights_q0.75_c_pair_sp2 <- paste(cci_conserved_Weights_q0.75$sp2_source, cci_conserved_Weights_q0.75$sp2_target, sep = "_")
  cci_conserved_Weights_q0.75_c_pair_sp3 <- paste(cci_conserved_Weights_q0.75$sp3_source, cci_conserved_Weights_q0.75$sp3_target, sep = "_")
  specie_conserved_CCC_sp1_q <- subset(specie_conserved_CCC_sp1, c_pair %in% cci_conserved_Weights_q0.75_c_pair_sp1)
  specie_conserved_CCC_sp2_q <- subset(specie_conserved_CCC_sp2, c_pair %in% cci_conserved_Weights_q0.75_c_pair_sp2)
  specie_conserved_CCC_sp3_q <- subset(specie_conserved_CCC_sp3, c_pair %in% cci_conserved_Weights_q0.75_c_pair_sp3)
  specie_conserved_CCC_sp1_q <- specie_conserved_CCC_sp1_q[!duplicated(specie_conserved_CCC_sp1_q), ]
  specie_conserved_CCC_sp2_q <- specie_conserved_CCC_sp2_q[!duplicated(specie_conserved_CCC_sp2_q), ]
  specie_conserved_CCC_sp3_q <- specie_conserved_CCC_sp3_q[!duplicated(specie_conserved_CCC_sp3_q), ]

  specie_conserved_CCC_sp1_q <- specie_conserved_CCC_sp1_q[, c("c_pair", "lr_p", "ligand_receptor_avg", "species")]
  specie_conserved_CCC_sp2_q <- specie_conserved_CCC_sp2_q[, c("c_pair", "lr_p", "ligand_receptor_avg", "species")]
  specie_conserved_CCC_sp3_q <- specie_conserved_CCC_sp3_q[, c("c_pair", "lr_p", "ligand_receptor_avg", "species")]
  colnames(specie_conserved_CCC_sp1_q) <- c("Cell_Pair", "LR", "ligand_receptor_avg", "species")
  colnames(specie_conserved_CCC_sp2_q) <- c("Cell_Pair", "LR", "ligand_receptor_avg", "species")
  colnames(specie_conserved_CCC_sp3_q) <- c("Cell_Pair", "LR", "ligand_receptor_avg", "species")
  merge_width_df <- Merge_data_for_heatmap(specie_conserved_CCC_sp1_q, specie_conserved_CCC_sp2_q)
  sp3_s_id <- seq(3, nrow(specie_conserved_CCC_sp1_q)*2, 3)
  for (i in seq_len(nrow(specie_conserved_CCC_sp3_q))) {
    first_two_rows <- merge_width_df[1:sp3_s_id[i]-1, ]
    remaining_rows <- merge_width_df[-(1:sp3_s_id[i]-1), ]
    new_row <- specie_conserved_CCC_sp3_q[i, ]
    merge_width_df <- rbind(first_two_rows, new_row, remaining_rows)
  }
  write.table(merge_width_df, "lr_avg_df_three_species.txt")

  merge_width_df1 <- merge_width_df
  merge_width_df1$LR <- rep(merge_width_df$LR[seq(1,nrow(merge_width_df), 3)], each = 3)
  merge_width_df_width <- merge_width_df1[,c("LR", "Cell_Pair", "ligand_receptor_avg")] %>%
    tidyr::pivot_wider(
      names_from = "LR",
      values_from = "ligand_receptor_avg"
    )
  make_rowname_lr2 <- merge_width_df_width$Cell_Pair
  merge_width_df_width <- merge_width_df_width[, -1]
  rownames(merge_width_df_width) <- make_rowname_lr2
  return(list(merge_width_df=merge_width_df_width, specie_conserved_CCC_sp1=specie_conserved_CCC_sp1, specie_conserved_CCC_sp2=specie_conserved_CCC_sp1, specie_conserved_CCC_sp3=specie_conserved_CCC_sp3))
}
get_average_expression <- function(specie_conserved_CCC, seurat_object, avg_group, species_set, assay_set) {
  specie_conserved_CCC$source <- sub(species_set, "", specie_conserved_CCC$source)
  specie_conserved_CCC$target <- sub(species_set, "", specie_conserved_CCC$target)
  unique_LR <- unique(c(specie_conserved_CCC$ligand, specie_conserved_CCC$receptor))
  unique_LR_avg <- as.data.frame(Seurat::AverageExpression(seurat_object, features = unique_LR, group.by = avg_group, assays = assay_set))
  colnames(unique_LR_avg) <- gsub(paste0(assay_set, "\\."), "", colnames(unique_LR_avg))
  colnames(unique_LR_avg) <- gsub("\\.", " ", colnames(unique_LR_avg))
  for (i in 1:nrow(specie_conserved_CCC)) {
    l_avg <- unique_LR_avg[match(specie_conserved_CCC$ligand[i], rownames(unique_LR_avg)), match(specie_conserved_CCC$source[i], colnames(unique_LR_avg))]
    r_avg <- unique_LR_avg[match(specie_conserved_CCC$receptor[i], rownames(unique_LR_avg)), match(specie_conserved_CCC$target[i], colnames(unique_LR_avg))]
    specie_conserved_CCC$ligand_source_avg[i] <- l_avg
    specie_conserved_CCC$receptor_target_avg[i] <- r_avg
    specie_conserved_CCC$ligand_receptor_avg[i] <- sqrt(l_avg * r_avg)
  }
  specie_conserved_CCC$species <- species_set
  return(specie_conserved_CCC)
}
Merge_data_for_heatmap <- function(LRpair_show_df1, LRpair_show_df2) {
  library(tidyverse)
  merge_LRpair_show <- data.frame()
  for (i in 1:max(nrow(LRpair_show_df1), nrow(LRpair_show_df2))) {
    if (i <= nrow(LRpair_show_df1)) {
      merge_LRpair_show <- rbind(merge_LRpair_show, LRpair_show_df1[i, ])
    }
    if (i <= nrow(LRpair_show_df2)) {
      merge_LRpair_show <- rbind(merge_LRpair_show, LRpair_show_df2[i, ])
    }
  }
  return(merge_LRpair_show)
}
#' Create a pheatmap for average expression of ligand-receptor pairs
#' @description make a pheatmap for average expression of ligand-receptor pairs using geometric mean of corresponding cell type
#'
#' @param width_df width dataframe, result of function make_pheatmap_LR_data
#' @param filename name to save the pheatmap, default is "pheatmap_LR.pdf"
#' @param color_pheatmap_set colors for the pheatmap, default is viridis::viridis(8)
#' @param species_colors vector, colors vectors for the species
#' @param species_name1 two character, like "mm", for species1
#' @param species_name2 two character, like "zf", for species2
#' @param species_name3 two character, like "ch", for species3
#' @param annotation_names_col logical, whether show the column title name
#' @param border_color, character, for the color of the border, like "white"
#' @param scale how to scale the data, one of "none", "row", "column"
#' @param cluster_rows logical, whether cluster the row
#' @param cluster_cols logical, whether cluster the column
#' @param legend logical, whether show the legend
#' @param legend_breaks break legend with this vector value
#' @param legend_labels labels for the break legend
#' @param show_rownames logical, whether show the rownames
#' @param show_colnames logical, whether show the colnames
#' @param fontsize numeric, size of the font
#' @param cellwidth numeric, width of the cell
#' @param cellheight numeric, heigth of the cell
#' @param width numeric, width of the heatmap
#' @param height numeric, height of the heatmap
#' @param ... other parameters in pheatmap
#'
#' @return pheatmap
#' @export
#'
#' @examples
#' merge_avg_width_df <- make_pheatmap_LR_data(cci_conserved_results = conserved_result_mm_zf_ch,
#' seurat_object_sp1 = Mm_seurat,
#' seurat_object_sp2 = Zf_seurat,
#' seurat_object_sp3 = ch_seurat,
#' avg_group = "cellname",
#' species_name1 = "mm",
#' species_name2 = "zf",
#' species_name3 = "ch",
#' cci_conserved_Weights = conserved_cci_result$cci_conserved_Weights)
#' pheatmap_LR_multi(width_df = merge_avg_width_df$merge_width_df, filename = "lr_avg_pheatmap.pdf")
pheatmap_LR_multi <- function(width_df,
                              filename = "pheatmap_LR.pdf",
                              color_pheatmap_set = NULL,
                              species_colors = NULL,
                              species_name1 = "mm",
                              species_name2 = "zf",
                              species_name3 = "ch",
                              annotation_names_col = TRUE,
                              border_color = "white",
                              scale = "none",
                              cluster_rows = FALSE,
                              cluster_cols = FALSE,
                              legend = TRUE,
                              legend_breaks = NA,
                              legend_labels = NA,
                              show_rownames = T,
                              show_colnames = T,
                              fontsize = 10,
                              cellwidth = 16,
                              cellheight = 12,
                              width = 14,
                              height = 14,
                              ...) {
  anno_row_df <- data.frame("source_target" = rownames(width_df), "species" = substr(rownames(width_df), start = 1, stop = 2))
  rownames(anno_row_df) <- anno_row_df$source_target
  anno_row_df <- anno_row_df[-1]
  LRpair_show_rename_row <- substr(rownames(width_df), start = 3, stop = nchar(rownames(width_df)))
  if (is.null(color_pheatmap_set)) {
    color_pheatmap = viridis::viridis(8)
  } else {
    color_pheatmap = color_pheatmap_set
  }
  if (is.null(species_colors)) {
    species_colors <-c(rgb(0/255,114/255,189/255), rgb(248/255,115/255,106/255), rgb(169/255,169/255,169/255))
  }
  names(species_colors) <- c(species_name1, species_name2, species_name3)
  species_list <- list(species = species_colors)
  pheatmap::pheatmap(width_df,
                     color = color_pheatmap,
                     annotation_row =  anno_row_df,
                     annotation_colors = species_list,
                     annotation_names_col = annotation_names_col,
                     labels_row = LRpair_show_rename_row,
                     border_color = border_color,
                     scale = scale,
                     cluster_rows = cluster_rows,
                     cluster_cols = cluster_cols,
                     legend = legend,
                     legend_breaks = legend_breaks,
                     legend_labels = legend_labels,
                     show_rownames = show_rownames,
                     show_colnames = show_colnames,
                     fontsize = fontsize,
                     cellwidth = cellwidth,
                     cellheight = cellheight,
                     filename = filename,
                     width = width,
                     height = height,
                     ...)
}
