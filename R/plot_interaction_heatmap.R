#' Heapmap of cell-cell interaction
#' @describeIn plot heapmap to show interaction between conserved cell types
#'
#' @param SingleCellSignalR_sp1_result dataframe, result of the SingCellSignalR, must contain "source", "target", "LRscore_scale"
#' @param SingleCellSignalR_sp2_result dataframe, result of the SingCellSignalR, must contain "source", "target", "LRscore_scale"
#' @param SingleCellSignalR_sp3_result dataframe, result of the SingCellSignalR, must contain "source", "target", "LRscore_scale"
#' @param value_name character, names of the sum weight, default is "scale_weight"
#' @param used_col default is c("source", "target", "LRscore_scale"), they should be in the column of result of SingCellSignalR
#' @param species_name1 character, species_name must be only two characters, like "mm"
#' @param species_name2 character, species_name must be only two characters, like "zf"
#' @param species_name3 character, species_name must be only two characters, like "ch"
#' @param sp1_conserved_celltype vector, a vector contained conserved cell type names in species1
#' @param sp2_conserved_celltype vector, a vector contained conserved cell type names in species2
#' @param sp3_conserved_celltype vector, a vector contained conserved cell type names in species3
#' @param celltype_colors vector, colors vector with names. names of celltype_colors must be paste0(species_name, celltype_name), like "mmCones"
#' @param brewer.pal_set character, pal in RColorBrewer
#' @param species_colors vector, colors vector without names. colors for species
#' @param ph_colors vector, colors for pheatmap
#' @param border_color, character, colors for border
#' @param scale, how to scale the data, one of "none", "row", "column"
#' @param cluster_rows logical, whether cluster the row
#' @param cluster_cols logical, whether cluster the column
#' @param legend logical, whether show the legend
#' @param show_rownames logical, whether show the rownames
#' @param show_colnames logical, whether show the colnames
#' @param fontsize numeric, size of the font
#' @param cellwidth numeric, width of the cell
#' @param cellheight numeric, heigth of the cell
#' @param angle_col numeric, angle of the colnames
#' @param annotation_names_col logical, whether show the name of column
#' @param annotation_names_row logical, whether show the name of row
#' @param filename the names used to save the heatmap
#' @param width numeric, width of the heatmap
#' @param height numeric, height of the heatmap
#' @param ... other parameters in pheatmap
#'
#' @return a list contain the data and the pheatmap object
#' @export
#'
#' @examples
#' mm_conserved_celltype <- c("Activated MG", "Cones", "GABAergic AC","Microglia", "Pericytes", "Resting MG", "RGC", "Rod BC", "Rods", "RPE", "V/E cells")
#' zf_conserved_celltype <- c("Activated MG",  "Cones", "GABAergic AC", "Microglia", "Pericytes", "Resting MG", "RGC","Cone BC", "Rods", "RPE", "V/E cells")
#' ch_conserved_celltype <- c("Activated MG", "Cones", "GABAergic AC", "Resting MG", "RGC", "Cone BC", "Rods")
#' plot_interaction_heatmap(SingleCellSignalR_sp1_result = SingleCellSignalR_mouse_result,
#'                          SingleCellSignalR_sp2_result = SingleCellSignalR_zebrafish_result,
#'                          SingleCellSignalR_sp3_result = SingleCellSignalR_chick_result,
#'                          species_name1 = "mm",
#'                          species_name2 = "zf",
#'                          species_name3 = "ch",
#'                          sp1_conserved_celltype = mm_conserved_celltype,
#'                          sp2_conserved_celltype = zf_conserved_celltype,
#'                          sp3_conserved_celltype = ch_conserved_celltype)
plot_interaction_heatmap <- function(SingleCellSignalR_sp1_result,
                                     SingleCellSignalR_sp2_result,
                                     SingleCellSignalR_sp3_result,
                                     value_name = "scale_weight",
                                     used_col = c("source", "target", "LRscore_scale"),
                                     species_name1 = "mm",
                                     species_name2 = "zf",
                                     species_name3 = "ch",
                                     sp1_conserved_celltype,
                                     sp2_conserved_celltype,
                                     sp3_conserved_celltype,
                                     celltype_colors = NULL,
                                     brewer.pal_set = "Set3",
                                     species_colors = NULL,
                                     ph_colors = viridis::viridis(8),
                                     border_color = "white",
                                     scale = "none",
                                     cluster_rows = FALSE,
                                     cluster_cols = FALSE,
                                     legend = TRUE,
                                     show_rownames = F,
                                     show_colnames = F,
                                     fontsize = 8,
                                     cellwidth = 20,
                                     cellheight = 24,
                                     angle_col = "45",
                                     annotation_names_col = TRUE,
                                     annotation_names_row = TRUE,
                                     filename = "CCI_pheatmap.pdf",
                                     width = 18,
                                     height = 10,
                                     ...
){
  all_weight_df_long_mm <- calculate_cell_pair_Weights(SingleCellSignalR_sp1_result[SingleCellSignalR_sp1_result$source %in% sp1_conserved_celltype & SingleCellSignalR_sp1_result$target %in% sp1_conserved_celltype,], used_col = used_col, value_name = value_name)
  all_weight_df_long_zf <- calculate_cell_pair_Weights(SingleCellSignalR_sp2_result[SingleCellSignalR_sp2_result$source %in% sp2_conserved_celltype & SingleCellSignalR_sp2_result$target %in% sp2_conserved_celltype,], used_col = used_col, value_name = value_name)
  all_weight_df_long_ch <- calculate_cell_pair_Weights(SingleCellSignalR_sp3_result[SingleCellSignalR_sp3_result$source %in% sp3_conserved_celltype & SingleCellSignalR_sp3_result$target %in% sp3_conserved_celltype,], used_col = used_col, value_name = value_name)
  subset_mm <- all_weight_df_long_mm %>% arrange(source) %>% arrange(target)
  subset_zf <- all_weight_df_long_zf %>% arrange(source) %>% arrange(target)
  subset_ch <- all_weight_df_long_ch %>% arrange(source) %>% arrange(target)
  subset_mm$species <- species_name1
  subset_zf$species <- species_name2
  subset_ch$species <- species_name3
  merge_subset <- data.frame()
  for (i in 1:max(nrow(subset_mm), nrow(subset_zf))) {
    if (i <= nrow(subset_mm)) {
      merge_subset <- rbind(merge_subset, subset_mm[i, ])
    }
    if (i <= nrow(subset_zf)) {
      merge_subset <- rbind(merge_subset, subset_zf[i, ])
    }
  }
  merge_subset$c_pair <- paste(merge_subset$source, merge_subset$target, sep = "_")
  subset_ch$c_pair <- paste(subset_ch$source, subset_ch$target, sep = "_")
  for (i in seq_len(nrow(subset_ch))) {
    id <- grep(as.character(subset_ch$c_pair[i]), merge_subset$c_pair)[length(grep(as.character(subset_ch$c_pair[i]), merge_subset$c_pair))]
    if (length(id)) {
      first_two_rows <- merge_subset[1:id, ]
      remaining_rows <- merge_subset[-(1:id), ]
      new_row <- subset_ch[i, ]
      merge_subset <- rbind(first_two_rows, new_row, remaining_rows)
    }else{
      new_row <- subset_ch[i, ]
      merge_subset <- rbind(merge_subset, new_row)
    }
  }
  merge_subset1 <- merge_subset %>% mutate(source = paste0(species, source))
  merge_subset_width_df <- merge_subset1[,c("source", "target", value_name)] %>%
    pivot_wider(
      names_from = "source",
      values_from = value_name
    )
  merge_subset_width_df[is.na(merge_subset_width_df)] <- 0
  make_rowname <- merge_subset_width_df$target
  merge_subset_width_df <- merge_subset_width_df[, -1]
  rownames(merge_subset_width_df) <- make_rowname
  annotate_col_df <- data.frame("source" = colnames(merge_subset_width_df))
  annotate_col_df$celltype <- substring(annotate_col_df$source, 3)
  annotate_col_df$species <- substr(annotate_col_df$source, 1, 2)
  make_rowname2 <- as.vector(unlist(annotate_col_df$source))
  annotate_col_df <- annotate_col_df[-1]
  rownames(annotate_col_df) <- make_rowname2
  colnames(annotate_col_df) <- c("source", "species")

  if (is.null(celltype_colors)) {
    celltype_colors <- colorRampPalette(RColorBrewer::brewer.pal(11, brewer.pal_set))(length(unique(annotate_col_df$source)))
    names(celltype_colors) <- unique(annotate_col_df$source)
  }
  if (is.null(species_colors)) {
    species_colors <-grDevices::colorRampPalette(c(rgb(0/255,114/255,189/255), rgb(248/255,115/255,106/255), rgb(169/255,169/255,169/255)))(length(unique(substr(colnames(CSCT), 1, 2))))
  }
  names(species_colors) <- c(species_name1, species_name2, species_name3)
  colors_list <- list(species = species_colors, source= celltype_colors)
  anno_row_df <- as.data.frame(rownames(merge_subset_width_df))
  colnames(anno_row_df) <- "target"
  rownames(anno_row_df) <- anno_row_df$target
  colors_row <- celltype_colors[anno_row_df$target]
  colors_row <- na.omit(colors_row)
  colors_list$target <- colors_row
  ph <- pheatmap::pheatmap(merge_subset_width_df,
                           color = ph_colors,
                           annotation_row = anno_row_df,
                           annotation_col = annotate_col_df,
                           annotation_colors = colors_list,
                           border_color = border_color,
                           scale = scale,
                           cluster_rows = cluster_rows,
                           cluster_cols = cluster_cols,
                           legend = legend,
                           show_rownames = show_rownames,
                           show_colnames = show_colnames,
                           fontsize = fontsize,
                           cellwidth = cellwidth,
                           cellheight = cellheight,
                           angle_col = angle_col,
                           annotation_names_col = annotation_names_col,
                           annotation_names_row = annotation_names_row,
                           filename = filename,
                           width = width,
                           height = height,
                           ...)
  return(list(ph_data=merge_subset_width_df, pheatmap=ph))
}
calculate_cell_pair_Weights <- function(data, used_col, value_name) {
  colnames(data) <- c("source", "target", "LRscore_scale")
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
