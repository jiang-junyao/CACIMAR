#' Create sankey plot
#' @description Build a sankey plot to show the cell-cell interaction profile of two species
#' @param links dataframe, a links file containing source, target, and weight
#' @param specie_name1 character, name for specie 1
#' @param specie_name2 character, name for specie 2
#' @param output_file character, save sankey plot with this file name
#' @param colors_file character, file path for the color file, it is NULL by default and it will assign colors automatically
#' @param brewer.pal_set if colors_file is null, the brewer.pal can be choose, like "Set1", "Set3", "Paired", ......
#' @param node_width numeric, width of the node
#' @param node_padding numeric, distance between the node
#' @param node_stroke_width numeric, width of the stroke around nodes
#' @param node_corner_radius numeric, the radius of the node
#' @param drag_x logical, TRUE indicates that the plot can be horizontally dragged
#' @param drag_y logical, TRUE indicates that the plot can be vertically dragged
#' @param units character, name for units
#' @param node_pos_x character, variable used for grouping nodes on the x-axis
#' @param align character, alignment of the nodes. One of 'right', 'left', 'justify', 'center', 'none'. If 'none', then the labels of the nodes are always to the right of the node
#' @param scale_node_breadths_by_string logical, Put nodes at positions relatively to string lengths - only work well currently with align='none'
#' @param show_node_values logical, Show values above nodes. Might require and increased node margin.
#' @param height numeric, height of the sankey plot for saving
#' @param width numeric, width of the sankey plot for saving
#' @param linkColor numeric, color of links
#' @param link_type character, one of 'bezier', 'l-bezier', 'trapezoid', 'path1' and 'path2'
#' @param curvature numeric, curvature parameter for bezier links - between 0 and 1
#' @param link_opacity numeric, opacity of links
#' @param link_gradient logical, add a gradient to the links
#' @param node_shadow logical, add a shadow to the nodes
#' @param node_label_margin numeric, distance between the node and font
#' @param zoom logical, value to enable (TRUE) or disable (FALSE) zooming
#' @param iterations numeric, number of iterations in the diagramm layout for computation of the depth (y-position) of each node. Note: this runs in the browser on the client so don't push it too high
#' @param x_scaling_factor numeric, scale the computed x position of the nodes by this value
#' @param font_size  numeric, size of the font of the plot
#' @param ... param of sankeyD3::sankeyNetwork()
#'
#' @return a sankey plot
#' @export
#'
#' @examples
#' create_sankey(links = all_weight_df_long[, c("Source2", "target", "scale_weight")],
#' output_file = "sankey_scale_weight.html",
#' specie_name1 = "Mm",
#' specie_name2 = "Zf",
#' colors_file = NULL)
#'
create_sankey <- function(links,
                          specie_name1 = "Mm",
                          specie_name2 = "Zf",
                          output_file = "sankey.html",
                          colors_file = NULL,
                          brewer.pal_set = "Set3",
                          node_width = 18,
                          node_padding = 15,
                          node_stroke_width = 0,
                          node_corner_radius = 0,
                          drag_x = TRUE,
                          drag_y = TRUE,
                          units = "TWh",
                          node_pos_x = "group",
                          align = "none",
                          scale_node_breadths_by_string = F,
                          show_node_values = FALSE,
                          height = 1400,
                          width = 1800,
                          linkColor = "#A0A0A0",
                          link_type = "bezier",
                          curvature = 0.5,
                          link_opacity = 5,
                          link_gradient = FALSE,
                          node_shadow = FALSE,
                          node_label_margin = 5,
                          zoom = TRUE,
                          font_family = NULL,
                          iterations = 0,
                          x_scaling_factor = 1.4,
                          font_size = 16,
                          ...
){
  # links
  colnames(links) <- c("Source", "target", "weight")
  # nodes
  library(tidyverse)
  nodes <- data.frame(name=c(as.character(links$Source), as.character(links$target)) %>% unique())
  # groups
  nodes$group <- c(
    rep(0, sum(grepl(specie_name1, nodes$name))),
    rep(1, (length(nodes$name) - sum(sum(grepl(specie_name1, nodes$name)), sum(grepl(specie_name2, nodes$name))))),
    rep(2, sum(grepl(specie_name2, nodes$name)))
  )

  # node name
  nodes$name_node <- sub(paste("^", specie_name1, sep = ""), "", nodes$name)
  nodes$name_node <- sub(paste("^", specie_name2, sep = ""), "", nodes$name_node)

  # node colors
  if (is.null(colors_file)) {
    colors <- colorRampPalette(RColorBrewer::brewer.pal(11, brewer.pal_set))(length(unique(nodes$name_node)))
    colors <- as.data.frame(colors)
    colnames(colors) <- "colors"
    colors$celltype <- unique(nodes$name_node)
    nodes$colors <- colors$colors[match(nodes$name_node, colors$celltype)]
  } else {
    colors <- colors_file
    nodes$colors <- colors$colors[match(nodes$name_node, colors$celltype)]
  }

  # ID start from 0
  links$IDsource <- match(links$Source, nodes$name) - 1
  links$IDtarget <- match(links$target, nodes$name) - 1
  colnames(links) <- c("source", "target", "weight", "IDsource", "IDtarget")

  #Plot
  p <- sankeyD3::sankeyNetwork(Links = links,
                               Nodes = nodes,
                               Source = "IDsource",
                               Target = "IDtarget",
                               Value = "weight",
                               NodeID = "name_node",
                               nodeWidth = node_width,
                               nodePadding = node_padding,
                               nodeStrokeWidth = node_stroke_width,
                               nodeCornerRadius = node_corner_radius,
                               NodeColor = "colors",
                               dragX = drag_x,
                               dragY = drag_y,
                               units = units,
                               NodePosX = node_pos_x,
                               align = align,
                               scaleNodeBreadthsByString = scale_node_breadths_by_string,
                               showNodeValues = show_node_values,
                               height = height,
                               width = width,
                               linkColor = linkColor,
                               linkType = link_type,
                               curvature = curvature,
                               linkOpacity = link_opacity,
                               linkGradient = link_gradient,
                               nodeShadow = node_shadow,
                               nodeLabelMargin = node_label_margin,
                               zoom = zoom,
                               iterations = iterations,
                               xScalingFactor = x_scaling_factor,
                               fontSize = font_size)
  # save
  sankeyD3::saveNetwork(p, output_file)
  print(p)
  webshot::webshot(output_file, paste(output_file, ".pdf", sep = ""))
  return(p)
}

#' Cell-cell interaction analysis with SingleCellSignalR algorithm
#' @description Perform cell-cell interaction with SingleCellSignalR algorithm contained in liana package.
#' @param seurat_obj seurat object, a seurat object with cell types in active idents
#' @param expre_cell numeric, filter genes expressed in more than expre_cell cells
#' @param select_DB character, databased used for analysis, more option can run liana::show_resources()
#' @param target_organism ncbi_taxid' or 'name' of the target organism. See ‘show_homologene' for available organisms via OmnipathR’s 'HomoloGene'
#' @param method character, method use for cell-cell interaction analysis. By default, we set it to SingleCellSignalR
#'
#' @return data frame of cell cell interaction analysis result
#' @export
#'
#' @examples
#'# Mouse analysis
#'SingleCellSignalR_mouse_result <- perform_CCI_analysis(seurat_obj=Mm_seurat, target_organism=10090)
#'
# Zebrafish analysis
#'SingleCellSignalR_zebrafish_result <- perform_CCI_analysis(seurat_obj=Zf_seurat, target_organism=7955)
perform_CCI_analysis <- function(seurat_obj, expre_cell = 10, select_DB = "Consensus", target_organism, method = 'sca') {
  seurat_obj <- seurat_obj[apply(GetAssayData(seurat_obj, slot = "counts"), 1, function(x) sum(x > 0) > expre_cell),]
  op_resource <- liana::select_resource(select_DB)[[1]]
  ortholog_resource <- liana::generate_homologs(op_resource = op_resource, target_organism = target_organism)
  liana_result <- liana::liana_wrap(seurat_obj,
                                    method = method,
                                    resource = 'custom',
                                    external_resource = ortholog_resource)
  return(liana_result)
}

#'  Sum weight of cell-cell interactions
#' @description Here calculate overall weight of cell-cell interactions within each pair of cell types
#' @param species1_cci data frame,  cell-cell interactions result with perform_CCI_analysis of species 1
#' @param species2_cci data frame,  cell-cell interactions result with perform_CCI_analysis of species 2
#' @param specie_name1 character, name for species 1, like "Mm"
#' @param specie_name2 character, name for species 2, like "Zf"
#'
#' @return data frame of weight result
#' @export
#'
#' @examples
#' all_weight_df_long <- calculate_Weights(species1_cci = SingleCellSignalR_mouse_result,
#' species2_cci = SingleCellSignalR_zebrafish_result)
#'
#' head(all_weight_df_long)
#'             source       target    weight species scale_weight             Source            Source2
#' 1     Activated MG Activated MG 152.07967      Mm  0.008256281     MmActivated MG     MmActivated MG
#' 2       Astrocytes Activated MG 230.97254      Mm  0.012539310       MmAstrocytes       MmAstrocytes
#' 3            Cones Activated MG  51.89929      Mm  0.002817570            MmCones            MmCones
#' 4     GABAergic AC Activated MG  99.92955      Mm  0.005425093     MmGABAergic AC     MmGABAergic AC
#' 5   Glycinergic AC Activated MG  78.95619      Mm  0.004286467   MmGlycinergic AC   MmGlycinergic AC
#' 6 Horizontal cells Activated MG 106.23687      Mm  0.005767513 MmHorizontal cells MmHorizontal cells
calculate_Weights <- function(
    species1_cci,
    species2_cci,
    specie_name1 ="Mm",
    specie_name2 = "Zf"){
  species1_all_weight_df_long <- calculateWeights(data = species1_cci,
                                                  specie_name = specie_name1)
  species2_all_weight_df_long <- calculateWeights(data = species2_cci,
                                                  specie_name = specie_name2)
  # bind df
  all_weight_df_long <- rbind(species1_all_weight_df_long, species2_all_weight_df_long)
  all_weight_df_long <- subset(all_weight_df_long, weight > 0)
  all_weight_df_long$Source <- paste0(all_weight_df_long$species, all_weight_df_long$source)
  # order for sankey plot
  all_weight_df_long$target <- as.character(all_weight_df_long$target)
  all_weight_df_long$Source2 <- all_weight_df_long$Source
  edit_target <- all_weight_df_long$Source[which(all_weight_df_long$species == specie_name2)]
  edit_source <- all_weight_df_long$target[which(all_weight_df_long$species == specie_name2)]
  all_weight_df_long$Source2[which(all_weight_df_long$species == specie_name2)] <- edit_source
  all_weight_df_long$target[which(all_weight_df_long$species == specie_name2)] <- edit_target
  return(all_weight_df_long)
}

calculateWeights <- function(data, specie_name) {
  ccc_weight <- reshape2::dcast(data[, c(2, 3, 11)],
                                source ~ target,
                                sum,
                                value.var = 'LRscore')
  rownames(ccc_weight) <- ccc_weight[, 1]
  all_weight_df_long <- reshape2::melt(ccc_weight,
                                       id.vars = "source",
                                       variable.name = "target",
                                       value.name = "weight")
  all_weight_df_long$species <- specie_name
  all_weight_df_long$scale_weight <- all_weight_df_long$weight/sum(all_weight_df_long$weight)

  return(all_weight_df_long)
}
