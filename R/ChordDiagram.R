#' Make a ChordDiagram
#' @description Construct a ChordDiagram to show the conservation score for intercellular interactions
#' @param net data frame, a data frame contain the source, target, and weight
#' @param Score_factor_size numeric, scale the weight score for showing
#' @param brewer_pal_used character, the brewer_pal in RColorBrewer, like "Set1",  "Set3", "Paired"
#' @param grid.col vector, color vectors with names. By default, it is NULL and it will asign colors automatically.
#' @param grid.border character, corlors for borders of grids. If it is NULL, the border color is same as grid color
#' @param is_link_colors_threshold logical, whether change the colors of links which below a threshold to a  specific color
#' @param link_colors_threshold numeric, the threshold used to change colors for links
#' @param link_colors_down_threshold_col character, specific corlor for links which below a threshold
#' @param filename character, save the plot with this filename
#' @param order_grid vector, a ordered vector for grid. By default, it is NULL
#' @param directional  Directions for the links. 1 means the direction is from the first column in df to the second column, -1 is the reverse, 0 is no direction, and 2 for two directional. The value can be a vector which has same length as number of rows in df
#' @param direction.type type for representing directions. Can be one or two values in "diffHeight" and "arrows"
#' @param diffHeight numeric, The difference of height between two 'roots' if directional is set to TRUE
#' @param annotationTrack annotation for the track. It can set to c("name","grid", "axis"), or one of them.
#' @param reduce numeric, if the ratio of the width of certain grid compared to the whole circle is less than this value, the grid is removed on the plot. Set it to value less than zero if you want to keep all tiny grid.
#' @param link.arr.type type for the arrow, "big.arrow" is set by default
#' @param link.border character, corlors for links
#' @param link.lwd numeric, width for link borders
#' @param link.sort logical, whether sort links on every sector based on the width of the links on it
#' @param link.decreasing logical, for link.sort
#' @param link.largest.ontop  logical, controls the order of adding links
#' @param link.overlap logical, whether the links that come or end in a same sector overlap
#' @param transparency numeric, transparency for link colors
#' @param cex numeric, font size for the text in annotationTrack
#' @param ylim_edit numeric, control the text height from the grid
#' @param facing control the direction of the text in annotationTrack, it can be one of c("inside", "outside", "reverse.clockwise", "clockwise", "downward", "bending", "bending.inside", "bending.outside")
#' @param niceFacing logical, adjusted the text to fit human eyes
#' @param col character, colors for the text in annotationTrack
#' @param adj offset for text, like c(0, 0.5). By default the text position adjustment is either horizontal or vertical in the canvas coordinate system.
#' @param start.degree numeric, rotation the ChordDiagram with this angle
#' @param picture_width numeric, width of the ChordDiagram for saving
#' @param picture_height numeric, height of the ChordDiagram for saving
#'
#' @return saving ChordDiagram in current directory
#' @export
#'
#' @examples
#' load(system.file("extdata", "CCC_conserved_sumary.rda", package = "CACIMAR"))
#' cci_data = CCC_conserved_sumary[, c("Source", "Target", "score_weight")]
#' ChordDiagram(net = cci_data, filename = "chordDiagram_cutoff0.75.pdf", link_colors_threshold = 0.75)
#'
ChordDiagram <- function(net,
                         Score_factor_size = 20,
                         brewer_pal_used = "Set3",
                         grid.col = NULL,
                         grid.border = "#F8F8F8",
                         is_link_colors_threshold = TRUE,
                         link_colors_threshold = 0.5,
                         link_colors_down_threshold_col = "grey",
                         filename = "chordDiagram.pdf",
                         order_grid = NULL,
                         directional = 1,
                         direction.type = c("arrows"),  # "diffHeight"
                         diffHeight = -0.03,
                         annotationTrack = "grid",
                         reduce = -1,
                         link.arr.type = "big.arrow",
                         link.border = "#F8F8F8",
                         link.lwd = 0.6,
                         link.sort = TRUE,
                         link.decreasing = TRUE,
                         link.largest.ontop = TRUE,
                         link.overlap = F,
                         transparency = 0,
                         cex = 1.3,
                         ylim_edit = 0,
                         facing = "clockwise",
                         niceFacing = TRUE,
                         col = "black",
                         adj = c(0, 0.5),
                         start.degree = 70,
                         picture_width = 7,
                         picture_height = 9
) {
  # link value
  # net <- as_tibble(net)
  colnames(net) <- c("source", "target", "Score")
  net$Score <- 2^(net$Score*Score_factor_size)

  if (is.null(grid.col)) {
    grid.col <- colorRampPalette(RColorBrewer::brewer.pal(11, brewer_pal_used))(length(c(unique(net$source), unique(net$target))))
    grid.col <- setNames(grid.col, c(unique(net$source), unique(net$target)))
  } else {
    grid.col <- grid.col
  }

  # grid colors
  grid.cols_df <- as.data.frame(grid.col)
  colnames(grid.cols_df) <- c("colors")
  grid.cols_df$celltype <- names(grid.col)

  # link colors
  net$colors <- grid.cols_df$colors[match(net$source, grid.cols_df$celltype)]
  net$link_colors <- net$colors

  if (is_link_colors_threshold) {
    Score <- as.vector(unlist(net$Score))
    median_value <- quantile(Score, probs = link_colors_threshold)
    net$link_colors[which(net$Score < median_value)] <- link_colors_down_threshold_col
  }

  pdf(file = filename, width = picture_width, height = picture_height)
  circlize::circos.par(start.degree = start.degree)
  circlize::chordDiagram(net,
                         grid.col = grid.col,
                         grid.border = grid.border,
                         col = net$link_colors,
                         order = order_grid,
                         directional = directional,
                         direction.type = direction.type,  # "diffHeight"
                         diffHeight = diffHeight,
                         annotationTrack = annotationTrack,
                         reduce = reduce,
                         link.arr.type = link.arr.type,
                         link.border = link.border,
                         link.lwd = link.lwd,
                         link.sort = link.sort,
                         link.decreasing = link.decreasing,
                         link.largest.ontop = link.largest.ontop,
                         link.overlap = link.overlap,
                         transparency = transparency,
                         preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(net)))))
  )
  # text
  for(si in circlize::get.all.sector.index()) {
    xlim = circlize::get.cell.meta.data("xcenter", sector.index = si, track.index = 1)
    ylim = circlize::get.cell.meta.data("ylim", sector.index = si, track.index = 1)
    circlize::circos.text(xlim,
                          ylim+ylim_edit,
                          si,
                          sector.index = si,
                          track.index = 1,
                          cex = cex,
                          facing = facing,
                          niceFacing = niceFacing,
                          col = col,
                          adj = adj)
  }
  dev.off()
}

