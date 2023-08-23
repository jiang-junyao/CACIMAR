#' Plot communication between cell types
#'
#' @description  This function takes a dataframe containing information about communication between cell types,
#' and plots a graph to visualize the communication pattern. The graph represents sender and
#' receiver cell types as nodes, and the strength of communication as the width of edges connecting them.
#'
#' @param graph_df The input dataframe containing communication data, see the example in tutorial.
#' @param sendercolumn single column name, The name of the column containing sender cell types.
#' @param receivercolumn single column name, The name of the column containing receiver cell types.
#' @param widthcolumn single column name, The name of the column containing the width of edges (communication strength).
#' @param colors Optional vector of colors for nodes and edges.
#' @param useLabels Logical value indicating whether to display labels for nodes.
#' @param nodeSize Size of the nodes.
#' @param vertex.label.color Color of node labels.
#' @param vertex.label.cex Size of node labels.
#' @param ...
#'
#' @return The plotted graph as an invisible object
#' @import igraph
#' @importFrom igraph graph_from_data_frame layout.circle
#' @export
#'
#' @examples
#' plot_graph_network(long_df, sendercolumn = "celltype1", receivercolumn = "celltype2", widthcolumn = "score", useLabels = T, colors = colors, nodeSize = 10, vertex.label.color = "black", vertex.label.cex = 1.5)
#'
Plot_Celltype.Communication <- function(graph_df, sendercolumn, receivercolumn, widthcolumn, colors = NULL, useLabels = TRUE, nodeSize = 5, vertex.label.color, edge.label.color, vertex.label.cex, ...) {
  graph_df$sender <- as.character(graph_df[[sendercolumn]])
  graph_df$receiver <- as.character(graph_df[[receivercolumn]])
  graph_df$width <- graph_df[[widthcolumn]]

  vertices <- data.frame(id = unique(graph_df$sender), size = nodeSize,
                         stringsAsFactors = FALSE, label = unique(graph_df$sender))

  if (!is.null(colors)) {
    vertices$label.color <- colors[vertices$id]
    vertices$color <- colors[vertices$id]
    graph_df$color <- colors[graph_df$sender]
  }

  plot_graph <- igraph::graph_from_data_frame(graph_df, vertices = vertices)

  layout <- layout.circle(plot_graph)

  if (useLabels) {
    plot(plot_graph, layout = layout, vertex.frame.color = NA, edge.arrow.mode = "-", edge.curved = 0.3,
         edge.label = NULL,
         edge.alpha = 0.7,
         vertex.label = vertices$label,
         vertex.label.color = vertex.label.color,
         vertex.label.cex = vertex.label.cex,
         ...)
  } else {
    plot(plot_graph, layout = layout, vertex.frame.color = NA, edge.arrow.mode = "-", vertex.label = NA,
         edge.curved = 0.3,
         edge.label = NULL,
         edge.alpha = 0.7,
         vertex.label = vertices$label,
         vertex.label.color = vertex.label.color,
         vertex.label.cex = vertex.label.cex,
         ...)
  }
  invisible(plot_graph)
}
