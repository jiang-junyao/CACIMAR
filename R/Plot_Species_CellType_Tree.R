#' @title Plot Phylogenetic Tree of Different Species Cell Types
#' @description Plot a tree of different species cell types using the conserved score of cell types
#'
#' @param dist_matrix Distance matrix
#' @param hcluster.method chariter, Hierarchical clustering method, possible values are "average" (default), "single", "complete"
#' @param species.vector vector, Species vector
#' @param layout.tree Tree layout ("rectangular", "dendrogram", "fan", "circular"), default is "rectangular"
#' @param offset Tiplab offset, default is 0.01
#' @param tiplab.size Tiplab size, default is 5
#' @param tippoint.shape Tippoint shape,  default is 22
#' @param tippoint.shape.size Tippoint size, default is 4
#' @param geom_nodepoint Nodepoint size, default is 0
#' @param col.value Color vector for species, default is c(rgb(102/255,46/255,115/255), rgb(31/255,153/255,139/255))
#'
#' @return a list, which contain the Tree data and the Tree plot
#' @export
#'
#' @examples load(system.file("extdata", "zf_mm_markers.rda", package = "CACIMAR"))
#' expression <- Identify_ConservedCellTypes(OrthG_Mm_Zf,Zf_marker,Mm_marker,'zf','mm')
#' dist_matrix <- expression[[2]]
#' species.vector <- substr(rownames(dist_matrix), 1, 2)
#' Species_CellType_Tree <- Plot_Species_CellType_Tree(dist_matrix = dist_matrix, species.vector = species.vector, hcluster.method = "average", geom_nodepoint = 0, layout.tree = "rectangular", offset = 0.02, tiplab.size = 5, tippoint.shape = 22, tippoint.shape.size = 4)
#' Species_CellType_Tree$Plot
Plot_Species_CellType_Tree <- function(dist_matrix, hcluster.method = c("average", "single", "complete"),
                                       species.vector,
                                       layout.tree = c('rectangular', 'dendrogram', 'fan', 'circular'),
                                       offset = 0.01, tiplab.size = 5, tippoint.shape = c(22, 21),
                                       tippoint.shape.size = 4, geom_nodepoint = 0,
                                       col.value = c(rgb(102/255,46/255,115/255), rgb(31/255,153/255,139/255))) {
  stopifnot(is.matrix(dist_matrix))
  stopifnot(length(hcluster.method) == 1 && hcluster.method %in% c("average", "single", "complete"))
  stopifnot(is.vector(species.vector))
  stopifnot(length(layout.tree) == 1 && layout.tree %in% c('rectangular', 'dendrogram', 'fan', 'circular'))
  stopifnot(is.numeric(offset) && offset >= 0)
  stopifnot(is.numeric(tiplab.size) && tiplab.size > 0)
  stopifnot(is.numeric(tippoint.shape) && length(tippoint.shape) == 1)
  stopifnot(is.numeric(tippoint.shape.size) && tippoint.shape.size > 0)
  stopifnot(is.numeric(geom_nodepoint) && geom_nodepoint >= 0)
  stopifnot(length(col.value) >= 2 && is.character(col.value))


  # Perform hierarchical clustering
  hclust_tree <- hclust(as.dist(1 - dist_matrix), method = hcluster.method)

  # Convert the clustering result into a phylogenetic tree object
  tree <- ape::as.phylo(hclust_tree)
  tree$species <- species.vector

  # Generate the plot of the cell type tree
  tree_plot <- Plot_CellType_Tree(tree, layout.tree, offset, tiplab.size, tippoint.shape,
                                  tippoint.shape.size, geom_nodepoint = geom_nodepoint, col.value = col.value)

  # Return the tree object and the plot in a list
  return(list("Tree" = tree, "Plot" = tree_plot))
}

Plot_CellType_Tree <- function(tree, layout.tree, geom_nodepoint, col.value, offset, tiplab.size, tippoint.shape, tippoint.shape.size){
  tree.p <- ggtree::ggtree(tree, layout=layout.tree, size=0.7, col="black")
  data <- fortify(tree)
  tree.df <- data.frame("tip.label" = tree$tip.label, "species" = tree$species)
  tregraph <- tree.p %<+% tree.df +
    geom_tiplab(aes(col = species), offset = offset, size = tiplab.size) +
    geom_tippoint(aes(col = species, fill = species), shape = tippoint.shape, size = tippoint.shape.size) +
    geom_nodepoint(size = geom_nodepoint) +
    theme_tree2() +
    xlim(NA, max(data$x) * 1.3) +
    scale_fill_manual(values = col.value) +
    scale_color_manual(values = col.value) +
    guides(shape = ggplot2::guide_legend(override.aes = list(shape = c(tippoint.shape, tippoint.shape)))) +
    theme(legend.text = element_text(size = 16),  # 设置图例文本大小为12
          legend.title = element_text(size = 16))  # 设置图例标题大小为14

  return(tregraph)
}
