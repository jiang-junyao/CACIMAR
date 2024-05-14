#' Phylogenetree of cell types across three species
#' @describeIn # construct a phylogenetree across three species with conserved score of cell types
#'
#' @param CSCT matrix, a matrix with conserved score of cell types. rownames and colnames are cell types that will appear in the tree.
#' @param species_colors vecotor, colors vectors for species, like c("red", "blue", "orange") for three species.
#' @param tree_method character, the mathod used to constuct the tree, like "Neighbor-Joining", "average", "ward.D", "ward.D2", "single", et al.
#' @param tree_layout default is "fan", support for layout in ggtree
#' @param tree_size numeric, size of the tree
#' @param sp_ann_offset numeric, distance of the heapmap annotation from the tree
#' @param sp_ann_width numeric, width of the annotation heapmap
#' @param sp_ann_font.size numeric, size of the annotation font
#' @param sp_ann_hjust numeric, hjust of the annotation
#' @param sp_ann_colnames_angle numeric, angle of the annotation
#' @param legend.position, one of the right, left, bottom, top, none
#' @param legend.t_size numeric, size of the legend text
#' @param legend.tl_size numeric, size of the legend title
#' @param tipl_offset numeric, distance of the tip label from the tree
#' @param tipl_size numeric, size of the tip label
#' @param nodesize numeric, size of the node
#' @param t_width numeric, save the piture with this width
#' @param t_height numeric, save the piture with this height
#'
#' @return the tree object
#' @export
#'
#' @examples
#' p <- phylogenetic_tree_circle(CSCT, tree_method = "Neighbor-Joining")
#' p
phylogenetic_tree_circle <- function(CSCT,
                              species_colors = NULL,
                              tree_method = c("Neighbor-Joining", "average", "ward.D", "ward.D2", "single", "complete",  "mcquitty", "median", "centroid"),
                              tree_layout="fan",
                              tree_size = 1.35,
                              sp_ann_offset = 0.2,
                              sp_ann_width=0.08,
                              sp_ann_font.size=0,
                              sp_ann_hjust = 0,
                              sp_ann_colnames_angle=0,
                              legend.position = "top",
                              legend.t_size = 16,
                              legend.tl_size = 16,
                              tipl_offset = 1.55,
                              tipl_size = 12,
                              nodesize = 1,
                              t_width = 29,
                              t_height = 29
) {
  library(ggplot2)
  library(ggtree)
  if (is.null(species_colors)) {
    species_colors <-grDevices::colorRampPalette(c(rgb(0/255,114/255,189/255), rgb(248/255,115/255,106/255), rgb(169/255,169/255,169/255)))(length(unique(substr(colnames(CSCT), 1, 2))))
  }
  tree <- ape::nj(as.dist(1 - CSCT))
  if (tree_method == "average") {
    tree <- stats::hclust(as.dist(1 - CSCT), method = tree_method)
    tree <- ape::as.phylo(tree)
  } else if (tree_method == "Neighbor-Joining") {
    tree <- ape::nj(as.dist(1 - CSCT))
  }
  pggtree <- ggtree::ggtree(tree, layout=tree_layout, ladderize = FALSE, branch.length = "none", size = tree_size)
  tree.df2 <- data.frame("species" = substr(tree$tip.label, 1, 2), row.names = tree$tip.label)
  pggtree2 <- ggtree::gheatmap(pggtree, tree.df2, offset = sp_ann_offset, width=sp_ann_width, font.size=sp_ann_font.size, hjust = sp_ann_hjust, colnames_angle=sp_ann_colnames_angle) +
    ggplot2::scale_fill_manual(values = species_colors, name = 'species') +
    theme(legend.position = legend.position, legend.text = element_text(size = legend.t_size), legend.title = element_text(size = legend.tl_size))
  tree_data <- fortify(tree)
  Tipd <- subset(tree_data, isTip == "TRUE")
  tregraph3 <- pggtree2 %<+% Tipd +
    geom_tiplab(aes(label = substring(label, 3)),
                offset = tipl_offset, size = tipl_size) +
    geom_nodepoint(size = nodesize)
  ggsave(paste(tree_method, "phylogenetic_tree.pdf", sep = "_"), tregraph3, width = t_width, height = t_height)
  return(tregraph3)
}
