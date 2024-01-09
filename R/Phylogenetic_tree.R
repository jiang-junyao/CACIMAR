#' Plot phylogenetic tree
#' @description Build a phylogenetic tree to identify the conservation of the cell types across species
#' @param SCT_matrix matrix, contain the conservation scores of cell types across species
#' @param tree_method character, method used for constructing phylogenetic tree, it can be "hierarchical clustering", or "classic Neighbor-Joining"
#' @param hcluster.method character, method used for clustering phylogenetic tree, it can be one of "average", "single", "complete", "mcquitty", "median", "centroid"
#' @param species.vector vector, length is equal to nrow(SCT_matrix), used for annotation bar
#' @param layout.tree character, method used for the layout of the tree. it can be one of 'rectangular', 'dendrogram', 'fan', 'circular'
#' @param tree_size numeric, size of the tree's branch
#' @param xlim_tree numeric, x limit for the picture, used for the width of the picture
#' @param fontface_for_tiplab logical, whether used fontface for tip label. If it is TRUE, then the tiplab_cols for tip label won't work
#' @param set_fontface_for_tiplab if fontface_for_tiplab is TRUE, the parameter will work. By default, it is set with c("bold.italic", "italic", "plain")
#' @param conserved_hm_celltype vector, contain the conserved cell types in heatmap, or celltypes defined by yourself
#' @param bar_width integer, width of annotation bar
#' @param Show_Pairwise_distances logical, whether save the plot of pairwise distances between input data and constructed tree or not
#' @param boot_times integer, the times for bootstrap test, usually 100, or you can set larger 1000
#' @param offset numeric, the distance of celltype annotation bar from tip
#' @param offset2 numeric, the distance of species annotation bar from tip
#' @param offset_tiplab numeric, the distance of labels from tip
#' @param fontface_tiplab charcter, the fontface used for tip label
#' @param tiplab.size numeric, the size of tip labels
#' @param fontface_tiplab_with_colors charcter, colors for tip label. if fontface_for_tiplab is FALSE, then can set one fontface for all tip labels, like "bold", "italic", "plain"
#' @param tiplab_cols vector, colors for tip group
#' @param tippoint.shape numeric, the shape used for tip point, such as 21, 17, 22, ......
#' @param tippoint.shape.size numeric, the size of the tip point shape
#' @param geom_nodepoint integer, the size of the node point
#' @param show_branch.length logical, whether show the branch length, default is TRUE
#' @param round_x numeric, the decimal places contained for round() function for branch length number
#' @param fill_brandlength character, colors for branch length
#' @param size_brandlength numeric, the size for the branch length label
#' @param annotation_colors_df dataframe, must contain two column, named "colors" and "celltype", default is null
#' @param brewer_pal_used character, the brewer_pal in RColorBrewer, like "Set1",  "Set3", "Paired". It set colors for tip label. It works when fontface_for_tiplab is FALSE
#' @param col.value chacterize vector, colors used for species bar annotation, must be the same length of the number of species
#' @param colors_labels chacterize vector, labels for the legend, the length of colors_breaks
#' @param show_colnames logical, whether show the colnames of annotation bar
#' @param colnames_angle  integer, the angle of the colnames of annotation bar
#' @param colnames_position character, the position of the colnames of annotation bar, can be top or bottom
#' @param font.size_colnames integer, the size of the colnames of annotation bar
#' @param colnames_offset_y numeric, the distance of the colnames from the annotation bar
#' @param custom_column_labels character, edit the colnames of annotation bar by yourself
#' @param bootstrap_value_size numeric, the size of the bootstrap values showed in the tree. If you don't want to showed the  bootstrap values in the tree, you can set the parameter to 0
#' @param bootstrap_value_col character, the color of the bootstrap values showed in the tree
#' @param legend.text_size integer, the size of the legend
#' @param width numeric, the width of the plot
#' @param height numeric, the height of the plot
#' @param plot.margin the plot margin, it should be something like margin(10, 10, 10, -10), the four number represent top, right, bottom, left, respectively
#' @param legend.position character, one of "right", "left", "top", "bottom", the position of the legend
#'
#' @return save plot in present file directory
#' @export
#'
#' @examples
#' library(CACIMAR, lib.loc = "/usr/local/lib/R/site-library")
#' load(system.file("extdata", "zf_mm_markers.rda", package = "CACIMAR"))
#' expression <- CACIMAR::Identify_ConservedCellTypes(OrthG_Mm_Zf,Zf_marker,Mm_marker,'zf','mm')
#' SCT_matrix <- expression[[2]]
#' # get the species names
#' species.vector <- substr(rownames(SCT_matrix), 1, 2)
#' # get the conserved celltypes in heatmap
#' SNT_h <- SCT_matrix[grep('mm',rownames(SCT_matrix)),as.numeric(grep('zf',colnames(SCT_matrix)))]
#' conserved_hm_celltypes <- get_conserved_hm_celltypes(SNT_h)
#' Plot_phylogenetic_tree(SCT_matrix = SCT_matrix,
#'                        species.vector = species.vector,
#'                        conserved_hm_celltype = conserved_hm_celltypes)
#'
Plot_phylogenetic_tree <- function(
    SCT_matrix,
    tree_method = "hierarchical clustering",
    hcluster.method = "average",
    species.vector,
    layout.tree = 'rectangular',
    tree_size = 1.4,
    xlim_tree = 0.65,
    fontface_for_tiplab = TRUE,
    set_fontface_for_tiplab = c("bold.italic", "italic", "plain"),
    conserved_hm_celltype = NULL,
    bar_width = .025,
    Show_Pairwise_distances = FALSE,
    boot_times = 100,
    offset = 0.04,
    offset2 = 0.022,
    offset_tiplab = 0.07,
    fontface_tiplab_with_colors = "bold",
    tiplab.size = 7,
    tiplab_cols = NULL,
    tippoint.shape = 21,
    tippoint.shape.size = 0,
    geom_nodepoint = 5,
    show_branch.length = TRUE,
    round_x = 2,
    fill_brandlength = 'lightgreen',
    size_brandlength = 6,
    annotation_colors_df = NULL,
    brewer_pal_used = "Set3",
    col.value = c(rgb(102/255,46/255,115/255), rgb(31/255,153/255,139/255)),
    colors_labels = NULL,
    show_colnames = FALSE,
    colnames_angle = 0,
    colnames_position = "top",
    font.size_colnames = 5,
    colnames_offset_y = 0.8,
    custom_column_labels = c("Species", "Celltypes"),
    bootstrap_value_size = 0,
    bootstrap_value_col = "black",
    legend.text_size = 18,
    width = 13,
    height = 22,
    plot.margin = margin(10, 10, 10, -10),
    legend.position = "top"
) {

  stopifnot(is.matrix(SCT_matrix))
  stopifnot(length(hcluster.method) == 1 && hcluster.method %in% c("average", "single", "complete", "mcquitty", "median", "centroid"))
  stopifnot(is.vector(species.vector))
  stopifnot(length(layout.tree) == 1 && layout.tree %in% c('rectangular', 'dendrogram', 'fan', 'circular'))
  stopifnot(is.numeric(xlim_tree))
  stopifnot(is.numeric(offset) && offset >= 0)
  stopifnot(is.numeric(offset2) && offset2 >= 0)
  stopifnot(is.numeric(tiplab.size) && tiplab.size > 0)
  stopifnot(is.numeric(tippoint.shape) && length(tippoint.shape) == 1)
  stopifnot(is.numeric(tippoint.shape.size) && tippoint.shape.size >= 0)
  stopifnot(is.numeric(geom_nodepoint) && geom_nodepoint >= 0)
  stopifnot(length(col.value) >= 2 && is.character(col.value))
  stopifnot(is.logical(show_colnames))
  stopifnot(is.numeric(colnames_angle))
  stopifnot(is.character(colnames_position) & colnames_position %in% c("top", "bottom"))
  stopifnot(is.numeric(font.size_colnames))
  stopifnot(is.numeric(colnames_offset_y))
  stopifnot(is.character(custom_column_labels))
  stopifnot(is.numeric(bootstrap_value_size))
  stopifnot(is.numeric(legend.text_size))
  stopifnot(is.character(bootstrap_value_col))
  stopifnot(is.numeric(width))
  stopifnot(is.numeric(height))

  dist <- as.dist(1 - SCT_matrix)
  if (tree_method == "hierarchical clustering") {
    tree <- stats::hclust(dist, method = hcluster.method)
    tree <- ape::as.phylo(tree)
  } else if (tree_method == "classic Neighbor-Joining") {
    tree <- ape::nj(dist)
  }

  tree$species <- species.vector
  if (tree_method == "classic Neighbor-Joining") {
    Pairwise_distances_main <- "classic Neighbor-Joining algorithm"
    cat("Creating tree with classic Neighbor-Joining algorithm.\n")
    cat("Neighbor-joining: Taking the two closest nodes of the tree and defining them as neighbors; repeat this process until all nodes are paired together.\n")
  } else if (tree_method == "hierarchical clustering" & hcluster.method == "average") {
    Pairwise_distances_main <- "UPGMA"
    cat("Creating tree with UPGMA method.\n")
    cat("UPGMA (Unweighted Pair Group Method with Arithmetic Mean): This is the simplest method for constructing trees, assuming the same evolutionary speed for all lineages.\n")
  } else {
    Pairwise_distances_main <- paste("hierarchical clustering", hcluster.method, sep = "_")
    cat(paste("Creating tree with hierarchical clustering", hcluster.method, sep = ":"))
  }

  if (Show_Pairwise_distances) {
    plot_Pairwise_distances(dist, tree, Pairwise_distances_main)
  }

  set.seed(1234)
  if (tree_method == "classic Neighbor-Joining") {
    myBoots <- ape::boot.phylo(tree, SCT_matrix, function(d) {
      dist <- as.dist(1 - d)
      dist[is.na(dist)] <- 1
      tree <- ape::nj(dist)
    },
    jumble = F,
    trees = T,
    B = boot_times,
    block = ceiling(ncol(SCT_matrix)/2),
    mc.cores = 1)

  } else if (tree_method == "hierarchical clustering") {
    myBoots <- ape::boot.phylo(tree, SCT_matrix, function(d, method = hcluster.method) {
      dist <- as.dist(1 - d)
      dist[is.na(dist)] <- 1
      tree <- stats::hclust(dist, method = method)
      ape::as.phylo(tree)
    },
    jumble = F,
    trees = T,
    B = boot_times,
    block = ceiling(ncol(SCT_matrix)/2),
    mc.cores = 1)

  }

  tree_plot <- Plot_Tree(tree,
                         layout.tree = layout.tree,
                         tree_size = tree_size,
                         xlim_tree = xlim_tree,
                         fontface_for_tiplab = fontface_for_tiplab,
                         set_fontface_for_tiplab = set_fontface_for_tiplab,
                         conserved_hm_celltype = conserved_hm_celltype,
                         bar_width = bar_width,
                         offset = offset,
                         offset2 = offset2,
                         offset_tiplab = offset_tiplab,
                         fontface_tiplab_with_colors = fontface_tiplab_with_colors,
                         tiplab.size = tiplab.size,
                         tiplab_cols = tiplab_cols,
                         tippoint.shape = tippoint.shape,
                         tippoint.shape.size = tippoint.shape.size,
                         geom_nodepoint = geom_nodepoint,
                         show_branch.length = show_branch.length,
                         round_x = round_x,
                         fill_brandlength = fill_brandlength,
                         size_brandlength = size_brandlength,
                         col.value = col.value,
                         annotation_colors_df = annotation_colors_df,
                         brewer_pal_used = brewer_pal_used,
                         colors_labels = colors_labels,
                         show_colnames = show_colnames,
                         colnames_angle = colnames_angle,
                         colnames_position = colnames_position,
                         font.size_colnames= colnames_position,
                         colnames_offset_y= colnames_offset_y,
                         custom_column_labels= custom_column_labels,
                         bootstrap_value_size = bootstrap_value_size,
                         bootstrap_value_col = bootstrap_value_col,
                         legend.text_size = legend.text_size,
                         bootstrap_values = myBoots$BP,
                         plot.margin = plot.margin,
                         legend.position = legend.position)
  ggsave(paste(Pairwise_distances_main, "Tree.pdf", sep = "_"), tree_plot$plot, width = width, height = height)
  return(list("Tree" = tree, "Bootstrap" = myBoots, "Plot" = tree_plot$plot, "conserved_table" = tree_plot$conserved_table))
}

plot_Pairwise_distances <- function(dist, tree, Pairwise_distances_main) {
  x <- as.vector(dist)
  y <- as.vector(as.dist(cophenetic(tree)))
  plot(x, y,
       xlab = "Original Pairwise Distances",
       ylab = "Pairwise Distances on the Tree",
       pch = 16, col = "black", cex = 1.5,
       xlim = c(min(x) * 0.95, max(x) * 1.05),
       ylim = c(min(y) * 0.95, max(y) * 1.05),
       cex.lab = 1.5, cex.axis = 1.2, cex.main = 1.5)
  abline(lm(y~x), col = "darkgrey", lwd = 1)
  box(lwd = 0.7, col = "black")
  legend("topleft", legend = "Regression Line", col = "darkgrey", lwd = 0.7, cex = 1.2)
  text(x = max(x), y = max(y) + 0.03 * diff(range(y)),
       labels = paste("Correlation: ", round(cor(x, y)^2, digits = 2)),
       pos = 2, col = "black", cex = 1.2)
  title(main = Pairwise_distances_main, line = 1, font.main = 2)
  axis(1, cex.axis = 1.2, lwd = 1, col.axis = "black")
  axis(2, cex.axis = 1.2, lwd = 1, col.axis = "black")
  par(family = "serif")
  dev.copy(pdf, paste(Pairwise_distances_main, "Pairwise_distances.pdf", sep = "_"))
  dev.off()
  cat("Correlation of the tree and distances:\n")
  cat(cor(x, y)^2, "\n")
}

get_conserved_hm_celltypes <- function(data) {
  conserved_hm_celltypes_df <- data.frame()
  for (i in 1:nrow(data)) {
    single_row <- data[i, ]
    max_col_id <- which.max(single_row)
    single_col <- data[, max_col_id]
    max_row_id <- which.max(single_col)

    if (max_row_id == i) {
      if (data[i, max_col_id] > quantile(data, 0.75)) {
        conserved_hm_celltypes_df <- rbind(conserved_hm_celltypes_df, c(rownames(data)[i], colnames(data)[max_col_id]))
      }
    }
  }
  conserved_hm_celltypes <- unlist(as.vector(conserved_hm_celltypes_df), use.names = FALSE)
  return(conserved_hm_celltypes)
}

Plot_Tree <- function(
    tree,
    layout.tree,
    tree_size,
    xlim_tree,
    fontface_for_tiplab,
    set_fontface_for_tiplab,
    conserved_hm_celltype,
    bar_width,
    geom_nodepoint,
    show_branch.length,
    round_x,
    fill_brandlength,
    size_brandlength,
    col.value, # species colors
    annotation_colors_df, # celltype colors df
    brewer_pal_used,
    colors_labels,
    show_colnames,
    colnames_angle,
    colnames_position,
    font.size_colnames,
    colnames_offset_y,
    custom_column_labels,
    bootstrap_value_size,
    bootstrap_value_col,
    legend.text_size,
    offset,
    offset2,
    offset_tiplab,
    tiplab.size,
    fontface_tiplab_with_colors,
    tiplab_cols,
    tippoint.shape,
    tippoint.shape.size,
    bootstrap_values,
    plot.margin,
    legend.position
) {
  require(ggtree)
  require(ggplot2)
  require(scales)

  tree.p <- ggtree::ggtree(tree, layout = layout.tree, size = tree_size, col = "black") + xlim(NA, xlim_tree) #text show
  data <- fortify(tree)

  tree.df <- data.frame("tip.label" = tree$tip.label)
  rownames(tree.df) <- tree$tip.label
  tree.df2 <- data.frame("species" = tree$species)
  rownames(tree.df2) <- tree$tip.label

  Tip1 <- filter(data, isTip == "TRUE")
  Tip_p <- as.data.frame(table(Tip1$parent))
  conserved_tree_node <- Tip_p$Var1[which(Tip_p$Freq == 2)]
  not_conserved_tree_node <- Tip_p$Var1[which(Tip_p$Freq == 1)]
  Tip1$Tip_group <- "not_conserved"
  conserved_t_id <- as.vector(sapply(conserved_tree_node, function(x) grep(x, Tip1$parent)))
  Tip1$Tip_group[conserved_t_id] <- "poorly_conserved"
  conserved_h_id <- as.vector(sapply(conserved_hm_celltype, function(x) grep(x, Tip1$label)))
  Tip1$Tip_group[conserved_h_id] <- "conserved"
  Tip1$Tip_group <- factor(Tip1$Tip_group, levels = c("conserved", "poorly_conserved", "not_conserved"))

  if (fontface_for_tiplab) {
    Tip1$fontface <- set_fontface_for_tiplab[3]
    Tip1$fontface[which(Tip1$Tip_group == "poorly_conserved")] <- set_fontface_for_tiplab[2]
    Tip1$fontface[which(Tip1$Tip_group == "conserved")] <- set_fontface_for_tiplab[1]
    Tip1$fontface <- factor(Tip1$fontface, levels = set_fontface_for_tiplab)
  }

  cat("The conservation of cell types are shown below:")
  print(table(Tip1$Tip_group))
  cat("\n\n")
  cat("conserved cell types are:\n")
  cat(Tip1$label[which(Tip1$Tip_group == "conserved")])
  cat("\n\n")
  cat("poorly conserved cell types are:\n")
  cat(Tip1$label[which(Tip1$Tip_group == "poorly_conserved")])
  cat("\n\n")
  cat("not conserved cell types are:\n")
  cat(Tip1$label[which(Tip1$Tip_group == "not_conserved")])
  cat("\n\n")

  if (!is.null(annotation_colors_df)) {
    colors_a <- c(annotation_colors_df$colors)
    colors_a <- setNames(colors_a, annotation_colors_df$celltype)
  } else {
    no_species_prefix_celltypes_all <- unique(substring(tree.df$tip.label, 3))
    species_prefix <- unique(substring(tree.df$tip.label, 1, 2))
    species_a_posible_celltypes <- paste0(species_prefix[1], no_species_prefix_celltypes_all)
    species_b_posible_celltypes <- paste0(species_prefix[2], no_species_prefix_celltypes_all)
    species_a_posible_celltypes_col <- colorRampPalette(RColorBrewer::brewer.pal(11, brewer_pal_used))(length(species_a_posible_celltypes))
    colors_a <- c(species_a_posible_celltypes_col, species_a_posible_celltypes_col)
    colors_a <- setNames(colors_a, c(species_a_posible_celltypes, species_b_posible_celltypes))
  }

  col.value <- setNames(col.value, rev(unique(tree.df2$species)))
  if (is.null(colors_labels)) {
    colors_labels_used <- unique(tree.df2$species)
  } else {
    colors_labels_used <- colors_labels
  }

  tregraph <- gheatmap(
    tree.p,
    tree.df,
    offset = offset,
    width = bar_width,
    colnames = show_colnames,
    colnames_position = colnames_position,
    font.size = font.size_colnames,
    colnames_angle = colnames_angle,
    colnames_offset_y = colnames_offset_y,
    custom_column_labels = custom_column_labels
  ) +
    scale_fill_manual(values = colors_a, breaks = NULL)
  tregraph1 <- tregraph + ggnewscale::new_scale_fill()
  tregraph2 <- gheatmap(
    tregraph1,
    tree.df2,
    offset = offset2,
    width = bar_width,
    colnames = show_colnames,
    colnames_position = colnames_position,
    font.size = font.size_colnames,
    colnames_angle = colnames_angle,
    colnames_offset_y = colnames_offset_y,
    custom_column_labels = custom_column_labels
  ) +
    scale_fill_manual(
      values = col.value,
      breaks = rev(unique(tree.df2$species)),
      labels = colors_labels_used
    )

  if (is_null(tiplab_cols) & length(unique(Tip1$Tip_group)) == 3) {
    tiplab_cols <- c("#003472","#4b5cc4", "#737373")
  } else  {
    tiplab_cols <-grDevices::colorRampPalette(c("#CC6666", "#9999CC", "#66CC99"))(length(unique(Tip1$Tip_group)))
  }

  if (fontface_for_tiplab) {
    tregraph2 <- tregraph2 + ggnewscale::new_scale_fill()
    tregraph3 <- tregraph2 %<+% Tip1 +
      geom_tiplab(aes(label = substring(label, 3)),
                  offset = offset_tiplab,
                  size = tiplab.size,
                  fontface = Tip1$fontface) +
      geom_tippoint(
        shape = tippoint.shape,
        size = tippoint.shape.size
      ) +
      geom_nodepoint(size = geom_nodepoint)
  } else {
    tregraph2 <- tregraph2 + ggnewscale::new_scale_fill()
    tregraph3 <- tregraph2 %<+% Tip1 +
      geom_tiplab(aes(label = substring(label, 3), colour = Tip_group),
                  offset = offset_tiplab,
                  size = tiplab.size, fontface = fontface_tiplab_with_colors) +
      scale_color_manual(values = tiplab_cols, breaks = NULL)+
      geom_tippoint(
        shape = tippoint.shape,
        size = tippoint.shape.size
      ) +
      geom_nodepoint(size = geom_nodepoint)
  }

  if (show_branch.length) {
    tregraph3 <- tregraph3 + geom_label(aes(x=branch, label= round(data$branch.length, round_x)), fill=fill_brandlength, size = size_brandlength)
  }

  internal_nodes <- data[data$isTip == FALSE, ]
  internal_nodes$bootstrap <- bootstrap_values
  tregraph3 <- tregraph3 +
    geom_text(
      data = internal_nodes,
      aes(x = x, label = bootstrap),
      vjust = -0.5,
      hjust = 1,
      size = bootstrap_value_size,
      color = bootstrap_value_col
    ) +
    theme(
      legend.text = element_text(size = legend.text_size),
      legend.title = element_blank()
    ) +
    theme(plot.margin = plot.margin)
  tregraph3 <- tregraph3 + theme(legend.position = legend.position,
                                 legend.spacing.x = unit(0.4, 'cm'),
                                 legend.key.size = unit(1.2, "lines"),
                                 legend.key.width = unit(1.5, "lines"),
                                 legend.text = element_text(margin = margin(r = 30)))
  return(list("plot" = tregraph3, "conserved_table" = Tip1))
}

