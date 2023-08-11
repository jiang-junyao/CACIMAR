# Run like this
# Identify_CellType_evolution(Zf_seurat, Mm_seurat, ConservedMarker, species1_name = "Zf", species2_name = "Mm", Used_species1_ID = "Used_zf_ID", Used_species2_ID = "Used_mm_ID", dist.method = "euclidean", hclust.method = "average", layout.tree = 'rectangular', geom_nodepoint = 0, col.value = c("#990000", "#660099"))

library(Seurat)
library(CACIMAR)
library(ggplot2)
library(ape)
library(ggtree)
## Identify_CellType_evolution
Identify_CellType_evolution <- function(species1_seurat, species2_seurat, ConservedMarker, species1_name = "Mm", species2_name = "Zf", Used_species1_ID, Used_species2_ID, dist.method = "euclidean", hclust.method = "average", layout.tree = "rectangular", geom_nodepoint = 0, col.value = c("#990000", "#660099")){
  print("Get mean expression")
  # process to get mean expression
  my_bind_mat <- process_mean_data(species1_seurat, species2_seurat, ConservedMarker, species1_name, species2_name, Used_species1_ID, Used_species2_ID)
  # Do the hclust analysis and plot
  p <- Plot_celltype_tree(my_bind_mat, dist.method, hclust.method, layout.tree, geom_nodepoint, col.value)
  print(p)
  return(p)
}

# The seurat object input should have been normalized 
# active.ident should be the celltype
process_mean_data <- function(species1_seurat, species2_seurat, ConservedMarker, species1_name, species2_name, Used_species1_ID, Used_species2_ID) {
  # Subset data for species 1
  species1_marker_seurat <- species1_seurat[ConservedMarker[[Used_species1_ID]], ]
  species1_marker_data <- GetAssayData(species1_marker_seurat, assay = "RNA", slot = "data")
  
  # Subset data for species 2
  species2_marker_seurat <- species2_seurat[ConservedMarker[[Used_species2_ID]], ]
  species2_marker_data <- GetAssayData(species2_marker_seurat, assay = "RNA", slot = "data")
  
  # Calculate mean expression for species 1
  species1_mean_mat <- aggregate(as.data.frame(t(as.matrix(species1_marker_data))), by = list(species1_marker_seurat@active.ident), FUN = "mean") # cell x gene
  species1_group_names <- species1_mean_mat$Group.1
  names(species1_group_names) <- NULL
  rownames(species1_mean_mat) <- species1_group_names
  species1_mean_mat <- species1_mean_mat[, -1]
  rownames(species1_mean_mat) <- paste(species1_name, rownames(species1_mean_mat), sep = "_")
  
  # Calculate mean expression for species 2
  species2_mean_mat <- aggregate(as.data.frame(t(as.matrix(species2_marker_data))), by = list(species2_marker_seurat@active.ident), FUN = "mean") # cell x gene
  species2_group_names <- species2_mean_mat$Group.1
  names(species2_group_names) <- NULL
  rownames(species2_mean_mat) <- species2_group_names
  species2_mean_mat <- species2_mean_mat[, -1]
  rownames(species2_mean_mat) <- paste(species2_name, rownames(species2_mean_mat), sep = "_")
  
  # Merge data for both species 
  colnames(species2_mean_mat) <- colnames(species1_mean_mat)
  combined_mat <- rbind(species1_mean_mat, species2_mean_mat)
  
  return(combined_mat)
}

Plot_celltype_tree <- function(expression_matrix, dist.method = c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"), hclust.method = c("average", "single", "complete"), layout.tree = c('rectangular', 'dendrogram', 'slanted', 'ellipse', 'roundrect', 'fan', 'circular', 'inward_circular', 'radial', 'equal_angle', 'daylight', 'ape'), geom_nodepoint = 0, col.value){
  tree <- Celltypes_hclust_tree(expression_matrix, dist.method, hclust.method)
  tree.plot <- Make_tree_plot(tree, layout.tree, geom_nodepoint, col.value)
  return(tree.plot)
}

Celltypes_hclust_tree <- function(expression_matrix, dist.method, hclust.method){
  distance_matrix <- dist(expression_matrix, method = dist.method)
  hc <- hclust(distance_matrix, method = hclust.method)
  tree <- ape::as.phylo(hc)
  # add species column to the tree
  first_parts <- vector("character", length(tree$tip.label))
  for (i in 1:length(tree$tip.label)) {
    split_parts <- strsplit(tree$tip.label[i], "_")  
    first_parts[i] <- split_parts[[1]][1] 
  }
  first_parts <- factor(first_parts, levels = unique(first_parts))
  tree$species <- first_parts
  return(tree)
}

Make_tree_plot <- function(tree, layout.tree, geom_nodepoint, col.value){
  tree.p <- ggtree::ggtree(tree, layout=layout.tree, size=0.8, col="black")
  data=fortify(tree)
  tree.df <- data.frame("tip.label"= tree$tip.label, "species" = tree$species)
  tregraph <- tree.p%<+%tree.df+
    geom_tiplab(aes(col=species), offset=0.1)+
    geom_tippoint(aes(col=species)) + 
    geom_nodepoint(color="black", alpha=1/4, size=geom_nodepoint) + 
    theme_tree2() + 
    xlim(NA, max(data$x)*1.3) +
    scale_fill_manual(values = col.value) +
    scale_color_manual(values = col.value) + 
    guides(shape = guide_legend(override.aes = list(shape = c(16, 16))))
  
  return(tregraph)
}
