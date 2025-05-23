% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pheatmap_LR1.R
\name{get_average_expression}
\alias{get_average_expression}
\title{Geometric mean of ligand and receptor}
\usage{
get_average_expression(
  specie_conserved_CCC,
  seurat_object,
  avg_group,
  species_set,
  assay_set
)
}
\arguments{
\item{specie_conserved_CCC}{dataframe, also result of the ligand receptor analysis(for example, ConservedCCI$Conserved_ligand_receptor$sp1_orthg_ccc_df). It must contain 4 column for ligand, receptor, source celltype, target celltype}

\item{seurat_object}{a seurat object contains the expression data for celltypes of source and target}

\item{avg_group}{a character, which variable used to average the expression, usally "celltype"}

\item{species_set}{species_set should only contain two characters, like "Mm" to label the species}

\item{assay_set}{Assay to get the expression, like "RNA", or "SCT"}
}
\description{
This function caculates the geometric mean of ligand and receptor. Average expression of ligand comes from source celltype, and verage expression of receptor comes from target celltype
}
\examples{
LRpair_show <- get_average_expression(CCC_conserved_sumary = ConservedCCI[[2]], specie_conserved_CCC = ConservedCCI$Conserved_ligand_receptor$sp1_orthg_ccc_df, seurat_object = Mm_seurat_object, avg_group = "celltype", species_set = "Mm")
}
