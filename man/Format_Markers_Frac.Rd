% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Identify_Markers.R
\name{Format_Markers_Frac}
\alias{Format_Markers_Frac}
\title{Format marker genes for plotting}
\usage{
Format_Markers_Frac(Marker_genes)
}
\arguments{
\item{Marker_genes}{data.frame, generated by \code{\link{Identify_Markers}}}
}
\description{
Order the gene expression in each cluster to make the heatmap
look better
}
\examples{
data("pbmc_small")
all.markers <- Identify_Markers(pbmc_small)
all.markers2 <- Format_Markers_Frac(all.markers)
}
