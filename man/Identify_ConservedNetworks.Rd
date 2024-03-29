% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CSNetworks.R
\name{Identify_ConservedNetworks}
\alias{Identify_ConservedNetworks}
\title{Identify conserved regulatory networks}
\usage{
Identify_ConservedNetworks(
  OrthG,
  Species1_GRN,
  Species2_GRN,
  Species_name1,
  Species_name2
)
}
\arguments{
\item{OrthG}{ortholog genes database}

\item{Species1_GRN}{gene regulatory network of species 1}

\item{Species2_GRN}{gene regulatory network of species 2}

\item{Species_name1}{character, indicating the species names of Species1_GRN}

\item{Species_name2}{character, indicating the species names of Species2_GRN}
}
\value{
list contains two df. First df contains details of conserved regulatory
network, second df contains NCS between module pairs
}
\description{
Use Score of Conserved network to identify conserved regulatory
network modules based on homologous genes databased and topology of networks
}
\examples{
load(system.file("extdata", "gene_network.rda", package = "CACIMAR"))
n1 <- Identify_ConservedNetworks(OrthG_Mm_Zf,mm_gene_network,zf_gene_network,'mm','zf')
}
