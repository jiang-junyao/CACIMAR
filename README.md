
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CCtMR

<!-- badges: start -->
<!-- badges: end -->

CCtMR is an R package to identify cross-species marker genes, cell types
and gene regulatory networks based on single cell sequencing.

## Installation

Install CCtMR from github, run:

``` r
# install.packages("devtools")
devtools::install_github("jiangjunyao123/CCtMR")
```

## Example

### Identify Cell types

First, using known marker genes to annotate each cluster. This method is
based on AUC (area under the receiver operating characteristic curve of
gene expression), and is very sensitive to the marker genes input.

``` r
Marker<-read.table('D:\\GIBH\\platform\\test data/Retinal_markersZf.txt',header = T)
rownames(Marker)<-Marker[,1];Marker<-Marker[,-1]
zfcelltype<-Identify_CellType(a,Marker)
```

### Identify markers

``` r
library(CCtMR)
setwd('D:\\GIBH\\platform\\test data')
a<-readRDS('D:\\GIBH\\platform\\test data/Zebrafishdata.rds')
### I only use 3 clusters here to reduce the running time
a<-subset(a,idents = c(1,2,3))
Marker1<-Identify_Markers(a,Spec1='Zf')
cc<-Seurat_Markers(a,Spec1 = 'Zf')
```

### Get cross-species markers

``` r
###Get_Wilcox_Markers_Cond ###???usage???
wilcox<-read.table('D:\\GIBH\\platform\\test data/mmP60RmmNMDA_mmP60mmLD_wilcoxMG_MarkerGenes.txt')
Cor<-read.table('D:\\GIBH\\platform\\test data/mmP60RmmNMDA_mmLD_pbmcSubC_MG_Bin50_R5_GeneCor.txt',header = T)
wilcox_result<-get_Wilcox_Markers_Cond(wilcox,Cor)
###Overlap_Markers_Cond
mmMarkers3_F3F0 <- read.delim("D:/GIBH/platform/test data/mmP60RmmNMDA_mmP60mmLD_P03_Markers3_F3F0.txt")
zfMarkers3_F3F0 <- read.delim("D:/GIBH/platform/test data/zfAdzfNMDA_zfAdzfLD_zfAdzfTR_P03_Markers3_F3F0.txt")
mmCelltype<-read.table('D:/GIBH/platform/test data/mmP60RmmNMDA_mmP60mmLD_Cell_Types.txt',header = T)
zfCelltype<-read.table('D:/GIBH/platform/test data/zfAdzfNMDA_zfAdzfLD_zfAdzfTR_Cell_Types.txt',header = T)
mmMarker<-Overlap_Markers_Cond(mmMarkers3_F3F0,mmCelltype,Spec2='mm')
zfMarker<-Overlap_Markers_Cond(zfMarkers3_F3F0,zfCelltype,Spec2='zf')
###Get_Used_OrthG
rownames(mmMarker)<-mmMarkers3_F3F0[,1]
rownames(zfMarker)<-zfMarkers3_F3F0[,1]
OrthG<-read.delim('D:/GIBH/platform/test data/RNA_genes_mmVSzf.txt')
ShMarker<-Get_Used_OrthG(OrthG,mmMarker,zfMarker,Species = c('mm','zf'))
refined_markers<-Refine_TwoSpecies(ShMarker,mmCelltype,zfCelltype,Species = c('mm','zf'))
```

### plot MarkersHeaptmap

``` r
Marker1_plot<-Format_Markers_Frac(Marker1)
plot_MarkersHeatmap(Marker1_plot)
```

### Cross-species celltype heamtmap

``` r
expression<-read.table('D:\\GIBH\\platform\\CellType_Comp\\CellType_Comp\\Data/mmP60RmmNMDA_chP10chNMDA_zfAdzfNMDA_Power01_SharedMarkers_Frac.txt')
celltypes<-read.delim("D:/GIBH/platform/CellType_Comp/CellType_Comp/Data/mmP60RmmNMDA_chP10chNMDA_zfAdzfNMDA_Cell_Types.txt")
a<-Heatmap_Cor(expression,celltypes,cluster_cols=T, cluster_rows=F)
Plot_tree(a)
```

### Cross-species regulatory networks
