
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CCtMR

<!-- badges: start -->
<!-- badges: end -->

CCtMR is an R package to identify cross-species marker genes, cell types
and gene regulatory networks based on scRNA-seq data.

## Installation

Install CCtMR from github, run:

``` r
# install.packages("devtools")
devtools::install_github("jiangjunyao123/CCtMR")
```

## Tutorial

### 1.Identify Cell types

Using known marker genes to annotate each cluster. This method is based
on AUC (area under the receiver operating characteristic curve of gene
expression), and is very sensitive to the marker genes input.

#### Inputs data

1.  Seurat object Seurat object should have clustering information in
    active.ident slot and meta.data slot

2.  Marker genes table Rownames of Marker genes table should be the same
    format as the rownames format of seurat object, and should contain
    CellType column

#### Example

``` r
Marker<-read.table('D:\\GIBH\\platform\\test data/Retinal_markersZf.txt',header = T)
head(Marker)
#>            EnsemblID  Symbol            CellType
#> 1 ENSDARG00000045904   nr2e3 Rods,Rodprogenitors
#> 2 ENSDARG00000019566 neurod1      RodProgenitors
#> 3 ENSDARG00000099572   hmgn2      RodProgenitors
#> 4 ENSDARG00000002193     rho                Rods
#> 5 ENSDARG00000100466     nrl                Rods
#> 6 ENSDARG00000011235    otx2                Rods
```

``` r
### I only use 3 clusters here to reduce the running time
seurat_object<-readRDS('D:\\GIBH\\platform\\test data/Zebrafishdata.rds')
seurat_object<-subset(a,idents = c(1,2,3))
rownames(Marker)<-Marker[,1];Marker<-Marker[,-1]
zfcelltype<-Identify_CellType(seurat_object,Marker)
```

### 2.Identify markers

``` r
library(CCtMR)
setwd('D:\\GIBH\\platform\\test data')
a<-readRDS('D:\\GIBH\\platform\\test data/Zebrafishdata.rds')
### I only use 3 clusters here to reduce the running time
a<-subset(a,idents = c(1,2,3))
Marker1<-Identify_Markers(a,Spec1='Zf')
cc<-Seurat_Markers(a,Spec1 = 'Zf')
```

### 3.Get cross-species markers

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

### 4.plot MarkersHeaptmap

``` r
Marker1_plot<-Format_Markers_Frac(Marker1)
plot_MarkersHeatmap(Marker1_plot)
```

### 5.Cross-species celltype heamtmap

``` r
expression<-read.table('D:\\GIBH\\platform\\CellType_Comp\\CellType_Comp\\Data/mmP60RmmNMDA_chP10chNMDA_zfAdzfNMDA_Power01_SharedMarkers_Frac.txt')
celltypes<-read.delim("D:/GIBH/platform/CellType_Comp/CellType_Comp/Data/mmP60RmmNMDA_chP10chNMDA_zfAdzfNMDA_Cell_Types.txt")
a<-Heatmap_Cor(expression,celltypes,cluster_cols=T, cluster_rows=F)
Plot_tree(a)
```

\#\#\#6. Cross-species regulatory networks
