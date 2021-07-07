
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CCtMR

<!-- badges: start -->
<!-- badges: end -->

The goal of CCtMR is to …

## Installation

You can install the released version of CCtMR from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("CCtMR")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("jiang-junyao/CCtMR")
```

## Example

### Identify markers

``` r
library(CCtMR)
setwd('D:\\GIBH\\platform\\test data')
a<-readRDS('Zebrafishdata.rds')
### I only use 3 clusters here to reduce the running time
a<-subset(a,idents = c(1,2,3))
Marker1<-Identify_Markers(a,Spec1='Zf')
```

### plot MarkersHeaptmap

``` r
Marker1_plot<-Format_Markers_Frac(Marker1)
plot_MarkersHeatmap(Marker1_plot)
```

### Get cross-species markers

``` r
###Get_Wilcox_Markers_Cond
wilcox<-read.table('D:\\GIBH\\platform\\test data/mmP60RmmNMDA_mmP60mmLD_wilcoxMG_MarkerGenes.txt')
Cor<-read.table('D:\\GIBH\\platform\\test data/mmP60RmmNMDA_mmLD_pbmcSubC_MG_Bin50_R5_GeneCor.txt',header = T)
wilcox_result<-Get_Wilcox_Markers_Cond(wilcox,Cor)
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
