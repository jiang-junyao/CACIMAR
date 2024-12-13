
<!-- README.md is generated from README.Rmd. Please edit that file -->

[![](https://img.shields.io/badge/Seurat-version4.2-red.svg)](https://satijalab.org/seurat/articles/get_started.html)

# CACIMAR: cross-species analysis of cell identities, markers and regulations using single-cell sequencing profiles

CACIMAR has the following 4 main functions:

- **Identify all cell-type specific markers in each species**

- **Identify species-specific or evolutionally conserved markers**

- **Identify conserved intracellular regulations (gene regulatory
  networks)**

- **Identify conserved intercellular interactions (cell-cell
  communications)**

<img src="Readme%20figure/Workflow.png"
style="width:70.0%;height:70.0%" />

## Installation

Install CACIMAR from github, run:

``` r
# install.packages("devtools")
devtools::install_github("jiang-junyao/CACIMAR")
```

## News

2024.10.28: remove warning message in ‘identify_conserved_marker’
function.

2024.10.24: Fixed bugs resulting from a package version update.

2024.08.02: Bugs in ‘Identify_CellType’ have been fixed.

### Version1.0.0: First version of CACIMAR was published!

## Qucik start (Run CACIMAR in 1 minutes)

For users who want to perform cross-species analysis quickly, we provide
a brief tutorial based on small single-cell RNA sequencing (scRNA-seq)
data to demonstrate the functionality and API of CACIMAR.

[Quick tutorial for identifying evolutionally conserved celltype and
markers based on scRNA-seq
data](https://jiang-junyao.github.io/CACIMAR/quick_celltype.html)

[Quick tutorial for identifying evolutionally conserved gene regulatory
networks](https://jiang-junyao.github.io/CACIMAR/quick_grn.html)

[Quick tutorial for identifying evolutionally conserved ligand-receptor
interaction and cell-cell
communications](https://jiang-junyao.github.io/CACIMAR/quick_interaction.html)

[Tutorial for build customized orthg database based on
MGI](https://jiang-junyao.github.io/CACIMAR/build_db.html)

## Citation

[CACIMAR: Cross-species Analysis of Cell Identities, Markers,
Regulations and Interactions Using Single-cell RNA Sequencing Data.
(preprint)](https://www.biorxiv.org/content/10.1101/2024.01.23.576964v1)

## Contact

If you find any bugs, please feel free to e-mail me (jyjiang@link.cuhk.edu.hk)
