# nnSVG-analyses

This repository contains code scripts to reproduce analyses and figures in our paper "nnSVG for the scalable identification of spatially variable genes using nearest-neighbor Gaussian processes".


## nnSVG

`nnSVG` is a method for computationally scalable identification of spatially variable genes (SVGs) in spatially-resolved transcriptomics data.

`nnSVG` is available as an R package from [Bioconductor](https://bioconductor.org/packages/nnSVG) and [GitHub](https://github.com/lmweber/nnSVG).

More details including installation instructions and a tutorial are available on the [Bioconductor](https://bioconductor.org/packages/nnSVG) and [GitHub](https://github.com/lmweber/nnSVG) pages.


## Contents

The code in this repository is organized into subdirectories in the [analyses/](analyses/) directory.


## Data

References to original datasets used in the analyses are provided in the manuscript.

Datasets can be downloaded as [SpatialExperiment](https://bioconductor.org/packages/SpatialExperiment) formatted data objects from the [STexampleData](https://bioconductor.org/packages/STexampleData) Bioconductor package within the preprocessing scripts.

R objects (`.RData` files) containing source data to reproduce individual figures from the manuscript within the analysis and evaluation scripts are available from [Figshare](https://doi.org/10.6084/m9.figshare.23561439.v2).


## Citation

A preprint describing `nnSVG` is available from bioRxiv:

- [Weber L.M. et al. (2022), "nnSVG: scalable identification of spatially variable genes using nearest-neighbor Gaussian processes", bioRxiv](https://www.biorxiv.org/content/10.1101/2022.05.16.492124v2)

