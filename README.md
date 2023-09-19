# multistateQTL R Package

multistateQTL is an R package for applying basic statistical tests, summarizing, and visualizing QTL summary statistics from multiple states (e.g., tissues, celltypes, environmental conditions). It works on the `QTLExperiment` (`QTLE`) object class (available [here](https://gitlab.svi.edu.au/biocellgen-public/qtlexperiment)), where rows represent features (e.g., genes, transcripts, genomic regions), columns represent states, and assays are the various summary statistics.

|                |               |
| -------------- | ------------- |
| Project Status | [![Project Status.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip) |
| Travis CI      | [![Pipeline Status](https://gitlab.svi.edu.au/biocellgen-public/multistateqtle/badges/master/pipeline.svg)](https://gitlab.svi.edu.au/biocellgen-public/multistateqtle/-/commits/master) |
| Test coverage  | [![Coverage Report](https://gitlab.svi.edu.au/biocellgen-public/multistateqtle/badges/master/coverage.svg)](https://gitlab.svi.edu.au/biocellgen-public/multistateqtle/-/commits/master)|
| Latest release  | [![Latest Release](https://gitlab.svi.edu.au/biocellgen-public/multistateqtle/-/badges/release.svg)](https://gitlab.svi.edu.au/biocellgen-public/multistateqtle/-/releases)|


## Installation and Usage

This package in stable but undergoing active development and currently only lives on GitHub. To install from GitHub, use devtools:

```
install.packages("devtools")
devtools::install_git("https://gitlab.svi.edu.au/biocellgen-public/multistateQTL", build_vignettes = TRUE)
```

We plan to submit QTLExperiment and multistateQTL to Bioconductor in the near future. Using the most recent version of R is strongly recommended (R 4.2.1 at the time of writing). 

There are several other packages from CRAN and Bioconductor that multistateQTL uses, so you will need to have these packages installed as well. The CRAN packages should install automatically when multistateQTL is installed, but you will need to install the Bioconductor packages manually.

Not all of the following are strictly necessary, but have been included here as they enhance the functionality of multistateQTL. The commands below should help with package installations.

### QTLExperiment

```{r install-qtlexperiment}
devtools::install_git("https://gitlab.svi.edu.au/biocellgen-public/qtlexperiment.git", build_vignettes = TRUE)
```

### CRAN

```{r install-cran}
install.packages(c("knitr", "dplyr", "collapse", "ggplot2", "circlize", "vroom", "mashr", "tidyr", "fitdistrplus"))
```

### Bioconductor 

```{r load-bioc}
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("SummarizedExperiment", "BiocGenerics", "S4Vectors", "ComplexHeatmap"))
```


## Getting started

Get started with multistateQTL by checking out the vignette. From inside an R session, load multistateQTL and then browse the vignettes:

```
library(multistateQTL)
browseVignettes("multistateQTL")
```

There is a detailed HTML document available that introduces the main features and functionality of multistateQTL.


## Planned package additions

Here is a running list of features that we plan to add to the multistateQTL package:

- Global and feature-wise FDR correction of p-values assay
- conditional FDR (cFDR)


## Disclaimer

The package is currently in a Beta state. The major functionality of the package is settled, but it is still under development so may change from time to time. Please do try it and contact me with bug reports, feedback, feature requests, questions and suggestions to improve the package.

Christina Brady Del Azodi, November 2022
