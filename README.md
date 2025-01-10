# multistateQTL R Package

multistateQTL is an R package for applying basic statistical tests, summarizing, and visualizing QTL summary statistics from multiple states (e.g., tissues, celltypes, environmental conditions). It works on the `QTLExperiment` (`QTLE`) object class (available in Bioconductor, [QTLExperiment](https://bioconductor.org/packages/release/bioc/html/QTLExperiment.html)), where rows represent features (e.g., genes, transcripts, genomic regions), columns represent states, and assays are the various summary statistics.

|                |               |
| -------------- | ------------- |
| Project Status | [![Project Status.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip) |


## Installation and Usage

To install from GitHub, use devtools:

```
install.packages("devtools")
devtools::install_git("https://github.com/dunstone-a/multistateQTL.git", build_vignettes = TRUE)
```

From Bioconductor, use:

```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("multistateQTL")
```

Using the most recent version of R is strongly recommended (R 4.4 at the time of writing). 

There are several other packages from CRAN and Bioconductor that multistateQTL uses, so you will need to have these packages installed as well. The CRAN packages should install automatically when multistateQTL is installed, but you will need to install the Bioconductor packages manually.

It is also necessary to install the package QTLExperiment. 

### QTLExperiment

```{r install-qtlexperiment}
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("QTLExperiment")
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

Written by Christina Brady Del Azodi, November 2022.
Revised by Amelia Dunstone, April 2024
