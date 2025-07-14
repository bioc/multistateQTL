## Changes in version 2.1.1 (2025-07-14)
* Moved data.table to imports instead of depends.
* Updated testing to avoid error with ggplot2 version 4.0.0.

## Changes in version 1.99.5 (2025-03-17)
* Added support for additional arguments in plotting functions to be passed to ComplexHeatmap.

## Changes in version 1.99.0 (2025-03-07)
* Patches to make multistateQTL compatible with major change to internal representations in QTLExperiment.

## Changes in version 1.3.2 (2025-01-28)
* Fix up a test with a rounding error
* Remove bug in getComplete where setting verbose = TRUE would lead to an error
* Update replaceNAs test

## Changes in version 1.3.1 (2025-01-14)
* Fix up a test that was failing due to machine tolerance

## Changes in version 1.1.0 (2024-10-23)
* Brought back in runSignificantFeatures function
* Added a new colour palette for heatmaps needing more than 10 colours
* More tests for plotting

## Changes in version 0.99.1 (2024-04-21)
* Small tweaks e.g. using seq_along()
* Changed some arguments to use camelCase instead of snake_case
* Improved vignette.

## Changes in version 0.99.0 (2024-03-22)
* Submitted to Bioconductor