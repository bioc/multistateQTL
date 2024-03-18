# Checks that parameter estimates work
# library(multistateQTL); library(testthat)
# source("setup.R"); source("test-estimate.R")

input_path <- system.file("extdata", package="multistateQTL")
state <- c("lung", "thyroid", "spleen", "blood")

input <- data.frame(
    state=state,
    path=paste0(input_path, "/GTEx_tx_", state, ".tsv"))

gtex <- sumstats2qtle(
    input,
    feature_id="molecular_trait_id",
    variant_id="rsid",
    betas="beta",
    errors="se",
    pvalues="pvalue",
    verbose=TRUE)


test_that("estimate parameters works for gtex data", {

    params <- qtleEstimate(gtex, thresh_sig = 0.05, thresh_null = 0.5)
    params

    expect_equal(round(params$cv.sig.shape, 4), 7.0706)
})


test_that("plot simulation params works", {

    params <- qtleEstimate(gtex, thresh_sig = 0.05, thresh_null = 0.5)
    params

    p1 <- plotSimulationParams(params=params)

    expect_equal(class(p1), c("gg", "ggplot"))
})
