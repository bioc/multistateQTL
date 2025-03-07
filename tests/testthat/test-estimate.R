# Checks that parameter estimates work
# library(multistateQTL); library(testthat)
# source("setup.R"); source("test-estimate.R")

test_that("estimate parameters works for gtex data", {
    
    # Parameters used for the simulation
    params <- list(
        betas.sig.shape = 2.160473,
        betas.sig.rate = 14.19137,
        cv.sig.shape = 16.6759,
        cv.sig.rate = 35.56338,
        betas.null.shape = 1.053135,
        betas.null.rate = 22.55316,
        cv.null.shape = 1.698207,
        cv.null.rate = 0.7553863)
    
    qtleSimulate

    params <- qtleEstimate(gtex, threshSig = 0.05, threshNull = 0.5)
    params

    expect_equal(round(params$cv.sig.shape, 2), 7.07)
})

