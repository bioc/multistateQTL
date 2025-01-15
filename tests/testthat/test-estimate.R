# Checks that parameter estimates work
# library(multistateQTL); library(testthat)
# source("setup.R"); source("test-estimate.R")

test_that("estimate parameters works for gtex data", {

    params <- qtleEstimate(gtex, threshSig = 0.05, threshNull = 0.5)
    params

    expect_equal(round(params$cv.sig.shape, 2), 7.07)
})
