# Checks that parameter estimates work
# library(multistateQTL); library(testthat)
# source("setup.R"); source("test-estimate.R")

test_that("estimate parameters works for gtex data", {

    params <- qtleEstimate(gtex, threshSig = 0.05, threshNull = 0.5)
    params

    expect_equal(round(params$cv.sig.shape, 2), 7.07)
    expect_equal(round(params$cv.sig.rate, 2), 7.72)
    expect_equal(round(params$cv.null.shape, 2), 1.90)
    expect_equal(round(params$cv.null.rate, 2), 0.90)
    expect_equal(round(params$betas.sig.shape, 2), 3.47)
    expect_equal(round(params$betas.sig.rate, 2), 16.27)
    expect_equal(round(params$betas.null.shape, 2), 1.36)
    expect_equal(round(params$betas.null.rate, 2), 11.03)
})

