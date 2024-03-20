# Checks that parameter estimates work
# library(multistateQTL); library(testthat)
# source("setup.R"); source("test-estimate.R")

test_that("estimate parameters works for gtex data", {

    params <- qtleEstimate(gtex, thresh_sig = 0.05, thresh_null = 0.5)
    params

    expect_equal(round(params$cv.sig.shape, 4), 7.0706)
})
