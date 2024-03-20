# Checks that QTLs can be categorised into types
# library(multistateQTL); library(testthat)
# source("setup.R"); source("test-runTestMetrics.R")

sim_top <- runPairwiseSharing(sim_top)
sim_top <- runTestMetrics(sim_top)

test_that("run test metrics works", {
    sim_top <- getTopHits(sim_sig, assay="lfsrs", mode="state")
    sim_top <- runPairwiseSharing(sim_top)
    sim_top <- runTestMetrics(sim_top)

    expect_equal(class(sim_top)[1], "QTLExperiment")
})

test_that("run test metrics errors work", {
    expect_error(runTestMetrics(mock))
})
# Add tests that the classifications are correct TODO


