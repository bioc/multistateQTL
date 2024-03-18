# Checks that QTLs can be categorised into types
# library(multistateQTL); library(testthat)
# source("setup.R"); source("test-runTestMetrics.R")

msqe <- mock
sim <- callSignificance(msqe)

sim <- qtleSimulate(nstates=10, nfeatures=100, ntests=1000,
                    global=0.2, multi=0.4, unique=0.2, k=2)
sim <- callSignificance(sim, mode="simple", assay="lfsrs",
                        thresh=0.0001, second.thresh=0.0002)

sim_sig <- getSignificant(sim)
sim_top <- getTopHits(sim_sig, assay="lfsrs", mode="state")
sim_top <- runPairwiseSharing(sim_top)
sim_top <- runTestMetrics(sim_top)

test_that("run test metrics works", {
    sim <- callSignificance(mock, mode = "simple", assay = "lfsrs")
    sim_sig <- getSignificant(sim)
    sim_top <- getTopHits(sim_sig, assay="lfsrs", mode="state")
    sim_top <- runPairwiseSharing(sim_top)
    sim_top <- runTestMetrics(sim_top)

    expect_equal(class(sim_top)[1], "QTLExperiment")
})

test_that("run test metrics errors work", {
    expect_error(runTestMetrics(mock))
})
# Add tests that the classifications are correct TODO


