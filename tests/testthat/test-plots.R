# Checks that plots work
# library(multistateQTL); library(testthat)
# source("setup.R"); source("test-plots.R")

sim <- qtleSimulate(nstates=10, nfeatures=100, ntests=1000,
                    global=0.2, multi=0.4, unique=0.2, k=2)
sim <- callSignificance(sim, mode="simple", assay="lfsrs",
                        thresh=0.0001, second.thresh=0.0002)


test_that("plot simulation params works", {

    params <- qtleEstimate(mock, thresh_sig = 0.05, thresh_null = 0.5)

    p1 <- plotSimulationParams(params=params)

    expect_equal(class(p1), c("gg", "ggplot"))
})

test_that("plot compare states works", {

    sim_sig <- getSignificant(sim)
    sim_top <- getTopHits(sim_sig, assay="lfsrs", mode="state")
    sim_top <- runPairwiseSharing(sim_top)
    sim_top <- runTestMetrics(sim_top)

    p1 <- plotCompareStates(sim_top, x="S01", y="S02")

    expect_equal(class(p1[[1]]), c("gg", "ggplot"))
    expect_equal(class(p1[[2]]), "table")
})

test_that("plot pairwise sharing works", {
    sim_sig <- getSignificant(sim)
    sim_top <- getTopHits(sim_sig, assay="lfsrs", mode="state")
    sim_top <- runPairwiseSharing(sim_top)

    plotPairwiseSharing(sim_top)
})

test_that("plot errors work", {

    sim_sig <- getSignificant(sim)
    sim_top <- getTopHits(sim_sig, assay="lfsrs", mode="state")
    sim_top <- runPairwiseSharing(sim_top)
    sim_top <- runTestMetrics(sim_top)
    sim_top_ms <- subset(sim_top, qtl_type_simple == "multistate")

    # plotQTLClusters needs a QTLExperiment object
    expect_error(plotQTLClusters(summary(mock)))

    # plotQTLClusters needs fill_by to be in assays()
    expect_error(
        plotQTLClusters(
            sim_top_ms,
            annotate_states = c("multistateGroup"),
                    annotate_tests = c("qtl_type", "mean_beta", "QTL"),
            fill_by = "beta"))

    # plotCompareStates
    expect_error(plotCompareStates(mock, x="S01", y="S02"))
})

test_that("plot QTL clusters works", {
    sim_sig <- getSignificant(sim)
    sim_top <- getTopHits(sim_sig, assay="lfsrs", mode="state")
    sim_top <- runPairwiseSharing(sim_top)
    sim_top <- runTestMetrics(sim_top)
    sim_top_ms <- subset(sim_top, qtl_type_simple == "multistate")

    p1 <- plotQTLClusters(sim_top_ms, annotate_states = c("multistateGroup"),
                annotate_tests = c("qtl_type", "mean_beta", "QTL"))

    expect_output(print(class(p1)), "Heatmap")

    # Test fill by "pvalues" TODO
})
