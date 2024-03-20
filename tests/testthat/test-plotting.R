# Checks that plots work
# library(multistateQTL); library(testthat)
# source("setup.R"); source("test-plots.R")

test_that("plot simulation params works", {

    params <- qtleEstimate(qtle, thresh_sig = 0.05, thresh_null = 0.5)

    p1 <- plotSimulationParams(params=params)

    expect_equal(class(p1), c("gg", "ggplot"))
})

test_that("plot compare states works", {
    sim_top <- runPairwiseSharing(sim_top)
    sim_top <- runTestMetrics(sim_top)

    p1 <- plotCompareStates(sim_top, x="S01", y="S02")

    expect_equal(class(p1[[1]]), c("gg", "ggplot"))
    expect_equal(class(p1[[2]]), "table")
})


test_that("plot pairwise sharing errors work", {
    expect_error(plotPairwiseSharing(summary(qtle)))

    expect_error(plotPairwiseSharing(sim))
})

test_that("plot pairwise sharing works", {
    sim_top <- runPairwiseSharing(sim_top)

    p1 <- plotPairwiseSharing(sim_top)

    expect_output(print(class(p1)), "Heatmap")
})


test_that("produce pairwise sharing plots with complex column annotations", {

    sim_top <- runPairwiseSharing(sim_top)

    p1 <- plotPairwiseSharing(sim_top)

    expect_output(print(class(p1)), "Heatmap")

    p2 <- plotPairwiseSharing(sim_top, annotate_cols = c("nSignificant", "multistateGroup"))

    expect_output(print(class(p2)), "Heatmap")
})


test_that("produce upset plots with complex row annotations", {
    p1 <- plotUpSet(sim_top, annotate_by = c("nSignificant", "multistateGroup"))
    expect_output(print(class(p1)), "Heatmap")
})

test_that("upset plot errors work", {

    # plotUpset needs significance assay
    expect_error(plotUpSet(qtle))
})


test_that("plot errors work", {
    sim_top <- runPairwiseSharing(sim_top)
    sim_top <- runTestMetrics(sim_top)
    sim_top_ms <- subset(sim_top, qtl_type_simple == "multistate")

    # plotQTLClusters needs a QTLExperiment object
    expect_error(plotQTLClusters(summary(qtle)))

    # plotQTLClusters needs fill_by to be in assays()
    expect_error(
        plotQTLClusters(
            sim_top_ms,
            annotate_states = c("multistateGroup"),
                    annotate_tests = c("qtl_type", "mean_beta", "QTL"),
            fill_by = "beta"))

    # plotCompareStates
    expect_error(plotCompareStates(qtle, x="S01", y="S02"))
})

test_that("plot QTL clusters works", {

    sim_top <- runPairwiseSharing(sim_top)
    sim_top <- runTestMetrics(sim_top)
    sim_top_ms <- subset(sim_top, qtl_type_simple == "multistate")

    p1 <- plotQTLClusters(sim_top_ms, annotate_states = c("multistateGroup"),
                annotate_tests = c("qtl_type", "mean_beta", "QTL"))

    expect_output(print(class(p1)), "Heatmap")
})






