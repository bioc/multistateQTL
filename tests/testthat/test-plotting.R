# Checks that plots work
# library(multistateQTL); library(testthat)
# source("setup.R"); source("test-plots.R")

test_that("plot simulation params works", {

    params <- qtleEstimate(qtle, threshSig = 0.05, threshNull = 0.5)

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

    p2 <- plotPairwiseSharing(sim_top, annotateColsBy = c("nSignificant", "multistateGroup"))

    expect_output(print(class(p2)), "Heatmap")
    
    sim <- qtleSimulate(
        nStates=12, nFeatures=100, nTests=1000,
        global=0.2, multi=0.4, unique=0.2, k=2)
    sim <- callSignificance(sim, mode="simple", assay="lfsrs",
        thresh=0.0001, secondThresh=0.0002)
    sim_sig <- getSignificant(sim)
    sim_top <- getTopHits(sim_sig, assay="lfsrs", mode="state")

    sim_top <- runPairwiseSharing(sim_top)
    
    # More than 10 colours
    p3 <- plotPairwiseSharing(sim_top, annotateColsBy = c("nSignificant", "state_id"))
    expect_output(print(class(p3)), "Heatmap")
})


test_that("produce upset plots with complex row annotations", {
    p1 <- plotUpSet(sim_top, annotateColsBy = c("nSignificant", "multistateGroup"))
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
            annotateColsBy = c("multistateGroup"),
                    annotate_tests = c("qtl_type", "mean_beta", "QTL"),
            fillBy = "beta"))

    # plotCompareStates
    expect_error(plotCompareStates(qtle, x="S01", y="S02"))
})

test_that("plot QTL clusters works", {

    sim_top <- runPairwiseSharing(sim_top)
    sim_top <- runTestMetrics(sim_top)
    sim_top_ms <- subset(sim_top, qtl_type_simple == "multistate")

    p1 <- plotQTLClusters(sim_top_ms, annotateColsBy = c("multistateGroup"),
                annotateRowsBy = c("qtl_type", "mean_beta", "QTL"))

    expect_output(print(class(p1)), "Heatmap")
})






