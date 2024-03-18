# Checks the complexheatmap plots
# library(multistateQTL); library(testthat)
# source("setup.R"); source("test-plotting-complexheatmaps.R")

sim <- qtleSimulate(nstates=10, nfeatures=100, ntests=1000,
                    global=0.2, multi=0.4, unique=0.2, k=2)
sim <- callSignificance(sim, assay="lfsrs", thresh=0.001)

sim_sig <- getSignificant(sim)
sim_top <- getTopHits(sim_sig, assay="lfsrs", mode="state")

test_that("plot pairwise sharing errors work", {
    expect_error(plotPairwiseSharing(summary(mock)))

    expect_error(plotPairwiseSharing(sim))
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
    expect_error(plotUpSet(mock))
})


