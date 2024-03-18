# Checks the state and test wise summarize functions
# library(multistateQTL); library(testthat)
# source("setup.R"); source("test-get-significant.R")

msqe <- qtle

test_that("getSignificant requires callSignificance", {
    expect_error(getSignificant(msqe))
})

test_that("getSignificant percentage inputs work", {
    percent <- 0.2
    nStates <- ncol(qtle)

    sim <- callSignificance(msqe)

    expect_equal(
        getSignificant(sim, n = 0.2),
        getSignificant(sim, n = 0.2*nStates))
})

test_that("getTopPerFeature selects the correct snp-feature pair", {
    sim <- callSignificance(msqe)
    sim_top <- getTopHits(sim, assay="lfsrs", mode="state", verbose = TRUE)
    expect_equal(
        min(pvalues(sim_top)[, 1]),
        min(pvalues(subset(msqe, feature_id == "geneC"))[, 1]))
})

test_that("getTopPerFeature works with global method", {
    sim <- callSignificance(msqe)
    sim_top <- getTopHits(sim, assay="lfsrs", mode="global", verbose = TRUE)
})


test_that("getSignificant in simple mode works with one and two thresholds", {

    test <- callSignificance(msqe, thresh=0.1, mode="simple")
    test <- getSignificant(test)

    expect_equal(rowSums(pvalues(test) <= 0.1), rowSums(assay(test, "significant")))

    test <- callSignificance(msqe, thresh=0.1, second.thresh=0.5, mode="simple")
    test <- getSignificant(test)

    expect_true(all(rowSums(assay(test, "significant")) >= rowSums(pvalues(test) <= 0.1)))
})

test_that("getSignificant argument matching works", {

    # Test that argument matching works for `mode`
    expect_equal(
        class(callSignificance(msqe, thresh = 0.1, mode = "simpl"))[1],
        "QTLExperiment")
    expect_equal(
        class(callSignificance(msqe, thresh = 0.1, mode = "glob"))[1],
        "QTLExperiment")
    expect_error(
        class(callSignificance(msqe, thresh = 0.1, mode = "featurewise")))

    # Test that argument matching works for `p.adjust.method``
    expect_equal(
        class(callSignificance(msqe, mode = "simple", p.adjust.method = "bonf"))[1],
        "QTLExperiment")
    expect_error(
        class(callSignificance(msqe, mode = "simple", p.adjust.method = "ho")))
})

test_that("getSignificant global and feature-wise modes works", {
    # Test that global-FDR works
    test_glob <- callSignificance(msqe, thresh=0.1, mode="global-FDR")
    test_glob <- getSignificant(test_glob)

    # Test that feature-wise-FDR works
    test_feat <- callSignificance(msqe, thresh=0.1, mode="feature-wise-FDR")
    test_feat <- getSignificant(test_feat, verbose = TRUE)

    expect_true(all(test_glob$significance_threshold <= test_feat$significance_threshold))

    test_bonf <- callSignificance(msqe, thresh=0.1, p.adjust.method="bonferroni")
    test_bonf <- getSignificant(test_bonf)
    test_fdr <- callSignificance(msqe, thresh=0.1)
    test_fdr <- getSignificant(test_fdr)
    expect_true(all(test_bonf$significance_threshold <= test_fdr$significance_threshold))

})

