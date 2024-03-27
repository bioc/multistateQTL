# Checks for proper functioning of the reducedDim methods.
# library(testthat); library(multistateQTL)
# source("test-simulations.R")

test_that("qtle can be simulated", {
    # Simulated QTLE where half of the QTLs are globally significant
    sim <- qtleSimulate(ntests=1000, nstates=6, global=0.5)

    expect_equal(class(sim)[1], "QTLExperiment")
    expect_equal(sum(rowData(sim)$QTL == "global")/nrow(sim), 0.5)

    sim <- qtleSimulate(nstates=10, nfeatures=100, ntests=1000,
        global=0.2, multi=0.4, unique=0.2, k=2)

    expect_equal(class(sim)[1], "QTLExperiment")
    expect_equal(
        as.vector(table(rowData(sim)$QTL)/nrow(sim)),
        c(0.2, 0.4, 0.2, 0.2)) # global, multistate, null, unique
})

test_that("subsetting based on multistate category works", {

    sim_unique <- subset(sim, QTL == "unique")
    expect_true(all(rowData(sim_unique)$QTL == "unique"))

    sim_multi <- subset(sim, QTL == "multistate")
    expect_true(all(rowData(sim_multi)$QTL == "multistate"))

})

test_that("performance metrics work", {
    perf_metrics <- simPerformance(sim)
   expect_equal(class(perf_metrics)[1], "list")
})

test_that("performance metrics errors work", {
    # Remove one state's rowData information
    rowData(sim) <- rowData(sim)[, 1:10]

    expect_error(simPerformance(sim))
})
