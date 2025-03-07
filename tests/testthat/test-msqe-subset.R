# Test that the S4 subset method works on all fields in the rowData
# library(multistateQTL); library(testthat)
# source("setup.R")


qtle_sig <- callSignificance(qtle)


test_that("Subsetting based on feature_ids works", {
    qtle_top <- getTopHits(qtle_sig, assay="lfsrs", mode="state", verbose = TRUE)

    expect_equal(
        dim(subset(qtle, feature_id == "geneC")),
        dim(qtle[feature_id(qtle) == "geneC",])
    )
})


test_that("Subsetting based on multistate category works", {
    
    sim_unique <- subset(sim, QTL == "unique")
    expect_true(all(rowData(sim_unique)$QTL == "unique"))
    
    sim_group1 <- subset(sim, multistateGroup == "Group1")
    expect_true(all(rowData(sim_group1)$multistateGroup == "Group1"))
    
})


test_that("subsetting based on qtl_type_simple works", {
    sim_top <- runPairwiseSharing(sim_top)
    sim_top <- runTestMetrics(sim_top)
    
    expect_identical(
        subset(sim_top, qtl_type_simple == "multistate"),
        sim_top[rowData(sim_top)$qtl_type_simple == "multistate", ]
    )
})

