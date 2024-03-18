# Checks that NAs can be removed or imputed.
# library(testthat); library(multiStateQTLExperiment)
# source("test-replaceNAs.R")

set.seed(1234)
sim <- qtleSimulate(
    nstates=10, nfeatures=100, ntests=1000,
    global=0.2, multi=0.4, unique=0.2, k=2)

# Add some NA values to the betas, errors and lfsrs
na_pattern <- sample(seq(1, ncol(sim)*nrow(sim)), 1000)

sim_na <- sim

assay(sim_na, "betas")[na_pattern] <- NA
assay(sim_na, "errors")[na_pattern] <- NA
assay(sim_na, "lfsrs")[na_pattern] <- NA

test_that("qtle can be subset to only rows with mostly complete entries", {

    # This needs a seed!!
    # Have data for at least half of the states
    expect_message(
        getComplete(sim_na, n=0.5, verbose=TRUE),
        "Removing 4 tests with NAs in >= 5 states...")

    sim_complete <- getComplete(sim_na, n=0.5, verbose=TRUE)

    # This is hard coded!!
    expect_equal(nrow(sim_complete), nrow(sim_na) - 4)

    # This is simplified
    expect_true(nrow(sim_complete) <= nrow(sim_na))
})

test_that("NAs can be replaced", {
    sim_comp <- getComplete(sim_na, n=0.5, verbose=TRUE)
    sim_final <- replaceNAs(sim_comp, verbose=TRUE)

    # Dimensions should be the same
    expect_equal(dim(sim_comp), dim(sim_final))

    # Shouldn't be any NA values
    expect_false(
        any(c(
            is.na(assay(sim_final, "betas")),
            is.na(assay(sim_final, "errors")),
            is.na(assay(sim_final, "lfsrs")))))
})

test_that("NAs in betas and pvalues are replaced with constants", {
    sim_comp <- getComplete(sim_na, n=0.5)
    sim_final <- replaceNAs(sim_comp)

    # Indices of NA values for first state
    locations <- is.na(betas(sim_comp)[,1])

    # betas get set to 0
    expect_true(all(betas(sim_final)[locations,1] == 0))

    # lfsrs get set to 1
    expect_true(all(lfsrs(sim_final)[locations,1] == 1))
})

test_that("NAs in errors can be replaced with the mean", {
    sim_comp <- getComplete(sim_na, n=0.5)
    sim_final <- replaceNAs(sim_comp)

    # Indices of NA values for first state
    locations <- is.na(errors(sim_comp)[,1])

    means <- apply(
        errors(sim_comp),
        MARGIN = 1,
        FUN = function(x) {mean(x, na.rm = TRUE)})

    # errors get set to the mean
    expect_true(all(errors(sim_final)[locations,1] == means[locations]))
})

test_that("NAs in errors can be replaced with the median", {
    sim_comp <- getComplete(sim_na, n=0.5)
    sim_final <- replaceNAs(
        sim_comp,
        methods=list(betas = 0, errors = "median", lfsrs = 1),
        verbose=TRUE)

    # Indices of NA values for first state
    locations <- is.na(errors(sim_comp)[,1])

    medians <- apply(
        errors(sim_comp),
        MARGIN = 1,
        FUN = function(x) {median(x, na.rm = TRUE)})

    # errors get set to the mean
    expect_true(all(errors(sim_final)[locations,1] == medians[locations]))
})


test_that("NAs can be replaced with specified constants", {
    sim_comp <- getComplete(sim_na, n=0.5)
    sim_final <- replaceNAs(
        sim_comp,
        methods = list(betas = -0.5, errors = "median", lfsrs = 0.5))

    # Indices of NA values for first state
    locations <- is.na(betas(sim_comp)[,1])

    # betas get set to a constant
    expect_true(all(betas(sim_final)[locations,1] == -0.5))

    # lfsrs get set to a constant
    expect_true(all(lfsrs(sim_final)[locations,1] == 0.5))
})

test_that("NA replacement error messages work", {
    sim_comp <- getComplete(sim_na, n=0.5)

    # Specify a method for every assay
    expect_error(
        replaceNAs(sim_comp, methods = list(errors = "median")))

    # Method to impute error assay doesn't exist
    expect_error(
        replaceNAs(sim_comp, methods = list(betas=0, errors="med", lfsrs=1)))


})
