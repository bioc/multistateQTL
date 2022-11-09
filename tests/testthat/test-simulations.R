# Checks for proper functioning of the reducedDim methods.
# library(testthat); library(multiStateQTLExperiment)
# source("test-simulations.R")

data <- readRDS("data/msqe_gtex.rds")
data <- data[sample(1:nrow(data), 1000), ]
msqe <- mockQTLE(nQTL = 10000)
msqe2 <- msqe

test_that("estimate parameters works", {
  params <- msqeEstimate(data, thresh = 0.1)
  params <- msqeEstimate(msqe, thresh = 0.1)
  

})
