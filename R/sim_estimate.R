#' @title Estimate parameters from real data for simulating multi-state QTL summary statistics 
#' 
#' @param data A `QTLExperiment` object or named list containing "betas" and
#'            "errors" matrices.
#' @param assay Assay containing test statistic information to use.
#' @param thresh_sig Max threshold (pval/lfsr) for calling tests as significant.
#' @param thresh_null Min threshold (pval/lfsr) for calling tests as null.  
#' @param verbose Logical.
#'
#' @details The simulation consists of user defined number of equal numbers of four different
#' types of effects: null, equal among conditions, present only in
#' first condition, independent across conditions
#'
#' @examples
#' default_params <- qtleParams()
#' qtleEstimate(mockQTLE())
#' 
#' @name qtleEstimate
#' @rdname qtle_simulations
#' 
#' @importFrom stats rnorm
#' @importFrom QTLExperiment mockQTLE
#'
#' @export
qtleEstimate <- function(data, assay="pvalues", 
                         thresh_sig=0.001, thresh_null=0.1, verbose=TRUE) {
  UseMethod("qtleEstimate")
}


#' @name qtleEstimate
#' @rdname qtle_simulations
#' @export
#' 
qtleEstimate.QTLExperiment <- function(data, assay="pvalues", 
                                       thresh_sig=0.001, thresh_null=0.1,
                                       verbose=TRUE) {
  
  if (assay %in% names(assays(data))) {
    data <- list(betas=assay(data, "betas"), 
                 errors=assay(data, "errors"),
                 test.statistics=assay(data, assay))
  } else{
    data <- list(betas=assay(data, "betas"), 
                 errors=assay(data, "errors"))
  }

  qtleEstimate(data, thresh_sig, thresh_null)
}


#' @name qtleEstimate
#' @rdname qtle_simulations
#' 
#' @export
#' @importFrom mashr mash_1by1
qtleEstimate.list <- function(data, assay="pvalues", thresh_sig=0.01, 
                              thresh_null=0.1, verbose=TRUE){
  
  if("test.statistics" %in% names(data)){
    betas <- abs(data[["betas"]])
    cvs <- data[["errors"]] / betas
    test.statistic <- data[["test.statistics"]]
  } else{
    x <- mash_1by1(mash_set_data(data[["betas"]], data[["errors"]]))
    betas <- abs(x$result$PosteriorMean)
    cvs <- x$result$PosteriorSD / betas
    test.statistic <- x$result$lfsr
  }
  
  index_sig <- which(test.statistic <= thresh_sig, arr.ind = TRUE)
  index_null <- which(test.statistic >= thresh_null, arr.ind = TRUE)

  sig_cvs <- cvs[index_sig]
  sig_betas <- betas[index_sig]
  
  if(length(sig_cvs) > 1e5){
    message("Down-sampling significant examples for faster estimation.")
    keep <- sample(1:length(sig_cvs), 1e5, replace = FALSE)
    sig_cvs <- sig_cvs[keep]
    sig_betas <- sig_betas[keep]
  }
  
  null_cvs <- cvs[index_null]
  null_betas <- betas[index_null]
  
  if(length(null_cvs) > 1e5){
    message("Down-sampling null examples for faster estimation.")
    keep <- sample(1:length(null_cvs), 1e5, replace = FALSE)
    null_cvs <- null_cvs[keep]
    null_betas <- null_betas[keep]
  }
  
  
  
  params <- .fit_gamma(sig_cvs, "cv.sig")
  params <- c(params, .fit_gamma(null_cvs, "cv.null"))
  params <- c(params, .fit_gamma(sig_betas, "betas.sig"))
  params <- c(params, .fit_gamma(null_betas, "betas.null"))
  
  return(params)
}

#' @importFrom fitdistrplus fitdist
.fit_gamma <- function(x, name="str") {
  
  l <- length(x)
  if(l <= 10){stop("Not enough data provided to accuraltey estimate parameters
  for ", name, ". Provide more data or change the significance threshold.")}
  
  if(l < 100){warning("Parameters for ", name, " are being estimated from ",
                      l, " examples, consider adding more data or changing the
                      significance threshold to get more robust estimates.")}
  

  fit <- fitdist(x, "gamma", method = "mge", gof = "CvM", lower = c(0, 0))

  if (fit$convergence > 0) {
    warning("Fitting ", name, "using the Goodness of Fit method failed,",
            " using the Method of Moments instead")
    fit <- fitdist(x, "gamma", method = "mme", lower = c(0, 0))
  }
  
  params <- list()
  params[[paste0(name, ".shape")]] = fit$estimate[["shape"]]
  params[[paste0(name, ".rate")]] = fit$estimate[["rate"]]
  
  return(params)
}


#' @title Default qtle simulation parameters
#' 
#' @export
qtleParams <- function() {
  list(betas.sig.shape = 6.020092,
       betas.sig.rate = 9.977374,
       cv.sig.shape = 3.387213,
       cv.sig.rate = 35.28081,
       betas.null.shape = 3.143392,
       betas.null.rate = 11.57387,
       cv.null.shape = 12.29935,
       cv.null.rate = 75.69514)
}


