#' Return multistateQTL object with NAs filled in
#'
#' @description
#' A convenience function for imputing or filling in NAs in a `QTLExperiment`
#' object.
#'
#'
#' @param object A `QTLExperiment` object
#' @param methods A named list with the method for each assay. Available methods
#'                are to replace with a given integer or with the row mean or
#'                median.
#' @param verbose logical. Whether to print progress messages.
#'
#' @importFrom SummarizedExperiment assays assays<-
#' @importFrom matrixStats rowMedians
#'
#' @export
#'
replaceNAs <- function(object,
                       methods=list(betas=0, errors="mean",
                                    pvalues=1, lfsrs=1), verbose=FALSE){
    if(!all(names(assays(object)) %in% names(methods))){
        stop("Provide replacement method for all assays...")
    }

    for(a in names(assays(object))){
        if(is.numeric(methods[[a]])){
            if (verbose) {message("Replacing NAs in ", a, " with ", methods[[a]], "...")}
            assays(object)[[a]][is.na(assays(object)[[a]])] <-  methods[[a]]

        } else if(methods[a] == "mean"){
            if (verbose) {message("Replacing NAs in ", a, " with the row mean...")}
            tmp <- assays(object)[[a]]
            k <- which(is.na(tmp), arr.ind=TRUE)
            assays(object)[[a]][k] <- rowMeans(tmp, na.rm=TRUE)[k[, 1]]

        } else if(methods[a] == "median"){
            if (verbose) {message("Replacing NAs in ", a, " with the row median...")}
            tmp <- assays(object)[[a]]
            k <- which(is.na(tmp), arr.ind=TRUE)
            assays(object)[[a]][k] <- rowMedians(tmp, na.rm=TRUE)[k[, 1]]

        } else{
            stop("Please specify valid replacement method for ", a)
        }
    }

    return(object)
}
