#' Return multistateQTL object with NAs filled in
#'
#' @description
#' A convenience function for imputing or filling in NAs in a `QTLExperiment`
#' object.
#'
#'
#' @param object A `QTLExperiment` object
#' @param methods A named list with the method for each assay. Available methods
#'   are to replace with a given integer or with the row mean or median.
#' @param verbose logical. Whether to print progress messages.
#' 
#' @return A `QTLExperiment` object with the same dimensions as the original object,
#'   but with the NA values replaced according to the input specifications.
#'  
#' @examples 
#' 
#' #' # Create a QTLExperiment object with NA values ------------------------------
#' qtle <- mockQTLE()
#'     
#' # Randomly remove 1000 elements from the betas matrix.
#' na_pattern <- sample(seq(1, ncol(qtle)*nrow(qtle)), 1000)
#' qtle_na <- qtle
#' assay(qtle_na, "betas")[na_pattern] <- NA
#' 
#' # There are some NA values in the "betas" assay
#' any(is.na(betas(qtle_na)))
#' 
#' qtle_complete <- replaceNAs(qtle_na)
#' 
#' # Now we don't have any NA values
#' any(is.na(betas(qtle_complete)))
#' 
#' ## Specify a specific method to impute NAs ----------------------------------
#' 
#' qtle_median <- replaceNAs(
#'     qtle_na, 
#'     methods=list(betas = 0, errors = "median", pvalues = 1), 
#'     verbose=TRUE)
#' 
#' @importFrom SummarizedExperiment assays assays<-
#' @importFrom matrixStats rowMedians
#'
#' @export
#'
replaceNAs <- function(object,
    methods=list(betas=0, errors="mean", pvalues=1, lfsrs=1), verbose=FALSE){
    
    if ( !is(object, "QTLExperiment") )
        stop("Object must be a QTLExperiment")
    
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
