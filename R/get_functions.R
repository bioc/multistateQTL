#' @title Filter QTLExperiment based on missing data
#'
#' @description
#' Method to filter \linkS4class{QTLExperiment} objects to remove tests with
#' greater than the permitted rate of missing values.
#'
#' @param qtle A `QTLExperiment` object
#' @param n Number (or percent if n < 1) of states requiring non-null values
#' @param verbose logical. Whether to print progress messages.
#'
#' @return a subset of the `QTLExperiment` object, with only tests with fewer
#' NAs than specified by n.
#' 
#' @examples
#' 
#' # Create a QTLExperiment object with NA values ------------------------------
#' sim <- qtleSimulate(
#'     nstates=10, nfeatures=100, ntests=1000,
#'     global=0.2, multi=0.4, unique=0.2, k=2)
#'     
#' # Randomly remove 1000 elements from the betas matrix.
#' na_pattern <- sample(seq(1, ncol(sim)*nrow(sim)), 1000)
#' sim_na <- sim
#' assay(sim_na, "betas")[na_pattern] <- NA
#' 
#' # Original object has more rows than the output of getComplete()
#' dim(sim_na)
#' 
#' sim_complete <- getComplete(sim_na)
#' dim(sim_complete)
#' 
#' 
#'
#' @importFrom QTLExperiment betas errors
#'
#' @name getComplete
#' @rdname getComplete
#' @export
#'
getComplete <- function(qtle, n = 1, verbose=FALSE){

    if(n < 1){
        n <- round(ncol(qtle) * n, 0)
    }

    keep <- rowSums(is.na(betas(qtle))) < n & rowSums(is.na(errors(qtle))) < n

    if(verbose) { message("Removing ", table(keep)[["FALSE"]],
                          " tests with NAs in >= ", n, " states...")}

    qtle <- qtle[which(keep), ]


    return(qtle)
}

#' Filter to only QTLs significant in at least one state
#'
#' @param qtle `QTLExperiment` object
#' @param n Number (or percent if n < 1) of states with significant association
#' @param assay The assay containing TRUE/FALSE significance calls for each QTL
#'   test.
#' @param verbose logical. Whether to print progress messages.
#' 
#' @return a subset of the `QTLExperiment` object, where all rows are
#'   significant in at least one state.
#'   
#' @examples
#' 
#' qtle <- mockQTLE()
#' 
#' qtle <- callSignificance(qtle)
#' dim(qtle)
#' qtle_sig <- getSignificant(qtle)
#' 
#' # There are fewer rows because we have removed tests which are not significant
#' # in any state.  
#' dim(qtle_sig)
#'
#' @importFrom SummarizedExperiment assay assays
#' 
#' @name getSignificant
#' @rdname getSignificant
#' @export
#'
getSignificant <- function(qtle, n=1,
                           assay = "significant",
                           verbose = FALSE){

    if( ! assay %in% names(assays(qtle)) ) {
        stop("First run callSignificance()...")
    }

    if(n < 1){
        n <- ncol(qtle) * n
    }

    if(verbose) { message("Removing QTL significant in < ", n, " states")}

    keep <- rowSums(assay(qtle, assay)) >= n
    qtle <- qtle[keep, ]

    if(verbose) { message("Number of remaining QTL: ", nrow(qtle)) }

    return(qtle)
}


#'
#' @title Filter QTLExperiment to keep only top hits
#'
#' @description
#' Method to return a subset of a \linkS4class{QTLExperiment} object containing
#' only the tests that are top hits. Top hits are defined as the test for each
#' feature with the most significant test statistic. Returns an array of the top 
#' QTL for each feature across all states
#'
#' @param qtle A `QTLExperiment` object
#' @param assay The assay containing the test statistic to minimize.
#' @param mode global/state to specify if the top hit per feature is desired
#'             from across all states or for each state.
#' @param assay_sig The assay containing TRUE/FALSE significance calls for each
#'                  QTL test.
#' @param verbose logical. Whether to print progress messages.
#'
#' @return A subset of the `QTLExperiment` object, with only tests that are the
#' top hits for each feature (`mode=global`) or for each feature for each
#' state (`mode=state`).
#' 
#' @examples
#' sumstats <- mockSummaryStats(nStates=10, nQTL=100, names=TRUE)
#' qtle <- QTLExperiment(
#'     assay=list(
#'     betas=sumstats$betas,
#'     errors=sumstats$errors,
#'     pvalues=sumstats$pvalues,
#'     lfsrs=sumstats$pvalues))
#' 
#' # Add 'significant' assay to object
#' qtle <- callSignificance(qtle)
#' 
#' # Filter to the top tests for each feature
#' qtle_glob <- getTopHits(qtle, assay="lfsrs", mode="global", verbose = TRUE)
#' # There are 3 rows corresponding to the three features.
#' table(feature_id(qtle_glob))
#' 
#' # At most one QTL is retained for each combination of feature_id and state_id
#' qtle_feat <- getTopHits(qtle, assay="lfsrs", mode="state", verbose = TRUE)
#' table(feature_id(qtle_feat))
#'
#' @importFrom dplyr group_by slice_min %>% filter
#' @importFrom matrixStats rowMins
#' @importFrom tidyr pivot_longer
#' @importFrom collapse fmutate
#'
#' @name getTopHits
#' @rdname getTopHits
#' @export
#'
getTopHits <- function(qtle, mode=c("global", "state"),
                       assay = "pvalues",
                       assay_sig = "significant",
                       verbose = FALSE){

    if(verbose) { message("Selecting top hits per feature from...") }
    keep <- c()

    if(mode == "global"){
        keep <- as.data.frame(list(feature_id = feature_id(qtle),
                                   id = rownames(qtle),
                                   value = rowMins(assay(qtle, assay)))) %>%
        group_by(feature_id) %>%
        slice_min(value, n = 1, with_ties = FALSE)

    } else if (mode == "state"){

        test_statistics <- as.data.frame(assay(qtle, assay))
    if (assay_sig %in% names(assays(qtle))) {
        test_statistics[!assay(qtle, assay_sig)] <- 1
    }
    keep <- test_statistics %>%
        fmutate(id=rownames(qtle),
                feature_id=feature_id(qtle)) %>%
        pivot_longer(-c(feature_id, id)) %>%
        group_by(feature_id, name) %>%
        slice_min(value, n = 1, with_ties = FALSE) %>%
        filter(value < 1)
    }

    keep <- unique(keep$id)

    if(verbose) { message("Total number of top associations: ", length(keep)) }

    return(qtle[rownames(qtle) %in% keep, ])
}



#' @importFrom data.table setDT .SD
.getTopPerFeature <- function(x, by, var, FUN){
    tmp <- setDT(list(x=x, by=by, var=var))
    tmp <- tmp[, .SD[FUN(x)], by = by]
    return(paste0(tmp$by, "|", tmp$var))
}
