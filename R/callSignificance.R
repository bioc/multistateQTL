#' Call associations as significant in each state.
#'
#' @details
#' This function adds a new assay to multistateQTL object with TRUE/FALSE
#' significance calls for each test for each state.
#'
#' @param object A \linkS4class{QTLExperiment} object.
#' @param thresh Significance threshold.
#' @param second.thresh Significance threshold for associations with significance in one state.
#' @param feature rowData column name with feature identifier.
#' @param assay assay name with significance score.
#' @param value Value to place in getSignificance
#' @param mode Method to determine significance threshold per state. Options are
#'   `simple`, `feature-wise-FDR`, and `global-FDR`.
#' @param p.adjust.method Method of multiple-test correction if mode != simple
#' @param ... arguments passed to \code{callSignificance}
#'
#' @export
#' @rdname callSignificance
#'
#' @return
#' A `QTLExperiment` object with a new assay called `significant` and with
#' a column called nSignificant added to the colData.
#'
#' @author
#' Christina B Azodi
#'
#'
#' @seealso
#' \code{\link{wilcox.test}}, on which this function is based.
#'
#' @examples
#' qtle <- mockQTLE()
#' 
#' assays(qtle)
#' qtle <- callSignificance(qtle)
#' 
#' # There is now an assay called 'significant'
#' assays(qtle)
#' 
#' # Use feature-wise FDR correction -------------------------------------------
#' qtle_feat <- callSignificance(qtle, thresh=0.1, mode="feature-wise-FDR")
#' 
#' @export
#' @name callSignificance
NULL

#' @import QTLExperiment
#' @importFrom SummarizedExperiment assay assay<-
#'
#' @rdname callSignificance
#'
#' @noRd
.callSignificance <- function(
    object, thresh = 0.05,
    second.thresh = thresh,
    feature = .feature_id,
    assay = "pvalues",
    mode = "simple",
    p.adjust.method = "fdr"){
    
    
    mode <- match.arg(mode, choices = c("simple", "feature-wise-FDR", "global-FDR"))
    p.adjust.method <- match.arg(
        p.adjust.method,
        choices = c("fdr", "holm", "hochberg", "hommel", "bonferroni",
            "BH", "BY", "none"))

    state_thresh <- .get_thresh_per_state(
        object, mode, thresh, p.adjust.method, assay)
    object$significance_threshold <- state_thresh[1, ]
    significant <- as.matrix(assay(object, assay)) <= state_thresh

    # If using less strict threshold for tests with at least one sig hit
    if(!is.null(second.thresh) & second.thresh != thresh){

        state_2thresh <- .get_thresh_per_state(
            object, mode, second.thresh, p.adjust.method, assay)
        object$significance_threshold2 <- state_2thresh[1, ]
        significant.second <- as.matrix(assay(object, assay)) <= state_2thresh

        first.sig <- which(rowSums(significant) > 0)
        significant[first.sig, ] <- significant.second[first.sig, ]
    }

    assay(object, "significant") <- significant
    object$nSignificant <- colSums(significant)

    return(object)

}


#' @export
#' @rdname callSignificance
#' @importFrom QTLExperiment QTLExperiment
setMethod("callSignificance", "QTLExperiment", 
    function(
        object,
        thresh = 0.05,
        second.thresh = thresh,
        feature = .feature_id,
        assay = "pvalues",
        mode = "simple",
        p.adjust.method="fdr", ...) {
    .callSignificance(object=object, thresh=thresh, second.thresh=second.thresh,
        feature=feature, assay=assay, mode=mode,
        p.adjust.method=p.adjust.method, ...)
})


###########################################################
# Internal functions
###########################################################


#' Get the threshold for each state
#'
#' @param object A \linkS4class{QTLExperiment} object.
#' @param mode Method to determine significance threshold per state
#' @param thresh Significance threshold
#' @param p.adjust.method Method of multiple-test correction if mode != simple
#' @param assay Name of assay with test statistic information to use.
#'
#' @importFrom SummarizedExperiment assay
#' @import QTLExperiment
#' 
#' @noRd
.get_thresh_per_state <- function(object, mode, thresh, p.adjust.method, assay){

    if (mode == "simple") {
        state_thresh <- rep(thresh, ncol(object))
        names(state_thresh) <- colnames(object)

    } else if (mode == "global-FDR"){
        state_thresh <- apply(assay(object, assay), 2, function(x)
            .get_corrected_thresh(
                x,
                thresh = thresh,
                p.adjust.method = p.adjust.method))
    } else if (mode == "feature-wise-FDR") {
        top <- getTopHits(object, mode="state", assay=assay)
        state_thresh <- apply(assay(top, assay), 2, function(x)
            .get_corrected_thresh(
                x,
                thresh = thresh,
                p.adjust.method = p.adjust.method))
    }

    state_thresh <- t(as.data.frame(state_thresh))
    state_thresh <- state_thresh[rep(seq_len(nrow(state_thresh)),
                                     each = nrow(object)), ]
    return(state_thresh)
}

#' Get corrected significance threshold for state
#' @param x One column from a \linkS4class{QTLExperiment} object.
#' @param thresh Significance threshold
#' @param p.adjust.method Method of multiple-test correction if mode != simple
#'
#' @importFrom stats p.adjust
#' @import QTLExperiment
#' 
#' @noRd
.get_corrected_thresh <- function(x, thresh, p.adjust.method){

    sig <- data.frame(list(raw = unname(x),
        adj = unname(p.adjust(x, method = p.adjust.method))))
    sig <- sig[sig$adj <= thresh, ]

    if (nrow(sig)==0) {emp_sig_thresh <- 0
    } else{ emp_sig_thresh <- max(sig$raw, na.rm = TRUE)}

    emp_sig_thresh

}

