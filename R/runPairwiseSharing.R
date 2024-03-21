#' @title Compute the proportion of (significant) QTL shared by pairs of
#' conditions
#'
#' @param qtle A `QTLExperiment` object.
#' @param assay The assay containing the metric used to determine sharing (i.e.
#'              the metric to be within a factor X to be considered shared).
#' @param assay_sig The assay containing significance information.
#' @param factor a number in [0,1] the factor within which effects are
#'               considered to be shared
#' @param FUN a function to be applied to the estimated effect sizes before
#'            assessing sharing. Default 'FUN=identity', 'FUN=abs' ignores
#'            the sign of the effects when assessing sharing.
#' @param ... Additional parameters to pass on to internal functions.
#'
#' @details For each pair of states, the effects that are significant
#' (as determined by `callSignificance`) in at least one of the two states are
#' identified. Then the fraction of those with an estimated effect size
#' (i.e. betas) within a factor `factor` of one another is computed and
#' returned.
#'
#' @return The `QTLExperiment` object with a matrix called pairwiseSharing added
#' to the metadata.
#'
#' @examples
#' m <- mockQTLE()
#' m <- callSignificance(m, assay="pvalues")
#' runPairwiseSharing(m) # sharing by magnitude (same sign)
#' runPairwiseSharing(m, factor=0) # sharing by sign
#' runPairwiseSharing(m, FUN=abs) # sharing by magnitude when sign is ignored
#'
#' @importFrom SummarizedExperiment assay
#' @importFrom S4Vectors metadata metadata<-
#'
#' @name runPairwiseSharing
#' @rdname runPairwiseSharing
#' @export
#'
runPairwiseSharing <- function(qtle,
    assay="betas",
    assay_sig="significant",
    factor=0.5,
    FUN=identity, ...){

    if( ! assay_sig %in% names(assays(qtle)) ) {
        stop("First run callSignificance()...")
    }

    nStates <- ncol(qtle)
    sigMat <- assay(qtle, assay_sig)
    effectMat <- assay(qtle, assay)
    S <- matrix(NA, nrow=nStates, ncol=nStates)

    for(i in 1:nStates){
        for(j in i:nStates){
            sig_i <- which(sigMat[, i])
            sig_j <- which(sigMat[, j])
            sig_union <- union(sig_i, sig_j)

            if(length(sig_union) > 0){
                ratio <- FUN(effectMat[sig_union, i]) / FUN(effectMat[sig_union, j] + 1e-8)
            } else{ ratio <- 0 }

            S[i,j] <- mean(ratio > factor & ratio < (1/factor))
        }
    }

    S[lower.tri(S, diag = FALSE)] <- t(S)[lower.tri(S, diag = FALSE)]
    colnames(S) <- row.names(S) <- colnames(qtle)
    metadata(qtle)$pairwiseSharing <- S

    return(qtle)
}



