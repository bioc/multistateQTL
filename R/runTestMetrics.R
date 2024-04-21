#' @title Classify multi-state QTL
#'
#' @description
#' Takes the results from `callSignificance()` and from the assay `betas` to
#' categorize each QTL test using two classification strategies:
#'
#' Strategy 1 (qtl_type):
#' (1) global-shared, (2) global-diverging, (3) multi-state-shared,
#' (4) multi-state-diverging, or (5) unique.
#'
#' Strategy 2 (qtl_type_simple):
#' (1) global, (2) multi-state, or (3) unique.
#'
#' @param qtle QTLExperiment qtle
#' @param assay Name of assay containing QTL effect size estimate (e.g. betas)
#' @param assaySig Name of assay with TRUE/FALSE significance calls
#' @param globalBuffer Number of states that can be not-significant and the
#'                      QTL will still be called as global, for example, if
#'                      globalBuffer=1, then a QTL will be considered global if
#'                      if is significant in all or all but 1 state.
#' @param ... arguments passed to \code{runTestMetrics}
#'
#' @details If a test is significant in more than one sign across different
#' states, returns TRUE in rowData(qtle)$diverging
#'
#' @return The `QTLExperiment` object with the following columns added to the
#' rowData: nSignificant, effect_sd, qtl_type, qtl_type_simple
#'
#' @examples
#' m <- mockQTLE()
#' m <- callSignificance(m)
#' m <- runTestMetrics(m)
#'
#' @importFrom SummarizedExperiment assay rowData<-
#' @importFrom stats sd
#'
#' @name runTestMetrics
#' @rdname runTestMetrics
#' @export
#'
runTestMetrics <- function(qtle, assay="betas",
    assaySig = "significant",
    globalBuffer = 0, ...) {
    
    if ( !is(qtle, "QTLExperiment") )
        stop("qtle must be a QTLExperiment")

    if( ! assaySig %in% names(assays(qtle)) ) {
        stop("Run callSignificance() or specify assay name that contains logical
         significance calls.")
    }

    significance_mat <- assay(qtle, assaySig)
    beta_mat <- assay(qtle, assay)
    beta_mat[!significance_mat] <- NA
    sign_mat <- sign(beta_mat)
    n <- ncol(qtle)
    nGlobal <- ncol(qtle) - globalBuffer
    nSignificant <- rowSums(significance_mat)

    if(n > 2){
        type <- ifelse(nSignificant >= nGlobal, "global",
            ifelse(nSignificant == 0, "not_sig",
                ifelse(nSignificant == 1, "unique", "multistate")))

    } else if(n==2){
        state1 <- state_id(qtle)[1]
        state2 <- state_id(qtle)[2]
        type <- apply(significance_mat, 1,
            FUN = function(x) ifelse(sum(x) == 2, "both",
                ifelse(sum(x) == 0, "not_sig",
                    ifelse(x[state1]==TRUE,
                        state1, state2))))
        both <- type[type=="both"]
    }

    effect_sd <- apply(beta_mat, 1, FUN = function(x) ifelse(sum(!is.na(x)) >= 2,
        sd(na.omit(x)), NA))

    effect_dir <- apply(sign_mat, 1,
        FUN = function(x) ifelse(sum(is.na(x)) >= n-1, "",
            ifelse(length(na.omit(unique(x)))==1,
                "_shared", "_diverging")))

    rowData(qtle)$nSignificant <- nSignificant
    rowData(qtle)$effect_sd <- effect_sd
    rowData(qtle)$qtl_type <- paste0(type, effect_dir)
    rowData(qtle)$qtl_type_simple <- type

    return(qtle)
}

