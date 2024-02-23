#' Function to summarize QTL summary statistics by feature ID
#'
#' This function pulls the top QTL for each feature (e.g., gene) and
#' outputs a **new** \linkS4class{QTLExperiment}. Note that the
#' output from this function should be assigned to a new variable as it will
#' overwrite if assigned to the same variable.
#'
#' @param object A \linkS4class{QTLExperiment} object with multiple
#'        QTL tests (i.e., rows) for at least one feature.
#' @param assay Assay containing T/F significance calls for each test.
#'
#' @importFrom collapse fmutate
#' @importFrom dplyr filter
#' @importFrom tidyr pivot_longer
#' @importFrom data.table setDT
#' @importFrom SummarizedExperiment assay rowData colData<-
#' @importFrom S4Vectors metadata metadata<-
#'
#' @export
#'
runSignificantFeatures <- function(object, assay = "significant") {

    eFeatures <- as.data.frame(assay(object, assay)) %>%
        fmutate(.feature_id = rowData(object)[[.feature_id]]) %>%
        pivot_longer(-.feature_id) %>% dplyr::filter(value == TRUE) %>% unique()

    eFeatures_list <- split(eFeatures$.feature_id, eFeatures$name)

    n_eFeatures <- lapply(eFeatures_list, FUN=function(x) {length(x)})

    colData(object)$nSignificantFeatures <- unlist(n_eFeatures[colnames(object)])
    metadata(object)$eFeatures <- eFeatures_list

    return(object)
}


