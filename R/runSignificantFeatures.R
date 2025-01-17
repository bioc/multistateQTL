#' Function to summarize QTL summary statistics by feature ID
#'
#' This function adds a summary of the features with significant QTL in each 
#' state to the metadata of the QTLExperiment object.
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
#' @return The `QTLExperiment` object with a summary of significant features in 
#' the metadata and a new column `nSignificantFeatures` in the colData.
#' 
#' @examples
#' 
#' qtle <- mockQTLE()
#' 
#' qtle <- callSignificance(qtle)
#' 
#' # There is an assay called 'significant'
#' assays(qtle)
#' 
#' # Obtain summary of significant features for each state
#' qtle <- runSignificantFeatures(qtle)
#' 
#' # There is a summary added to the metadata of the object
#' metadata(qtle)
#' 
#' @name runSignificantFeatures
#' @rdname runSignificantFeatures
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