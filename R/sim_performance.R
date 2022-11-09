#' @title Performance metrics for multistateQTL simulations 
#' 
#' @param qtle A `multistateQTL` object.
#' @param assay Name of the `multistateQTL` assay containing the significance calls 
#'
#' @examples
#' sim <- qtleSimulate()
#' sim <- callSignificance(sim, assay="lfsrs", thresh=0.1)
#' simPerformance(sim)
#'
#' @importFrom stats rnorm
#' @importFrom QTLExperiment mockQTLE
#' 
#' @name simPerformance
#' @rdname qtle_simulations
#' 
#' @export

simPerformance <- function(qtle, assay="significant") {
  
  state_ids <- state_id(qtle)
  if(!all(state_ids %in% names(rowData(qtle)))){
    stop("The qtle rowData is missing the TRUE/FALSE calls for some state(s)")
  }
  
  # Macro performance metrics
  simulated <- as.vector(as.matrix(rowData(qtle)[, state_ids]))
  called <- as.vector(assay(qtle, assay))
  
  cm = as.matrix(table(simulated = simulated, called = called))
  n = sum(cm) # number of instances
  diag = diag(cm) # number of correctly classified instances per class 
  rowsums = apply(cm, 1, sum) # number of instances per class
  colsums = apply(cm, 2, sum) # number of predictions per class

  accuracy = sum(diag) / n 
  precision = diag / colsums 
  recall = diag / rowsums 
  f1 = 2 * precision * recall / (precision + recall) 

  return(list(accuracy = accuracy, precision=precision[["TRUE"]], 
              recall=recall[["TRUE"]], f1=f1[["TRUE"]],
              cm = cm))
}