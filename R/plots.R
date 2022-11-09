#' @title Heatmap of pairwise QTL sharing with state-level annotations
#'
#' @description 
#' Methods to plot at heatmap of the pairwise sharing of QTL as calculated
#' by `runPairwiseSharing`.
#' 
#' @param object A \code{QTLExperiment} object
#' @param slot Name of slot in metadata list with Pairwise Sharing matrix. 
#' @param annotate_rows character or array of characters specifying the 
#'        column(s) in colData to be plotted as row annotations.
#' @param annotate_cols character or array of characters specifying the 
#'        column(s) in colData to be plotted as column annotations.
#' @param cell_annotate Logical to annoate cells with their values.
#' @param col_range Optional range for the color legend
#' @param name character specifying the column in colData to use to label rows
#'        and columns. Default is colnames(qtle).
#' @param dist.method Distance method used for hierarchical clustering. Valid 
#'        values are the supported methods in dist() function.
#' @param size numeric scalar giving default font size for plotting theme.
#' @param ... further arguments passed to \code{\link[Rtsne]{Rtsne}}
#'
#' 
#' @return Returns a \code{ComplexHeatmap} object.
#' @name plotPairwiseSharing
#' @rdname plotting
#' 
#' @importFrom SummarizedExperiment colData
#' @importFrom grid gpar
#' @importFrom S4Vectors metadata
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom grid grid.text
#' @importFrom methods is
#' 
#' @export
#'
plotPairwiseSharing <- function(object, slot = "pairwiseSharing",
                                annotate_rows = NULL, annotate_cols = NULL,
                                cell_annotate = FALSE, col_range = NULL,
                                name = "colnames",
                                dist.method = "pearson",
                                size = 8, ...) {

  if ( !is(object, "QTLExperiment") )
    stop("Object must be a QTLExperiment")

  if (!(slot %in% names(metadata(object)))) {
    stop("Run runPairwiseSharing or provide metadata slot with data to plot.")
  }

  mat <- metadata(object)[[slot]]
  
  if (all(mat == 0)) { stop("No pairwise sharing detected...")}
  if (any(is.na(mat))) { stop("NaNs in pairwise sharing results...")}
  
  if(name != "colnames"){
    row.names(mat) <- colData(object)[, name]
    colnames(mat) <- colData(object)[, name]
  }
  
  if(!is.null(col_range)){
    mat.cols <- .resolve_complexheatmap_colors(col_range)
  } else{
    mat.cols <- .resolve_complexheatmap_colors(array(mat))
  }

  if(!is.null(annotate_rows)){
    anns <- as.data.frame(colData(object)[, annotate_rows])
    colnames(anns) <- annotate_rows
    ann_rows <- .resolve_annotations(anns)

  }else {ann_rows <- NULL}

  if(!is.null(annotate_cols)){
    anns <- as.data.frame(colData(object)[, annotate_cols])
    colnames(anns) <- annotate_cols
    ann_cols <- .resolve_annotations(anns, FUN=HeatmapAnnotation)

  }else {ann_cols <- NULL}
  
  if(cell_annotate) {
    cell_annotate <- function(j, i, x, y, width, height, fill) {
      grid.text(sprintf("%.2f", mat[i, j]), x, y, 
                gp = gpar(fontsize = size-2))
      } 
  } else{
    cell_annotate <- NULL
  }
  
  Heatmap(mat, name = slot, col = mat.cols,
          clustering_distance_rows = dist.method,
          clustering_distance_columns = dist.method,
          top_annotation = ann_cols, right_annotation = ann_rows,
          cell_fun = cell_annotate,
          column_names_gp = gpar(fontsize = size),
          row_names_gp = gpar(fontsize = size))

}


#' @title Upset plot of QTL sharing between states with state-level annotations
#'
#' @description 
#' Convenient method to plot an UpSet plot showing the number of QTL that are 
#' significant in sets of states. 
#' 
#' @param object an \code{QTLExperiment} object
#' @param assay Name of assay to use to assess significance.
#' @param min_shared minimum number of shared QTL for set to be included. 
#' @param min_degree minimum degree of sharing for set to be included.
#' @param max_degree maximum degree of sharing for set to be included.
#' @param annotate_by character or array of characters specifying the 
#'        column(s) in colData to be plotted as annotations.
#' @param name character specifying the column in colData to use to label rows
#'        and columns. Default is colnames(qtle).
#' @param comb_order characters specifying how sets should be ordered. Options
#'        include the set size (set_size), combination size (comb_size), degree (deg). 
#' @param set_order Array specifying order of states. 
#' @param ... further arguments passed to \code{\link[Rtsne]{Rtsne}}
#'
#' @return Returns a \code{ComplexHeatmap} object.
#' @name plotUpSet
#' @rdname plotting
#'
#' @importFrom ComplexHeatmap make_comb_mat UpSet comb_size comb_degree set_size
#' 
plotUpSet <- function(object, 
                      assay = "significant",
                      name = "colnames",
                      min_shared=10,
                      min_degree=2,
                      max_degree=NULL,
                      annotate_by = NULL,
                      comb_order="comb_size",
                      set_order = order(ss)){

  if( ! assay %in% names(assays(object)) ) {
    stop("Run callSignificance() first...")
  }
  
  sigMatrix <- as.matrix(assay(object, assay))
  sigMatrix <- sigMatrix[(rowSums(sigMatrix[ , ]) != 0), ]
  if(name != "colnames"){
    colnames(sigMatrix) <- colData(object)[, name]
  }

  combMatrix <- make_comb_mat(sigMatrix)
  combMatrix <- combMatrix[comb_size(combMatrix) >= min_shared]
  combMatrix <- combMatrix[comb_degree(combMatrix) >= min_degree]

  if(!is.null(max_degree)){
    combMatrix <- combMatrix[comb_degree(combMatrix) <= max_degree]
  }

  if(!is.null(annotate_by)){
    if(!is(annotate_by, "HeatmapAnnotation")){
      anns <- as.data.frame(colData(object)[, annotate_by])
      colnames(anns) <- annotate_by
      rowAnns <- .resolve_annotations(anns, FUN=rowAnnotation)
    } else{ rowAnns <- annotate_by}
  }else {rowAnns <- NULL}

  ss <- -set_size(combMatrix)

  if (!is.null(comb_order) & comb_order != "set_size") {
    cs <- comb_size(combMatrix)
    cd <- comb_degree(combMatrix)
    if(comb_order == "comb_size") {
      comb_order <- order(-cs)
    } else if (comb_order == "comb_degree") {
      comb_order <- order(cd, -cs)
    } else { comb_order <- order(cd) }
  }


  UpSet(combMatrix, set_order = set_order,
        comb_order = comb_order) + rowAnns

}


#' @title Heatmap of QTL across states
#'
#' @description 
#' Convenience method for plotting values from any assay specified by fill_by 
#' across states and tests. 
#' 
#' 
#' @param object an \code{QTLExperiment} object.
#' @param fill_by name of assay to use for main heatmap object.
#' @param FUN Function to be applied to fill_by assay before plotting (e.g. 
#'            identity, abs, log10).
#' @param min_sig minimum number of significant states for QTL to be included.
#' @param annotate_states character or array of characters specifying the 
#'        column(s) in colData to be plotted to annotate states.
#' @param annotate_tests character or array of characters specifying the 
#'        column(s) in rowData to be plotted to annotate QTL tests.
#' @param state_id colData column to use to label states.
#' @param column_order Null for clustering or array to overwrite state order
#' @param row_order Null for clustering or array to overwrite test order
#' @param show_row_names Logical to plot row (i.e. test) names.
#' @param row_km Set k for k-means clustering of tests. 
#' @param col_km Set k for k-means clustering of states
#'
#' @return Returns a \code{ComplexHeatmap} object.
#' @name plotQTLClusters
#' @rdname plotting
#' 
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom methods is
#' 
#' @export
#'
plotQTLClusters <- function(object,
                            fill_by="betas",
                            FUN=identity,
                            min_sig = 1,
                            annotate_states=NULL,
                            annotate_tests=NULL,
                            show_row_names=FALSE,
                            state_id = "state_id",
                            column_order = NULL,
                            row_order = NULL,
                            row_km = 0,
                            col_km = 0) {

  if (!is(object, "QTLExperiment")) {
    stop("Object must be a QTLExperiment")
  }

  if ( !(fill_by %in% names(assays(object)))) {
    stop("Specify which assay to fill and cluster by")
  }

  mat <- FUN(as.matrix(assays(object)[[fill_by]]))
  if(fill_by %in% c("p", "pvalues", "lfsrs")){
    mat <- -log10(mat)
  }
  
  if(state_id != "state_id"){
    colnames(mat) <- colData(object)[, state_id]
  }
  
  mat.cols <- .resolve_complexheatmap_colors(array(mat))

  if(!is.null(annotate_states)){
    if(!is(annotate_states, "HeatmapAnnotation")){
      anns <- as.data.frame(colData(object)[, annotate_states])
      colnames(anns) <- annotate_states
      ann_state <- .resolve_annotations(anns, FUN=HeatmapAnnotation)
    } else{ 
      ann_state <- annotate_states }
  }else {ann_state <- NULL}

  if(!is.null(annotate_tests)){
    if(!is(annotate_states, "HeatmapAnnotation")){
      anns <- as.data.frame(rowData(object)[, annotate_tests])
      colnames(anns) <- annotate_tests
      ann_tests <- .resolve_annotations(anns)
    } else{ann_tests <- annotate_tests}
  }else {ann_tests <- NULL}

  if(is.null(row_order)){row_order <- row.names(mat)}
  if(is.null(column_order)){column_order <- names(mat)}
  
  Heatmap(mat, name = fill_by, col = mat.cols, 
          show_row_names = show_row_names,
          top_annotation = ann_state, 
          row_km = row_km, column_km = col_km,
          right_annotation = ann_tests)
}



#' @title Compare QTL between two states
#'
#' @description 
#' Convenience method for comparing the assay value, specified by assay, 
#' between two states. 
#' 
#' @param object an \code{QTLExperiment} object.
#' @param x Name of state for x-axis
#' @param y Name of state for y-axis
#' @param assay name of assay to plot.
#' @param significance_assay name of assay with TRUE/FALSE significance calls.
#' @param FUN Function to be applied to fill_by assay before plotting (e.g. 
#'            identity, abs, log10).
#' @param alpha Transparency.
#' @param col_both Color for tests significant in both states.
#' @param col_diverging Color for tests significant in both states, with 
#'        diverging effect sizes.
#' @param col_neither Color for null tests.
#' @param col_x Color for tests significant in the x-axis state only.
#' @param col_y Color for tests significant in the y-axis state only.
#' 
#' @return Returns a list containing the counts for each color category 
#'         and the plot object.
#'
#' @name plotCompareStates
#' @rdname plotting

#' @importFrom ggplot2 geom_abline geom_point aes .data xlim ylim
#' 
#' @export

plotCompareStates <- function(object, x, y,
                              assay="betas", FUN=identity,
                              significance_assay = "significant",
                              alpha=0.2,
                              col_both="#4477AA",
                              col_diverging="#EE6677",
                              col_neither="gray50",
                              col_x="#CCBB44",
                              col_y="#AA3377") {
  
  if( ! significance_assay %in% names(assays(object)) ) {
    stop("Run callSignificance() or specify assay name that contains logical 
         significance calls.")
  }
  
  object <- object[, c(x, y)]
  
  object <- runTestMetrics(object)
  
  to_plot <- FUN(as.data.frame(assay(object, assay)))
  to_plot[, "qtl_type"] <- rowData(object)[["qtl_type"]]
  to_plot$name <- rownames(to_plot)
  diverging <- row.names(to_plot[grepl("_diverging", to_plot[["qtl_type"]]), ])
  
  cols <- list(both_shared=col_both, 
               both_diverging = col_diverging, 
               not_sig=col_neither)
  cols[[x]] <- col_x
  cols[[y]] <- col_y

  min_val <- min(c(to_plot[[x]], to_plot[[y]]))
  max_val <- max(c(to_plot[[x]], to_plot[[y]]))

  plot <- ggplot(to_plot, aes(x = .data[[x]], y = .data[[y]], color=qtl_type)) + 
    geom_point(size = 1, alpha=alpha) +
    scale_color_manual(values=cols, na.value = "#000000") + 
    xlim(c(min_val, max_val)) + ylim(c(min_val, max_val)) +
    geom_abline(linetype = "dashed", slope = 1) 
  
  counts <- table(to_plot[, "qtl_type"])
  
  return(list(plot=plot, counts=counts))
 
}


#' @title Distribution of estimated simulation parameters
#' 
#' @param params Estimated simulation parameters from qtleEstimate.
#' @param n Number of random values to sample for plots.
#' 
#' @return A ggplot2 object
#' 
#' @name plotSimulationParams
#' @rdname plotting
#' 
#' @importFrom ggplot2 ggplot geom_density facet_grid theme_classic aes
#' @importFrom stats rgamma
#' @export
#' 
plotSimulationParams <- function(params, n=1e5){
  
  demo_data <- as.data.frame(list(statistic=c(rep("beta", n*2), rep("CV", n*2)),
                                  qtl_type=rep(c(rep("significant", n), 
                                             rep("not significant", n)), 2),
                                  value=c(rgamma(n, params$betas.sig.shape, 
                                                  params$betas.sig.rate),
                                           rgamma(n, params$betas.null.shape, 
                                                  params$betas.null.rate),
                                           rgamma(n, params$cv.sig.shape, 
                                                  params$cv.sig.rate),
                                           rgamma(n, params$cv.null.shape, 
                                                  params$cv.null.rate))))
  
  ggplot(demo_data, aes(x=value, fill=qtl_type, color=qtl_type)) +
    geom_density(alpha=0.2) + 
    facet_grid(.~statistic, scales="free") + 
    theme_classic()
  
}

