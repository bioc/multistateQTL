#' @title Heatmap of pairwise QTL sharing with state-level annotations
#'
#' @description
#' Methods to plot a heatmap of the pairwise sharing of QTL as calculated
#' by `runPairwiseSharing`.
#'
#' @param object A \code{QTLExperiment} object
#' @param slot Name of slot in metadata list with Pairwise Sharing matrix.
#' @param annotateRowsBy character or array of characters specifying the
#'        column(s) in colData to be plotted as row annotations.
#' @param annotateColsBy character or array of characters specifying the
#'        column(s) in colData to be plotted as column annotations.
#' @param annotateCells Logical to annotate cells with their values.
#' @param colourRange Optional range for the color legend
#' @param name character specifying the column in colData to use to label rows
#'        and columns. Default is colnames(qtle).
#' @param distMethod Distance method used for hierarchical clustering. Valid
#'        values are the supported methods in dist() function.
#' @param size numeric scalar giving default font size for plotting theme.
#' @param ... further arguments passed to \code{\link[Rtsne]{Rtsne}}
#'
#' @return Returns a \code{ComplexHeatmap} object.
#' 
#' @examples 
#' sim <- qtleSimulate(
#'     nStates=10, nFeatures=100, nTests=1000,
#'     global=0.2, multi=0.4, unique=0.2, k=2)
#' sim <- callSignificance(sim, mode="simple", assay="lfsrs", 
#'     thresh=0.0001, secondThresh=0.0002)
#' sim_sig <- getSignificant(sim)
#' sim_top <- getTopHits(sim_sig, assay="lfsrs", mode="state")
#' sim_top <- runPairwiseSharing(sim_top)
#' 
#' plotPairwiseSharing(sim_top)
#' 
#' # Plot with complex column annotations
#' plotPairwiseSharing(sim_top, annotateColsBy = c("nSignificant", "multistateGroup"))
#' 
#' 
#' @name plotPairwiseSharing
#' @rdname plotPairwiseSharing
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
    annotateRowsBy = NULL, annotateColsBy = NULL,
    annotateCells = FALSE, colourRange = NULL,
    name = "colnames",
    distMethod = "pearson",
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

    if(!is.null(colourRange)){
        mat.cols <- .resolve_complexheatmap_colors(colourRange)
    } else{
        mat.cols <- .resolve_complexheatmap_colors(array(mat))
    }

    if(!is.null(annotateRowsBy)){
        anns <- as.data.frame(colData(object)[, annotateRowsBy])
        colnames(anns) <- annotateRowsBy
        ann_rows <- .resolve_annotations(anns)

    }else {ann_rows <- NULL}

    if(!is.null(annotateColsBy)){
        anns <- as.data.frame(colData(object)[, annotateColsBy])
        colnames(anns) <- annotateColsBy
        ann_cols <- .resolve_annotations(anns, FUN=HeatmapAnnotation)

    }else {ann_cols <- NULL}

    if(annotateCells) {
        annotateCells <- function(j, i, x, y, width, height, fill) {
            grid.text(sprintf("%.2f", mat[i, j]), x, y,
                gp = gpar(fontsize = size-2))
        }
    } else{
        annotateCells <- NULL
    }

    Heatmap(mat, name = slot, col = mat.cols,
        clustering_distance_rows = distMethod,
        clustering_distance_columns = distMethod,
        top_annotation = ann_cols, right_annotation = ann_rows,
        cell_fun = annotateCells,
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
#' @param minShared minimum number of shared QTL for set to be included.
#' @param minDegree minimum degree of sharing for set to be included.
#' @param maxDegree maximum degree of sharing for set to be included.
#' @param annotateColsBy character or array of characters specifying the
#'        column(s) in colData to be plotted as annotations.
#' @param name character specifying the column in colData to use to label rows
#'        and columns. Default is colnames(qtle).
#' @param comb_order characters specifying how sets should be ordered. Options
#'        include the set size (set_size), combination size (comb_size), degree (deg).
#' @param set_order Array specifying order of states.
#' @param ... further arguments passed to \code{\link[Rtsne]{Rtsne}}
#'
#' @return Returns a \code{ComplexHeatmap} object.
#' 
#' @examples 
#' sim <- qtleSimulate(
#'     nStates=10, nFeatures=100, nTests=1000,
#'     global=0.2, multi=0.4, unique=0.2, k=2)
#' sim <- callSignificance(sim, mode="simple", assay="lfsrs", 
#'     thresh=0.0001, secondThresh=0.0002)
#' sim_sig <- getSignificant(sim)
#' sim_top <- getTopHits(sim_sig, assay="lfsrs", mode="state")
#' sim_top <- runPairwiseSharing(sim_top)
#'
#' plotUpSet(sim_top)
#' 
#' # Upset plot with complex row annotations
#' plotUpSet(sim_top, annotateColsBy = c("nSignificant", "multistateGroup"))

#' @name plotUpSet
#' @rdname plotUpSet
#'
#' @importFrom ComplexHeatmap make_comb_mat UpSet comb_size comb_degree set_size
#'
#' @export
#'
plotUpSet <- function(object,
    assay = "significant",
    name = "colnames",
    minShared=10,
    minDegree=2,
    maxDegree=NULL,
    annotateColsBy = NULL,
    comb_order="comb_size",
    set_order = order(ss), 
    ...){
    
    if ( !is(object, "QTLExperiment") )
        stop("Object must be a QTLExperiment")

    if( ! assay %in% names(assays(object)) ) {
        stop("Run callSignificance() first...")
    }

    sigMatrix <- as.matrix(assay(object, assay))
    sigMatrix <- sigMatrix[(rowSums(sigMatrix[ , ]) != 0), ]
    if(name != "colnames"){
        colnames(sigMatrix) <- colData(object)[, name]
    }

    combMatrix <- make_comb_mat(sigMatrix)
    combMatrix <- combMatrix[comb_size(combMatrix) >= minShared]
    combMatrix <- combMatrix[comb_degree(combMatrix) >= minDegree]

    if(!is.null(maxDegree)){
        combMatrix <- combMatrix[comb_degree(combMatrix) <= maxDegree]
    }

    if(!is.null(annotateColsBy)){
        if(!is(annotateColsBy, "HeatmapAnnotation")){
            anns <- as.data.frame(colData(object)[, annotateColsBy])
            colnames(anns) <- annotateColsBy
            rowAnns <- .resolve_annotations(anns, FUN=rowAnnotation)
        } else{ rowAnns <- annotateColsBy}
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
#' Convenience method for plotting values from any assay specified by fillBy
#' across states and tests.
#'
#'
#' @param object an \code{QTLExperiment} object.
#' @param fillBy name of assay to use for main heatmap object.
#' @param FUN Function to be applied to fillBy assay before plotting (e.g.
#'            identity, abs, log10).
#' @param minSig minimum number of significant states for QTL to be included.
#' @param annotateColsBy character or array of characters specifying the
#'        column(s) in colData to be plotted to annotate states.
#' @param annotateRowsBy character or array of characters specifying the
#'        column(s) in rowData to be plotted to annotate QTL tests.
#' @param state_id colData column to use to label states.
#' @param columnOrder Null for clustering or array to overwrite state order
#' @param rowOrder Null for clustering or array to overwrite test order
#' @param show_row_names Logical to plot row (i.e. test) names.
#' @param row_km Set k for k-means clustering of tests.
#' @param column_km Set k for k-means clustering of states
#'
#' @return Returns a \code{ComplexHeatmap} object.
#' 
#' @examples
#' sim <- qtleSimulate(
#'     nStates=10, nFeatures=100, nTests=1000,
#'     global=0.2, multi=0.4, unique=0.2, k=2)
#' sim <- callSignificance(sim, mode="simple", assay="lfsrs", 
#'     thresh=0.0001, secondThresh=0.0002)
#'     
#' sim_sig <- getSignificant(sim)
#' sim_top <- getTopHits(sim_sig, assay="lfsrs", mode="state")
#' sim_top <- runTestMetrics(sim_top)
#' sim_top <- runPairwiseSharing(sim_top)
#' sim_top_ms <- subset(sim_top, qtl_type_simple == "multistate")
#' 
#' plotQTLClusters(sim_top)
#' 
#' # Plot with annotations for multi state group
#' plotQTLClusters(sim_top_ms, annotateColsBy = c("multistateGroup"),
#'     annotateRowsBy = c("qtl_type", "mean_beta", "QTL"))
#' 
#' @name plotQTLClusters
#' @rdname plotQTLClusters
#'
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom methods is
#'
#' @export
#'
plotQTLClusters <- function(object,
    fillBy="betas",
    FUN=identity,
    minSig = 1,
    annotateColsBy=NULL,
    annotateRowsBy=NULL,
    show_row_names=FALSE,
    state_id = "state_id",
    columnOrder = NULL,
    rowOrder = NULL,
    row_km = 0,
    column_km = 0) {

    if (!is(object, "QTLExperiment")) {
        stop("Object must be a QTLExperiment")
    }

    if ( !(fillBy %in% names(assays(object)))) {
        stop("Specify which assay to fill and cluster by")
    }

    mat <- FUN(as.matrix(assays(object)[[fillBy]]))
    if(fillBy %in% c("p", "pvalues", "lfsrs")){
        mat <- -log10(mat)
    }

    if(state_id != "state_id"){
        colnames(mat) <- colData(object)[, state_id]
    }

    mat.cols <- .resolve_complexheatmap_colors(array(mat))

    if(!is.null(annotateColsBy)){
        if(!is(annotateColsBy, "HeatmapAnnotation")){
            anns <- as.data.frame(colData(object)[, annotateColsBy])
            colnames(anns) <- annotateColsBy
            ann_state <- .resolve_annotations(anns, FUN=HeatmapAnnotation)
        } else{
            ann_state <- annotateColsBy }
    }else {ann_state <- NULL}

    if(!is.null(annotateRowsBy)){
        if(!is(annotateColsBy, "HeatmapAnnotation")){
            anns <- as.data.frame(rowData(object)[, annotateRowsBy])
            colnames(anns) <- annotateRowsBy
            ann_tests <- .resolve_annotations(anns)
        } else{ann_tests <- annotateRowsBy}
    }else {ann_tests <- NULL}

    if(is.null(rowOrder)){rowOrder <- row.names(mat)}
    if(is.null(columnOrder)){columnOrder <- names(mat)}

    Heatmap(mat, name = fillBy, col = mat.cols,
        show_row_names = show_row_names,
        top_annotation = ann_state,
        row_km = row_km, column_km = column_km,
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
#' @param assaySig name of assay with TRUE/FALSE significance calls.
#' @param FUN Function to be applied to fillBy assay before plotting (e.g.
#'            identity, abs, log10).
#' @param alpha Transparency.
#' @param colBoth Color for tests significant in both states.
#' @param colDiverging Color for tests significant in both states, with
#'        diverging effect sizes.
#' @param colNeither Color for null tests.
#' @param colX Color for tests significant in the x-axis state only.
#' @param colY Color for tests significant in the y-axis state only.
#'
#' @return Returns a list containing the counts for each color category
#'         and the plot object.
#'         
#' @examples 
#' sim <- qtleSimulate(
#'     nStates=10, nFeatures=100, nTests=1000,
#'     global=0.2, multi=0.4, unique=0.2, k=2)
#' sim <- callSignificance(sim, mode="simple", assay="lfsrs", 
#'     thresh=0.0001, secondThresh=0.0002)
#' sim_sig <- getSignificant(sim)
#' sim_top <- getTopHits(sim_sig, assay="lfsrs", mode="state")
#' sim_top <- runPairwiseSharing(sim_top)
#' sim_top <- runTestMetrics(sim_top)
#' plotCompareStates(sim_top, x="S01", y="S02")
#'
#' @name plotCompareStates
#' @rdname plotCompareStates

#' @importFrom ggplot2 geom_abline geom_point aes .data xlim ylim
#'
#' @export

plotCompareStates <- function(object, x, y,
    assay="betas", FUN=identity,
    assaySig = "significant",
    alpha=0.2,
    colBoth="#4477AA",
    colDiverging="#EE6677",
    colNeither="gray50",
    colX="#CCBB44",
    colY="#AA3377") {

    if ( !is(object, "QTLExperiment") )
        stop("Object must be a QTLExperiment")
    
    if( ! assaySig %in% names(assays(object)) ) {
        stop("Run callSignificance() or specify assay name that contains logical
         significance calls.")
    }

    object <- object[, c(x, y)]

    object <- runTestMetrics(object)

    to_plot <- FUN(as.data.frame(assay(object, assay)))
    to_plot[, "qtl_type"] <- rowData(object)[["qtl_type"]]
    to_plot$name <- rownames(to_plot)
    diverging <- row.names(to_plot[grepl("_diverging", to_plot[["qtl_type"]]), ])

    cols <- list(both_shared=colBoth,
        both_diverging = colDiverging,
        not_sig=colNeither)
    cols[[x]] <- colX
    cols[[y]] <- colY

    min_val <- min(c(to_plot[[x]], to_plot[[y]]))
    max_val <- max(c(to_plot[[x]], to_plot[[y]]))

    plot <- ggplot(to_plot, aes(x = .data[[x]], y = .data[[y]], color=qtl_type)) +
        geom_point(size = 1, alpha=alpha) +
        scale_color_manual(values=cols, na.value = "#000000") +
        xlim(c(min_val, max_val)) + ylim(c(min_val, max_val)) +
        geom_abline(linetype = "dashed", slope = 1) +
        theme_classic()

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
#' @examples
#' qtle <- mockQTLE() 
#' params <- qtleEstimate(qtle, threshSig = 0.05, threshNull = 0.5)
#' plotSimulationParams(params=params)
#' 
#'
#' @name plotSimulationParams
#' @rdname plotSimulationParams
#'
#' @importFrom ggplot2 ggplot geom_density facet_grid theme_classic aes
#' @importFrom stats rgamma
#' @export
#'
plotSimulationParams <- function(params, n=1e5){
    
    if ( !is(params, "list") )
        stop("params must be a list with names as described in ?qtleParams")

    demo_data <- as.data.frame(
        list(
            statistic=c(rep("beta", n*2), rep("CV", n*2)),
            qtl_type=rep(c(rep("significant", n), rep("not significant", n)), 2),
            value=c(
                rgamma(n, params$betas.sig.shape,
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

