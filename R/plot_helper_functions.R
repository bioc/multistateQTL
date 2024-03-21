#' Color palettes for ComplexHeatmaps and their annotations. If more than 8
#' continuous or 4 qualitative palettes are requested they will start repeating.
#'
#' @param color_by array with values to color by
#' @param i index of color scheme to select
#'
#' @importFrom circlize colorRamp2
#' @importFrom viridis viridis
#' @importFrom grDevices hcl.colors
#' @importFrom stats setNames na.omit
#' 
#' @noRd
.resolve_complexheatmap_colors <- function(color_by, i=1) {

    if (is.numeric(color_by)) {
        color_by <- na.omit(color_by)
        i <- i %% 8
        min <- min(color_by)
        max <- max(color_by)

        if(all(color_by <= 1) & all(color_by >= 0)){
            n <- 9
            palette <- c("YlGnBu", "YlOrRd", "RdPu", "GnBu",
                "Greys", "YlOrBr", "YlGn", "BuPu")[i]

        } else if(any(color_by < 0) & any(color_by > 0)){
            n <- 11
            min <- -1 * max(abs(min), abs(max))
            max <- max(abs(min), abs(max))
            palette <- c("RdBu", "PRGn", "BrBG", "PiYg",
                "PuOr", "RdGy", "Spectral", "RdYlBu")[i]

        } else {
            n <- 9
            palette <- c("YlGnBu", "YlOrRd", "RdPu", "GnBu",
                "Greys", "YlOrBr", "YlGn", "BuPu")[i]
        }
        palette
        colors <- colorRamp2(seq(max, min, length = n),
            hcl.colors(n, palette))
    }
    else {
        levels <- unique(color_by)
        nLevels <- length(levels)

        if (nLevels <= 10){
            paulTolPallets <- list(
                bright = c("#66CCEE", "#228833", "#CCBB44",
                           "#AA3377", "#4477AA", "#EE6677",
                           "#EE7733", "#997700", "#BBCC33",
                           "#BBBBBB"),
                muted = c("#332288", "#CC6677", "#999933",
                          "#44AA99", "#882255", "#DDCC77",
                          "#117733", "#88CCEE", "#AA4499",
                          "#BBBBBB"),
                vibrant = c("#009988", "#EE7733", "#33BBEE",
                            "#EE3377", "#CC3311", "#0077BB",
                            "#DDAA33", "#BBCC33", "#AA4499",
                            "#BBBBBB"),
                light = c("#77AADD", "#BBCC33", "#EE8866",
                         "#44BB99", "#FFAABB", "#99DDFF",
                         "#BBBBBB", "#AAAA00", "#EEDD88",
                         "#AA44DD"),
                medium = c("#EECC66", "#EE99AA", "#6699CC",
                          "#CCDDAA", "#997700", "#994455",
                          "#004488", "#225522", "#DDDDDD",
                          "#555555"))
            colors <- setNames(paulTolPallets[[(i %% 5) + 1]][1:nLevels], levels)
        } else {
            colors <- setNames(viridis(nLevels), levels)
        }
    }

    return(colors)
}


#' Color palettes for ggplot and other (non-ComplexHeatmap) plotting functions.
#'
#' @param plot_out plot object
#' @param color_by array with values to color by
#' @param color_by_name name of color_by array
#' @param fill Logical indicating if the palette is for fill or color (e.g.
#'        scale_fill_brewer vs scale_color_brewer).
#'
#' @importFrom ggplot2 scale_fill_brewer scale_color_brewer scale_fill_manual scale_color_manual
#' @importFrom viridis scale_fill_viridis scale_color_viridis
#'
#' @noRd
.resolve_plot_colors <- function(plot_out, color_by, color_by_name,
    fill = FALSE) {
    
    if ( fill ) {
        if ( is.numeric(color_by) ) {
            if (all(color_by <= 1) & all(color_by >= 0)){
                plot_out <- plot_out +
                    scale_fill_viridis(name = color_by_name, option="mako")
            } else if (all(color_by <= 1) & all(color_by >= -1)){
                plot_out <- plot_out +
                    scale_fill_brewer(name = color_by_name, palette = "RdBu")
            } else{
                plot_out <- plot_out + scale_fill_viridis(name = color_by_name)
            }

        } else {
            nlevs_color_by <- nlevels(as.factor(color_by))
            if (nlevs_color_by <= 20) {
                plot_out <- plot_out + scale_fill_manual(
                    values = c("#66CCEE", "#228833", "#CCBB44", "#AA3377", "#4477AA",
                               "#EE6677", "#EE7733", "#997700", "#BBCC33", "#BBBBBB",
                               "#77AADD", "#BBCC33", "#EE8866", "#44BB99", "#FFAABB",
                               "#99DDFF", "#AAAA00", "#EEDD88", "#AA44DD", "#DDDDDD"),
                    name = color_by_name)
            } else {
                plot_out <- plot_out +
                    scale_fill_viridis(name = color_by_name, discrete = TRUE)
            }
        }
    } else {
        if ( is.numeric(color_by) ) {
            if (all(color_by <= 1) & all(color_by >= 0)){
                plot_out <- plot_out +
                    scale_color_viridis(name = color_by_name, option="mako")
            } else if (all(color_by <= 1) & all(color_by >= -1)){
                plot_out <- plot_out +
                    scale_color_brewer(name = color_by_name, palette = "RdBu")
            } else{
                plot_out <- plot_out + scale_color_viridis(name = color_by_name)
            }
        } else {
            nlevs_color_by <- nlevels(as.factor(color_by))
            if (nlevs_color_by <= 10) {
                plot_out <- plot_out + scale_color_manual(
                    values = c("#66CCEE", "#228833", "#CCBB44", "#AA3377", "#4477AA",
                               "#EE6677", "#EE7733", "#997700", "#BBCC33", "#BBBBBB",
                               "#77AADD", "#BBCC33", "#EE8866", "#44BB99", "#FFAABB",
                               "#99DDFF", "#AAAA00", "#EEDD88", "#AA44DD", "#DDDDDD"),
                    name = color_by_name)
            } else {
                plot_out <- plot_out + scale_color_viridis(name = color_by_name,
                                                           discrete = TRUE)
            }
        }
    }
    plot_out
}


#' Get ComplexHeatmap annotations for rows and columns
#'
#' @param anns data.frame with annotations as columns
#' @param FUN function to use to generate ComplexHeatmap annotation object.
#'        HeatmapAnnotation for columns and rowAnnotation for rows.
#'
#' @importFrom ComplexHeatmap HeatmapAnnotation rowAnnotation
#'
#' @noRd
.resolve_annotations <- function(anns, FUN=rowAnnotation) {

    annColors <- list()
    for(i in seq(1, ncol(anns))) {
        aName <- colnames(anns)[i]
        annColors[[aName]] <- .resolve_complexheatmap_colors(anns[, aName], i=i+1)
    }

    annotationObject <- do.call(FUN, list(df=anns, col = annColors))
    return (annotationObject)
}


