#' HSV encoded plots
#'
#' Using an HSV encoding of feature values, this functions can be used to color
#' code expression profiles of multiple features and visualize spatially.
#'
#' Using RGB encoding, we can show up to 3 features at the same time in the
#' "red", "green" and "blue" color channels. Whenever two or three features overlap,
#' the color will be a mix of the three channels, e.g. 50% green and 50% red will give a yellow color.
#' This strategy is very effective when looking at features values with significant
#' overlap but is limited to show maximum three features.
#'
#' If we want to show more than three features in the same plot, this
#' function provides a strategy to do this as long as the overlap between features
#' is relatively low. First, a color is assigned to each of N features by cutting
#' the hue (H) into N values with an even interval. The feature values (e.g. gene expression)
#' are then scaled to a 0-1 range which is encoded in the Value channel (V).
#' For each spot, the color with the highest V is selected meaning that only the
#' feature with the highest value will be shown in the plot. This strategy works well
#' for features with no or very little overlap but gets cluttered when to many
#' features are included.
#'
#' @param object Seurat object
#' @param features
#' \itemize{
#'     \item An \code{Assay} feature (e.g. a gene name - "MS4A1")
#'     \item A column name from meta.data (e.g. mitochondrial percentage - "percent.mito")
#' }
#' @param indices Numeric vector specifying sample indices to include in plot. Default is to show all samples.
#' @param spots Vector of spots to plot (default is all spots)
#' @param min.cutoff,max.cutoff Vector of minimum and maximum cutoff values for each feature,
#'  may specify quantile in the form of 'q##' where '##' is the quantile (eg, 'q1', 'q10')
#' @param slot Which slot to pull expression data from?
#' @param blend Scale and blend expression values to visualize coexpression of two features (This options will override other coloring parameters)
#' @param pt.size Adjust point size for plotting
#' @param shape.by If NULL, all points are circles (default). You can specify any
#' cell attribute (that can be pulled with FetchData) allowing for both
#' different colors and different shapes on cells
#' @param delim delimiter passed to \code{\link{GetCoords}} if adjusted ST coordinates are missing in the meta data
#' @param return.plot.list should the plots be returned as a list? By default, the plots are arranged into a grid
#' @param grid.ncol Number of columns for display when combining plots
#' @param verbose Print messages
#' @param ... Extra parameters passed on to \code{\link{STPlot}}
#'
#' @inheritParams STPlot
#' @importFrom cowplot plot_grid
#' @importFrom scales rescale
#' @importFrom ggplot2 ggplot theme theme_void
#' @importFrom zeallot %<-%
#'
#' @return A ggplot object
#' @export

HSVFeaturePlot <- function (
  object,
  features,
  indices = NULL,
  spots = NULL,
  min.cutoff = NA,
  max.cutoff = NA,
  slot = "data",
  pt.size = 1,
  shape.by = NULL,
  cols = NULL,
  dark.theme = TRUE,
  ncol = NULL,
  delim = NULL,
  return.plot.list = FALSE,
  grid.ncol = NULL,
  verbose = FALSE,
  theme = theme_void(),
  ...
) {

  # Collect data
  spots <- spots %||% colnames(x = object)
  data <- FetchData(object = object, vars = c(features), cells = spots, slot = slot)
  data.type <- unique(sapply(data, class))

  # Stop if feature classes are not numeric/integer
  if (!all(data.type %in% c("numeric", "integer"))) {
    stop("Classes features of class 'integer' or 'numeric' are allowed ... ")
  }

  # Check that group.by variable is present in meta.data slot, otherwise assume that there's only one sample present in Seurat object
  if ("sample" %in% colnames(object[[]])) {
    group.by <- "sample"
    data[,  group.by] <- object[[group.by, drop = TRUE]]
  } else {
    group.by <- NULL
  }

  # Add shape column if specified
  if (!is.null(x = shape.by)) {
    if (!shape.by %in% colnames(object[[]])) {
      stop(paste0("Shaping variable (shape.by) ", shape.by, " not found in meta.data slot"), call. = F)
    }
    data[, shape.by] <- as.character(object[[shape.by, drop = TRUE]])
  }

  # Obtain array coordinates
  image.type <- "empty"
  c(data, xlim, ylim, image.type) %<-% obtain.array.coords(object, data, image.type, delim)

  # Raise error if features are not present in Seurat object
  if (ncol(x = data) < 3) {
    stop("None of the requested features were found: ",
         paste(features, collapse = ", "),
         " in slot ",
         slot,
         call. = FALSE)
  }

  if (unique(sapply(data[, features], class)) == "numeric") {
    data <- feature.scaler(data, features, min.cutoff, max.cutoff, spots)
  }

  # Generate HSV encoded colors
  if (verbose) cat(paste0("Defining Hue for ", length(x = features), " features ... \n"))
  hue_breaks <- seq(0, 1, length.out = length(x = features) + 1)[1:length(x = features)]
  hsv.matrix <- t(matrix(c(hue_breaks, rep(1, length(hue_breaks )), rep(1, length(hue_breaks))), ncol = 3))
  rownames(hsv.matrix) <- c("h", "s", "v")
  ann.cols <- hsv2hex(hsv.matrix)

  # Select highest V
  d <- list()
  if (verbose) cat("Converting values to HSV ... \n")
  for (i in 1:length(features)) {
    ftr <- features[i]
    s <- data.frame(data,
                    h = hue_breaks[i],
                    s = 1,
                    v = scales::rescale(data[, ftr]))
    s$cols <- hsv2hex(t(s[, c("h", "s", "v")]))
    d <- c(d, list(s))
  }
  red.cols <- data.frame()
  if (verbose) cat("Selecting HSV colors for each spot ... \n")
  for (i in 1:nrow(data)) {
    n <- which.max(unlist(lapply(d, function(x) {
      x[i, "v"]
    })))

    red.cols <- rbind(red.cols, d[[n]][i, ])
  }

  if (verbose) cat("Plotting features:",
                   ifelse(length(features) == 1, features,  paste0(paste(features[1:(length(features) - 1)], collapse = ", "), " and ", features[length(features)])))
  # Normal visualization -------------------------------------------------------------------------------------
  plot <- STPlot(data, data.type, group.by, shape.by, NULL, pt.size,
                 palette = NULL, cols = NULL, rev.cols = F, ncol, spot.colors = red.cols$cols,
                 center.zero = F, center.tissue = F, plot.title = "",
                 xlim, ylim, FALSE, theme = theme, ...)

  if (dark.theme) {
    plot <- plot + dark_theme()
  }
  plot <- plot +
    geom_point(data = data.frame(x = rep(-1, length(features)), y = rep(-1, length(features)), features), aes(x, y, colour = features)) +
    scale_color_manual(values = setNames(ann.cols, features))

  suppressWarnings({print(plot)})
}
