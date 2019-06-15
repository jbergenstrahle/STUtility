#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Dimensional reduction plots
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Dimensional reduction plot on ST coordinates
#'
#' Graphs the output of a dimensional reduction technique on a 2D scatter plot where each point is a
#' spot and it's positioned based on the ST array coordinates in a grid pattern. Spots are colored
#' by the dimensional reduction vectors specified by the user.
#'
#' Function built upon the DimPlot() function from Seurat (https://github.com/satijalab/seurat/blob/master/R/utilities.R)
#'
#' @param object Seurat object
#' @param dims Dimensions to plot, must numeric vectoir specifying number of dimensions to plot
#' @param spots Vector of spots to plot (default is all spots)
#' @param pt.size Adjust point size for plotting
#' @param reduction Which dimensionality reduction to use. If not specified, first searches for umap, then tsne, then pca
#' @param group.by Name of a metadata column to facet plot by (deault is sampleID)
#' @param shape.by If NULL, all points are circles (default). You can specify any
#' cell attribute (that can be pulled with FetchData) allowing for both
#' different colors and different shapes on cells
#' @param palette color palette to use for heatmap, see \code{\link{palette.select}}
#' @param rev.cols logical specifying whether colorscale should be reversed
#' @param combine Combine plots into a single gg object; note that if TRUE; themeing will not work when plotting multiple features
#' @param ncol Number of columns for display of samples for each dimensionality reduction vector
#' @param delim delimiter passed to \code{\link{GetCoords}} if adjusted ST coordinates are missing in the meta data
#' @param return.plot.list should the plots be returned as a list? By default, the plots are arranged into a grid
#' @param grid.ncol Number of columns for display when combining plots
#' @param ... Extra parameters passed on to \code{\link{STPlot}}
#'
#' @return A ggplot object
#'
#' @importFrom rlang !!
#' @importFrom ggplot2 ggplot facet_wrap vars sym
#' @importFrom viridis magma
#'
#' @export

ST.DimPlot <- function(
  object,
  dims = c(1, 2),
  spots = NULL,
  pt.size = 1,
  reduction = NULL,
  group.by = NULL,
  shape.by = NULL,
  palette = "MaYl",
  rev.cols = F,
  dark.theme = F,
  combine = TRUE,
  ncol = NULL,
  delim = NULL,
  return.plot.list = F,
  grid.ncol = NULL,
  ...
) {
  reduction <- reduction %||% {
    default.reductions <- c('umap', 'tsne', 'pca')
    object.reductions <- FilterObjects(object = object, classes.keep = 'DimReduc')
    reduc.use <- min(which(x = default.reductions %in% object.reductions))
    default.reductions[reduc.use]
  }
  spots <- spots %||% colnames(x = object)
  data <- Embeddings(object = object[[reduction]])[spots, dims, drop = FALSE]
  data <- as.data.frame(x = data)
  dims <- paste0(Key(object = object[[reduction]]), dims)

  if (!is.null(x = group.by)) {
    data[,  group.by] <- object[[group.by, drop = TRUE]]
  } else if ("sample" %in% colnames(object[[]])) {
    if (length(unique(object[["sample"]]))) {
      warning("sample column found in meta data but not specified as the group.by variable ...\n  setting sample as group.by variable ...")
      group.by <- "sample"
      data[, group.by] <- object[[group.by, drop = TRUE]]
    }
  }

  if (!is.null(x = shape.by)) {
    data[, shape.by] <- object[[shape.by, drop = TRUE]]
  }

  # Obtain array coordinates
  if (all(c("adj_x", "adj_y") %in% colnames(object[[]]))) {
    data <- cbind(data, setNames(object[[c("adj_x", "adj_y")]], nm = c("x", "y")))
  } else {
    if(is.null(delim)) {
      stop("adjusted coordinates are not present in meta data and delimiter is missing ...")
    }
    coords <- GetCoords(colnames(object), delim)
    data <- cbind(data, coords[, c("x", "y")])
  }

  # Create plots
  plots <- lapply(X = dims, FUN = function(d) {
    plot <- STPlot(data,
                   group.by,
                   d,
                   pt.size,
                   palette,
                   rev.cols,
                   ncol,
                   ...)

    if (dark.theme) {
      plot <- plot + dark_theme()
    }
    return(plot)
  })

  if (return.plot.list) {
    return(plots)
  } else {
    cowplot::plot_grid(plotlist = plots, ncol = grid.ncol)
  }
}


#' Graphs ST spots colored by continuous variable, e.g. dimensional reduction vector
#'
#' @param data data.frame containing (x, y) coordinates, a group vector and a continuous variable vector
#' @param group.by specifies column to facet the plots by, e.g. sampleID
#' @param variable name of continuous variable
#' @param pt.size point size of each ST spot
#' @param palette color palette used for spatial heatmap
#' @param rev.cols logical specifying whether colorscale should be reversed
#' @param ... parameters passed to geom_point()
#' @export

STPlot <- function(
  data,
  group.by,
  variable,
  pt.size = 1,
  palette = "MaYl",
  rev.cols = F,
  ncol = NULL,
  ...
) {
  cols <- palette.select(palette)(3)
  if (rev.cols) {
    cols <- rev(cols)
  }
  p <- ggplot() +
      geom_point(data = data, mapping = aes_string(x = "x", y = "64 - y", color = variable), size = pt.size, ...) +
      scale_x_continuous(limits = c(0, 67)) +
      scale_y_continuous(limits = c(0, 64)) +
      theme_void() +
      facet_wrap(as.formula(paste("~", group.by)), ncol = ncol) +
      labs(title = variable, color = "") +
      scale_color_gradientn(colours = cols)
}


#' Generates a dark theme for STPlot

dark_theme <- function() {
  theme(plot.background = element_rect(fill = "black"),
        plot.title = element_text(colour = "white"),
        strip.text = element_text(colour = 'white'),
        legend.title = element_text(colour = "white"),
        legend.text = element_text(colour = "white"))
}
