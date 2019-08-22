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
#' @param dims Dimensions to plot, a numeric vector specifying number of dimensions to plot
#' @param spots Vector of spots to plot (default is all spots)
#' @param plot.type Specify the type of plot to use (default: "spots"). Available options are; "spots", "smooth" and "tesselation"
#' @param blend Scale and blend expression values to visualize coexpression of two features (This options will override other coloring parameters)
#' @param min.cutoff,max.cutoff Vector of minimum and maximum cutoff values for each feature,
#'  may specify quantile in the form of 'q##' where '##' is the quantile (eg, 'q1', 'q10')
#' @param pt.size Adjust point size for plotting
#' @param reduction Which dimensionality reduction to use. If not specified, first searches for umap, then tsne, then pca
#' @param shape.by If NULL, all points are circles (default). You can specify any
#' cell attribute (that can be pulled with FetchData) allowing for both
#' different colors and different shapes on cells
#' @param delim delimiter passed to \code{\link{GetCoords}} if adjusted ST coordinates are missing in the meta data
#' @param return.plot.list should the plots be returned as a list? By default, the plots are arranged into a grid
#' @param grid.ncol Number of columns for display when combining plots
#' @param verbose Print messages
#'
#' @param ... Extra parameters passed on to \code{\link{STPlot}}
#'
#' @inheritParams STPlot
#'
#' @return A ggplot object
#'
#' @importFrom rlang !!
#' @importFrom ggplot2 ggplot facet_wrap vars sym
#' @importFrom viridis magma
#' @importFrom Seurat FetchData Embeddings
#'
#' @export

ST.DimPlot <- function (
  object,
  dims = c(1, 2),
  spots = NULL,
  plot.type = "spots",
  blend = FALSE,
  min.cutoff = NA,
  max.cutoff = NA,
  pt.size = 1,
  reduction = NULL,
  shape.by = NULL,
  palette = "MaYl",
  rev.cols = F,
  dark.theme = F,
  ncol = NULL,
  delim = NULL,
  return.plot.list = F,
  grid.ncol = NULL,
  center.zero = T,
  channels.use = NULL,
  center.tissue = FALSE,
  verbose = FALSE,
  ...
) {

  reduction <- reduction %||% {
    default.reductions <- c('umap', 'tsne', 'pca', 'ica')
    object.reductions <- FilterObjects(object = object, classes.keep = 'DimReduc')
    reduc.use <- min(which(x = default.reductions %in% object.reductions))
    default.reductions[reduc.use]
  }

  # prepare data
  signs <- sign(dims); dims <- abs(dims)
  spots <- spots %||% colnames(x = object)
  data <- Embeddings(object = object[[reduction]])[spots, dims, drop = FALSE]
  if (verbose) cat(paste0("Selected ", length(spots), " spots"))
  data <- as.data.frame(x = t(t(data)*signs))
  dims <- paste0(Key(object = object[[reduction]]), dims)

  # Select colorscale if palette is NULL
  palette.info <- palette.select(info = T)
  palette <- palette %||% {
    palette <- subset(palette.info, category == "div")$palette[1]
  }

  # Check that the number of dimensions are 2 or three if blending is active
  if (blend & !length(x = dims) %in% c(2, 3)) {
    stop(paste0("Blending dim plots only works with two or three dimensions. \n",
                "Number of dimensions provided: ", length(x = dims)), call. = F)
  }

  # Check that group.by variable is present in meta.data slot, otherwise assume that there's only one sample present in Seurat object
  if ("sample" %in% colnames(object[[]])) {
    group.by <- "sample"
    data[,  group.by] <- object[[group.by, drop = TRUE]]
  } else {
    group.by <- NULL
  }

  # Extract shape.by column from meta data if applicable
  if (!is.null(x = shape.by)) {
    if (!shape.by %in% colnames(object[[]])) {
      stop(paste0("Shaping variable (shape.by) ", shape.by, " not found in meta.data slot"), call. = F)
    }
    data[, shape.by] <- as.character(object[[shape.by, drop = TRUE]])
  }

  # Obtain array coordinates
  image.type <- NULL
  if (all(c("warped_x", "warped_y") %in% colnames(object[[]]))) {
    data <- cbind(data, setNames(object[[c("warped_x", "warped_y")]], nm = c("x", "y")))
    xlim <- c(0, max(unlist(lapply(object@tools$dims, function(x) as.numeric(x[2]))))); ylim <- c(0, max(unlist(lapply(object@tools$dims, function(x) as.numeric(x[3])))))
    image.type <- "processed"
  } else if ("raw" %in% names(object@tools)) {
    data <- cbind(data, setNames(object[[c("pixel_x", "pixel_y")]], nm = c("x", "y")))
    xlim <- c(0, max(unlist(lapply(object@tools$dims, function(x) as.numeric(x[2]))))); ylim <- c(0, max(unlist(lapply(object@tools$dims, function(x) as.numeric(x[3])))))
  } else if (all(c("ads_x", "ads_y") %in% colnames(object[[]]))) {
    data <- cbind(data, setNames(object[[c("ads_x", "ads_y")]], nm = c("x", "y")))
    xlim <- ylim <- NULL
  } else {
    if(is.null(delim)) {
      stop("adjusted coordinates are not present in meta data and delimiter is missing ...")
    }
    coords <- GetCoords(colnames(object), delim)
    data <- cbind(data, coords[, c("x", "y")])
    xlim <- ylim <- NULL
  }

  # Scale data values
  data <- feature.scaler(data, dims, min.cutoff, max.cutoff, spots)

  # blend colors or plot each dimension separately
  if (blend) {
    colored.data <- apply(data[, 1:(ncol(data) - 3)], 2, rescale)
    channels.use <- channels.use %||% c("red", "green", "blue")[1:ncol(colored.data)]

    if (verbose) cat(paste0("Blending colors for dimensions ",
                            paste0(ifelse(length(dims) == 2, paste0(dims[1], " and ", dims[2]), paste0(dims[1], dims[2], " and ", dims[2]))),
                            ": \n", paste(paste(dims, channels.use, sep = ":"), collapse = "\n")))

    spot.colors <- ColorBlender(colored.data, channels.use)
    data <- data[, (ncol(data) - 2):ncol(data)]
    plot <- STPlot(data, data.type = "numeric", group.by, shape.by, NULL,
                   pt.size, palette, rev.cols, ncol, spot.colors, center.zero, center.tissue, xlim, ylim,
                   plot.title = paste(paste(dims, channels.use, sep = ":"), collapse = ", "), ...)
    if (dark.theme) {
      plot <- plot + dark_theme()
    }
    return(plot)
  } else {
    spot.colors <- NULL
    if (verbose) cat("Plotting dimensions:",
                     ifelse(length(dims) == 1, dims,  paste0(paste(dims[1:(length(dims) - 1)], collapse = ", "), " and ", dims[length(dims)])))

    # Create plots
    if (plot.type == "spots") {
      # Normal visualization -------------------------------------------------------------------------------------
      plots <- lapply(X = dims, FUN = function(d) {
        plot <- STPlot(data, data.type = "numeric", group.by, shape.by, d, pt.size,
                       palette, rev.cols, ncol, spot.colors, center.zero, center.tissue, NULL, xlim, ylim, ...)

        if (dark.theme) {
          plot <- plot + dark_theme()
        }
        return(plot)
      })

      # Draw plots
      if (return.plot.list) {
        return(plots)
      } else {
        plot_grid(plotlist = plots, ncol = grid.ncol)
      }

    } else if (plot.type == "smooth") {
      # Smooth visualization -------------------------------------------------------------------------------------
      plots <- lapply(X = dims, FUN = function(d) {
        plot <- SmoothPlot(data, image.type, data.type = "numeric", group.by, d,
                           palette, rev.cols, ncol, center.zero, xlim, ylim, ...)
        return(plot)
      })

      # Draw plots
      grid.ncol <- grid.ncol %||% round(sqrt(length(x = plots)))
      grid.nrow <- ceiling(length(x = plots)/grid.ncol)

      stack <- c()
      for (i in 1:grid.nrow) {
        i <- i - 1
        stack <- c(stack, image_append(Reduce(c, plots[(i*grid.ncol + 1):(i*grid.ncol + grid.ncol)])))
      }

      final_img <- image_append(Reduce(c, stack), stack = T)
      print(final_img)
    }
  }
}



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Feature plots
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#'  Visualize 'features' on an ST array grid
#'
#' Colors spots on an an ST array grid according to a 'feature'
#' (i.e. gene expression (raw counts or scaled) and features available in the meta data slot)
#'
#' Function built upon the FeaturePlot() function from Seurat (https://github.com/satijalab/seurat/blob/master/R/visualization.R)
#'
#' @param object Seurat object
#' @param features
#' \itemize{
#'     \item An \code{Assay} feature (e.g. a gene name - "MS4A1")
#'     \item A column name from meta.data (e.g. mitochondrial percentage - "percent.mito")
#' }
#' @param indices Numeric vector specifying sample indices to include in plot. Default is to show all samples.
#' @param spots Vector of spots to plot (default is all spots)
#' @param plot.type Specify the type of plot to use (default: "spots"). Available options are; "spots", "smooth" and "tesselation"
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
#'
#' @return A ggplot object
#' @export

ST.FeaturePlot <- function (
  object,
  features,
  split.labels = FALSE,
  indices = NULL,
  spots = NULL,
  plot.type = "spots",
  min.cutoff = NA,
  max.cutoff = NA,
  slot = "data",
  blend = FALSE,
  pt.size = 1,
  shape.by = NULL,
  palette = NULL,
  rev.cols = FALSE,
  dark.theme = FALSE,
  ncol = NULL,
  delim = NULL,
  return.plot.list = FALSE,
  grid.ncol = NULL,
  center.zero = FALSE,
  channels.use = NULL,
  center.tissue = FALSE,
  theme = theme_void(),
  verbose = FALSE,
  ...
) {
  spots <- spots %||% colnames(x = object)
  data <- FetchData(object = object, vars = c(features), cells = spots, slot = slot)

  data.type <- unique(sapply(data, class))

  if (length(data.type) > 1 & !all(data.type %in% c("numeric", "integer"))) {
    stop("Mixed classes (", paste(unique(sapply(data, class)), collapse = ", "), ") are not allowed in features ... ")
  }

  if (!blend && length(x = features) %in% c(2, 3) & !all(data.type %in% c("numeric", "integer"))) {
    stop("Blending feature plots only works with two or three numeric features")
  }

  # Select colorscale
  palette.info <- palette.select(info = T)
  palette <- palette %||% {
    palette <- subset(palette.info, category == "seq")$palette[1]
  }

  # Check that group.by variable is present in meta.data slot, otherwise assume that there's only one sample present in Seurat object
  if ("sample" %in% colnames(object[[]])) {
    group.by <- "sample"
    data[,  group.by] <- object[[group.by, drop = TRUE]]
  } else {
    group.by <- NULL
  }

  if (!is.null(x = shape.by)) {
    if (!shape.by %in% colnames(object[[]])) {
      stop(paste0("Shaping variable (shape.by) ", shape.by, " not found in meta.data slot"), call. = F)
    }
    data[, shape.by] <- as.character(object[[shape.by, drop = TRUE]])
  }

  # Obtain array coordinates
  image.type <- NULL
  if (all(c("warped_x", "warped_y") %in% colnames(object[[]]))) {
    data <- cbind(data, setNames(object[[c("warped_x", "warped_y")]], nm = c("x", "y")))
    xlim <- c(0, max(unlist(lapply(object@tools$dims, function(x) as.numeric(x[2]))))); ylim <- c(0, max(unlist(lapply(object@tools$dims, function(x) as.numeric(x[3])))))
    image.type <- "processed"
  } else if ("raw" %in% names(object@tools)) {
    data <- cbind(data, setNames(object[[c("pixel_x", "pixel_y")]], nm = c("x", "y")))
    xlim <- c(0, max(unlist(lapply(object@tools$dims, function(x) as.numeric(x[2]))))); ylim <- c(0, max(unlist(lapply(object@tools$dims, function(x) as.numeric(x[3])))))
  } else if (all(c("ads_x", "ads_y") %in% colnames(object[[]]))) {
    data <- cbind(data, setNames(object[[c("ads_x", "ads_y")]], nm = c("x", "y")))
    xlim <- ylim <- NULL
  } else {
    if(is.null(delim)) {
      stop("adjusted coordinates are not present in meta data and delimiter is missing ...")
    }
    coords <- GetCoords(colnames(object), delim)
    data <- cbind(data, coords[, c("x", "y")])
    xlim <- ylim <- NULL
  }

  if (ncol(x = data) < 3) {
    stop("None of the requested features were found: ",
         paste(features, collapse = ", "),
         " in slot ",
         slot,
         call. = FALSE)
  }

  if (class(unique(sapply(data, class))) == "numeric") {
    data <- feature.scaler(data, min.cutoff, max.cutoff, spots)
  }

  # Subset by index
  if (!all(as.character(indices) %in% data[, group.by])) stop(paste0("Index out of range. "), call. = FALSE)
  if (!is.null(indices)) data <- data[data[, group.by] %in% as.character(indices), ]

  if (blend) {
    colored.data <- apply(data[, 1:(ncol(data) - 3)], 2, rescale)
    channels.use <- channels.use %||% c("red", "green", "blue")[1:ncol(colored.data)]

    if (verbose) cat(paste0("Blending colors for features ",
                            paste0(ifelse(length(features) == 2, paste0(features[1], " and ", features[2]), paste0(features[1], features[2], " and ", features[2]))),
                            ": \n", paste(paste(features, channels.use, sep = ":"), collapse = "\n")))

    spot.colors <- ColorBlender(colored.data, channels.use)
    data <- data[, (ncol(data) - 3):ncol(data)]
    plot <- STPlot(data, data.type, group.by, shape.by, NULL, pt.size,
                   palette, rev.cols, ncol, spot.colors, center.zero, center.tissue,
                   plot.title = paste(paste(features, channels.use, sep = ":"), collapse = ", "),
                   xlim, ylim, FALSE, theme, ...)
    if (dark.theme) {
      plot <- plot + dark_theme()
    }
    return(plot)
  } else {
    spot.colors <- NULL
    if (verbose) cat("Plotting features:",
                     ifelse(length(features) == 1, features,  paste0(paste(features[1:(length(features) - 1)], collapse = ", "), " and ", features[length(features)])))
    # Create plots
    if (plot.type == "spots") {
      # Normal visualization -------------------------------------------------------------------------------------
      plots <- lapply(X = features, FUN = function(d) {
        plot <- STPlot(data, data.type, group.by, shape.by, d, pt.size,
                       palette, rev.cols, ncol, spot.colors, center.zero,
                       center.tissue, NULL, xlim, ylim, split.labels, theme, ...)

        if (dark.theme) {
          plot <- plot + dark_theme()
        }
        return(plot)
      })

      # Draw plots
      if (return.plot.list) {
        return(plots)
      } else {
        plot_grid(plotlist = plots, ncol = grid.ncol)
      }

    } else if (plot.type == "smooth") {
      if (data.type %in% c("factor", "character")) stop(paste0("Smoothing has not yet been implemented for categorical variables"), call. = FALSE)
      # Smooth visualization -------------------------------------------------------------------------------------
      plots <- lapply(X = features, FUN = function(d) {
        plot <- SmoothPlot(data, image.type, data.type, group.by, d,
                       palette, rev.cols, ncol, center.zero, xlim, ylim, ...)
        return(plot)
      })

      # Draw plots
      grid.ncol <- grid.ncol %||% round(sqrt(length(x = plots)))
      grid.nrow <- ceiling(length(x = plots)/grid.ncol)

      stack <- c()
      for (i in 1:grid.nrow) {
        i <- i - 1
        stack <- c(stack, image_append(Reduce(c, plots[(i*grid.ncol + 1):(i*grid.ncol + grid.ncol)])))
      }

      final_img <- image_append(Reduce(c, stack), stack = T)
      print(final_img)
    }
  }
}


#' Graphs ST spots colored by continuous variable, e.g. dimensional reduction vector
#'
#' @param data data.frame containing (x, y) coordinates, a group vector and a continuous variable vector
#' @param data.type type of data, e.g. numeric or integer
#' @param group.by specifies column to facet the plots by, e.g. sample
#' @param shape.by specifies column to shape points by, e.g. morphological region
#' @param variable name of continuous variable
#' @param pt.size point size of each ST spot
#' @param palette color palette used for spatial heatmap
#' @param rev.cols logical specifying whether colorscale should be reversed
#' @param ncol number of columns in \code{facet_wrap}
#' @param spot.colors character vector woth color names that overrides default coloring with ggplot2
#' @param center.zero should the colorscale be centered around 0? Set to TRUE for scaled data
#' @param center.tissue Adjust coordinates so that the center of the tissue is in the middle of the array along the y-axis
#' @param plot.title Add title to plot
#' @param xlim,ylim Set x/y-axis limits
#' @param split.labels Only active if the features are specified by character vectors.
#' The plot will be split into one plot for each group label, highlighting the labelled spots.
#' @param theme Object of class "theme" used to change the background theme
#' @param ... parameters passed to geom_point()
#'
#' @importFrom ggplot2 geom_point aes_string scale_x_continuous scale_y_continuous theme_void theme_void labs scale_color_gradient2 scale_color_gradientn scale_color_manual
#'
#' @export

STPlot <- function (
  data,
  data.type = NULL,
  group.by,
  shape.by,
  variable,
  pt.size = 1,
  palette = "MaYl",
  rev.cols = F,
  ncol = NULL,
  spot.colors = NULL,
  center.zero = TRUE,
  center.tissue = F,
  plot.title = NULL,
  xlim = NULL,
  ylim = NULL,
  split.labels = FALSE,
  theme = theme_void(),
  ...
) {

  xlim <- xlim %||% c(0, 67); ylim <- ylim %||% c(0, 64)

  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }

  if (is.null(spot.colors)) {
    if (class(data[, variable]) == "factor") {
      label.colors <- gg_color_hue(length(levels(data[, variable])))
      names(label.colors) <- levels(data[, variable])
    } else if (class(data[, variable]) == "character") {
      label.colors <- gg_color_hue(length(unique(data[, variable])))
      names(label.colors) <- unique(data[, variable])
    }
  }

  # Stop if split.labels is activated and there are more than 1 samples
  if (any(data.type %in% c("character", "factor")) & split.labels) {
    if (length(unique(as.character(data[, group.by]))) > 1) stop(paste0("Splitting of group labels only work for one sample. Please set a single sample index with indices. "), call. = FALSE)
    plot.title <- paste0("Sample ", unique(as.character(data[, group.by])), ": ", variable)
    new.data <- data.frame()
    # Order by decreasing size
    if (class(data[, variable]) != "factor") data[, variable] <- factor(data[, variable], levels = unique(data[, variable]))
    levels.keep <- levels(data[, variable])

    for (lbl in levels.keep) {
      dt <- data
      dt[, variable] <- ifelse(dt[, variable] == lbl, lbl, "-")
      dt[, group.by] <- lbl
      new.data <- rbind(new.data, dt)
    }

    data <- new.data
    levels.keep <- c("-", levels.keep)
    data[, variable] <- factor(data[, variable], levels = levels.keep)
    label.colors <- c("-" = "lightgray", label.colors)
  }

  # Center tissue along y-axis
  if (center.tissue) {
    if (!is.null(group.by)) {
      data <- do.call(rbind, lapply(split(data, data[, group.by]), function(d) {
        d[, "y"] <- d[, "y"] - median(d[, "y"]) + ylim[2]/2
        return(d)
      }))
    } else {
      data[, "y"] <- data[, "y"] - median(data[, "y"]) + ylim[2]/2
    }
  }

  # Obtain colors from selected palette
  cols <- palette.select(palette)(3)
  if (palette == "heat") cols <- palette.select(palette)(4)
  if (rev.cols) {
    cols <- rev(cols)
  }

  # Create new plot
  p <- ggplot()
  # Make sure that levels are correct
  data[, group.by] <- factor(data[, group.by], levels = unique(data[, group.by]))
  if (length(spot.colors) > 0) {

    # Add shape aesthetic and blend colors if blend is active
    if (!is.null(shape.by)) {
      p <- p + geom_point(data = data, mapping = aes_string(x = "x", y = paste0(ylim[2], " - y"), shape = shape.by), color = spot.colors, size = pt.size, ...)
    } else {
      p <- p + geom_point(data = data, mapping = aes_string(x = "x", y = paste0(ylim[2], " - y")), color = spot.colors, size = pt.size, ...)
    }

  } else {

    # Add shape aesthetic only
    if (!is.null(shape.by)) {
      p <- p + geom_point(data = data, mapping = aes_string(x = "x", y = paste0(ylim[2], " - y"), color = paste0("`", variable, "`"), shape = shape.by), size = pt.size, ...)
    } else {
      p <- p + geom_point(data = data, mapping = aes_string(x = "x", y = paste0(ylim[2], " - y"), color = paste0("`", variable, "`")), size = pt.size, ...)
    }

  }

  # Add ST array dimensions and plot title
  p <- p +
    scale_x_continuous(limits = xlim) +
    scale_y_continuous(limits = ylim) +
    #theme_void() +
    labs(title = ifelse(!is.null(plot.title), plot.title, variable), color = "")

  # Set theme
  p <- p + theme

  # Facet plots by group variable
  if (!is.null(group.by)) {
    p <- p +
      facet_wrap(as.formula(paste("~", group.by)), ncol = ncol)
  }

  # Center colorscale at 0
  if (center.zero & !any(data.type %in% c("character", "factor"))) {
    p <- p +
      scale_color_gradient2(low = cols[1], mid = cols[2], high = cols[3], midpoint = 0)
  } else if (any(data.type %in% c("character", "factor"))) {
    p <- p +
      labs(color = variable) +
      scale_color_manual(values = label.colors)
  } else {
    p <- p +
      scale_color_gradientn(colours = cols)
  }

  return(p)
}


#' Graphs a smooth interpolation heatmap colored by continuous variable, e.g. dimensional reduction vector
#'
#' @param image.type Specifies the image is "processed", otherwise NULL
#' @param darken Adds black color to the bottom of the colorscale to increase the contrast of the colorscale
#' @param whiten Adds white color to the top of the colorscale to increase the contrast of the colorscale
#'
#' @inheritParams STPlot
#'
#' @importFrom ggplot2 ggplot aes geom_raster scale_x_continuous scale_y_continuous theme_void guides scale_fill_gradient2 labs scale_fill_gradientn ggsave
#' @importFrom magick image_read image_border image_annotate image_composite
#' @importFrom imager as.cimg
#' @importFrom grDevices as.raster

SmoothPlot <- function (
  data,
  image.type,
  data.type = NULL,
  group.by,
  variable,
  palette = "MaYl",
  rev.cols = F,
  ncol = NULL,
  center.zero = TRUE,
  xlim = NULL,
  ylim = NULL,
  darken = FALSE,
  whiten = FALSE,
  ...
) {

  xlim <- xlim %||% c(0, 67); ylim <- ylim %||% c(0, 64)

  image.masks <- NULL
  if (image.type == "processed") {
    image.masks <- object@tools$processed.masks
  } else if ("masked.masks" %in% names(object@tools)) {
    image.masks <- object@tools$masked.masks
  }
  samplenames <- names(object@tools$raw)

  # Set colors
  ncolors <- ifelse(palette == "heat", 4, 5)
  cols <- palette.select(palette)(ncolors)
  if (rev.cols) {
    cols <- rev(cols)
  }

  if (darken) cols <- c("black", cols); if (whiten) cols <- c(cols, "white")

  val.limits <- range(data[, variable])

  # Create legend
  lg <- g_legend(data, variable, center.zero, cols, val.limits)

  # Subset only based on one value's expression
  p.list <- lapply(1:length(unique(data[, group.by])), function(i) {
    xdim <- object@tools$xdim; ydim <- round(as.numeric(object@tools$dims[[i]][2])/(as.numeric(object@tools$dims[[i]][2])/object@tools$xdim))
    data.subset <- subset(data, sample == i)
    #data.subset <- data.subset[data.subset[, variable] != 0, ]
    s.xy <- as.numeric(object@tools$dims[[1]][, 2:3])/c(xdim, ydim)

    data.subset[, c("x", "y")] <- data.subset[, c("x", "y")]/s.xy
    x <- data.subset[, "x"]; y <- data.subset[, "y"]
    min.x <- min(x); min.y <- min(y); max.x <- max(x); max.y <- max(y)
    tissue.width <- max.x - min.x; tissue.height <- max.y - min.y;

    # Run interpolation
    s =  interp(data.subset[, "x"], data.subset[, "y"], data.subset[, variable], nx = tissue.width, ny = tissue.height, extrap = FALSE, linear = FALSE, xo = 1:xdim, yo = 1:ydim)

    z <- t(s$z)
    x <- 1:ncol(z)
    y <- 1:nrow(z)
    gg <- data.frame(x = rep(x, each = nrow(z)), y = rep(y, times = ncol(z)), val = as.numeric(z))
    gg$val[gg$val < val.limits[1]] <- val.limits[1]; gg$val[gg$val > val.limits[2]] <- val.limits[2]


    p <- ggplot(gg, aes(x, 64 - y, fill = val)) +
      geom_raster(interpolate = TRUE) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      theme_void() +
      guides(fill = FALSE)

    # Center colorscale at 0
    if (center.zero) {
      p <- p +
        scale_fill_gradient2(low = cols[1], mid = cols[median(seq_along(cols))], high = cols[length(cols)], midpoint = 0, na.value = "#000000", limits = val.limits)
    } else if (any(data.type %in% c("character", "factor"))) {
      p <- p +
        labs(fill = variable)
    } else {
      p <- p +
        scale_fill_gradientn(colours = cols, na.value = "#000000", limits = val.limits)
    }

    tmp.file <- tempfile(pattern = "", fileext = ".png")

    png(width = xdim, height = ydim, file = tmp.file)
    par(mar = c(0, 0, 0, 0))
    plot(p)
    dev.off()

    p <- image_read(tmp.file)

    msk <- as.cimg(image.masks[[i]])
    masked.plot <- as.raster(magick2cimg(p)*(msk/255))
    masked.plot <- as.raster(image_annotate(image_read(masked.plot), text = samplenames[i], size = round(xdim/10), color = "#FFFFFF"))
    return(masked.plot)
  })

  # Draw on new device
  ncol <- ncol %||% ceiling(sqrt(length(p.list)))
  nrow <- ceiling(length(p.list)/ncol)

  tmp.file <- tempfile(pattern = "", fileext = ".png")

  png(width = xdim*ncol, height = ydim*nrow, file = tmp.file)
  par(mfrow = c(nrow, ncol), mar = c(0, 0, 0, 0), bg = "#000000")
  for (p in p.list) {
    plot(p)
  }
  dev.off()

  im <- image_read(tmp.file)
  im <- image_border(im, "#000000", paste(xdim/10, xdim/10, sep = "x"))
  im <- image_annotate(im, text = variable, size = round(xdim/10), color = "#FFFFFF")

  # Draw legend
  tmp.file <- tempfile(pattern = "", fileext = ".png")
  ggsave(plot = lg, width = 2.8/5, height = 7.8/5, filename = tmp.file, dpi = 150, units = "in")
  lgim <- image_read(tmp.file)

  im <- image_composite(image = im, composite_image = lgim, offset = paste0("+", xdim*ncol, "+", (ydim*nrow)/2))

  return(im)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Dimensional reduction plots on HE images
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#'  Visualize dimensionality reduction vectors on an ST array grid overlayed on top of HE image
#'
#' Colors spots on an an ST array grid according to a dimension
#' (i.e. gene expression (raw counts or scaled) and features available in the meta data slot)
#'
#' @param sample.index Index specifying the sample that you want to use for plotting
#' @param type Image type to plot on
#' @param slot Which slot to pull expression data from?
#' @param ... Extra parameters passed on to \code{\link{ST.ImagePlot}}
#'
#' @inheritParams ST.ImagePlot
#' @inheritParams ST.DimPlot
#' @importFrom cowplot plot_grid
#' @importFrom scales rescale
#'
#' @return A ggplot object
#' @export

DimOverlay <- function(
  object,
  dims = c(1:2),
  sample.index = 1,
  type = NULL,
  min.cutoff = NA,
  max.cutoff = NA,
  blend = FALSE,
  pt.size = 1,
  reduction = NULL,
  shape.by = NULL,
  palette = NULL,
  rev.cols = FALSE,
  dark.theme = FALSE,
  delim = NULL,
  return.plot.list = FALSE,
  grid.ncol = NULL,
  center.zero = FALSE,
  channels.use = NULL,
  verbose = FALSE,
  ...
) {

  # Check length of sample index
  if (length(sample.index) > 1) stop(paste0("Only one sample index can be selected."), call. = FALSE)

  reduction <- reduction %||% {
    default.reductions <- c('umap', 'tsne', 'pca', 'ica')
    object.reductions <- FilterObjects(object = object, classes.keep = 'DimReduc')
    reduc.use <- min(which(x = default.reductions %in% object.reductions))
    default.reductions[reduc.use]
  }

  type <- type %||% {
    choices <- c("processed", "masked", "raw", "processed.masks", "masked.masks")
    choices[min(as.integer(na.omit(match(names(object@tools), choices))))]
  }

  # Check that selected image type is present in Seurat object
  msgs <- c("raw" = "LoadImages()", "masked" = "MaskImages()", "processed" = "WarpImages()", "masked.masks" = "MaskImages()", "processed.masks" = "WarpImages()")
  if (!type %in% names(object@tools)) stop(paste0("You need to run ", msgs[type], " before using FeatureOverlay() on '", type, "' images"), call. = FALSE)

  # Check that image pointer is alive)
  if (!sample.index %in% names(object@tools[[type]])) {
    stop(paste0("sample.index ", sample.index, " does not match any of the images present in the Seurat object or is out of range"), call. = T)
  }
  image <- as.raster(image_annotate(image_read(object@tools[[type]][[sample.index]]), text = paste(sample.index), size = round(object@tools$xdim/10)))
  imdims <- object@tools$dims[[sample.index]]

  group.var = "sample"
  if (group.var %in% colnames(object[[]])) {
    sample.index <- ifelse(class(sample.index) == "numeric", unique(object[[group.var, drop = T]])[sample.index], sample.index)
    # Select spots matching sample.index
    spots <- colnames(object)[object[[group.var, drop = T]] == sample.index]
    if (verbose) cat(paste0("Selected ", length(spots), " spots matching index ", sample.index))
  } else {
    # Assuming that there's only one sample in the Seurat object
    spots <- NULL
    spots <- spots %||% colnames(x = object)
    if (verbose) cat("Selecting all spots")
  }

  signs <- sign(dims); dims <- abs(dims)
  data <- Embeddings(object = object[[reduction]])[spots, dims, drop = FALSE]
  data <- as.data.frame(x = t(t(data)*signs))
  dims <- paste0(Key(object = object[[reduction]]), dims)

  # Select colorscale
  palette.info <- palette.select(info = T)
  palette <- palette %||% {
    palette <- subset(palette.info, category == "seq")$palette[1]
  }

  # Check that the number of dimensions are 2 or three if blending is active
  if (blend & !length(x = dims) %in% c(2, 3)) {
    stop(paste0("Blending dim plots only works with two or three dimensions. \n",
                "Number of dimensions provided: ", length(x = dims)), call. = F)
  }

  # Obtain array coordinates
  px.ids <- ifelse(rep(type %in% c("raw", "masked", "masked.masks"), 2), c("pixel_x", "pixel_y"), c("warped_x", "warped_y"))

  if (all(px.ids %in% colnames(object[[]]))) {
    data <- cbind(data, setNames(object[[px.ids]][spots, ], nm = c("x", "y")))
  } else {
    stop(paste0(paste(px.ids, collapse = " and "), " coordinates are not present in meta data."), call. = FALSE)
  }

  data <- feature.scaler(data, dims, min.cutoff, max.cutoff, spots)
  data[, group.var] <- sample.index

  if (blend) {
    colored.data <- apply(data[, 1:(ncol(data) - 2)], 2, rescale)
    channels.use <- channels.use %||% c("red", "green", "blue")[1:ncol(colored.data)]

    if (verbose) cat(paste0("Blending colors from features ", paste(paste(dims, channels.use, sep = ":"), collapse = ", ")))

    spot.colors <- ColorBlender(colored.data, channels.use)
    data <- data[, (ncol(data) - 1):ncol(data)]
    plot <- ST.ImagePlot(data, data.type = "numeric", group.var, shape.by, variable, image, imdims, pt.size, palette,
                         rev.cols, ncol = NULL, spot.colors, center.zero,
                         plot.title = paste(paste(features, channels.use, sep = ":"), collapse = ", "), FALSE, ...)
    if (dark.theme) {
      plot <- plot + dark_theme()
    }
    return(plot)
  } else {
    spot.colors <- NULL

    if (verbose) cat("Plotting features:",
                     ifelse(length(dims) == 1, dims,  paste0(paste(dims[1:(length(dims) - 1)], collapse = ", "), " and ", dims[length(dims)])))

    # Create plots
    plots <- lapply(X = dims, FUN = function(d) {
      plot <- ST.ImagePlot(data, data.type = "numeric", group.var, shape.by, d, image, imdims, pt.size, palette,
                           rev.cols, ncol = NULL, spot.colors, center.zero, plot.title = NULL, FALSE, ...)

      if (dark.theme) {
        plot <- plot + dark_theme()
      }
      return(plot)
    })

    if (return.plot.list) {
      return(plots)
    } else {
      plot_grid(plotlist = plots, ncol = grid.ncol)
    }
  }
}




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Feature plots on HE images
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#'  Visualize 'features' on an ST array grid overlayed on top of HE image
#'
#' Colors spots on an an ST array grid according to a 'feature'
#' (i.e. gene expression (raw counts or scaled) and features available in the meta data slot)
#'
#' @param sample.index Index specifying the sample that you want to use for plotting
#' @param type Image type to plot on
#' @param slot Which slot to pull expression data from?
#' @param ... Extra parameters passed on to \code{\link{ST.ImagePlot}}
#'
#' @inheritParams ST.ImagePlot
#' @inheritParams ST.FeaturePlot
#' @importFrom cowplot plot_grid
#' @importFrom scales rescale
#'
#' @return A ggplot object
#' @export

FeatureOverlay <- function(
  object,
  features,
  sample.index = 1,
  type = NULL,
  min.cutoff = NA,
  max.cutoff = NA,
  slot = "data",
  blend = FALSE,
  pt.size = 2,
  shape.by = NULL,
  palette = NULL,
  rev.cols = FALSE,
  dark.theme = FALSE,
  delim = NULL,
  return.plot.list = FALSE,
  grid.ncol = NULL,
  center.zero = FALSE,
  channels.use = NULL,
  verbose = FALSE,
  split.labels = FALSE,
  ...
) {

  # Check length of sample index
  if (length(sample.index) > 1) stop(paste0("Only one sample index can be selected."), call. = FALSE)

  type <- type %||% {
    choices <- c("processed", "masked", "raw", "processed.masks", "masked.masks")
    choices[min(as.integer(na.omit(match(names(object@tools), choices))))]
  }

  # Check that selected image type is present in Seurat object
  msgs <- c("raw" = "LoadImages()", "masked" = "MaskImages()", "processed" = "WarpImages()", "masked.masks" = "MaskImages()", "processed.masks" = "WarpImages()")
  if (!type %in% names(object@tools)) stop(paste0("You need to run ", msgs[type], " before using FeatureOverlay() on '", type, "' images"), call. = FALSE)

  # Check that image pointer is alive)
  if (!sample.index %in% names(object@tools[[type]])) {
    stop(paste0("sample.index ", sample.index, " does not match any of the images present in the Seurat object or is out of range"), call. = T)
  }
  image <- as.raster(image_annotate(image_read(object@tools[[type]][[sample.index]]), text = paste(sample.index), size = round(object@tools$xdim/10)))
  imdims <- object@tools$dims[[sample.index]]

  group.var = "sample"
  if (group.var %in% colnames(object[[]])) {
    sample.index <- ifelse(class(sample.index) == "numeric", unique(object[[group.var, drop = T]])[sample.index], sample.index)
    # Select spots matching sample.index
    spots <- colnames(object)[object[[group.var, drop = T]] == sample.index]
    if (verbose) cat(paste0("Selected ", length(spots), " spots matching index ", sample.index))
  } else {
    # Assuming that there's only one sample in the Seurat object
    spots <- NULL
    spots <- spots %||% colnames(x = object)
    if (verbose) cat("Selecting all spots")
  }

  data <- FetchData(object = object, vars = c(features), cells = spots, slot = slot)
  data.type <- unique(sapply(data, class))

  if ((blend & !length(x = features) %in% c(2, 3)) & !all(data.type %in% c("numeric", "integer"))) {
    stop(paste0("Blending feature plots only works with two or three features of class numeric/integer. \n",
                "Number of features provided: ", length(x = features), "\n",
                "feature class: ", data.type), call. = F)
  }

  # Select colorscale
  palette.info <- palette.select(info = T)
  palette <- palette %||% {
    palette <- subset(palette.info, category == "seq")$palette[1]
  }

  # Obtain array coordinates
  px.ids <- ifelse(rep(type %in% c("raw", "masked", "masked.masks"), 2), c("pixel_x", "pixel_y"), c("warped_x", "warped_y"))

  if (all(px.ids %in% colnames(object[[]]))) {
    data <- cbind(data, setNames(object[[px.ids]][spots, ], nm = c("x", "y")))
  } else {
    stop(paste0(paste(px.ids, collapse = " and "), " coordinates are not present in meta data."), call. = FALSE)
  }

  if (ncol(x = data) < 3) {
    stop("None of the requested features were found: ",
         paste(features, collapse = ", "),
         " in slot ",
         slot,
         call. = FALSE)
  }

  if (class(unique(sapply(data, class))) == "numeric") {
    data <- feature.scaler(data, min.cutoff, max.cutoff, spots)
  }

  # Add index column
  data[, group.var] <- sample.index

  if (blend) {
    colored.data <- apply(data[, 1:(ncol(data) - 3)], 2, rescale)
    channels.use <- channels.use %||% c("red", "green", "blue")[1:ncol(colored.data)]

    if (verbose) cat(paste0("Blending colors from features ", paste(paste(features, channels.use, sep = ":"), collapse = ", ")))

    spot.colors <- ColorBlender(colored.data, channels.use)
    data <- data[, (ncol(data) - 2):ncol(data)]
    plot <- ST.ImagePlot(data, data.type, group.var, shape.by, variable, image, imdims, pt.size, palette,
                         rev.cols, ncol = NULL, spot.colors, center.zero,
                         plot.title = paste(paste(features, channels.use, sep = ":"), collapse = ", "), split.labels)#, ...)
    if (dark.theme) {
      plot <- plot + dark_theme()
    }
    return(plot)
  } else {
    spot.colors <- NULL

    if (verbose) cat("Plotting features:",
                     ifelse(length(features) == 1, features,  paste0(paste(features[1:(length(features) - 1)], collapse = ", "), " and ", features[length(features)])))

    # Create plots
    plots <- lapply(X = features, FUN = function(d) {
      plot <- ST.ImagePlot(data, data.type, group.var, shape.by, d, image, imdims, pt.size, palette,
                           rev.cols, ncol = NULL, spot.colors, center.zero, NULL, split.labels)#, ...)

      if (dark.theme) {
        plot <- plot + dark_theme()
      }
      return(plot)
    })

    if (return.plot.list) {
      return(plots)
    } else {
      plot_grid(plotlist = plots, ncol = grid.ncol)
    }
  }
}


#' Graphs ST spots colored by continuous variable, e.g. dimensional reduction vector
#'
#' @importFrom ggplot2 geom_point aes_string scale_x_continuous scale_y_continuous theme_void theme_void labs scale_color_gradient2 scale_color_gradientn annotation_custom scale_color_manual
#' @importFrom magick image_info
#' @importFrom grid rasterGrob unit
#' @importFrom grDevices as.raster
#'
#' @param image image of class "raster" to use as background for plotting
#' @param dims Dimensions of original image
#'
#' @inheritParams STPlot
#'
#' @export

ST.ImagePlot <- function (
  data,
  data.type,
  group.by,
  shape.by,
  variable,
  image,
  dims,
  pt.size = 2,
  palette = "MaYl",
  rev.cols = F,
  ncol = NULL,
  spot.colors = NULL,
  center.zero = T,
  plot.title = NULL,
  split.labels = FALSE,
  ...
) {

  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }

  if (is.null(spot.colors)) {
    if (class(data[, variable]) == "factor") {
      label.colors <- gg_color_hue(length(levels(data[, variable])))
      names(label.colors) <- levels(data[, variable])
    } else if ((class(data[, variable]) == "character")) {
      label.colors <- gg_color_hue(length(unique(data[, variable])))
      names(label.colors) <- unique(data[, variable])
    }
  }

  # Stop if split.labels is activated and there are more than 1 samples
  if (any(data.type %in% c("character", "factor")) & split.labels) {
    plot.title <- paste0("Sample ", unique(as.character(data[, group.by])), ": ", variable)
    new.data <- data.frame()
    # Order by decreasing size
    if (class(data[, variable]) != "factor") data[, variable] <- factor(data[, variable], levels = unique(data[, variable]))
    levels.keep <- levels(data[, variable])

    for (lbl in levels.keep) {
      dt <- data
      dt[, variable] <- ifelse(dt[, variable] == lbl, lbl, "-")
      dt[, group.by] <- lbl
      new.data <- rbind(new.data, dt)
    }

    data <- new.data
    levels.keep <- c("-", levels.keep)
    data[, variable] <- factor(data[, variable], levels = levels.keep)
    label.colors <- c("-" = "lightgray", label.colors)
  }

  # Obtain colors from selected palette
  cols <- palette.select(palette)(3)
  if (palette == "heat") cols <- palette.select(palette)(4)
  if (rev.cols) {
    cols <- rev(cols)
  }

  # Obtain image dimensions
  x_dim <- as.numeric(dims[2])
  y_dim <- as.numeric(dims[3])

  # Draw image
  g <- rasterGrob(image, width = unit(1, "npc"), height = unit(1, "npc"), interpolate = TRUE)

  # Create new plot
  plots <- lapply(unique(data[, group.by]), function(v) {
    p <- ggplot() +
      annotation_custom(g, -Inf, Inf, -Inf, Inf)

    if (length(spot.colors) > 0) {

      # Add shape aesthetic and blend colors if blend is active
      if (!is.null(shape.by)) {
        p <- p + geom_point(data = data[data[, group.by] == v, ], mapping = aes_string(x = "x", y = paste0(y_dim, " - y"), shape = shape.by), color = spot.colors, size = pt.size)#, ...)
      } else {
        p <- p + geom_point(data = data[data[, group.by] == v, ], mapping = aes_string(x = "x", y = paste0(y_dim, " - y")), color = spot.colors, size = pt.size)#, ...)
      }

    } else {

      # Add shape aesthetic only
      if (!is.null(shape.by)) {
        p <- p + geom_point(data = data[data[, group.by] == v, ], mapping = aes_string(x = "x", y = paste0(y_dim, " - y"), color = paste0("`", variable, "`"), shape = shape.by), size = pt.size)#, ...)
      } else {
        p <- p + geom_point(data = data[data[, group.by] == v, ], mapping = aes_string(x = "x", y = paste0(y_dim, " - y"), color = paste0("`", variable, "`")), size = pt.size)#, ...)
      }
    }

    # Add ST array dimensions
    p <- p +
      scale_x_continuous(limits = c(0, x_dim), expand = c(0, 0)) +
      scale_y_continuous(limits = c(0, y_dim), expand = c(0, 0)) +
      theme_void() +
      labs(title = ifelse(!is.null(plot.title), plot.title, variable), color = "")

    # Center colorscale at 0
    if (center.zero & !any(data.type %in% c("character", "factor"))) {
      p <- p +
        scale_color_gradient2(low = cols[1], mid = cols[2], high = cols[3], midpoint = 0)
    } else if (any(data.type %in% c("character", "factor"))) {
      p <- p +
        labs(color = variable) +
        scale_color_manual(values = label.colors)
    } else {
      p <- p +
        scale_color_gradientn(colours = cols)
    }

    return(p)
  })

  final.plot <- plot_grid(plotlist = plots)
}

#' Apply DimOverlay to multiple samples
#'
#' @param object Seurat object
#' @param sampleids Names of samples to plot
#' @param method Display method
#' @param ncols Number of columns in output image
#' @param ... Parameters passed to DimOverlay
#' @inheritParams DimOverlay
#'
#' @export
#'

MultiDimOverlay <- function(
  object,
  sampleids,
  method = "viewer",
  ncols = NULL,
  dims = c(1:2),
  type = NULL,
  min.cutoff = NA,
  max.cutoff = NA,
  blend = FALSE,
  pt.size = 1,
  reduction = NULL,
  shape.by = NULL,
  palette = NULL,
  rev.cols = FALSE,
  dark.theme = FALSE,
  center.zero = FALSE,
  channels.use = NULL,
  verbose = FALSE,
  ...
) {

  ncols <- ncols %||% round(sqrt(length(x = sampleids)))
  nrows <- round(length(x = sampleids)/ncols)

  p.list <- lapply(sampleids, function(i) {
    DimOverlay(object, dims = dims, sample.index = i, type = type, min.cutoff = min.cutoff, max.cutoff = max.cutoff, blend = blend, pt.size = pt.size, reduction = reduction, shape.by = shape.by, palette = palette, rev.cols = rev.cols, dark.theme = dark.theme, delim = NULL, return.plot.list = FALSE, grid.ncol = NULL, center.zero = center.zero, channels.use = channels.use, verbose = verbose, ... = ...)
  })

  tmp.file <- tempfile(pattern = "", fileext = ".png")

  colf <- ceiling(sqrt(length(x = dims)))
  colr <- round(length(x = dims)/colf)

  png(width = object@tools$xdim*ncols*colf, height = object@tools$xdim*nrows*colr, file = tmp.file)
  par(mar = c(0, 0, 0, 0))
  plot(cowplot::plot_grid(plotlist = p.list, ncol = ncols))
  dev.off()

  p <- image_read(tmp.file)

  if (method == "viewer") {
    print(p)
    unlink(tmp.file)
  } else if (method == "raster") {
    par(mar = c(0, 0, 0, 0))
    plot(as.raster(p))
    unlink(tmp.file)
  } else {
    stop(paste0("Invalid method ", method), call. = FALSE)
  }
}


#' Apply FeatureOverlay to multiple samples
#'
#' @param object Seurat object
#' @param sampleids Names of samples to plot
#' @param method Display method
#' @param ncols Number of columns in output image
#' @param ... Parameters passed to DimOverlay
#' @inheritParams FeatureOverlay
#'
#' @export
#'

MultiFeatureOverlay <- function(
  object,
  sampleids,
  method = "viewer",
  ncols = NULL,
  features,
  type = NULL,
  min.cutoff = NA,
  max.cutoff = NA,
  slot = "data",
  blend = FALSE,
  pt.size = 2,
  shape.by = NULL,
  palette = NULL,
  rev.cols = FALSE,
  dark.theme = FALSE,
  center.zero = FALSE,
  channels.use = NULL,
  verbose = FALSE,
  ...
) {

  ncols <- ncols %||% round(sqrt(length(x = sampleids)))
  nrows <- round(length(x = sampleids)/ncols)

  p.list <- lapply(sampleids, function(s) {
    FeatureOverlay(object, features = features, sample.index = s, type = type,
                   min.cutoff = min.cutoff, max.cutoff = max.cutoff, slot = slot,
                   blend = blend, pt.size = pt.size, shape.by = shape.by,
                   palette = palette, rev.cols = rev.cols, dark.theme = dark.theme,
                   delim = NULL, return.plot.list = FALSE, grid.ncol = NULL,
                   center.zero = center.zero, channels.use = channels.use,
                   verbose = verbose, ... = ...)
  })

  tmp.file <- tempfile(pattern = "", fileext = ".png")

  colf <- ceiling(sqrt(length(x = features)))
  colr <- round(length(x = features)/colf)

  png(width = object@tools$xdim*ncols*colf, height = object@tools$xdim*nrows*colr, file = tmp.file)
  par(mar = c(0, 0, 0, 0))
  plot(cowplot::plot_grid(plotlist = p.list, ncol = ncols))
  dev.off()

  p <- image_read(tmp.file)

  if (method == "viewer") {
    print(p)
    unlink(tmp.file)
  } else if (method == "raster") {
    par(mar = c(0, 0, 0, 0))
    plot(p)
    unlink(tmp.file)
  } else {
    stop(paste0("Invalid method ", method), call. = FALSE)
  }
}



#' Squeeze 2 or 3 column feature data into the unit cube and converts into RGB space
#'
#' @param data data.frame containing feature values and coordinates
#' @param channels.use Select channels to use for blending. Default is red, green and blue but the order can be shuffled.
#' For 2 features, the default is red and green. (options: "red", "green" and "blue")

ColorBlender <- function(
  data,
  channels.use = NULL
) {
  rgb.order <- setNames(1:3, c("red", "green", "blue"))
  if (!length(channels.use) == ncol(data)) {
    stop(paste0("channels.use must be same length as number of features or dimensions"))
  } else if (!all(channels.use %in% names(rgb.order))) {
    stop("Invalid color names in channels.use. Valid options are: 'red', 'green' and 'blue'")
  } else if (sum(duplicated(channels.use))){
    stop("Duplicate color names are not allowed in channels.use")
  }
  col.order <- rgb.order[channels.use]

  if (ncol(data) == 2) {
    data <- cbind(data, rep(0, nrow(data)))
    col.order <- c(col.order, setdiff(1:3, col.order))
    data <- data[, col.order]
  } else if (ncol(data) == 3) {
    data <- data[, col.order]
  }
  color.codes <- rgb(data)
}


#' Find the quantile of a data
#'
#' Converts a quantile in character form to a number regarding some data
#' String form for a quantile is represented as a number prefixed with 'q'
#' For example, 10th quantile is 'q10' while 2nd quantile is 'q2'
#'
#' Will only take a quantile of non-zero data values
#'
#' @param cutoff The cutoff to turn into a quantile
#' @param data The data to find the quantile of
#'
#' @references \url{https://github.com/satijalab/seurat/blob/master/R/visualization.R}
#'
#' @return The numerical representation of the quantile
#'
#' @importFrom stats quantile

SetQuantile <- function(
  cutoff,
  data
) {
  if (grepl(pattern = '^q[0-9]{1,2}$', x = as.character(x = cutoff), perl = TRUE)) {
    this.quantile <- as.numeric(x = sub(
      pattern = 'q',
      replacement = '',
      x = as.character(x = cutoff)
    )) / 100
    data <- unlist(x = data)
    data <- data[data > 0]
    cutoff <- quantile(x = data, probs = this.quantile)
  }
  return(as.numeric(x = cutoff))
}

#' Function used to scale numerical features
#'
#' @param data data.frame containing x, y coordinates and columns with numerical features
#' @param min.cutoff,max.cutoff Vector of minimum and maximum cutoff values for each feature,
#'  may specify quantile in the form of 'q##' where '##' is the quantile (eg, 'q1', 'q10')

feature.scaler <- function (
  data,
  features,
  min.cutoff,
  max.cutoff,
  spots
) {
  min.cutoff <- setNames(mapply(
    FUN = function(cutoff, feature) {
      return(ifelse(
        test = is.na(x = cutoff),
        yes = min(data[, feature]),
        no = cutoff
      ))
    },
    cutoff = min.cutoff,
    feature = features
  ), nm = features)
  max.cutoff <- setNames(mapply(
    FUN = function(cutoff, feature) {
      return(ifelse(
        test = is.na(x = cutoff),
        yes = max(data[, feature]),
        no = cutoff
      ))
    },
    cutoff = max.cutoff,
    feature = features
  ), nm = features)
  check.lengths <- unique(x = vapply(
    X = list(features, min.cutoff, max.cutoff),
    FUN = length,
    FUN.VALUE = numeric(length = 1)
  ))
  if (length(x = check.lengths) != 1) {
    stop("There must be the same number of minimum and maximum cuttoffs as there are features")
  }

  data[, features] <- sapply(X = features,
                                           FUN = function(index) {
                                             data.feature <- as.vector(x = data[, index])
                                             min.use <- SetQuantile(cutoff = min.cutoff[index], data.feature)
                                             max.use <- SetQuantile(cutoff = max.cutoff[index], data.feature)
                                             data.feature[data.feature < min.use] <- min.use
                                             data.feature[data.feature > max.use] <- max.use
                                             return(data.feature)
                                           })
  return(data)
}

#' Generates a dark theme for STPlot
#' @importFrom ggplot2 element_rect element_text theme

dark_theme <- function() {
  theme(plot.background = element_rect(fill = "black"),
        plot.title = element_text(colour = "white"),
        strip.text = element_text(colour = 'white'),
        legend.title = element_text(colour = "white"),
        legend.text = element_text(colour = "white"))
}

#' Store rgb values after ST.Dimplot
#' OBS! could be merged with ST.Dimplot?
#'
#' @param plotObj A ST.DimPlot output object
#'
#' @importFrom plotly ggplotly
#' @importFrom ggplot2 ggplot_build
#'
#' @export

ST.rgbDimPlot <- function(
  plotObj){
  ggBuild <- ggplot_build(gg)
  rgbTransform <- as.data.frame(t(col2rgb(ggBuild$data[[1]]$colour)))
  plotObj$layers[[1]]$mapping$rgbs <- paste("R:", round(rgbTransform$red/2.55),
                                       "% G:", round(rgbTransform$green/2.55),
                                       "% B:", round(rgbTransform$blue/2.55), "%",
                                       sep="")

  plot <- ggplotly(gg , tooltip="rgbs")
  return(plot)

}


#' Select spots by rgb values after ST.DimPlot
#'
#' @param object S4 object
#' @param dimPlot A ST.DimPlot output object
#' @param label1 Optional name of the 1st group of capture-spots
#' @param label2 Optional name of the 2nd group of capture-spots
#'
#' @export

labelRGB <- function(object,
                      dimPlot,
                      label1 = "label1",
                      label2 = NULL,
                      red1=c(0,255),
                      green1=c(0,255),
                      blue1=c(0,255),
                      red2=c(0,255),
                      green2=c(0,255),
                      blue2=c(0,255)){
  # --------------------- Add this part to parameter in ST.Dimplot
  ggBuild <- ggplot_build(dimPlot)
  rgbTransform <- as.data.frame(t(col2rgb(ggBuild$data[[1]]$colour)))
  ggBuild$layers[[1]]$mapping$rgbs <- paste("R:", round(rgbTransform$red/2.55),
                                            "% G:", round(rgbTransform$green/2.55),
                                            "% B:", round(rgbTransform$blue/2.55), "%",
                                            sep="")

  object$rgb.red <- round(rgbTransform$red/2.55)
  object$rgb.green <- round(rgbTransform$green/2.55)
  object$rgb.blue <- round(rgbTransform$blue/2.55)
  # -------------------------------------------------------------
  object$rgbInfo <- 0
  object$rgbInfo[which(object$rgb.red > red1[1] & object$rgb.red < red1[2]
                       & object$rgb.green > green1[1] & object$rgb.green < green1[2]
                       & object$rgb.blue > blue1[1] & object$rgb.blue < blue1[2])] <- label1

  if(!is.null(label2)){
    object$rgbInfo[which(object$rgb.red > red2[1] & object$rgb.red < red2[2]
                         & object$rgb.green > green2[1] & object$rgb.green < green2[2]
                         & object$rgb.blue > blue2[1] & object$rgb.blue < blue2[2])] <- label2
  }

  object$rgbInfo  <- as.factor(object$rgbInfo)

  #seuratobject <- AddMetaData(object, metadata=object$rgbInfo)
  Idents(object = object) <- object$rgbInfo #skriver over?

  return(object)
}


#' Create legend
#'
#' Creates a colorscale legend for a dataset in grob format
#'
#' @param data data.frame with values that should be mapped onto colorscale
#' @param variable character string specifying tha name of the column (variable) containing values
#' @param center.zero Specifies whether or not the colorscale should be centered at zero
#' @param cols Character vector with color ids to be used in colorscale
#' @param val.limits Specifies the limits for values in colorscale
#'
#' @importFrom ggplot2 ggplot geom_point aes_string scale_fill_gradient2 labs scale_fill_gradientn element_rect element_text ggplot_gtable ggplot_build theme

g_legend <- function (
  data,
  variable,
  center.zero,
  cols,
  val.limits
) {
  lg <- ggplot() + geom_point(data = data, aes_string("x", "y", fill = paste0("`", variable, "`")), shape = 22, alpha = 0)
  if (center.zero) {
    lg <- lg +
      scale_fill_gradient2(low = cols[1], mid = cols[median(seq_along(cols))], high = cols[length(cols)], midpoint = 0, na.value = "#000000", limits = val.limits)
  } else if (any(data.type %in% c("character", "factor"))) {
    lg <- lg +
      labs(fill = variable)
  } else {
    lg <- lg +
      scale_fill_gradientn(colours = cols, na.value = "#000000", limits = val.limits)
  }
  lg <- lg + theme(legend.background = element_rect(fill = "#000000"), legend.key = element_rect(fill = "#FFFFFF"), legend.text = element_text(color = "#FFFFFF"), legend.title = element_text(colour = "#FFFFFF"))
  tmp <- ggplot_gtable(ggplot_build(lg))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

