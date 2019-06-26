#' Function used to read HE images in jpeg format
#'
#' @param object Seurat object
#' @param image.paths Paths to HE images. This is only required if image paths are missing in the Seurat object.
#' @param labels Optional sample labels
#' @param xdim Sets the pixel width for scaling, e.g. "400"
#'
#' @importFrom magick image_read
#'
#' @export

LoadImages <- function(
  object,
  image.paths = NULL,
  labels = NULL,
  xdim = "200"
) {

  if (!is.null(image.paths)) {
    if (!length(x = unique(object[["sample"]]) == length(x = image.paths))) {
      stop(paste0("Number of images (",
                  length(image.paths),
                  ") does not match the number of samples (",
                  length(x = unique(object[["sample"]]), ")")))
    } else {
      object@tools <- list(imgs = image.paths)
    }
  }

  if (!"imgs" %in% names(object@tools)) {
    stop(paste0("Image paths are not present in Seurat object. Provide image.paths manually."))
  }

  imgs <- c()
  for (path in object@tools$imgs) {
    imgs <- c(imgs, image_read(path))
  }

  if (!is.null(labels)) {
    if (length(labels) != length(imgs)) {
      stop("Number of labels does not match the number of images")
    }
    imgs <- setNames(lapply(seq_along(imgs), function(i) {
      image_annotate(image_scale(imgs[[i]], xdim), text = labels[i], size = round(as.numeric(xdim)/20))
    }), nm = labels)
  } else {
    imgs <- lapply(seq_along(imgs), function(i) {
      image_scale(imgs[[i]], xdim)
    })
  }

  object@tools$pointers <- imgs
  object@tools$rasters <- lapply(imgs, as.raster)

  return(object)
}

#' Function used to plot HE images obtained with \code{\link{LoadImages}}
#'
#' @param object Seurat object
#' @param index Image index
#' @inheritParams LoadImages
#'
#' @importFrom magick image_append image_annotate image_scale
#'
#' @export

ImagePlot <- function(
  object,
  index = NULL
) {

  if (!"imgs" %in% names(object@tools)) {
    stop(paste0("Image paths are not present in Seurat object. Run LoadImages before plotting."))
  }

  images <- object@tools$pointers
  check_pointers <- any(sapply(lapply(images, function(images) {
    try(image_info(images))
  }), class) == "try-error")

  if (check_pointers) {
    warning(paste0("Image pointer is dead. You cannot save or cache image objects between R sessions. \n",
                   "Rerun ImageRead if you want to set the image sizes manually. \n",
                   "Setting image size to '200' pixels "), call. = F)
    images <- c()
    for (path in object@tools$imgs) {
      images <- c(images, image_read(path))
    }
    if (!is.null(names(images))) {
      images <- setNames(lapply(seq_along(images), function(i) {
        image_annotate(image_scale(images[[i]], xdim), text = labels[i], size = round(as.numeric(xdim)/20))
      }), nm = names(images))
    } else {
      images <- lapply(seq_along(imgs), function(i) {
        image_scale(images[[i]], xdim)
      })
    }
    object@tools$pointers <- images
  }

  if (is.null(index)) {
    ncols <- round(sqrt(length(x = images)))
    nrows <- round(length(x = images)/ncols)
    stack <- c()
    for (i in 1:nrows) {
      i <- i - 1
      stack <- c(stack, image_append(Reduce(c, images[(i*ncols + 1):(i*ncols + ncols)])))
    }
    print(image_append(Reduce(c, stack), stack = T))
  } else {
    if (!index %in% 1:length(images)) stop("Image index out of bounds", call. = F)
    print(images[[index]])
  }

  return(object)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Feature plots on HE images
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#'  Visualize 'features' on an ST array grid overlayed on top of HE image
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
#' @param image HE image of class "magick-image" (see \code{\link{ImageRead}})
#' @param spots Vector of spots to plot (default is all spots)
#' @param min.cutoff,max.cutoff Vector of minimum and maximum cutoff values for each feature,
#'  may specify quantile in the form of 'q##' where '##' is the quantile (eg, 'q1', 'q10')
#' @param slot Which slot to pull expression data from?
#' @param blend Scale and blend expression values to visualize coexpression of two features (This options will override other coloring parameters)
#' @param pt.size Adjust point size for plotting
#' @param shape.by If NULL, all points are circles (default). You can specify any
#' cell attribute (that can be pulled with FetchData) allowing for both
#' different colors and different shapes on cells
#' @param return.plot.list should the plots be returned as a list? By default, the plots are arranged into a grid
#' @param grid.ncol Number of columns for display when combining plots
#' @param ... Extra parameters passed on to \code{\link{ST.ImagePlot}}
#'
#' @inheritParams ST.ImagePlot
#' @importFrom cowplot plot_grid
#' @importFrom scales rescale
#'
#' @return A ggplot object
#' @export

ST.ImageFeaturePlot <- function(
  object,
  features,
  image,
  spots,
  min.cutoff = NA,
  max.cutoff = NA,
  slot = "scale.data",
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
  ...
) {
  data <- FetchData(object = object, vars = c(features), cells = spots, slot = slot)
  data <- as.data.frame(lapply(data, function(x) {
    new_x <- ifelse(test = sapply(x, function(n) {class(n) == "factor"}), yes = as.character(x), no = x)
    return(new_x)
  }))

  data.type <- unique(sapply(data, class))

  if (!blend && length(x = features) %in% c(2, 3) & !all(data.type %in% c("numeric", "integer"))) {
    stop("Blending feature plots only works with two or three numeric features")
  }

  # Select colorscale
  palette.info <- palette.select(info = T)
  palette <- palette %||% {
    palette <- subset(palette.info, category == "seq")$palette[1]
  }

  # Obtain array coordinates
  if (all(c("pixel_x", "pixel_y") %in% colnames(object[[]]))) {
    data <- cbind(data, setNames(object[[c("pixel_x", "pixel_y")]][spots, ], nm = c("x", "y")))
  } else {
    stop("pixel coordinates are not present in meta data.", call. = FALSE)
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

  if (blend) {
    colored.data <- apply(data[, 1:(ncol(data) - 2)], 2, rescale)
    channels.use <- channels.use %||% c("red", "green", "blue")[1:ncol(colored.data)]
    spot.colors <- ColorBlender(colored.data, channels.use)
    data <- data[, (ncol(data) - 1):ncol(data)]
    plot <- ST.ImagePlot(data, data.type, shape.by, variable, image, pt.size, palette,
      rev.cols, ncol, spot.colors, center.zero,
      plot.title = paste(paste(features, channels.use, sep = ":"), collapse = ", "), ...)
    if (dark.theme) {
      plot <- plot + dark_theme()
    }
    return(plot)
  } else {
    spot.colors <- NULL
    # Create plots
    plots <- lapply(X = features, FUN = function(d) {
      plot <- ST.ImagePlot(data, data.type, shape.by, d, image, pt.size, palette,
                           rev.cols, ncol, spot.colors, center.zero,
                           plot.title = paste(paste(features, channels.use, sep = ":"), collapse = ", "), ...)

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
#' @param data data.frame containing (x, y) coordinates, a group vector and a continuous variable vector
#' @param data.type type of data, e.g. numeric or integer
#' @param shape.by specifies column to shape points by, e.g. morphological region
#' @param variable name of continuous variable
#' @param pt.size point size of each ST spot
#' @param palette color palette used for spatial heatmap
#' @param rev.cols logical specifying whether colorscale should be reversed
#' @param ncol number of columns in \code{facet_wrap}
#' @param spot.colors character vector woth color names that overrides default coloring with ggplot2
#' @param center.zero should the colorscale be centered around 0? Set to TRUE for scaled data
#' @param plot.title  Add title to plot
#' @param ... parameters passed to geom_point()
#'
#' @importFrom ggplot2 geom_point aes_string scale_x_continuous scale_y_continuous theme_void theme_void labs scale_color_gradient2 scale_color_gradientn annotation_custom
#' @importFrom magick image_info
#' @importFrom grid rasterGrob
#' @importFrom grDevices as.raster
#'
#' @export

ST.ImagePlot <- function(
  data,
  data.type,
  shape.by,
  variable,
  image,
  pt.size = 1,
  palette = "MaYl",
  rev.cols = F,
  ncol = NULL,
  spot.colors = NULL,
  center.zero = T,
  plot.title = NULL,
  ...
) {
  # Obtain colors from selected palette
  cols <- palette.select(palette)(3)
  if (rev.cols) {
    cols <- rev(cols)
  }

  # Obtain image dimensions
  dims <- image_info(image)
  x_dim <- as.numeric(dims[2])
  y_dim <- as.numeric(dims[3])

  # Draw image
  g <- rasterGrob(as.raster(image), width = unit(1, "npc"), height = unit(1, "npc"), interpolate = TRUE)

  # Create new plot
  p <- ggplot() +
    annotation_custom(g, -Inf, Inf, -Inf, Inf)

  if (length(spot.colors) > 0) {

    # Add shape aesthetic and blend colors if blend is active
    if (!is.null(shape.by)) {
      p <- p + geom_point(data = data, mapping = aes_string(x = "x", y = paste0(y_dim, " - y"), shape = shape.by), color = spot.colors, size = pt.size, ...)
    } else {
      p <- p + geom_point(data = data, mapping = aes_string(x = "x", y = paste0(y_dim, " - y")), color = spot.colors, size = pt.size, ...)
    }

  } else {

    # Add shape aesthetic only
    if (!is.null(shape.by)) {
      p <- p + geom_point(data = data, mapping = aes_string(x = "x", y = paste0(y_dim, " - y"), color = variable, shape = shape.by), size = pt.size, ...)
    } else {
      p <- p + geom_point(data = data, mapping = aes_string(x = "x", y = paste0(y_dim, " - y"), color = variable), size = pt.size, ...)
    }
  }

  # Add ST array dimensions
  p <- p +
    scale_x_continuous(limits = c(0, x_dim), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, y_dim), expand = c(0, 0)) +
    theme_void() +
    labs(title = ifelse(!is.null(plot.title), plot.title, variable), color = "")

  # Center colorscale at 0
  if (center.zero) {
    p <- p +
      scale_color_gradient2(low = cols[1], mid = cols[2], high = cols[3], midpoint = 0)
  } else if (data.type %in% c("character", "factor")) {
    p <- p +
      labs(color = variable)
  } else {
    p <- p +
      scale_color_gradientn(colours = cols)
  }
  return(p)
}
