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
#' @param blend Scale and blend expression values to visualize coexpression of two features
#' @param min.cutoff,max.cutoff Vector of minimum and maximum cutoff values for each feature,
#'  may specify quantile in the form of 'q##' where '##' is the quantile (eg, 'q1', 'q10')
#' @param pt.size Adjust point size for plotting
#' @param reduction Which dimensionality reduction to use. If not specified, first searches for umap, then tsne, then pca
#' @param group.by Name of a metadata column to facet plot by (deault is sampleID)
#' @param shape.by If NULL, all points are circles (default). You can specify any
#' cell attribute (that can be pulled with FetchData) allowing for both
#' different colors and different shapes on cells
#' @param delim delimiter passed to \code{\link{GetCoords}} if adjusted ST coordinates are missing in the meta data
#' @param return.plot.list should the plots be returned as a list? By default, the plots are arranged into a grid
#' @param grid.ncol Number of columns for display when combining plots
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
#'
#' @export

ST.DimPlot <- function(
  object,
  dims = c(1, 2),
  spots = NULL,
  blend = FALSE,
  min.cutoff = NA,
  max.cutoff = NA,
  pt.size = 1,
  reduction = NULL,
  group.by = NULL,
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

  # Select colorscale
  palette.info <- palette.select(info = T)
  palette <- palette %||% {
    palette <- subset(palette.info, category == "div")$palette[1]
  }

  if (!blend && length(x = dims) %in% c(2, 3)) {
    stop("Blending feature plots only works with two or three dimensions")
  }

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
    data[, shape.by] <- as.character(object[[shape.by, drop = TRUE]])
    features <- colnames(x = data)[1:(ncol(x = data) - 4)]
  } else {
    features <- colnames(x = data)[1:(ncol(x = data) - 3)]
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

  data <- feature.scaler(data, features, min.cutoff, max.cutoff, spots)

  if (blend) {
    colored.data <- apply(data[, 1:(ncol(data) - 3)], 2, rescale)
    channels.use <- channels.use %||% c("red", "green", "blue")[1:ncol(colored.data)]
    spot.colors <- ColorBlender(data, channels.use)
    data <- data[, (ncol(data) - 2):ncol(data)]
    plot <- STPlot(data,
                   data.type = "numeric",
                   group.by,
                   shape.by,
                   d,
                   pt.size,
                   palette,
                   rev.cols,
                   ncol,
                   spot.colors,
                   center.zero,
                   plot.title = paste(paste(dims, channels.use, sep = ":"), collapse = ", "))
    return(plot)
  } else {
    spot.colors <- NULL
    # Create plots
    plots <- lapply(X = dims, FUN = function(d) {
      plot <- STPlot(data,
                     data.type,
                     group.by,
                     shape.by,
                     d,
                     pt.size,
                     palette,
                     rev.cols,
                     ncol,
                     spot.colors,
                     center.zero)

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
#' @param spots Vector of spots to plot (default is all spots)
#' @param min.cutoff,max.cutoff Vector of minimum and maximum cutoff values for each feature,
#'  may specify quantile in the form of 'q##' where '##' is the quantile (eg, 'q1', 'q10')
#' @param slot Which slot to pull expression data from?
#' @param blend Scale and blend expression values to visualize coexpression of two features
#' @param pt.size Adjust point size for plotting
#' @param group.by Name of a metadata column to facet plot by (deault is sampleID)
#' @param shape.by If NULL, all points are circles (default). You can specify any
#' cell attribute (that can be pulled with FetchData) allowing for both
#' different colors and different shapes on cells
#' @param delim delimiter passed to \code{\link{GetCoords}} if adjusted ST coordinates are missing in the meta data
#' @param return.plot.list should the plots be returned as a list? By default, the plots are arranged into a grid
#' @param grid.ncol Number of columns for display when combining plots
#' @param ... Extra parameters passed on to \code{\link{STPlot}}
#'
#' @inheritParams STPlot
#' @importFrom cowplot plot_grid
#' @importFrom scales rescale
#' @importFrom ggplot2 ggplot theme
#'
#' @return A ggplot object
#' @export

ST.FeaturePlot <- function(
  object,
  features,
  spots = NULL,
  min.cutoff = NA,
  max.cutoff = NA,
  slot = "scale.data",
  blend = FALSE,
  pt.size = 1,
  group.by = NULL,
  shape.by = NULL,
  palette = NULL,
  rev.cols = F,
  dark.theme = F,
  ncol = NULL,
  delim = NULL,
  return.plot.list = F,
  grid.ncol = NULL,
  center.zero = F,
  channels.use = NULL,
  ...
) {
  spots <- spots %||% colnames(x = object)
  data <- FetchData(object = object, vars = c(features), cells = spots, slot = slot)
  data <- as.data.frame(lapply(data, function(x) {
    new_x <- ifelse(test = sapply(x, function(n) {class(n) == "factor"}), yes = as.character(x), no = x)
    return(new_x)
  }))

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
    stopifnot(shape.by %in% colnames(object[[]]))
    data[, shape.by] <- object[[shape.by, drop = TRUE]]
  }

  # Obtain array coordinates
  if (all(c("adj_x", "adj_y") %in% colnames(object[[]]))) {
    data <- cbind(data, setNames(object[[c("adj_x", "adj_y")]], nm = c("x", "y")))
  } else {
    if(is.null(delim)) {
      stop("adjusted coordinates are not present in meta data and delimiter is missing ...", call. = FALSE)
    }
    coords <- GetCoords(colnames(object), delim)
    data <- cbind(data, coords[, c("x", "y")])
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
    colored.data <- apply(data[, 1:(ncol(data) - 3)], 2, rescale)
    channels.use <- channels.use %||% c("red", "green", "blue")[1:ncol(colored.data)]
    spot.colors <- ColorBlender(data, channels.use)
    data <- data[, (ncol(data) - 2):ncol(data)]
    plot <- STPlot(data,
                   data.type,
                   group.by,
                   shape.by,
                   d,
                   pt.size,
                   palette,
                   rev.cols,
                   ncol,
                   spot.colors,
                   center.zero,
                   plot.title = paste(paste(features, channels.use, sep = ":"), collapse = ", "))
    return(plot)
  } else {
    spot.colors <- NULL
    # Create plots
    plots <- lapply(X = features, FUN = function(d) {
      plot <- STPlot(data,
                     data.type,
                     group.by,
                     shape.by,
                     d,
                     pt.size,
                     palette,
                     rev.cols,
                     ncol,
                     spot.colors,
                     center.zero)

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

  if (!length(channels.use) == ncol(colored.data)) {
    stop(paste0("channels.use must be same length as number of features"))
  } else if (!all(channels.use %in% names(rgb.order))) {
    stop("Invalid color names")
  } else if (sum(duplicated(channels.use))){
    stop("Duplicate color names not allowed")
  }
  col.order <- rgb.order[channels.use]

  if (ncol(colored.data) == 2) {
    colored.data <- cbind(colored.data, rep(0, nrow(colored.data)))
    col.order <- c(col.order, setdiff(1:3, col.order))
    colored.data <- colored.data[, col.order]
  } else if (ncol(colored.data) == 3) {
    colored.data <- colored.data[, col.order]
  }
  color.codes <- rgb(colored.data)
}


#' Graphs ST spots colored by continuous variable, e.g. dimensional reduction vector
#'
#' @param data data.frame containing (x, y) coordinates, a group vector and a continuous variable vector
#' @param data.type type of data, e.g. numeric or integer
#' @param group.by specifies column to facet the plots by, e.g. sampleID
#' @param shape.by specifies column to shape points by, e.g. morphological region
#' @param variable name of continuous variable
#' @param pt.size point size of each ST spot
#' @param palette color palette used for spatial heatmap
#' @param rev.cols logical specifying whether colorscale should be reversed
#' @param ncol number of columns in \code{facet_wrap}
#' @param spot.colors character vector woth color names that overrides default coloring with ggplot2
#' @param center.zero should the colorscale be centered around 0? Set to TRUE for scaled data
#' @param ... parameters passed to geom_point()
#'
#' @importFrom ggplot2 geom_point aes_string scale_x_continuous scale_y_continuous theme_void theme_void labs scale_color_gradient2 scale_color_gradientn
#'
#' @export

STPlot <- function(
  data,
  data.type,
  group.by,
  shape.by,
  variable,
  pt.size = 1,
  palette = "MaYl",
  rev.cols = F,
  ncol = NULL,
  spot.colors = NULL,
  center.zero = T,
  plot.title = NULL,
  ...
) {

  cols <- palette.select(palette)(3)
  if (rev.cols) {
    cols <- rev(cols)
  }

  p <- ggplot()
  if (length(spot.colors) > 0) {
    if (!is.null(shape.by)) {
      p <- p + geom_point(data = data, mapping = aes_string(x = "x", y = "64 - y", shape = shape.by), color = spot.colors, size = pt.size)
    } else {
      p <- p + geom_point(data = data, mapping = aes_string(x = "x", y = "64 - y"), color = spot.colors, size = pt.size)
    }
  } else {
    if (!is.null(shape.by)) {
      p <- p + geom_point(data = data, mapping = aes_string(x = "x", y = "64 - y", color = variable, shape = shape.by), size = pt.size)
    } else {
      p <- p + geom_point(data = data, mapping = aes_string(x = "x", y = "64 - y", color = variable), size = pt.size)
    }
  }

  p <- p +
      scale_x_continuous(limits = c(0, 67)) +
      scale_y_continuous(limits = c(0, 64)) +
      theme_void() +
      facet_wrap(as.formula(paste("~", group.by)), ncol = ncol) +
      labs(title = ifelse(!is.null(plot.title), plot.title, variable), color = "")
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

feature.scaler <- function(
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

dark_theme <- function() {
  theme(plot.background = element_rect(fill = "black"),
        plot.title = element_text(colour = "white"),
        strip.text = element_text(colour = 'white'),
        legend.title = element_text(colour = "white"),
        legend.text = element_text(colour = "white"))
}

#' QC histograms on the number of unqiue transcripts and genes per spot
#'
#' @param se S4 object
#' @export

qc.hist <- function(se, metric="transcripts"){

  #OBS OBS ADD NON.SEURAT ALTERNATIVE
  par(mfrow=c(1,2))
  hist(cm$nCount_RNA, main="Unique transcripts per spot", xlab="nr of transcripts", breaks="FD")
  hist(cm$nFeature_RNA, main="Unique genes per spot", xlab="nr of genes", breaks="FD")

}

#' Correlation plots between all samples
#'
#' @param se S4 object
#' @export

qc.corr <- function(se, metric="transcripts"){

  #OBS OBS ADD NON.SEURAT ALTERNATIVE


}
