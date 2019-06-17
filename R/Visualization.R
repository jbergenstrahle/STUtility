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
#' @param blend.threshold The color cutoff from weak signal to strong signal; ranges from 0 to 1.
#' @param pt.size Adjust point size for plotting
#' @param group.by Name of a metadata column to facet plot by (deault is sampleID)
#' @param shape.by If NULL, all points are circles (default). You can specify any
#' cell attribute (that can be pulled with FetchData) allowing for both
#' different colors and different shapes on cells
#' @param palette color palette to use for heatmap (only relevant for numeric variables, e.g. gene expression), see \code{\link{palette.select}}
#' @param rev.cols logical specifying whether colorscale should be reversed
#' @param combine Combine plots into a single gg object; note that if TRUE; themeing will not work when plotting multiple features
#' @param ncol Number of columns for display of samples for each dimensionality reduction vector
#' @param delim delimiter passed to \code{\link{GetCoords}} if adjusted ST coordinates are missing in the meta data
#' @param return.plot.list should the plots be returned as a list? By default, the plots are arranged into a grid
#' @param grid.ncol Number of columns for display when combining plots
#' @param ... Extra parameters passed on to \code{\link{STPlot}}
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
  blend.threshold = 0.5,
  pt.size = 1,
  group.by = NULL,
  shape.by = NULL,
  palette = "MaYl",
  rev.cols = F,
  dark.theme = F,
  ncol = NULL,
  delim = NULL,
  return.plot.list = F,
  grid.ncol = NULL,
  ...
) {
  spots <- spots %||% colnames(x = object)
  data <- FetchData(object = object, vars = c(features), cells = spots, slot = slot)
  data <- as.data.frame(lapply(data, function(x) {
    new_x <- ifelse(test = sapply(x, function(n) {class(n) == "factor"}), yes = as.character(x), no = x)
    new_x <- ifelse(test = sapply(x, function(n) {class(n) == "integer"}), yes = as.numeric(x), no = x)
    return(new_x)
  }))

  if (length(unique(sapply(data, class))) > 1) {
    stop("Mixed classes (", paste(unique(sapply(data, class)), collapse = ", "), ") are not allowed in features ... ")
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
    data <- feature.scaler(data, min.cutoff, max.cutoff)
  }

  # Create plots
  plots <- lapply(X = features, FUN = function(d) {
    plot <- STPlot(data,
                   group.by,
                   d,
                   pt.size,
                   palette,
                   rev.cols,
                   ncol)

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
#' Graphs ST spots colored by numeric variable, e.g. dimensional reduction vector or a character
#'
#' @param data data.frame containing (x, y) coordinates, a group vector and feature vectors of class
#' numeric/character
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
=======
  p <- ggplot() +
    geom_point(data = data, mapping = aes_string(x = "x", y = "64 - y", color = variable), size = pt.size, ...) +
    scale_x_continuous(limits = c(0, 67)) +
    scale_y_continuous(limits = c(0, 64)) +
    theme_void() +
    facet_wrap(as.formula(paste("~", group.by)), ncol = ncol)
  if (class(data[, variable, drop = T]) == "numeric") {
    cols <- palette.select(palette)(3)
    if (rev.cols) {
      cols <- rev(cols)
    }
    p <- p +
      labs(title = variable, color = "value") +
      scale_color_gradientn(colours = cols)
  } else if (class(data[, variable, drop = T]) %in% c("character", "factor")) {
    p <- p +
      labs(color = variable)
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
  min.cutoff,
  max.cutoff
) {
  features <- colnames(x = data)[1:(ncol(x = data) - 3)]
  min.cutoff <- mapply(
    FUN = function(cutoff, feature) {
      return(ifelse(
        test = is.na(x = cutoff),
        yes = min(data[, feature]),
        no = cutoff
      ))
    },
    cutoff = min.cutoff,
    feature = features
  )
  max.cutoff <- mapply(
    FUN = function(cutoff, feature) {
      return(ifelse(
        test = is.na(x = cutoff),
        yes = max(data[, feature]),
        no = cutoff
      ))
    },
    cutoff = max.cutoff,
    feature = features
  )
  check.lengths <- unique(x = vapply(
    X = list(features, min.cutoff, max.cutoff),
    FUN = length,
    FUN.VALUE = numeric(length = 1)
  ))
  if (length(x = check.lengths) != 1) {
    stop("There must be the same number of minimum and maximum cuttoffs as there are features")
  }

  data[, 1:(ncol(x = data) - 2)] <- sapply(X = 1:(ncol(x = data) - 2),
                                           FUN = function(index) {
                                             data.feature <- as.vector(x = data[, index])
                                             min.use <- SetQuantile(cutoff = min.cutoff[index], data.feature)
                                             max.use <- SetQuantile(cutoff = max.cutoff[index], data.feature)
                                             data.feature[data.feature < min.use] <- min.use
                                             data.feature[data.feature > max.use] <- max.use
                                             return(data.feature)
                                           })
  colnames(x = data)[1:(ncol(x = data) - 3)] <- features
  rownames(x = data) <- spots
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
