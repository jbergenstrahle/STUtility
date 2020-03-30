#' @include generics.R Staffli.R Visualization_utilities.R
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Dimensional reduction plots
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Dimensional reduction plot on ST coordinates
#'
#' Graphs the selected vectors of a dimensional reduction technique on a 2D grid of spots.
#'
#' @section Blending values:
#' The blend option can be useful if you wish to visualize multiple dimensions simultaneuosly and works for two or three dimensionality
#' reduction vectors. Each of the selected vectors are rescaled from 0 to 1 and are used as RGB color channels to produce mixed color for each
#' spot. This can be particularly useful when looking at overlapping value vectors. For example, if you are looking at two overlapping value vectors
#' "A" and "B" and use the blend option, the "A" values will be encoded in the "red" channel and the "B" values in the "green" channel. If a spot is
#' purely "A" the color will be red and if it is purely "B" it will green. Any mixture of "A" and "B" will produce a color between red and green
#' where a 50/50 mixture gives yellow color. The amplitude if the values will also determine the brightness of the color.
#'
#' @section Colorscale:
#' Note that the dimensionality reduction outputs from method s.a. PCA, ICA and UMAP are typically centered at 0 and it is tehrefore
#' appropriate to use a divergent colorscale for this type of output, e.g. the "RdBu" color palette in RColorBrewer. However,
#' if you are plotting factor from the Non-negative Matrix Factorization (NMF) method, you will only have positive value and the
#' `center.zero` argument should tehrefore be set to FALSE.
#'
#' @param object Seurat object
#' @param dims Dimensions to plot [default: 1, 2]
#' @param spots Character vector with spot IDs to plot [default: all spots]
#' @param indices Indices to subset data by
#' @param plot.type Specify the type of plot to use [default: "spots"]. Available options are; "spots" (a "smooth" options will be added soon)
#' @param blend Scale and blend expression values to visualize coexpression of two features (this options will override other coloring parameters).
#' See 'Blending values' below for a more thourough description.
#' @param min.cutoff,max.cutoff Vector of minimum and maximum cutoff values for each feature, may specify quantile in the form of 'q##' where '##'
#' is the quantile (eg, 'q1', 'q10'). This can be useful if you have outlier values that skew the colorscale in the plot. For example, if you specify
#' 'q1', you will trim of values below the 1st percentile. [default: no cuttoffs]
#' @param pt.size Adjust point size for plotting [default: 1]
#' @param pt.alpha Adjust point opacity for plotting [default: 1]
#' @param pt.border Should a border be drawn around the spots? [default: TRUE]
#' @param reduction Which dimensionality reduction to use. If not specified, first searches for "umap", then "tsne", then "pca"
#' @param shape.by You can specify any spot attribute (that can be pulled with FetchData) allowing for both different colors
#' and different shapes on spots
#' @param grid.ncol Number of columns for display when combining plots. This option will only have an effect on the sample level structure.
#' @param channels.use Color channels to use for blending. Has to be a character vector of length 2 or 3 with "red", "green" and "blue"
#' color names specified [default: c("red", "green", "blue)]
#' @param sb.size Size of scalebar [default: 2.5]
#' @param show.sb Should a scalebar be drawn? [default: TRUE]
#' @param value.scale Defines how the dimensionality reduction values should be mapped to the colorbar. `value.scale = "samplewise"` will
#' scale each feature independantly while `value.scale = "all"` will use the same value range for all vectors
#' @param verbose Print messages
#'
#' @param ... Extra parameters passed on to \code{\link{t}}
#'
#' @inheritParams draw_scalebar
#' @inheritParams STPlot
#'
#' @return A ggplot object
#'
#' @importFrom rlang !!
#' @importFrom ggplot2 ggplot facet_wrap vars sym
#' @importFrom viridis magma
#' @importFrom Seurat FetchData Embeddings
#' @importFrom zeallot %<-%
#'
#' @examples
#'
#' se <- RunPCA(se)
#'
#' # Plot the first 5 dimensions
#' ST.DimPlot(se, dims = 1:5, reduction = "pca")
#'
#' # Plot the first 5 dimensions and put each sample in a separate row
#' ST.DimPlot(se, dims = 1:5, reduction = "pca", grid.ncol = 1)
#'
#' # Blend values for dimensions 1 and 2
#' ST.DimPlot(se, dims = 1:2, reduction = "pca", blend = T)
#'
#' # Plot the first 5 dimensions and trim off 1st percentile values
#' ST.DimPlot(se, dims = 1:5, reduction = "pca", min.cutoff = 'q1')
#'
#' @export
#'
#' @seealso \code{\link{ST.FeaturePlot}} for how to plot features,
#' \code{\link{FeatureOverlay}} and \code{\link{DimOverlay}} for how to overlay plots on the
#' HE images.
#'

ST.DimPlot <- function (
  object,
  dims = c(1, 2),
  spots = NULL,
  indices = NULL,
  plot.type = "spots",
  blend = FALSE,
  min.cutoff = NA,
  max.cutoff = NA,
  pt.size = 1,
  pt.alpha = 1,
  pt.border = FALSE,
  reduction = NULL,
  shape.by = NULL,
  palette = "MaYl",
  cols = NULL,
  dark.theme = FALSE,
  ncol = NULL,
  grid.ncol = NULL,
  center.zero = TRUE,
  channels.use = NULL,
  center.tissue = FALSE,
  verbose = FALSE,
  sb.size = 2.5,
  show.sb = TRUE,
  value.scale = c("samplewise", "all"),
  custom.theme = NULL,
  ...
) {

  # Check to see if Staffli object is present
  if (!"Staffli" %in% names(object@tools)) stop("Staffli object is missing from Seurat object. Cannot plot without coordinates", call. = FALSE)
  st.object <- object@tools$Staffli

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

  # Add group column to data
  data[,  "sample"] <- st.object[[spots, "sample", drop = TRUE]]

  # Extract shape.by column from meta data if applicable
  if (!is.null(x = shape.by)) {
    if (!shape.by %in% colnames(object[[]])) {
      stop(paste0("Shaping variable (shape.by) ", shape.by, " not found in meta.data slot"), call. = F)
    }
    data[, shape.by] <- as.character(object[[shape.by, drop = TRUE]])
  }

  # Obtain array coordinates
  image.type <- "empty"
  c(data, image.type) %<-% obtain.array.coords(st.object, data, image.type, spots)

  # Scale data values
  data <- feature.scaler(data, dims, min.cutoff, max.cutoff)

  # Subset by index
  if (!is.null(indices)) {
    if (!all(as.character(indices) %in% data[, "sample"])) stop(paste0("Index out of range. "), call. = FALSE)
    data <- data[data[, "sample"] %in% as.character(indices), ]
  } else {
    indices <- unique(data[,  "sample"]) %>% as.numeric()
  }

  # Fetch dims
  if (image.type != "empty") {
    dims.list <- lapply(st.object@dims, function(x) {x[2:3] %>% as.numeric()})
  } else {
    dims.list <- st.object@limits
  }

  # Subset dims by indices
  if (!is.null(indices)) dims.list <- dims.list[indices]

  # Prepare data for scalebar
  # --------------------------------------------------------------------
  pxum <- prep.sb(st.object, data, data.type, indices, FALSE, dims, dims.list, show.sb)
  # --------------------------------------------------------------------

  # Set feature scale limits
  value.scale <- match.arg(value.scale, c("samplewise", "all"))
  if (value.scale == "all") {
    limits <- c(min(data[, dims]), max(data[, dims]))
  } else if (value.scale == "samplewise") {
    limits <- NULL
  }

  # blend colors or plot each dimension separately
  if (blend) {
    colored.data <- apply(data[, 1:(ncol(data) - ifelse(!is.null(shape.by), 4, 3))], 2, scales::rescale)
    channels.use <- channels.use %||% c("red", "green", "blue")[1:ncol(colored.data)]

    if (verbose) cat(paste0("Blending colors for dimensions ",
                            paste0(ifelse(length(dims) == 2, paste0(dims[1], " and ", dims[2]), paste0(dims[1], dims[2], " and ", dims[2]))),
                            ": \n", paste(paste(dims, channels.use, sep = ":"), collapse = "\n")))

    spot.colors <- ColorBlender(colored.data, channels.use)
    data <- data[, (ncol(data) - ifelse(!is.null(shape.by), 3, 2)):ncol(data)]
    p <- STPlot(data, data.type = "numeric", shape.by, NULL, pt.size, pt.alpha, pt.border,
                   palette, cols, ncol, spot.colors, center.zero, center.tissue,
                   plot.title = paste(paste(dims, channels.use, sep = ":"), collapse = ", "),
                   dims.list, split.labels = FALSE, pxum = pxum, sb.size = sb.size,
                   dark.theme, NULL, custom.theme, ...)
    if (dark.theme) {
      p <- p + dark_theme()
    }
    return(p)
  } else {
    if (plot.type == "spots") {
      spot.colors <- NULL
      if (verbose) cat("Plotting dimensions:",
                       ifelse(length(dims) == 1, dims,  paste0(paste(dims[1:(length(dims) - 1)], collapse = ", "), " and ", dims[length(dims)])))

      plots <- lapply(X = dims, FUN = function(d) {
        plot <- STPlot(data, data.type = "numeric", shape.by, d, pt.size, pt.alpha, pt.border,
                       palette, cols, ncol, spot.colors, center.zero, center.tissue,
                       d, dims.list, FALSE, pxum = pxum, sb.size = sb.size,
                       dark.theme, limits, custom.theme, ...)
        if (dark.theme) {
          plot <- plot + dark_theme()
        }
        return(plot)
      })

      # Draw plots
      p <- plot_grid(plotlist = plots, ncol = grid.ncol)
      if (dark.theme) {
        p <- p + dark_theme()
      }
      return(p)
    } else if (plot.type == "smooth") {
      stop("Smooth options not yet implemented for dimred output ...")
      plots <- lapply(X = dims, FUN = function(d) {
        plot <- SmoothPlot(st.object, data, image.type, data.type = "numeric", d,
                           palette, cols, ncol, center.zero, dark.theme, highlight.edges, ...)
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

#' Feature plot on ST coordinates
#'
#' Graphs selected features on a 2D grid of spots, for example raw gene counts, normalized gene counts or cluster labels.
#'
#' @details Note that you can only graph one specific class of features at the same time. This function does not support mixing of
#' numeric and character vectors for example.
#'
#' @section Splitting categorical features:
#' If you are plotting a categorical feature, e.g.cluster labels, you have the option to split each label into facets using \code{split.labels=TRUE}.
#' This is very useful if you have many different labels which can make it difficult to distinguish the different colors. However, splitting only
#' works for one sample at the time which has to be specified by the \code{indices} argument.
#'
#' @section Blending values:
#' The blend option can be useful if you wish to visualize multiple features simultaneuosly and works for two or three vectors. Each of the
#' selected vectors are rescaled from 0 to 1 and are used as RGB color channels to produce mixed color for each
#' spot. This can be particularly useful when looking at overlapping value vectors. For example, if you are looking at two overlapping features
#' "A" and "B" and use the blend option, the "A" values will be encoded in the "red" channel and the "B" values in the "green" channel. If a spot is
#' purely "A" the color will be red and if it is purely "B" it will green. Any mixture of "A" and "B" will produce a color between red and green
#' where a 50/50 mixture gives yellow color. The amplitude if the values will also determine the brightness of the color.
#'
#' @param object Seurat object
#' @param features
#' \itemize{
#'     \item An \code{Assay} feature (e.g. a gene name - "MS4A1")
#'     \item A column name from meta.data (e.g. mitochondrial percentage - "percent.mito")
#' }
#' @param indices Integer vector specifying sample indices to include in the plot [default: show all samples]
#' @param spots Character vector with spot IDs to plot [default: all spots]
#' @param plot.type Specify the type of plot to use [default: "spots"]. Available options are; "spots" (a "smooth" options will be added soon)
#' @param blend Scale and blend expression values to visualize coexpression of two features (this options will override other coloring parameters).
#' See 'Blending values' below for a more thourough description.
#' @param min.cutoff,max.cutoff Vector of minimum and maximum cutoff values for each feature, may specify quantile in the form of 'q##' where '##'
#' is the quantile (eg, 'q1', 'q10'). This can be useful if you have outlier values that skew the colorscale in the plot. For example, if you specify
#' 'q1', you will trim of values below the 1st percentile. [default: no cuttoffs]
#' @param slot Which slot to pull the data from? [default: 'data']
#' @param pt.size Adjust point size for plotting [default: 1]
#' @param pt.alpha Adjust point opacity for plotting [default: 1]
#' @param shape.by You can specify any spot attribute (that can be pulled with FetchData) allowing for both different colors
#' and different shapes on spots
#' @param grid.ncol Number of columns for display when combining plots. This option will only have an effect on the sample level structure.
#' @param channels.use Color channels to use for blending. Has to be a character vector of length 2 or 3 with "red", "green" and "blue"
#' color names specified [default: c("red", "green", "blue)]
#' @param sb.size Size of scalebar [default: 2.5]
#' @param show.sb Should a scalebar be drawn? [default: TRUE]
#' @param value.scale Defines how the dimensionality reduction values should be mapped to the colorbar. `value.scale = "samplewise"` will
#' scale each feature independantly while `value.scale = "all"` will use the same value range for all vectors
#' @param verbose Print messages
#'
#' @param ... Extra parameters passed on to \code{\link{t}}
#'
#' @inheritParams STPlot
#' @inheritParams SmoothPlot
#' @importFrom cowplot plot_grid
#' @importFrom ggplot2 ggplot theme theme_void
#' @importFrom zeallot %<-%
#'
#' @return A ggplot object
#'
#' @examples
#'
#' # Plot the number of unique genes and the number of UMIs per spot
#' ST.featurePlot(se, features = c("nFeature_RNA", "nCount_RNA"))
#'
#' # Plot selected genes
#' ST.featurePlot(se, features = c("Cck", "Dcn"))
#'
#' # Plot normalized values
#' se <- SCTransform(se)
#' ST.featurePlot(se, features = c("Cck", "Dcn"))
#'
#' # Change to scaled data
#' ST.featurePlot(se, features = c("Cck", "Dcn"), slot = "scale.data", center.zero = TRUE)
#'
#' # Cluster spots and plot cluster labels
#' se <- se %>% RunPCA()
#'    FindNeighbors(dims = 1:10, reduction = "pca") %>%
#'    FindClusters()
#' ST.featurePlot(se, features = "seurat_clusters")
#' # Split cluster labels into facets
#' ST.featurePlot(se, features = "seurat_clusters", split.labels = TRUE)
#'
#' @export
#'
#' @seealso \code{\link{ST.DimPlot}} for how to plot dimensionality reduction output,
#' \code{\link{FeatureOverlay}} and \code{\link{DimOverlay}} for how to overlay plots on the
#' HE images.
#'

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
  pt.alpha = 1,
  pt.border = FALSE,
  shape.by = NULL,
  palette = NULL,
  cols = NULL,
  dark.theme = FALSE,
  highlight.edges = TRUE,
  ncol = NULL,
  grid.ncol = NULL,
  center.zero = FALSE,
  channels.use = NULL,
  center.tissue = FALSE,
  verbose = FALSE,
  sb.size = 2.5,
  show.sb = TRUE,
  value.scale = c("samplewise", "all"),
  custom.theme = NULL,
  ...
) {
  # Check to see if Staffli object is present
  if (!"Staffli" %in% names(object@tools)) stop("Staffli object is missing from Seurat object. Cannot plot without coordinates", call. = FALSE)
  st.object <- object@tools$Staffli

  # Collect data
  spots <- spots %||% colnames(x = object)
  data <- FetchData(object = object, vars = c(features), cells = spots, slot = slot)
  data.type <- unique(sapply(data, class))

  # Stop if feature classes are mixed
  if (length(data.type) > 1 & !all(data.type %in% c("numeric", "integer"))) {
    stop("Mixed classes (", paste(unique(sapply(data, class)), collapse = ", "), ") are not allowed in features ... ")
  }

  # If blend option is set, stop if the number of features are not equal to 2 or 3
  if (blend) {
    if (plot.type == "smooth") stop(paste0("Smoothing option not available when blending ..."), call. = FALSE)
    if (!(length(x = features) %in% c(2, 3)) | !all(data.type %in% c("numeric", "integer"))) {
      stop("Blending feature plots only works with two or three numeric features")
    }
  }

  # Select colorscale
  palette.info <- palette.select(info = T)
  palette <- palette %||% {
    palette <- subset(palette.info, category == "seq")$palette[1]
  }

  # Add group column to data
  data[,  "sample"] <- st.object[[spots, "sample", drop = TRUE]]

  # Add shape column if specified
  if (!is.null(x = shape.by)) {
    if (!shape.by %in% colnames(object[[]])) {
      stop(paste0("Shaping variable (shape.by) ", shape.by, " not found in meta.data slot"), call. = F)
    }
    data[, shape.by] <- as.character(object[[shape.by, drop = TRUE]])
  }

  # Obtain array coordinates
  image.type <- "empty"
  c(data, image.type) %<-% obtain.array.coords(st.object, data, image.type, spots)

  # Raise error if features are not present in Seurat object
  if (ncol(x = data) < 4) {
    stop("None of the requested features were found: ",
         paste(features, collapse = ", "),
         " in slot ",
         slot,
         call. = FALSE)
  }

  if (all(data.type %in% c("numeric", "integer"))) {
    data <- feature.scaler(data, features = features, min.cutoff, max.cutoff)
  }

  # Subset by index
  if (!is.null(indices)) {
    if (!all(as.character(indices) %in% data[, "sample"])) stop(paste0("Index out of range. "), call. = FALSE)
    data <- data[data[, "sample"] %in% as.character(indices), ]
  } else {
    indices <- unique(data[,  "sample"]) %>% as.numeric()
  }

  # Fetch dims
  if (image.type != "empty") {
    dims <- lapply(st.object@dims, function(x) {x[2:3] %>% as.numeric()})
  } else {
    dims <- st.object@limits
  }

  # Subset dims by indices
  if (!is.null(indices)) dims <- dims[indices]

  # Prepare data for scalebar
  # --------------------------------------------------------------------
  pxum <- prep.sb(st.object, data, data.type, indices, split.labels, features, dims, show.sb)
  # --------------------------------------------------------------------

  # Set feature scale limits
  value.scale <- match.arg(value.scale, c("samplewise", "all"))
  if (value.scale == "all" & all(data.type %in% c("numeric", "integer"))) {
    limits <- c(min(data[, features]), max(data[, features]))
  } else if (value.scale == "samplewise") {
    limits <- NULL
  }

  if (blend) {
    colored.data <- apply(data[, 1:(ncol(data) - ifelse(!is.null(shape.by), 4, 3))], 2, scales::rescale)
    channels.use <- channels.use %||% c("red", "green", "blue")[1:ncol(colored.data)]

    if (verbose) cat(paste0("Blending colors for features ",
                            paste0(ifelse(length(features) == 2, paste0(features[1], " and ", features[2]), paste0(features[1], features[2], " and ", features[2]))),
                            ": \n", paste(paste(features, channels.use, sep = ":"), collapse = "\n")))

    spot.colors <- ColorBlender(colored.data, channels.use)
    data <- data[, (ncol(data) - ifelse(!is.null(shape.by), 4, 3)):ncol(data)]
    plot <- STPlot(data, data.type, shape.by, NULL, pt.size, pt.alpha, pt.border,
                   palette, cols, ncol, spot.colors, center.zero, center.tissue,
                   plot.title = paste(paste(features, channels.use, sep = ":"), collapse = ", "),
                   dims, FALSE, pxum, sb.size, dark.theme, NULL, custom.theme, ...)
    if (dark.theme) {
      plot <- plot + dark_theme()
    }
    plot
  } else {
    spot.colors <- NULL
    if (verbose) cat("Plotting features:",
                     ifelse(length(features) == 1, features,  paste0(paste(features[1:(length(features) - 1)], collapse = ", "), " and ", features[length(features)])))

    # Create plots
    if (plot.type == "spots") {
      plots <- lapply(X = features, FUN = function(ftr) {
            plot <- STPlot(data, data.type, shape.by, ftr, pt.size, pt.alpha, pt.border,
                           palette, cols, ncol, spot.colors, center.zero,
                           center.tissue, ftr, dims, split.labels, pxum,
                           sb.size, dark.theme, limits, custom.theme, ...)

            if (dark.theme) {
              plot <- plot + dark_theme()
            }
            return(plot)
          })
      p <- plot_grid(plotlist = plots, ncol = grid.ncol)
      if (dark.theme) {
          p <- p + dark_theme()
      }
      p
    } else if (plot.type == "smooth") {
      if (any(data.type %in% c("factor", "character"))) {
        stop(paste0("Smoothing has not yet been implemented for categorical variables"), call. = FALSE)
      }
      plots <- lapply(X = features, function(ftr) {
        plot <- SmoothPlot(st.object, data, image.type, data.type, ftr,
                           palette, cols, ncol, center.zero, dark.theme,
                           highlight.edges, ...)
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
      par(mar = c(0, 0, 0, 0), bg = ifelse(dark.theme, "black", "white"))
      plot(final_img %>% as.raster())
    }
  }
}


#' Graphs ST spots colored by continuous or categorical features
#'
#' @param data Object of class 'data.frame' containing at least (x, y) coordinates, a "sample" vector with labels for each sample
#' and one column with the feature values. Can also include an additional column for shapes.
#' @param data.type String specifying the class of the features in data to be plotted
#' @param shape.by String specifying the column where the shaping label is stored
#' @param variable Name of feature column
#' @param pt.size Point size of each ST spot [default: 1]
#' @param pt.alpha Opacity of each ST spot [default: 1]
#' @param pt.border Should a border be drawn around the spots? [default: TRUE]
#' @param palette Color palette used for spatial heatmap (see \code{palette.select(info = T)} for available options).
#' Disabled if a color vector is provided (see \code{cols} below).
#' @param cols A vector of colors to use for colorscale, e.g. \code{cols = c("blue", "white", "red")} will
#' create a gradient color scale going from blue to white to red. This options will deactivate the \code{palette}
#' option.
#' @param ncol Number of columns to arrange the samples into. This can for example be useful to adjust if you want to visualize the samples
#' in just in one row or one column.
#' @param spot.colors Character vector with color names that overrides default coloring with ggplot2
#' @param center.zero Specifies whther or not the colorscale should be centered around 0. For some values, such as Principal Component vectors,
#' the distribution of values is centered at 0 and in that case it can be appropriate to use a divergent colorscale with a predefined value for 0.
#' If this parameter is set to TRUE, the ggplot2 function \code{scale_color_gradient2} will be used to control the coloring instead of
#' \code{scale_color_gradientn}. If center.zero is set to FALSE, the colorscale will simply map the values in equally spaced intervals which could skew
#' the interpretaion of the output plot.
#' @param center.tissue Adjust coordinates so that the center of the tissue is in the middle of the array along the y-axis. This can be useful if your
#' samples have been placed in very different parts of the capture area and you want to center the plots in the middle. This is however unnecessary if
#' you have already aligned the sample data (see \code{\link{AlignImages}}, \code{\link{WarpImages}} and \code{\link{ManualAlignImages}})
#' @param plot.title String specifying the title of the plot
#' @param dims List of dimensions for x and y scales. If you have mixed datasets from different arrays (platforms) with different resolution,
#' this list of dimensions will be used to specify the limits along the x- and y-axis of the array for each sample.
#' @param split.labels Only works if the features are specified by character vectors.
#' The plot will be split into one plot for each group label, highlighting the labelled spots.
#' @param custom.theme Object of class 'theme' used to change the background theme (see for example \url{https://ggplot2.tidyverse.org/reference/theme.html})
#' @param sb.size Defines the size of the scalebar
#' @param limits Sets the range of the colorbar values
#' @param ... Parameters passed to geom_point()
#'
#' @inheritParams draw_scalebar
#'
#' @importFrom ggplot2 geom_point aes_string scale_x_continuous scale_y_continuous theme_void theme_void labs scale_color_gradient2 scale_color_gradientn scale_color_manual
#'
#' @export

STPlot <- function (
  data,
  data.type = NULL,
  shape.by,
  variable,
  pt.size = 1,
  pt.alpha = 1,
  pt.border = FALSE,
  palette = "MaYl",
  cols = NULL,
  ncol = NULL,
  spot.colors = NULL,
  center.zero = TRUE,
  center.tissue = F,
  plot.title = NULL,
  dims = NULL,
  split.labels = FALSE,
  pxum = NULL,
  sb.size = 2.5,
  dark.theme = FALSE,
  limits = NULL,
  custom.theme = NULL,
  ...
) {

  # Remove NA values from data
  data <- na.omit(data)

  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }

  if (is.null(spot.colors)) {
    if (class(data[, variable]) == "factor") {
      if (!is.null(cols)) {
        stopifnot(length(cols) >= length(unique(data[, variable])))
        if (is.null(names(cols))) {
          label.colors <- cols[1:length(unique(data[, variable]))]
        } else {
          label.colors <- cols[levels(data[, variable])]
        }
      } else {
        label.colors <- gg_color_hue(length(levels(data[, variable])))
      }
      names(label.colors) <- levels(data[, variable])
    } else if (class(data[, variable]) == "character") {
      if (!is.null(cols)) {
        stopifnot(length(cols) >= length(unique(data[, variable])))
        if (is.null(names(cols))) {
          label.colors <- cols[1:length(unique(data[, variable]))]
        } else {
          label.colors <- cols[unique(data[, variable])]
        }
      } else {
        label.colors <- gg_color_hue(length(unique(data[, variable])))
      }
      names(label.colors) <- unique(data[, variable])
    }
  }

  # Stop if split.labels is activated and there are more than 1 samples
  if (any(data.type %in% c("character", "factor")) & split.labels) {
    if (length(unique(as.character(data[, "sample"]))) > 1) stop(paste0("Splitting of group labels only work for one sample. Please set a single sample index with the 'indices' argument. "), call. = FALSE)
    plot.title <- paste0("Sample ", unique(as.character(data[, "sample"])), ": ", variable)
    new.data <- data.frame()
    # Order by decreasing size
    if (class(data[, variable]) != "factor") data[, variable] <- factor(data[, variable], levels = unique(data[, variable]))
    levels.keep <- levels(data[, variable])

    for (lbl in levels.keep) {
      dt <- data
      dt[, variable] <- ifelse(dt[, variable] == lbl, lbl, "-")
      dt[, "sample"] <- lbl
      new.data <- rbind(new.data, dt)
    }

    data <- new.data
    levels.keep <- c("-", levels.keep)
    data[, variable] <- factor(data[, variable], levels = levels.keep)
    label.colors <- c("-" = "lightgray", label.colors)
  }

  # Center tissue along y-axis
  if (center.tissue) {
    data.split <- split(data, data[, "sample"])
    data <- do.call(rbind, lapply(seq_along(data.split), function(i) {
      d <- data.split[[i]]
      d[, "y"] <- d[, "y"] - median(d[, "y"]) + dims[[i]][2]/2
      return(d)
    }))
  }

  # Obtain colors from selected palette or from provided cols
  cols <- cols %||% {
    palette.select(palette)
  }

  # Define scales for each facet
  if (!split.labels) {
    limits_x <- lapply(seq_along(dims), function(i) {
      d <- dims[[i]]
      scale_override(i, scale_x_continuous(limits = c(0, d[1])))
    })
    limits_y <- lapply(seq_along(dims), function(i) {
      d <- dims[[i]]
      scale_override(i, scale_y_continuous(limits = c(0, d[2])))
    })
    limits_override <- c(limits_x, limits_y)

    # Invert y-axis by sample
    split.data <- lapply(paste0(seq_along(unique(data[, "sample"]))), function(s) {
      data[data[, "sample"] == s, ]
    })
    split.data <- lapply(seq_along(split.data), function(i) {
      spld <- split.data[[i]]
      ydim <- dims[[i]][2]
      spld$y <- ydim - spld$y
      return(spld)
    })
    data <- do.call(rbind, split.data)
  } else {
    lims <- unlist(dims)
    data$y <- lims[2] - data$y
  }

  # Create new plot
  p <- ggplot()
  # Make sure that levels are correct
  if (!split.labels) {
    data[, "sample"] <- factor(data[, "sample"], levels = unique(data[, "sample"]))
  } else {
    data[, "sample"] <- factor(data[, "sample"], levels = levels.keep)
  }
  if (length(spot.colors) > 0) {

    # Add shape aesthetic and blend colors if blend is active
    if (!is.null(shape.by)) {
      p <- p + geom_point(data = data, mapping = aes_string(x = "x", y = "y", shape = shape.by),
                          shape = 21, stroke = ifelse(pt.border, 0.2, 0),
                          fill = spot.colors, size = pt.size, alpha = pt.alpha, ...)
    } else {
      p <- p + geom_point(data = data, mapping = aes_string(x = "x", y = "y"),
                          shape = 21, stroke = ifelse(pt.border, 0.2, 0),
                          fill = spot.colors, size = pt.size, alpha = pt.alpha, ...)
    }

  } else {

    # Add shape aesthetic only
    if (!is.null(shape.by)) {
      p <- p + geom_point(data = data, mapping = aes_string(x = "x", y = "y", fill = paste0("`", variable, "`"), shape = shape.by),
                          shape = 21, stroke = ifelse(pt.border, 0.2, 0),
                          size = pt.size, alpha = pt.alpha, ...)
    } else {
      p <- p + geom_point(data = data, mapping = aes_string(x = "x", y = "y", fill = paste0("`", variable, "`")),
                          shape = 21, stroke = ifelse(pt.border, 0.2, 0),
                          size = pt.size, alpha = pt.alpha, ...)
    }

  }

  # Add ST array dimensions and plot title
  p <- p +
    labs(title = ifelse(!is.null(plot.title), plot.title, ""), fill = ifelse(all(data.type %in% c("numeric", "integer")), "value", "label"))

  # Set theme
  p <- p + theme_void()
  p <- p + theme(plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "npc"))

  if (!is.null(custom.theme)) {
    p <- p + custom.theme
  }

  # Facet plots by group variable
  if (!split.labels) {
    p <- p +
      facet_wrap_custom(~sample, scales = "free", ncol = ncol, scale_overrides = limits_override)
  } else {
    p <- p + facet_wrap(~sample, ncol = ncol) +
      scale_x_continuous(limits = c(0, lims[1])) +
      scale_y_continuous(limits = c(0, lims[2]))
  }

  ## Set the scale bar
  if (!is.null(pxum)) {
    p <- draw_scalebar(p, pxum = pxum, sb.size = sb.size, dark.theme = dark.theme)
  }

  # Center colorscale at 0
  if (is.null(spot.colors)) {
    if (center.zero & !any(data.type %in% c("character", "factor"))) {
      if (!is.null(limits)) {
        max_val <- max(limits)
      } else {
        max_val <- max(abs(data[, variable]))
      }
      limits <- c(-max_val, max_val)
      p <- p +
        scale_fill_gradientn(colours = cols, limits = limits)
        #scale_color_gradient2(low = cols[1], mid = cols[2], high = cols[3], midpoint = 0, limits = limits)
    } else if (any(data.type %in% c("character", "factor"))) {
      p <- p +
        labs(color = variable) +
        scale_fill_manual(values = label.colors)
    } else {
      p <- p +
        scale_fill_gradientn(colours = cols, limits = limits)
    }
  }

  return(p)
}


#' Graphs a smooth interpolation heatmap colored by continuous variable, e.g. dimensional reduction vector
#'
#' @param st.object A Staffli object
#' @param image.type Specifies the image is "processed", otherwise NULL
#' @param highlight.edges SHould edges be highlighted? [default: TRUE]
#'
#' @inheritParams STPlot
#'
#' @importFrom ggplot2 ggplot aes geom_raster scale_x_continuous scale_y_continuous theme_void guides scale_fill_gradient2 labs scale_fill_gradientn ggsave
#' @importFrom magick image_read image_border image_annotate image_composite
#' @importFrom imager as.cimg
#' @importFrom grDevices as.raster

SmoothPlot <- function (
  st.object,
  data,
  image.type,
  data.type = NULL,
  variable,
  palette = "MaYl",
  cols = NULL,
  ncol = NULL,
  center.zero = TRUE,
  dark.theme = FALSE,
  highlight.edges = TRUE,
  ...
) {

  image.masks <- NULL
  if (image.type == "processed") {
    msk.type <- paste0(image.type, ".masks")
    image.masks <- st.object[msk.type]
  } else if ("masked.masks" %in% names(object@tools)) {
    msk.type <- paste0(image.type, ".masks")
    image.masks <- st.object[msk.type]
  }
  samplenames <- names(st.object@samplenames)

  # Obtain colors from selected palette or from provided cols
  cols <- cols %||% {
    palette.select(palette)
  }

  val.limits <- range(data[, variable])

  # Create legend
  lg <- g_legend(data, "numeric", variable, center.zero, cols, val.limits, dark.theme = dark.theme)

  if (image.type %in% c('processed', 'masked')) {
    masks <- lapply(st.object[msk.type], as.cimg)
  } else {
    masks <- NULL
  }

  # Subset only based on one value's expression
  edges.list <- list()
  p.list <- lapply(1:length(unique(data[, "sample"])), function(i) {
    data_subset <- subset(data, sample == i)
    dims <- st.object@rasterlists$processed.masks[[i]] %>% dim()
    if (image.type %in% c('raw', 'masked', 'processed')) {
      extents <- st.object@dims[[i]][2:3] %>% as.numeric()
      data_subset[, c("x", "y")] <- data_subset[, c("x", "y")]/(extents[1]/dims[2])
    } else {
      extents <- st.object@limits[[i]]
      data_subset[, c("x", "y")] <- data_subset[, c("x", "y")]/(extents[1]/dims[2])
    }

    ow <- owin(xrange = c(0, dims[2]), yrange = c(0, dims[1]))
    p <- ppp(x = data_subset[, "x"], y = data_subset[, "y"], window = ow, marks = data_subset[, variable])
    suppressWarnings({s <- Smooth(p, 2, dimyx = dims)})
    m <-  as.matrix(s)
    m[m < 0] <- 0
    m <- m/max(m)

    if (image.type %in% c('processed', 'masked')) {
      msk <- masks[[i]]
      if (highlight.edges) {
        edges.list[[i]] <<- imgradient(msk, "xy") %>% enorm()
      }
      msk <- msk[, , , 1] %>% as.cimg() %>% threshold()
      m <- t(m) %>% as.cimg()
      masked.m <- m*msk
      p <- masked.m
    } else {
      p <- m %>% as.cimg()
    }
  })

  # Draw on new device
  ncol <- ncol %||% ceiling(sqrt(length(p.list)))
  nrow <- ceiling(length(p.list)/ncol)

  cscale <- scales::gradient_n_pal(cols, seq(val.limits[1], 1, length.out = length(x = cols)))

  ims <- list()
  for (i in seq_along(p.list)) {
    p <- p.list[[i]]
    tmp.file2 <- tempfile(pattern = "", fileext = ".png")
    png(width = dim(p)[1], height = dim(p)[2], file = tmp.file2)
    par(mar = c(0, 0, 0, 0), bg = ifelse(dark.theme, "#000000", "#FFFFFF"))
    plot(p, rescale = FALSE, colourscale = cscale)
    dev.off()
    im <- image_read(tmp.file2)
    ims[[i]] <- magick2cimg(im)
  }

  if (!dark.theme & !is.null(masks)) {
    ims <- lapply(seq_along(ims), function(i) {
      im <- ims[[i]]
      msk <- masks[[i]]
      im.masked <- im + !msk
      im.masked[im.masked > 1] <- 1
      ims[[i]] <- im.masked
    })
  }

  if (dark.theme) {
    ims <- lapply(seq_along(ims), function(i) {
      im <- ims[[i]]
      msk <- masks[[i]] %>% threshold()
      im.masked <- im*msk
      ims[[i]] <- im.masked
    })
  }

  if (highlight.edges) {
    ims <- lapply(seq_along(ims), function(i) {
      im <- ims[[i]]
      im <- im + edges.list[[i]]
      im[im > 1] <- 1
      return(im)
    })
  }

  ims <- lapply(ims, function(im) {
    im <- im %>% as.raster() %>% image_read()
    im <- image_border(im, ifelse(dark.theme, "#000000", "#FFFFFF"), paste(st.object@xdim/10, st.object@xdim/10, sep = "x"))
    im <- image_annotate(im, text = i, size = round(st.object@xdim/10), color = ifelse(dark.theme, "#FFFFFF", "#000000"))
  })

  tmp.file <- tempfile(pattern = "", fileext = ".png")
  png(width = dim(p)[1]*ncol, height = dim(p)[1]*nrow, file = tmp.file)
  par(mfrow = c(nrow, ncol), mar = c(0, 0, 0, 0), bg = ifelse(dark.theme, "#000000", "#FFFFFF"))
  for (im in ims) plot(im)
  dev.off()

  im <- image_read(tmp.file)
  im <- image_border(im, ifelse(dark.theme, "#000000", "#FFFFFF"), paste(st.object@xdim/7, st.object@xdim/10, sep = "x"))
  im <- image_annotate(im, text = variable, size = round(st.object@xdim/20), color = ifelse(dark.theme, "#FFFFFF", "#000000"))

  # Draw legend
  tmp.file <- tempfile(pattern = "", fileext = ".png")
  #ggsave(plot = lg, width = 2.8/5, height = 7.8/5, filename = tmp.file, dpi = 150, units = "in")
  #lgim <- image_read(tmp.file)

  grobHeight <- function(x) {
    grid::convertHeight(sum(x$heights), "in", TRUE)
  }

  grobWidth <- function(x) {
    grid::convertWidth(sum(x$widths), "in", TRUE)
  }

  ggsave(plot = lg, width = grobWidth(lg), height = grobHeight(lg), filename = tmp.file)
  iminf <- image_info(im)[2:3] %>% as.numeric()
  lgim <- image_read(tmp.file) %>% image_scale(paste0(iminf[2]/7))
  iminf.lgm <- image_info(lgim)[2:3] %>% as.numeric()

  im <- image_composite(image = im, composite_image = lgim, offset = paste0("+", st.object@xdim*ncol, "+", (iminf[2])/2 - (iminf.lgm[2])/2))

  return(im)
}

# TODO: fix "all" scaling with blend option

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Dimensional reduction plots on HE images
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Dimensional reduction plot on ST coordinates on top of HE image
#'
#' @param sample.index Index specifying the sample that you want to use for plotting
#' @param spots Character vector with spot IDs to plot [default: all spots]
#' @param type Image type to plot on. Here you can specify any of the images available in your Seurat object. To get this list you can
#' run the \code{\link{rasterlists}} function on your Seurat object. If the type is not specified, the images will be prioritized in the following
#' order if they are available; "processed", "masked" and "raw".
#' @param sample.label Should the sample label be included in the image? [default: TRUE]
#' @param show.sb Should the size bar be displayed? [default: TRUE]
#' @param value.scale Defines how the feature values should be mapped to the colorbar. If `value.scale = "samplewise"`, each feature will be
#' scaled independently and if `value.scale = "all"` the features will all have the same value reange.
#' @param ... Extra parameters passed on to \code{\link{ST.ImagePlot}}
#'
#' @inheritParams ST.ImagePlot
#' @inheritParams ST.DimPlot
#' @importFrom cowplot plot_grid
#'
#' @return A ggplot object
#'

spatial_dim_plot <- function (
  object,
  dims = 1:2,
  sample.index = 1,
  spots = NULL,
  type = NULL,
  min.cutoff = NA,
  max.cutoff = NA,
  blend = FALSE,
  pt.size = 1,
  pt.alpha = 1,
  pt.border = FALSE,
  reduction = NULL,
  shape.by = NULL,
  palette = "MaYl",
  cols = NULL,
  grid.ncol = NULL,
  center.zero = TRUE,
  channels.use = NULL,
  verbose = FALSE,
  dark.theme = FALSE,
  sample.label = TRUE,
  show.sb = TRUE,
  value.scale = c("samplewise", "all"),
  custom.theme = NULL,
  ...
) {

  # Check to see if Staffli object is present
  if (!"Staffli" %in% names(object@tools)) stop("Staffli object is missing from Seurat object. Cannot plot without coordinates", call. = FALSE)
  st.object <- object@tools$Staffli

  # Collect data
  spots <- spots %||% colnames(x = object)

  # Check length of sample index
  if (length(sample.index) > 1) stop(paste0("Only one sample index can be selected."), call. = FALSE)

  # Check if reduction is available
  if (!is.null(reduction)) {
    if (!reduction %in% names(object@reductions)) stop(paste0("Reduction ", reduction, " not found in Seurat object."))
  }
  reduction <- reduction %||% {
    default.reductions <- c('umap', 'tsne', 'pca', 'ica')
    object.reductions <- FilterObjects(object = object, classes.keep = 'DimReduc')
    reduc.use <- min(which(x = default.reductions %in% object.reductions))
    default.reductions[reduc.use]
  }

  # Check type
  type <- type %||% {
    if (is.null(rasterlists(st.object))) stop("There are no images present in the Seurat object. Run LoadImages() first.", call. = FALSE)
    choices <- c("processed", "masked", "raw", "processed.masks", "masked.masks")
    match.arg(choices, rasterlists(st.object), several.ok = T)[1]
  }

  # Check that selected image type is present in Seurat object
  msgs <- c("raw" = "LoadImages()", "masked" = "MaskImages()", "processed" = "WarpImages()", "masked.masks" = "MaskImages()", "processed.masks" = "WarpImages()")
  if (!type %in% names(msgs)) stop(paste0(type, " not a valid type"), call. = FALSE)
  if (!type %in% rasterlists(st.object)) stop(paste0("You need to run ", msgs[type], " before using DimOverlay() on '", type, "' images"), call. = FALSE)

  # Check that sample.index is OK
  if (!sample.index %in% names(st.object)) {
    stop(paste0("sample.index ", sample.index, " does not match any of the images present in the Seurat object or is out of range"), call. = T)
  }

  # Collect image
  image <- st.object[type][[sample.index]]
  if (dark.theme & type %in% c("masked", "processed")) {
    image[image == "#FFFFFF"]  <- "#000000"
  }
  if (sample.label) {
    image <- as.raster(image_annotate(image_read(image), text = paste(sample.index), color = ifelse(dark.theme, "#FFFFFF", "#000000"), size = round(st.object@xdim/10)))
  }
  imdims <- st.object@dims[[sample.index]][2:3] %>% as.numeric()

  # Select spots matching sample index
  sample.index <- ifelse(class(sample.index) == "numeric", unique(st.object[[, "sample", drop = T]])[sample.index], sample.index)
  # Select spots matching sample.index
  spots <- intersect(colnames(object)[st.object[[, "sample", drop = T]] == sample.index], spots)
  if (length(spots) == 0) stop(paste0("All selected spots are missing from sample ", sample.index, " ... \n"), call. = FALSE)
  if (verbose) cat(paste0("Selected ", length(spots), " spots matching index ", sample.index))

  # Collect dim-red data
  signs <- sign(dims); dims <- abs(dims)
  data <- Embeddings(object = object[[reduction]])[spots, dims, drop = FALSE]
  data <- as.data.frame(x = t(t(data)*signs))
  dims <- paste0(Key(object = object[[reduction]]), dims)

  # Select colorscale
  palette.info <- palette.select(info = T)
  palette <- palette %||% {
    palette <- subset(palette.info, category == "div")$palette[1]
  }

  # Check that the number of dimensions are 2 or three if blending is active
  if (blend & !length(x = dims) %in% c(2, 3)) {
    stop(paste0("Blending dim plots only works with two or three dimensions. \n",
                "Number of dimensions provided: ", length(x = dims)), call. = F)
  }

  # Obtain array coordinates
  px.ids <- ifelse(rep(type %in% c("raw", "masked", "masked.masks"), 2), c("pixel_x", "pixel_y"), c("warped_x", "warped_y"))

  if (all(px.ids %in% colnames(st.object[[]]))) {
    data <- cbind(data, setNames(st.object[[, px.ids]][spots, ], nm = c("x", "y")))
  } else {
    stop(paste0(paste(px.ids, collapse = " and "), " coordinates are not present in meta data."), call. = FALSE)
  }

  data <- feature.scaler(data, dims, min.cutoff, max.cutoff)
  data[, "sample"] <- sample.index

  # Set scalebar input
  if (show.sb) {
    pixels.per.um <- st.object@pixels.per.um[sample.index]
  } else {
    pixels.per.um <- NULL
  }

  # Set feature scale limits
  if (is.list(value.scale) & length(value.scale) == length(dims)) {
    limits <- setNames(value.scale, dims)
  } else {
    value.scale <- match.arg(value.scale, c("samplewise", "all"))
    if (value.scale == "all") {
      limits <- c(min(data[, dims]), max(data[, dims]))
    } else if (value.scale == "samplewise") {
      limits <- NULL
    }
  }

  if (blend) {
    colored.data <- apply(data[, 1:(ncol(data) - 3)], 2, scales::rescale)
    channels.use <- channels.use %||% c("red", "green", "blue")[1:ncol(colored.data)]

    if (verbose) cat(paste0("Blending colors from dimensions ", paste(paste(dims, channels.use, sep = ":"), collapse = ", ")))

    spot.colors <- ColorBlender(colored.data, channels.use)
    data <- data[, (ncol(data) - 2):ncol(data)]
    plot <- ST.ImagePlot(data, data.type = "numeric", shape.by, variable, image, imdims, pt.size, pt.alpha, pt.border,
                         palette, cols, ncol = NULL, spot.colors, center.zero,
                         plot.title = paste(paste(dims, channels.use, sep = ":"), collapse = ", "), FALSE,
                         dark.theme, pixels.per.um = pixels.per.um, NULL, custom.theme, ...)
    return(plot)
  } else {
    spot.colors <- NULL

    if (verbose) cat("Plotting dimensions:",
                     ifelse(length(dims) == 1, dims,  paste0(paste(dims[1:(length(dims) - 1)], collapse = ", "), " and ", dims[length(dims)])))

    # Create plots
    plots <- lapply(X = dims, FUN = function(d) {
      plot <- ST.ImagePlot(data, data.type = "numeric", shape.by, d, image, imdims, pt.size, pt.alpha, pt.border, palette, cols,
                           ncol = NULL, spot.colors, center.zero, plot.title = d, FALSE, dark.theme, pixels.per.um = pixels.per.um,
                           limits[[d]], custom.theme, ...)

      return(plot)
    })
    p <- plot_grid(plotlist = plots, ncol = grid.ncol)
    if (dark.theme) {
      p <- p + dark_theme()
    }
    return(p)
  }
}



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Feature plots on HE images
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#'  Visualize 'features' on one selected HE image
#'
#' Colors spots on an an ST array grid according to a 'feature'
#' (i.e. gene expression (raw counts or scaled) and features available in the meta data slot).
#' NOTE that this function only draws a plot for one sample at the time.
#'
#' @param sample.index Index specifying the sample that you want to use for plotting
#' @param spots Character vector with spot IDs to plot [default: all spots]
#' @param type Image type to plot on. Here you can specify any of the images available in your Seurat object. To get this list you can
#' run the \code{\link{rasterlists}} function on your Seurat object. If the type is not specified, the images will be prioritized in the following
#' order if they are available; "processed", "masked" and "raw".
#' @param slot Which slot to pull expression data from? [dafault: 'data']
#' @param sample.label Should the sample label be included in the image? [default: TRUE]
#' @param show.sb Should the size bar be displayed? [default: TRUE]
#' @param value.scale Defines how the feature values should be mapped to the colorbar. If `value.scale = "samplewise"`, each feature will be
#' scaled independently and if `value.scale = "all"` the features will all have the same value reange.
#' @param ... Extra parameters passed on to \code{\link{ST.ImagePlot}}
#'
#' @inheritParams ST.ImagePlot
#' @inheritParams ST.FeaturePlot
#' @importFrom cowplot plot_grid
#'
#' @return A ggplot object
#'

spatial_feature_plot <- function (
  object,
  features,
  sample.index = 1,
  spots = NULL,
  type = NULL,
  min.cutoff = NA,
  max.cutoff = NA,
  slot = "data",
  blend = FALSE,
  pt.size = 2,
  pt.alpha = 1,
  pt.border = FALSE,
  shape.by = NULL,
  palette = NULL,
  cols = NULL,
  ncol = NULL,
  grid.ncol = NULL,
  center.zero = FALSE,
  channels.use = NULL,
  verbose = FALSE,
  split.labels = FALSE,
  dark.theme = FALSE,
  sample.label = TRUE,
  show.sb = TRUE,
  value.scale = c("samplewise", "all"),
  custom.theme = NULL,
  ...
) {

  # Check to see if Staffli object is present
  if (!"Staffli" %in% names(object@tools)) stop("Staffli object is missing from Seurat object. Cannot plot without coordinates", call. = FALSE)
  st.object <- object@tools$Staffli

  # Obtain spots
  spots <- spots %||% colnames(object)

  # Check length of sample index
  if (length(sample.index) > 1) stop(paste0("Only one sample index can be selected."), call. = FALSE)

  type <- type %||% {
    if (is.null(rasterlists(st.object))) stop("There are no images present in the Seurat object. Run LoadImages() first.", call. = FALSE)
    choices <- c("processed", "masked", "raw", "processed.masks", "masked.masks")
    match.arg(choices, rasterlists(st.object), several.ok = T)[1]
  }

  # Check that selected image type is present in Seurat object
  msgs <- c("raw" = "LoadImages()", "masked" = "MaskImages()", "processed" = "WarpImages()", "masked.masks" = "MaskImages()", "processed.masks" = "WarpImages()")
  if (!type %in% names(msgs)) stop(paste0(type, " not a valid type"), call. = FALSE)
  if (!type %in% rasterlists(st.object)) stop(paste0("You need to run ", msgs[type], " before using DimOverlay() on '", type, "' images"), call. = FALSE)

  # Check that sample.index is OK
  if (!sample.index %in% names(st.object)) {
    stop(paste0("sample.index ", sample.index, " does not match any of the images present in the Seurat object or is out of range"), call. = T)
  }

  # Collect image
  image <- st.object[type][[sample.index]]
  if (dark.theme & type %in% c("masked", "processed")) {
    image[image == "#FFFFFF"]  <- "#000000"
  }
  if (sample.label) {
    image <- as.raster(image_annotate(image_read(image), text = paste(sample.index), color = ifelse(dark.theme, "#FFFFFF", "#000000"), size = round(st.object@xdim/10)))
  }
  imdims <- st.object@dims[[sample.index]][2:3] %>% as.numeric()

  # Select spots matching sample index
  sample.index <- ifelse(class(sample.index) == "numeric", unique(st.object[[, "sample", drop = T]])[sample.index], sample.index)
  spots <- intersect(colnames(object)[st.object[[, "sample", drop = T]] == sample.index], spots)
  if (length(spots) == 0) stop(paste0("All selected spots are missing from sample ", sample.index, " ... \n"), call. = FALSE)
  if (verbose) cat(paste0("Selected ", length(spots), " spots matching index ", sample.index))

  data <- FetchData(object = object, vars = c(features), cells = spots, slot = slot)
  data.type <- unique(sapply(data, class))

  # Select colorscale
  palette.info <- palette.select(info = T)
  palette <- palette %||% {
    palette <- subset(palette.info, category == "seq")$palette[1]
  }

  # Check that the number of dimensions are 2 or three if blending is active
  if (blend & !length(x = features) %in% c(2, 3)) {
    stop(paste0("Blending dim plots only works with two or three features. \n",
                "Number of dimensions provided: ", length(x = features)), call. = F)
  }

  # Obtain array coordinates
  px.ids <- ifelse(rep(type %in% c("raw", "masked", "masked.masks"), 2), c("pixel_x", "pixel_y"), c("warped_x", "warped_y"))

  if (all(px.ids %in% colnames(st.object[[]]))) {
    data <- cbind(data, setNames(st.object[[, px.ids]][spots, ], nm = c("x", "y")))
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

  if (all(data.type %in% c("numeric", "integer"))) {
    data <- feature.scaler(data, features, min.cutoff, max.cutoff)
  }

  # Add index column
  data[, "sample"] <- sample.index

  # Set scalebar input
  if (show.sb) {
    pixels.per.um <- st.object@pixels.per.um[sample.index]
  } else {
    pixels.per.um <- NULL
  }

  # Set feature scale limits
  if (is.list(value.scale) & length(value.scale) == length(features)) {
    limits <- setNames(value.scale, features)
  } else {
    value.scale <- match.arg(value.scale, c("samplewise", "all"))
    if (value.scale == "all") {
      limits <- c(min(data[, features]), max(data[, features]))
    } else if (value.scale == "samplewise") {
      limits <- NULL
    }
  }

  if (blend) {
    colored.data <- apply(data[, 1:(ncol(data) - 3)], 2, scales::rescale)
    channels.use <- channels.use %||% c("red", "green", "blue")[1:ncol(colored.data)]

    if (verbose) cat(paste0("Blending colors from features ", paste(paste(features, channels.use, sep = ":"), collapse = ", ")))

    spot.colors <- ColorBlender(colored.data, channels.use)
    data <- data[, (ncol(data) - 2):ncol(data)]
    plot <- ST.ImagePlot(data, data.type, shape.by, variable, image, imdims, pt.size, pt.alpha, pt.border,
                         palette, cols, ncol = NULL, spot.colors, center.zero,
                         plot.title = paste(paste(features, channels.use, sep = ":"), collapse = ", "),
                         split.labels, dark.theme, pixels.per.um = pixels.per.um, NULL, custom.theme, ...)
    return(plot)
  } else {
    spot.colors <- NULL

    if (verbose) cat("Plotting features:",
                     ifelse(length(features) == 1, features,  paste0(paste(features[1:(length(features) - 1)], collapse = ", "), " and ", features[length(features)])))

    # Create plots
    plots <- lapply(X = features, FUN = function(d) {
      plot <- ST.ImagePlot(data, data.type, shape.by, d, image, imdims, pt.size, pt.alpha, pt.border, palette, cols,
                           ncol = ncol, spot.colors, center.zero, d, split.labels, dark.theme, pixels.per.um = pixels.per.um,
                           limits[[d]], custom.theme, ...)

      return(plot)
    })

    p <- plot_grid(plotlist = plots, ncol = grid.ncol)
    if (dark.theme) {
      p <- p + dark_theme()
    }
    return(p)
  }
}


#' Graphs ST spots colored by continuous variable, e.g. dimensional reduction vector
#'
#' @param pixels.per.um Defines the number of pixels per micrometer to draw the scale bar
#' @param limits Sets the limits of the colorbar
#'
#' @importFrom ggplot2 geom_point aes_string scale_x_continuous scale_y_continuous theme_void theme_void labs scale_color_gradient2 scale_color_gradientn annotation_custom scale_color_manual
#' @importFrom magick image_info
#' @importFrom grid rasterGrob unit
#' @importFrom grDevices as.raster
#'
#' @param image Image of class "raster" to use as background for plotting
#' @param dims List of dimensions for original images. This list has to contain one element for each sample and each element
#' should be a vector of length 2 specifying the dimensions of the original HE image.
#'
#' @inheritParams STPlot
#'
#' @export

ST.ImagePlot <- function (
  data,
  data.type,
  shape.by,
  variable,
  image,
  dims,
  pt.size = 2,
  pt.alpha = 1,
  pt.border = FALSE,
  palette = "MaYl",
  cols = NULL,
  ncol = NULL,
  spot.colors = NULL,
  center.zero = TRUE,
  plot.title = NULL,
  split.labels = FALSE,
  dark.theme = FALSE,
  pixels.per.um = NULL,
  limits = NULL,
  custom.theme = NULL,
  ...
) {

  # Remove NA values from data
  data <- na.omit(data)

  # Define function to generate ggplot2 default colors
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }

  # Run checks to obtain spot colors
  if (is.null(spot.colors)) {
    if (class(data[, variable]) == "factor") {
      if (!is.null(cols)) {
        if (!length(cols) >= length(unique(data[, variable]))) stop("Not enough colors provided ... \n")
        if (is.null(names(cols))) {
          label.colors <- cols[1:length(unique(data[, variable]))]
        } else {
          if (!all(unique(as.character(data[, variable])) %in% names(cols))) stop("cols names must match variables. Please check the names of the cols parameter ...")
          label.colors <- cols[unique(as.character(data[, variable]))]
          names(label.colors) <- unique(as.character(data[, variable]))
        }
      } else {
        label.colors <- gg_color_hue(length(levels(data[, variable])))
      }
      if (is.null(names(label.colors))) {
        names(label.colors) <- levels(data[, variable])
      }
    } else if (class(data[, variable]) == "character") {
      if (!is.null(cols)) {
        stopifnot(length(cols) >= length(unique(data[, variable])))
        if (is.null(names(cols))) {
          label.colors <- cols[1:length(unique(data[, variable]))]
        } else {
          if (!all(unique(data[, variable]) %in% names(cols))) stop("cols names must match variables. Please check the names of the cols parameter ...")
          label.colors <- cols[unique(data[, variable])]
        }
      } else {
        label.colors <- gg_color_hue(length(unique(data[, variable])))
      }
      names(label.colors) <- unique(data[, variable])
    }
  }

  # Stop if split.labels is activated and there are more than 1 samples
  if (any(data.type %in% c("character", "factor")) & split.labels) {
    plot.title <- paste0("Sample ", unique(as.character(data[, "sample"])), ": ", variable)
    new.data <- data.frame()
    # Order by decreasing size
    if (class(data[, variable]) != "factor") data[, variable] <- factor(data[, variable], levels = unique(data[, variable]))
    levels.keep <- levels(data[, variable])

    for (lbl in levels.keep) {
      dt <- data
      dt[, variable] <- ifelse(dt[, variable] == lbl, lbl, "-")
      dt[, "sample"] <- lbl
      new.data <- rbind(new.data, dt)
    }

    data <- new.data
    levels.keep <- c("-", levels.keep)
    data[, variable] <- factor(data[, variable], levels = levels.keep)
    label.colors <- c("-" = "lightgray", label.colors)
  }

  # Obtain colors from selected palette or from provided cols
  cols <- cols %||% {
    palette.select(palette)
  }

  # Define limits
  c(x_dim, y_dim) %<-% dims

  # Draw image
  g <- rasterGrob(image, width = unit(1, "npc"), height = unit(1, "npc"), interpolate = TRUE)

  # Create new plot
  plots <- lapply(unique(data[, "sample"]), function(v) {
    p <- ggplot() +
      annotation_custom(g, -Inf, Inf, -Inf, Inf)

    if (length(spot.colors) > 0) {

      # Add shape aesthetic and blend colors if blend is active
      if (!is.null(shape.by)) {
        p <- p + geom_point(data = data[data[, "sample"] == v, ], mapping = aes_string(x = "x", y = paste0(y_dim, " - y"), shape = shape.by),
                            shape = 21, stroke = ifelse(pt.border, 0.2, 0),
                            fill = spot.colors, size = pt.size, alpha = pt.alpha, ...)
      } else {
        p <- p + geom_point(data = data[data[, "sample"] == v, ], mapping = aes_string(x = "x", y = paste0(y_dim, " - y")),
                            shape = 21, stroke = ifelse(pt.border, 0.2, 0),
                            fill = spot.colors, size = pt.size, alpha = pt.alpha, ...)
      }

    } else {

      # Add shape aesthetic only
      if (!is.null(shape.by)) {
        p <- p + geom_point(data = data[data[, "sample"] == v, ], mapping = aes_string(x = "x", y = paste0(y_dim, " - y"), fill = paste0("`", variable, "`"), shape = shape.by),
                            shape = 21, stroke = ifelse(pt.border, 0.2, 0),
                            size = pt.size, alpha = pt.alpha, ...)
      } else {
        p <- p + geom_point(data = data[data[, "sample"] == v, ], mapping = aes_string(x = "x", y = paste0(y_dim, " - y"), fill = paste0("`", variable, "`")),
                            shape = 21, stroke = ifelse(pt.border, 0.2, 0),
                            size = pt.size, alpha = pt.alpha, ...)
      }
    }

    # Add ST array dimensions
    p <- p +
      scale_x_continuous(limits = c(0, x_dim), expand = c(0, 0)) +
      scale_y_continuous(limits = c(0, y_dim), expand = c(0, 0)) +
      theme_void() +
      theme(plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "npc"))

    # Add ST array dimensions and plot title
    p <- p +
      labs(title = ifelse(!is.null(plot.title), plot.title, ""), fill = ifelse(all(data.type %in% c("numeric", "integer")), "value", "label"))

    ## Set the scale bar
    if (!is.null(pixels.per.um)) {
      hewidth <- dims[1]
      sb500 <- pixels.per.um*500
      p <- draw_scalebar(p, x = 7*hewidth/9, xend = 7*hewidth/9 + sb500, y = dims[2] - dims[2]/8, dark.theme = dark.theme)
    }

    # Center colorscale at 0
    if (is.null(spot.colors)) {
      #print(!any(data.type %in% c("character", "factor")))
      if (center.zero & !any(data.type %in% c("character", "factor"))) {
        if (!is.null(limits)) {
          max_val <- max(limits)
        } else {
          max_val <- max(abs(data[, variable]))
        }
        limits <- c(-max_val, max_val)
        p <- p +
          scale_fill_gradientn(colours = cols, limits = limits)
        #scale_color_gradient2(low = cols[1], mid = cols[2], high = cols[3], midpoint = 0, limits = limits)
      } else if (any(data.type %in% c("character", "factor"))) {
        p <- p +
          labs(fill = variable) +
          scale_fill_manual(values = label.colors)
      } else {
        p <- p +
          scale_fill_gradientn(colours = cols, limits = limits)
      }
    }

    if (dark.theme) {
      p <- p + dark_theme()
    }

    if (!is.null(custom.theme)) {
      p <- p + custom.theme
    }

    return(p)
  })

  final.plot <- plot_grid(plotlist = plots, ncol = ncol)
}


#' Overlay dimensionality reduction vectors on HE images
#'
#' Graphs the selected vectors of a dimensional reduction technique on a 2D grid of spots overlaid on top of an HE images.
#' Draws sample 1 as default, but can take multiple samples as well.
#'
#' @details It is typically difficult to explore details in the HE image when diplaying multiple samples side by side,
#' so we recommend to draw the plots for one sample at the time. If you have higher resolution images,
#' it could also take significant time to draw the plots.
#'
#' @section Blending values:
#' The blend option can be useful if you wish to visualize multiple dimensionality reduction simultaneuosly and works for two or three value vectors.
#' Each of the selected vectors are rescaled from 0 to 1 and are used as RGB color channels to produce mixed color for each
#' spot. This can be particularly useful when looking at overlapping value vectors. For example, if you are looking at two overlapping value vectors
#' "A" and "B" and use the blend option, the "A" values will be encoded in the "red" channel and the "B" values in the "green" channel. If a spot is
#' purely "A" the color will be red and if it is purely "B" it will green. Any mixture of "A" and "B" will produce a color between red and green
#' where a 50/50 mixture gives yellow color.
#'
#' @section Arrange plots:
#'
#' The `ncols.dims` argument will determine how each subplot called using
#' \code{\link{DimOverlay}} is arranged and will by default put all dims in 1 row, i.e.
#' `ncols.dims = length(dims)`. The `ncols.samples` argument will determine how these subplots
#' are arranged and will by default use 1 column, meaning that each subplot is put in its own row.
#' The output layout matrix would then be `length(samples)*length(dims)`
#'
#' @param object Seurat object
#' @param sampleids Integer vector specifying sample indices to include in the plot [default: 1]
#' @param ncols.dims Number of columns passed to \code{\link{DimOverlay}}. For example,
#' if you are plotting 4 dims, `ncols.dims = 2` will arrange the \code{\link{DimOverlay}}
#' plots into a 2x2 grid [default: `length(dims)`]. (see \emph{Arrange plots*} for a detailed description)
#' @param ncols.samples Number of columns in the layout grid for the samples. For example,
#' if you are plotting 4 samples, `ncols.samples = 2` will arrange the plots into a 2x2 grid [default: `1`].
#' (see \emph{Arrange plots} for a detailed description)
#' @param ... Parameters passed to other methods
#'
#' @inheritParams spatial_dim_plot
#'
#' @examples
#'
#' # Load images and run PCA
#' se <- LoadImages(se) %>%
#'    RunPCA()
#'
#' # Plot the first 2 dimensions on the first two tissue sections
#' DimOverlay(se, dims = 1:2, reduction = "pca", sampleids = 1:2)
#'
#' # Blend values for dimensions 1 and 2 on the first two tissue sections
#' DimOverlay(se, dims = 1:2, reduction = "pca", sampleids = 1:2, blend = T)
#'
#' # Plot the first 2 dimensions and trim off 1st percentile values on the first two tissue sections
#' DimOverlay(se, dims = 1:2, reduction = "pca", sampleids = 1:2, min.cutoff = 'q1')
#'
#' # Mask images and plot the first 2 dimensions on the masked images for samples 1 and 2
#' se <- MaskImages(se)
#' DimOverlay(se, dims = 1:2, reduction = "pca", sampleids = 1:2, type = "masked")
#'
#' @export
#'
#' @seealso \code{\link{ST.FeaturePlot}} and \code{\link{ST.DimPlot}} for how to plot features
#' without the HE image and \code{\link{FeatureOverlay}} for how to overlay feature plots on the HE images.
#'

DimOverlay <- function (
  object,
  dims = c(1:2),
  reduction = NULL,
  sampleids = 1,
  spots = NULL,
  ncols.dims = NULL,
  ncols.samples = NULL,
  type = NULL,
  min.cutoff = NA,
  max.cutoff = NA,
  blend = FALSE,
  pt.size = 1,
  pt.alpha = 1,
  pt.border = TRUE,
  shape.by = NULL,
  palette = "MaYl",
  cols = NULL,
  center.zero = TRUE,
  channels.use = NULL,
  dark.theme = FALSE,
  sample.label = TRUE,
  show.sb = TRUE,
  value.scale = c("samplewise", "all"),
  custom.theme = NULL,
  verbose = FALSE,
  ...
) {

  # Check to see if Staffli object is present
  if (!"Staffli" %in% names(object@tools)) stop("Staffli object is missing from Seurat object. Cannot plot without coordinates", call. = FALSE)
  st.object <- object@tools$Staffli

  # Select spots
  Staffli_meta <- subset(st.object[[]], sample %in% paste0(sampleids))
  selected.spots <- rownames(Staffli_meta)
  spots <- spots %||% intersect(colnames(object), selected.spots)
  if (length(spots) == 0) stop(paste0("None of the selected spots are present in samples ", paste(sampleids, collapse = ", "), " ... \n"), call. = FALSE)

  # Check that spots are present in all sampleids samples
  Staffli_meta_subset <- Staffli_meta[spots, ]
  remaining_samples <- unique(Staffli_meta_subset$sample)[which(unique(Staffli_meta_subset$sample) %in% sampleids)]
  if (length(x = remaining_samples) != length(x = sampleids)) warning(paste0("The selected spots are not present in all samples ", paste(sampleids, collapse = ", "), " ... \n",
                                                                             "Subsetting data to include samples ", paste(remaining_samples, collapse = ", "), "... \n"), call. = FALSE)

  ncols.dims <- ncols.dims %||% length(x = dims)
  ncols.samples <- ncols.samples %||% 1

  # Set scale if provided
  if (is.numeric(value.scale) & length(value.scale) == 2) {
    value.scale.list <- rep(list(value.scale), length(dims))
  } else {
    value.scale <- match.arg(value.scale, c("samplewise", "all"))
    if (value.scale == "all") {
      reduction <- reduction %||% {
        default.reductions <- c('umap', 'tsne', 'pca', 'ica')
        object.reductions <- FilterObjects(object = object, classes.keep = 'DimReduc')
        reduc.use <- min(which(x = default.reductions %in% object.reductions))
        default.reductions[reduc.use]
      }
      signs <- sign(dims); dims <- abs(dims)
      data <- Embeddings(object = object[[reduction]])[spots, dims, drop = FALSE]
      data <- as.data.frame(x = t(t(data)*signs))
      value.scale.list <- lapply(data, range)
    } else {
      value.scale.list <- "samplewise"
    }
  }
  p.list <- lapply(remaining_samples, function(i) {
    spatial_dim_plot(object, dims = dims, sample.index = i, spots = spots, type = type, min.cutoff = min.cutoff,
               max.cutoff = max.cutoff, blend = blend, pt.size = pt.size, pt.alpha = pt.alpha, pt.border = pt.border,
               reduction = reduction, shape.by = shape.by, palette = palette,
               cols = cols, grid.ncol = ncols.dims,
               center.zero = center.zero, channels.use = channels.use, dark.theme = dark.theme,
               sample.label = sample.label, show.sb = show.sb,
               value.scale = value.scale.list, custom.theme = custom.theme, verbose = verbose, ... = ...)
  })
  p <- cowplot::plot_grid(plotlist = p.list, ncol = ncols.samples)
  if (dark.theme) p <- p + dark_theme()
  return(p)
}


#' Overlay features on HE images
#'
#' Graphs the selected features on a 2D grid of spots overlaid on top of an HE images.
#' If numerical features are selected, e.g. gene expression values, a "spatial heatmap" will be drawn and
#' if a categorical variable is selected (such as cluster labels) each group of spots will be assigned a
#' color. Categorical and numerical values cannot be combined.
#'
#' NOTE that this function draws sample 1 as default, but can take multiple samples as well using the `sampleids argument`.
#'
#' @details It is typically difficult to explore details in the HE image when diplaying multiple samples side by side,
#' so we recommend to draw the plots for one sample at the time. If you have higher resolution images,
#' it could also take significant time to draw the plots.
#'
#' @section Blending values:
#' The blend option can be useful if you wish to visualize multiple numerical features simultaneuosly and works for two or three feature value vectors.
#' Each of the selected vectors are rescaled from 0 to 1 and are used as RGB color channels to produce mixed color for each
#' spot. This can be particularly useful when looking at overlapping value vectors. For example, if you are looking at two overlapping value vectors
#' "A" and "B" and use the blend option, the "A" values will be encoded in the "red" channel and the "B" values in the "green" channel. If a spot is
#' purely "A" the color will be red and if it is purely "B" it will green. Any mixture of "A" and "B" will produce a color between red and green
#' where a 50/50 mixture gives yellow color.
#'
#' @section Arrange plots:
#'
#' The `ncols.features` argument will determine how each subplot called using
#' \code{\link{DimOverlay}} is arranged and will by default put all dims in 1 row, i.e.
#' `ncols.features = length(features)`. The `ncols.samples` argument will determine how these subplots
#' are arranged and will by default use 1 column, meaning that each subplot is put in its own row.
#' The output layout matrix would then have the dimensions `length(samples)xlength(features)`
#'
#' @section Splitting categorical features:
#' If you are plotting a categorical feature, e.g.cluster labels, you have the option to split each label into facets using \code{split.labels=TRUE}.
#' This is very useful if you have many different labels which can make it difficult to distinguish the different colors.
#'
#' @section Arrange plots:
#'
#' The `ncols.features` argument will determine how each subplot is arranged and will by default put all features in 1 row, i.e.
#' `ncols.features = length(features)`. The `ncols.samples` argument will determine how these subplots
#' are arranged and will by default use 1 column, meaning that each subplot is put in its own row.
#' The output layout matrix would then have the dimensions `length(samples)xlength(features)`
#'
#' @param object Seurat object
#' @param sampleids Names of samples to plot
#' @param ncols.features Number of columns passed to \code{\link{FeatureOverlay}}. For example,
#' if you are plotting 4 features, `ncols.features = 2` will arrange the \code{\link{FeatureOverlay}}
#' plots into a 2x2 grid [default: `length(features)`]. (see \emph{Arrange plots*} for a detailed description)
#' @param ncols.samples Number of columns in the layout grid for the samples. For example,
#' if you are plotting 4 samples, `ncols.samples = 2` will arrange the plots obtained
#' from \code{\link{FeatureOverlay}} plots into a 2x2 grid [default: `1`].
#' (see \emph{Arrange plots*} for a detailed description)
#' @param split.feature.ncol Sets the number of columns on if split.labels is active
#' @param ... Parameters passed to DimOverlay
#' @inheritParams spatial_feature_plot
#'
#' @examples
#'
#' # Load images
#' se <- LoadImages(se)
#'
#' # Overlay the number of unique genes and the number of
#' # UMIs per spot on the first two tissue sections
#' FeatureOverlay(se, features = c("nFeature_RNA", "nCount_RNA"), sampleids = 1:2)
#'
#' # Plot selected genes on the first two tissue sections
#' FeatureOverlay(se, features = c("Cck", "Dcn"), sampleids = 1:2)
#'
#' # Plot normalized values on the first two tissue sections
#' se <- SCTransform(se)
#' FeatureOverlay(se, features = c("Cck", "Dcn"), sampleids = 1:2)
#'
#' # Change to scaled data
#' FeatureOverlay(se, features = c("Cck", "Dcn"), sampleids = 1:2,
#'                slot = "scale.data", center.zero = TRUE)
#'
#' # Mask images and plot plot the slected genes on the masked images for samples 1 and 2
#' se <- MaskImages(se)
#' FeatureOverlay(se, features = c("Cck", "Dcn"), sampleids = 1:2, type = "masked")
#'
#' @export
#'
#' @seealso \code{\link{ST.FeaturePlot}} and \code{\link{ST.DimPlot}} for how to plot features
#' without the HE image and \code{\link{DimOverlay}} for how to overlay dimensionality reduction output on the
#' HE images.
#'

FeatureOverlay <- function (
  object,
  features,
  sampleids = 1,
  spots = NULL,
  ncols.features = NULL,
  ncols.samples = NULL,
  type = NULL,
  min.cutoff = NA,
  max.cutoff = NA,
  slot = "data",
  blend = FALSE,
  pt.size = 2,
  pt.alpha = 1,
  pt.border = TRUE,
  shape.by = NULL,
  palette = NULL,
  cols = NULL,
  split.labels = FALSE,
  center.zero = FALSE,
  channels.use = NULL,
  dark.theme = FALSE,
  sample.label = TRUE,
  show.sb = TRUE,
  value.scale = c("samplewise", "all"),
  custom.theme = NULL,
  split.feature.ncol = NULL,
  verbose = FALSE,
  ...
) {

  # Check to see if Staffli object is present
  if (!"Staffli" %in% names(object@tools)) stop("Staffli object is missing from Seurat object. Cannot plot without coordinates", call. = FALSE)
  st.object <- object@tools$Staffli

  # Select spots
  Staffli_meta <- subset(st.object[[]], sample %in% paste0(sampleids))
  selected.spots <- rownames(Staffli_meta)
  spots <- spots %||% intersect(colnames(object), selected.spots)
  if (length(spots) == 0) stop(paste0("None of the selected spots are present in samples ", paste(sampleids, collapse = ", "), " ... \n"), call. = FALSE)

  # Check that spots are present in all sampleids samples
  Staffli_meta_subset <- Staffli_meta[spots, ]
  remaining_samples <- unique(Staffli_meta_subset$sample)[which(unique(Staffli_meta_subset$sample) %in% sampleids)]
  if (length(x = remaining_samples) != length(x = sampleids)) warning(paste0("The selected spots are not present in all samples ", paste(sampleids, collapse = ", "), " ... \n",
                                                                             "Subsetting data to include samples ", paste(remaining_samples, collapse = ", "), "... \n"), call. = FALSE)

  ncols.features <- ncols.features %||% length(x = features)
  ncols.samples <- ncols.samples %||% 1

  # Set scale if provided
  if (is.numeric(value.scale) & length(value.scale) == 2) {
    value.scale.list <- rep(list(value.scale), length(features))
  } else {
    value.scale <- match.arg(value.scale, c("samplewise", "all"))
    if (value.scale == "all") {
      data <- FetchData(object = object, vars = c(features), cells = spots, slot = slot)
      value.scale.list <- lapply(data, range)
    } else {
      value.scale.list <- "samplewise"
    }
  }

  p.list <- lapply(remaining_samples, function(s) {
    spatial_feature_plot(object, features = features, sample.index = s, spots = spots, type = type,
                   min.cutoff = min.cutoff, max.cutoff = max.cutoff, slot = slot,
                   blend = blend, pt.size = pt.size, pt.alpha = pt.alpha, pt.border = pt.border, shape.by = shape.by,
                   palette = palette, cols = cols, ncol = split.feature.ncol,
                   grid.ncol = ncols.features, center.zero = center.zero,
                   channels.use = channels.use, split.labels = split.labels, dark.theme = dark.theme,
                   sample.label = sample.label, show.sb = show.sb, value.scale = value.scale.list,
                   custom.theme = custom.theme, verbose = verbose, ... = ...)
  })

  p <- cowplot::plot_grid(plotlist = p.list, ncol = ncols.samples)
  if (dark.theme) p <- p + dark_theme()
  return(p)
}

