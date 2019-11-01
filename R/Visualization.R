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
#' @param object Seurat object
#' @param dims Dimensions to plot [default: 1, 2]
#' @param spots Character vector with spot IDs to plot [default: all spots]
#' @param plot.type Specify the type of plot to use [default: "spots"]. Available options are; "spots" (a "smooth" options will be added soon)
#' @param blend Scale and blend expression values to visualize coexpression of two features (this options will override other coloring parameters).
#' See 'Blending values' below for a more thourough description.
#' @param min.cutoff,max.cutoff Vector of minimum and maximum cutoff values for each feature, may specify quantile in the form of 'q##' where '##'
#' is the quantile (eg, 'q1', 'q10'). This can be useful if you have outlier values that skew the colorscale in the plot. For example, if you specify
#' 'q1', you will trim of values below the 1st percentile. [default: no cuttoffs]
#' @param pt.size Adjust point size for plotting [default: 1]
#' @param pt.alpha Adjust point opacity for plotting [default: 1]
#' @param reduction Which dimensionality reduction to use. If not specified, first searches for "umap", then "tsne", then "pca"
#' @param shape.by You can specify any spot attribute (that can be pulled with FetchData) allowing for both different colors
#' and different shapes on spots
#' @param grid.ncol Number of columns for display when combining plots. This option will only have an effect on the sample level structure.
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
#' HE images and \code{\link{MultiFeatureOverlay}} and \code{\link{MultiDimOverlay}} for how
#' to overlay plots on the HE images in multiple samples.
#'

ST.DimPlot <- function (
  object,
  dims = c(1, 2),
  spots = NULL,
  plot.type = "spots",
  blend = FALSE,
  min.cutoff = NA,
  max.cutoff = NA,
  pt.size = 1,
  pt.alpha = 1,
  reduction = NULL,
  shape.by = NULL,
  palette = "MaYl",
  cols = NULL,
  rev.cols = FALSE,
  dark.theme = FALSE,
  ncol = NULL,
  grid.ncol = NULL,
  center.zero = TRUE,
  channels.use = NULL,
  center.tissue = FALSE,
  verbose = FALSE,
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
  data[,  "sample"] <- st.object[[, "sample", drop = TRUE]]

  # Extract shape.by column from meta data if applicable
  if (!is.null(x = shape.by)) {
    if (!shape.by %in% colnames(object[[]])) {
      stop(paste0("Shaping variable (shape.by) ", shape.by, " not found in meta.data slot"), call. = F)
    }
    data[, shape.by] <- as.character(object[[shape.by, drop = TRUE]])
  }

  # Obtain array coordinates
  image.type <- "empty"
  c(data, image.type) %<-% obtain.array.coords(st.object, data, image.type)

  # Scale data values
  data <- feature.scaler(data, dims, min.cutoff, max.cutoff, spots)

  # blend colors or plot each dimension separately
  if (blend) {
    colored.data <- apply(data[, 1:(ncol(data) - ifelse(!is.null(shape.by), 4, 3))], 2, rescale)
    channels.use <- channels.use %||% c("red", "green", "blue")[1:ncol(colored.data)]

    if (verbose) cat(paste0("Blending colors for dimensions ",
                            paste0(ifelse(length(dims) == 2, paste0(dims[1], " and ", dims[2]), paste0(dims[1], dims[2], " and ", dims[2]))),
                            ": \n", paste(paste(dims, channels.use, sep = ":"), collapse = "\n")))
    if (image.type != "empty") {
      dims.list <- lapply(iminfo(st.object), function(x) {x[2:3] %>% as.numeric()})
    } else {
      dims.list <- st.object@limits
    }

    spot.colors <- ColorBlender(colored.data, channels.use)
    data <- data[, (ncol(data) - ifelse(!is.null(shape.by), 3, 2)):ncol(data)]
    plot <- STPlot(data, data.type = "numeric", shape.by, NULL, pt.size, pt.alpha,
                   palette, cols, rev.cols, ncol, spot.colors, center.zero, center.tissue,
                   plot.title = paste(paste(dims, channels.use, sep = ":"), collapse = ", "), dims.list, ...)
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
      if (image.type != "empty") {
        dims.list <- lapply(iminfo(st.object), function(x) {x[2:3] %>% as.numeric()})
      } else {
        dims.list <- st.object@limits
      }
      plots <- lapply(X = dims, FUN = function(d) {
        plot <- STPlot(data, data.type = "numeric", shape.by, d, pt.size, pt.alpha,
                       palette, cols, rev.cols, ncol, spot.colors, center.zero, center.tissue, NULL, dims.list, ...)

        if (dark.theme) {
          plot <- plot + dark_theme()
        }
        return(plot)
      })

      # Draw plots
      plot_grid(plotlist = plots, ncol = grid.ncol)

    } else if (plot.type == "smooth") {
      # Smooth visualization -------------------------------------------------------------------------------------
      stop("Smooth options not yet implemented")
      plots <- lapply(X = dims, FUN = function(d) {
        plot <- SmoothPlot(data, image.type, data.type = "numeric", d,
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
#' @param verbose Print messages
#'
#' @param ... Extra parameters passed on to \code{\link{STPlot}}
#'
#' @inheritParams STPlot
#' @importFrom cowplot plot_grid
#' @importFrom scales rescale
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
#' HE images and \code{\link{MultiFeatureOverlay}} and \code{\link{MultiDimOverlay}} for how
#' to overlay plots on the HE images in multiple samples.
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
  shape.by = NULL,
  palette = NULL,
  cols = NULL,
  rev.cols = FALSE,
  dark.theme = FALSE,
  ncol = NULL,
  grid.ncol = NULL,
  center.zero = FALSE,
  channels.use = NULL,
  center.tissue = FALSE,
  theme = theme_void(),
  verbose = FALSE,
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
  data[,  "sample"] <- st.object[[, "sample", drop = TRUE]]

  # Add shape column if specified
  if (!is.null(x = shape.by)) {
    if (!shape.by %in% colnames(object[[]])) {
      stop(paste0("Shaping variable (shape.by) ", shape.by, " not found in meta.data slot"), call. = F)
    }
    data[, shape.by] <- as.character(object[[shape.by, drop = TRUE]])
  }

  # Obtain array coordinates
  image.type <- "empty"
  c(data, image.type) %<-% obtain.array.coords(st.object, data, image.type)

  # Raise error if features are not present in Seurat object
  if (ncol(x = data) < 4) {
    stop("None of the requested features were found: ",
         paste(features, collapse = ", "),
         " in slot ",
         slot,
         call. = FALSE)
  }

  if (data.type %in% c("numeric", "integer")) {
    data <- feature.scaler(data, features = features, min.cutoff, max.cutoff, spots)
  }

  # Subset by index
  if (!is.null(indices)) {
    if (!all(as.character(indices) %in% data[, "sample"])) stop(paste0("Index out of range. "), call. = FALSE)
    data <- data[data[, "sample"] %in% as.character(indices), ]
  }


  if (blend) {
    colored.data <- apply(data[, 1:(ncol(data) - ifelse(!is.null(shape.by), 4, 3))], 2, rescale)
    channels.use <- channels.use %||% c("red", "green", "blue")[1:ncol(colored.data)]

    if (verbose) cat(paste0("Blending colors for features ",
                            paste0(ifelse(length(features) == 2, paste0(features[1], " and ", features[2]), paste0(features[1], features[2], " and ", features[2]))),
                            ": \n", paste(paste(features, channels.use, sep = ":"), collapse = "\n")))
    if (image.type != "empty") {
      dims <- lapply(iminfo(st.object), function(x) {x[2:3] %>% as.numeric()})
    } else {
      dims <- st.object@limits
    }

    if (!is.null(indices)) dims <- dims[indices]

    spot.colors <- ColorBlender(colored.data, channels.use)
    data <- data[, (ncol(data) - ifelse(!is.null(shape.by), 4, 3)):ncol(data)]
    plot <- STPlot(data, data.type, shape.by, NULL, pt.size, pt.alpha,
                   palette, cols, rev.cols, ncol, spot.colors, center.zero, center.tissue,
                   plot.title = paste(paste(features, channels.use, sep = ":"), collapse = ", "),
                   dims, FALSE, theme, ...)
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
      if (image.type != "empty") {
        dims <- lapply(iminfo(st.object), function(x) {x[2:3] %>% as.numeric()})
      } else {
        dims <- st.object@limits
      }

      if (!is.null(indices)) dims <- dims[indices]

      plots <- lapply(X = features, FUN = function(ftr) {
        plot <- STPlot(data, data.type, shape.by, ftr, pt.size, pt.alpha,
                       palette, cols, rev.cols, ncol, spot.colors, center.zero,
                       center.tissue, NULL, dims, split.labels, theme, ...)

        if (dark.theme) {
          plot <- plot + dark_theme()
        }
        return(plot)
      })

      plot_grid(plotlist = plots, ncol = grid.ncol)

    } else if (plot.type == "smooth") {
      stop("Smooth options not yet implemented")
      if (data.type %in% c("factor", "character")) stop(paste0("Smoothing has not yet been implemented for categorical variables"), call. = FALSE)
      # Smooth visualization -------------------------------------------------------------------------------------
      plots <- lapply(X = features, FUN = function(d) {
        plot <- SmoothPlot(data, image.type, data.type, d,
                       palette, cols, rev.cols, ncol, center.zero, xlim, ylim, ...)
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


#' Graphs ST spots colored by continuous or categorical features
#'
#' @param data Object of class 'data.frame' containing at least (x, y) coordinates, a "sample" vector with labels for each sample
#' and one column with the feature values. Can also include an additional column for shapes.
#' @param data.type String specifying the class of the features in data to be plotted
#' @param shape.by String specifying the column where the shaping label is stored
#' @param variable Name of feature column
#' @param pt.size Point size of each ST spot [default: 1]
#' @param pt.alpha Opacity of each ST spot [default: 1]
#' @param palette Color palette used for spatial heatmap (see \code{palette.select(info = T)} for available options).
#' Disabled if a color vector is provided (see \code{cols} below).
#' @param cols A vector of colors to use for colorscale, e.g. \code{cols = c("blue", "white", "red")} will
#' create a gradient color scale going from blue to white to red. This options will deactivate the \code{palette}
#' option.
#' @param rev.cols Logical specifying whether colorscale should be reversed [default: FALSE]
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
#' @param theme Object of class 'theme' used to change the background theme (see for example \url{https://ggplot2.tidyverse.org/reference/theme.html})
#' @param ... Parameters passed to geom_point()
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
  palette = "MaYl",
  cols = NULL,
  rev.cols = F,
  ncol = NULL,
  spot.colors = NULL,
  center.zero = TRUE,
  center.tissue = F,
  plot.title = NULL,
  dims = NULL,
  split.labels = FALSE,
  theme = theme_void(),
  ...
) {

  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }

  if (is.null(spot.colors)) {
    if (class(data[, variable]) == "factor") {
      if (!is.null(cols)) {
        stopifnot(length(cols) >= length(unique(data[, variable])))
        label.colors <- cols[1:length(unique(data[, variable]))]
      } else {
        label.colors <- gg_color_hue(length(levels(data[, variable])))
      }
      names(label.colors) <- levels(data[, variable])
    } else if (class(data[, variable]) == "character") {
      if (!is.null(cols)) {
        stopifnot(length(cols) >= length(unique(data[, variable])))
        label.colors <- cols[1:length(unique(data[, variable]))]
      } else {
        label.colors <- gg_color_hue(length(levels(data[, variable])))
      }
      names(label.colors) <- unique(data[, variable])
    }
  }

  # Stop if split.labels is activated and there are more than 1 samples
  if (any(data.type %in% c("character", "factor")) & split.labels) {
    if (length(unique(as.character(data[, "sample"]))) > 1) stop(paste0("Splitting of group labels only work for one sample. Please set a single sample index with indices. "), call. = FALSE)
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
    ifelse(rep(palette == "heat", 3), palette.select(palette)(4), palette.select(palette)(3))
  }

  if (rev.cols) {
    cols <- rev(cols)
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
    split.data <- split(data, data[, "sample"])
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
  data[, "sample"] <- factor(data[, "sample"], levels = unique(data[, "sample"]))
  if (length(spot.colors) > 0) {

    # Add shape aesthetic and blend colors if blend is active
    if (!is.null(shape.by)) {
      p <- p + geom_point(data = data, mapping = aes_string(x = "x", y = "y", shape = shape.by), color = spot.colors, size = pt.size, alpha = pt.alpha, ...)
    } else {
      p <- p + geom_point(data = data, mapping = aes_string(x = "x", y = "y"), color = spot.colors, size = pt.size, alpha = pt.alpha, ...)
    }

  } else {

    # Add shape aesthetic only
    if (!is.null(shape.by)) {
      p <- p + geom_point(data = data, mapping = aes_string(x = "x", y = "y", color = paste0("`", variable, "`"), shape = shape.by), size = pt.size, alpha = pt.alpha, ...)
    } else {
      p <- p + geom_point(data = data, mapping = aes_string(x = "x", y = "y", color = paste0("`", variable, "`")), size = pt.size, alpha = pt.alpha, ...)
    }

  }

  # Add ST array dimensions and plot title
  p <- p +
    labs(title = ifelse(!is.null(plot.title), plot.title, variable), color = "")

  # Set theme
  p <- p + theme_void() + theme

  # Facet plots by group variable
  if (!split.labels) {
    p <- p +
      facet_wrap_custom(~sample, scales = "free", ncol = ncol, scale_overrides = limits_override)
  } else {
    p <- p + facet_wrap(~sample, ncol = ncol) +
      scale_x_continuous(limits = c(0, lims[1])) +
      scale_y_continuous(limits = c(0, lims[2]))
  }

  # Center colorscale at 0
  if (is.null(spot.colors)) {
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
  variable,
  palette = "MaYl",
  cols = NULL,
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
  # Obtain colors from selected palette or from provided cols
  cols <- cols %||% {
    ifelse(rep(palette == "heat", 3), palette.select(palette)(4), palette.select(palette)(5))
  }
  if (rev.cols) {
    cols <- rev(cols)
  }

  if (darken) cols <- c("black", cols); if (whiten) cols <- c(cols, "white")

  val.limits <- range(data[, variable])

  # Create legend
  lg <- g_legend(data, variable, center.zero, cols, val.limits)

  # Subset only based on one value's expression
  p.list <- lapply(1:length(unique(data[, "sample"])), function(i) {
    data.subset <- subset(data, sample == i)
    if ("xdim" %in% names(object@tools)) {
      xdim <- 400; ydim <- round(as.numeric(object@tools$dims[[i]][2])/(as.numeric(object@tools$dims[[i]][2])/xdim))
      #data.subset <- data.subset[data.subset[, variable] != 0, ]
      s.xy <- as.numeric(object@tools$dims[[1]][, 2:3])/c(xdim, ydim)
      data.subset[, c("x", "y")] <- data.subset[, c("x", "y")]/s.xy
    } else {
      xdim <- 67; ydim <- 64
    }

    x <- data.subset[, "x"]; y <- data.subset[, "y"]
    min.x <- min(x); min.y <- min(y); max.x <- max(x); max.y <- max(y)
    tissue.width <- max.x - min.x; tissue.height <- max.y - min.y

    # Run interpolation
    s =  akima::interp(data.subset[, "x"], data.subset[, "y"], data.subset[, variable], nx = tissue.width, ny = tissue.height, extrap = FALSE, linear = FALSE, xo = 1:xdim, yo = 1:ydim)

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

    if (image.type != "empty") {
      msk <- as.cimg(image.masks[[i]])
      masked.plot <- as.raster(magick2cimg(p)*(msk/255))
      masked.plot <- as.raster(image_annotate(image_read(masked.plot), text = samplenames[i], size = round(xdim/10), color = "#FFFFFF"))
      return(masked.plot)
    } else {
      return(as.raster(p))
    }
  })

  xdim <- ydim <- 400

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

#' Dimensional reduction plot on ST coordinates on top of HE image
#'
#' Graphs the selected vectors of a dimensional reduction technique on a 2D grid of spots overlaid on top of an HE images.
#'
#' @details Note that this function only allows you to plot one sample at the time. It is typically difficult to explore details in the HE image
#' when diplaying multiple samples side by side, so we decided to limit the choice to plot one sample at the time. If you have higher resolution images,
#' it could also take significant time to draw the plots. If you still wish to show multiple samples together you can use the \code{\link{MultiFeatureOverlay}}
#' function.
#'
#' @section Blending values:
#' The blend option can be useful if you wish to visualize multiple dimensions simultaneuosly and works for two or three dimensionality
#' reduction vectors. Each of the selected vectors are rescaled from 0 to 1 and are used as RGB color channels to produce mixed color for each
#' spot. This can be particularly useful when looking at overlapping value vectors. For example, if you are looking at two overlapping value vectors
#' "A" and "B" and use the blend option, the "A" values will be encoded in the "red" channel and the "B" values in the "green" channel. If a spot is
#' purely "A" the color will be red and if it is purely "B" it will green. Any mixture of "A" and "B" will produce a color between red and green
#' where a 50/50 mixture gives yellow color. The amplitude if the values will also determine the brightness of the color.
#'
#'
#' @param sample.index Index specifying the sample that you want to use for plotting
#' @param type Image type to plot on. Here you can specify any of the images available in your Seurat object. To get this list you can
#' run the \code{\link{rasterlists}} function on your Seurat object. If the type is not specified, the images will be prioritized in the following
#' order if they are available; "processed", "masked" and "raw".
#' @param slot Which slot to pull expression data from? [default: 'data']
#' @param ... Extra parameters passed on to \code{\link{ST.ImagePlot}}
#'
#' @inheritParams ST.ImagePlot
#' @inheritParams ST.DimPlot
#' @importFrom cowplot plot_grid
#' @importFrom scales rescale
#'
#' @return A ggplot object
#'
#' @examples
#'
#' # Load images and run PCA
#' se <- LoadImages(se) %>%
#'    RunPCA()
#'
#' # Plot the first 2 dimensions
#' DimOverlay(se, dims = 1:2, reduction = "pca", sample.index = 1)
#'
#' # Blend values for dimensions 1 and 2
#' DimOverlay(se, dims = 1:2, reduction = "pca", sample.index = 1, blend = T)
#'
#' # Plot the first 2 dimensions and trim off 1st percentile values
#' DimOverlay(se, dims = 1:2, reduction = "pca", sample.index = 1, min.cutoff = 'q1')
#'
#' # Mask images and plot plot the first 2 dimensions on the masked image
#' se <- MaskImages(se)
#' DimOverlay(se, dims = 1:2, reduction = "pca", sample.index = 1, type = "masked")
#'
#' @export
#'
#' @seealso \code{\link{ST.FeaturePlot}} and \code{\link{ST.DimPlot}} for how to plot features
#' without the HE image, \code{\link{FeatureOverlay}} for how to overlay feature plots on the
#' HE images and \code{\link{MultiFeatureOverlay}} and \code{\link{MultiDimOverlay}} for how
#' to overlay plots on the HE images in multiple samples.
#'

DimOverlay <- function (
  object,
  dims = 1:2,
  sample.index = 1,
  type = NULL,
  min.cutoff = NA,
  max.cutoff = NA,
  blend = FALSE,
  pt.size = 1,
  pt.alpha = 1,
  reduction = NULL,
  shape.by = NULL,
  palette = NULL,
  cols = NULL,
  rev.cols = FALSE,
  grid.ncol = NULL,
  center.zero = TRUE,
  channels.use = NULL,
  verbose = FALSE,
  dark.theme = FALSE,
  ...
) {

  # Check to see if Staffli object is present
  if (!"Staffli" %in% names(object@tools)) stop("Staffli object is missing from Seurat object. Cannot plot without coordinates", call. = FALSE)
  st.object <- object@tools$Staffli

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
  image <- as.raster(image_annotate(image_read(image), text = paste(sample.index), color = ifelse(dark.theme, "#FFFFFF", "#000000"), size = round(st.object@xdim/10)))
  imdims <- iminfo(st.object)[[sample.index]][2:3] %>% as.numeric()

  # Select spots matching sample index
  sample.index <- ifelse(class(sample.index) == "numeric", unique(st.object[[, "sample", drop = T]])[sample.index], sample.index)
  # Select spots matching sample.index
  spots <- colnames(object)[st.object[[, "sample", drop = T]] == sample.index]
  if (verbose) cat(paste0("Selected ", length(spots), " spots matching index ", sample.index))

  # Collect dim-red data
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

  if (all(px.ids %in% colnames(st.object[[]]))) {
    data <- cbind(data, setNames(st.object[[, px.ids]][spots, ], nm = c("x", "y")))
  } else {
    stop(paste0(paste(px.ids, collapse = " and "), " coordinates are not present in meta data."), call. = FALSE)
  }

  data <- feature.scaler(data, dims, min.cutoff, max.cutoff, spots)
  data[, "sample"] <- sample.index

  if (blend) {
    colored.data <- apply(data[, 1:(ncol(data) - 3)], 2, rescale)
    channels.use <- channels.use %||% c("red", "green", "blue")[1:ncol(colored.data)]

    if (verbose) cat(paste0("Blending colors from features ", paste(paste(dims, channels.use, sep = ":"), collapse = ", ")))

    spot.colors <- ColorBlender(colored.data, channels.use)
    data <- data[, (ncol(data) - 2):ncol(data)]
    plot <- ST.ImagePlot(data, data.type = "numeric", shape.by, variable, image, imdims, pt.size, pt.alpha,
                         palette, cols, rev.cols, ncol = NULL, spot.colors, center.zero,
                         plot.title = paste(paste(dims, channels.use, sep = ":"), collapse = ", "), FALSE, dark.theme, ...)
    return(plot)
  } else {
    spot.colors <- NULL

    if (verbose) cat("Plotting features:",
                     ifelse(length(dims) == 1, dims,  paste0(paste(dims[1:(length(dims) - 1)], collapse = ", "), " and ", dims[length(dims)])))

    # Create plots
    plots <- lapply(X = dims, FUN = function(d) {
      plot <- ST.ImagePlot(data, data.type = "numeric", shape.by, d, image, imdims, pt.size, pt.alpha, palette, cols,
                           rev.cols, ncol = NULL, spot.colors, center.zero, plot.title = NULL, FALSE, dark.theme, ...)

      return(plot)
    })
    plot_grid(plotlist = plots, ncol = grid.ncol)
  }
}


# TODO: add cols options to overlay functions

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Feature plots on HE images
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#'  Visualize 'features' on an ST array grid overlayed on top of HE image
#'
#' Colors spots on an an ST array grid according to a 'feature'
#' (i.e. gene expression (raw counts or scaled) and features available in the meta data slot)
#'
#' @param sample.index Index specifying the sample that you want to use for plotting
#' @param type Image type to plot on. Here you can specify any of the images available in your Seurat object. To get this list you can
#' run the \code{\link{rasterlists}} function on your Seurat object. If the type is not specified, the images will be prioritized in the following
#' order if they are available; "processed", "masked" and "raw".
#' @param slot Which slot to pull expression data from? [dafault: 'data']
#' @param ... Extra parameters passed on to \code{\link{ST.ImagePlot}}
#'
#' @inheritParams ST.ImagePlot
#' @inheritParams ST.FeaturePlot
#' @importFrom cowplot plot_grid
#' @importFrom scales rescale
#'
#' @return A ggplot object
#'
#' @examples
#'
#' # Load images
#' se <- LoadImages(se)
#'
#' # Overlay the number of unique genes and the number of UMIs per spot on sample 1 HE image
#' FeatureOverlay(se, features = c("nFeature_RNA", "nCount_RNA"), sample.index = 1)
#'
#' # Plot selected genes
#' FeatureOverlay(se, features = c("Cck", "Dcn"), sample.index = 1)
#'
#' # Plot normalized values
#' se <- SCTransform(se)
#' FeatureOverlay(se, features = c("Cck", "Dcn"), sample.index = 1)
#'
#' # Change to scaled data
#' FeatureOverlay(se, features = c("Cck", "Dcn"), sample.index = 1, slot = "scale.data", center.zero = TRUE)
#'
#' # Mask images and plot plot the slected genes on the masked image
#' se <- MaskImages(se)
#' FeatureOverlay(se, features = c("Cck", "Dcn"), sample.index = 1, type = "masked")
#'
#' @export
#'
#' @seealso \code{\link{ST.FeaturePlot}} and \code{\link{ST.DimPlot}} for how to plot features
#' without the HE image, \code{\link{DimOverlay}} for how to overlay dimensionality reduction output on the
#' HE images and \code{\link{MultiFeatureOverlay}} and \code{\link{MultiDimOverlay}} for how
#' to overlay plots on the HE images in multiple samples.
#'

FeatureOverlay <- function (
  object,
  features,
  sample.index = 1,
  type = NULL,
  min.cutoff = NA,
  max.cutoff = NA,
  slot = "data",
  blend = FALSE,
  pt.size = 2,
  pt.alpha = 1,
  shape.by = NULL,
  palette = NULL,
  cols = NULL,
  rev.cols = FALSE,
  grid.ncol = NULL,
  center.zero = FALSE,
  channels.use = NULL,
  verbose = FALSE,
  split.labels = FALSE,
  dark.theme = FALSE,
  ...
) {

  # Check to see if Staffli object is present
  if (!"Staffli" %in% names(object@tools)) stop("Staffli object is missing from Seurat object. Cannot plot without coordinates", call. = FALSE)
  st.object <- object@tools$Staffli

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
  image <- as.raster(image_annotate(image_read(image), text = paste(sample.index), color = ifelse(dark.theme, "#FFFFFF", "#000000"), size = round(st.object@xdim/10)))
  imdims <- iminfo(st.object)[[sample.index]][2:3] %>% as.numeric()

  # Select spots matching sample index
  sample.index <- ifelse(class(sample.index) == "numeric", unique(st.object[[, "sample", drop = T]])[sample.index], sample.index)
  # Select spots matching sample.index
  spots <- colnames(object)[st.object[[, "sample", drop = T]] == sample.index]
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

  if (data.type %in% c("numeric", "integer")) {
    data <- feature.scaler(data, features, min.cutoff, max.cutoff, spots)
  }

  # Add index column
  data[, "sample"] <- sample.index

  if (blend) {
    colored.data <- apply(data[, 1:(ncol(data) - 3)], 2, rescale)
    channels.use <- channels.use %||% c("red", "green", "blue")[1:ncol(colored.data)]

    if (verbose) cat(paste0("Blending colors from features ", paste(paste(features, channels.use, sep = ":"), collapse = ", ")))

    spot.colors <- ColorBlender(colored.data, channels.use)
    data <- data[, (ncol(data) - 2):ncol(data)]
    plot <- ST.ImagePlot(data, data.type, shape.by, variable, image, imdims, pt.size, pt.alpha,
                         palette, cols, rev.cols, ncol = NULL, spot.colors, center.zero,
                         plot.title = paste(paste(features, channels.use, sep = ":"), collapse = ", "), split.labels, dark.theme, ...)
    return(plot)
  } else {
    spot.colors <- NULL

    if (verbose) cat("Plotting features:",
                     ifelse(length(features) == 1, features,  paste0(paste(features[1:(length(features) - 1)], collapse = ", "), " and ", features[length(features)])))

    # Create plots
    plots <- lapply(X = features, FUN = function(d) {
      plot <- ST.ImagePlot(data, data.type, shape.by, d, image, imdims, pt.size, pt.alpha, palette, cols,
                           rev.cols, ncol = NULL, spot.colors, center.zero, NULL, split.labels, dark.theme, ...)

      return(plot)
    })
    plot_grid(plotlist = plots, ncol = grid.ncol)
  }
}


#' Graphs ST spots colored by continuous variable, e.g. dimensional reduction vector
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
  palette = "MaYl",
  cols = NULL,
  rev.cols = F,
  ncol = NULL,
  spot.colors = NULL,
  center.zero = T,
  plot.title = NULL,
  split.labels = FALSE,
  dark.theme = FALSE,
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
    ifelse(rep(palette == "heat", 3), palette.select(palette)(4), palette.select(palette)(3))
  }

  if (rev.cols) {
    cols <- rev(cols)
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
        p <- p + geom_point(data = data[data[, "sample"] == v, ], mapping = aes_string(x = "x", y = paste0(y_dim, " - y"), shape = shape.by), color = spot.colors, size = pt.size, alpha = pt.alpha, ...)
      } else {
        p <- p + geom_point(data = data[data[, "sample"] == v, ], mapping = aes_string(x = "x", y = paste0(y_dim, " - y")), color = spot.colors, size = pt.size, alpha = pt.alpha, ...)
      }

    } else {

      # Add shape aesthetic only
      if (!is.null(shape.by)) {
        p <- p + geom_point(data = data[data[, "sample"] == v, ], mapping = aes_string(x = "x", y = paste0(y_dim, " - y"), color = paste0("`", variable, "`"), shape = shape.by), size = pt.size, alpha = pt.alpha, ...)
      } else {
        p <- p + geom_point(data = data[data[, "sample"] == v, ], mapping = aes_string(x = "x", y = paste0(y_dim, " - y"), color = paste0("`", variable, "`")), size = pt.size, alpha = pt.alpha, ...)
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

    if (dark.theme) {
      p <- p + dark_theme()
    }

    return(p)
  })

  final.plot <- plot_grid(plotlist = plots)
}

#' Apply DimOverlay to multiple samples
#'
#' @param object Seurat object
#' @param sampleids Integer vector specifying sample indices to include in the plot
#' @param method Display method; "viewer" or "raster" [default: "viewer"]. The option "viewer" only works within an
#' RStudio session and makes use of the Viewer pane. If you want to include images in an rmarkdown document or export
#' it to a file you have to set \code{method = "raster"}.
#' @param ncols Number of columns in output image. This will override the default which is to arrange the plots in
#' \code{ncols = ceiling(sqrt(#samples)); nrows = ceiling(#samples/ncols)}
#' @param ... Parameters passed to DimOverlay
#'
#' @inheritParams DimOverlay
#'
#' @examples
#'
#' # Load images and run PCA
#' se <- LoadImages(se) %>%
#'    RunPCA()
#'
#' # Plot the first 2 dimensions on the first two samples
#' MultiDimOverlay(se, dims = 1:2, reduction = "pca", sampleids = 1:2)
#'
#' # Blend values for dimensions 1 and 2 on the first two samples
#' MultiDimOverlay(se, dims = 1:2, reduction = "pca", sampleids = 1:2, blend = T)
#'
#' # Plot the first 2 dimensions and trim off 1st percentile values on the first two samples
#' MultiDimOverlay(se, dims = 1:2, reduction = "pca", sampleids = 1:2, min.cutoff = 'q1')
#'
#' # Mask images and plot the first 2 dimensions on the masked images for samples 1 and 2
#' se <- MaskImages(se)
#' MultiDimOverlay(se, dims = 1:2, reduction = "pca", sampleids = 1:2, type = "masked")
#'
#' @export
#'
#' @seealso \code{\link{ST.FeaturePlot}} and \code{\link{ST.DimPlot}} for how to plot features
#' without the HE image, \code{\link{FeatureOverlay}} and \code{\link{DimOverlay}} for how to
#' overlay feature plots on the HE images and \code{\link{MultiFeatureOverlay}} and for how
#' to overlay feature plots on the HE images in multiple samples.
#'

MultiDimOverlay <- function (
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
  pt.alpha = 1,
  reduction = NULL,
  shape.by = NULL,
  palette = NULL,
  cols = NULL,
  rev.cols = FALSE,
  center.zero = FALSE,
  channels.use = NULL,
  verbose = FALSE,
  dark.theme = FALSE,
  ...
) {

  # Check to see if Staffli object is present
  if (!"Staffli" %in% names(object@tools)) stop("Staffli object is missing from Seurat object. Cannot plot without coordinates", call. = FALSE)
  st.object <- object@tools$Staffli

  ncols <- ncols %||% ceiling(sqrt(length(x = sampleids)))
  nrows <- ceiling(length(x = sampleids)/ncols)

  p.list <- lapply(sampleids, function(i) {
    DimOverlay(object, dims = dims, sample.index = i, type = type, min.cutoff = min.cutoff,
               max.cutoff = max.cutoff, blend = blend, pt.size = pt.size, pt.alpha,
               reduction = reduction, shape.by = shape.by, palette = palette,
               cols = cols, rev.cols = rev.cols, grid.ncol = NULL,
               center.zero = center.zero, channels.use = channels.use, verbose = verbose, dark.theme = dark.theme, ... = ...)
  })

  tmp.file <- tempfile(pattern = "", fileext = ".png")

  colf <- ceiling(sqrt(length(x = dims)))
  colr <- round(length(x = dims)/colf)

  png(width = st.object@xdim*ncols*colf, height = st.object@xdim*nrows*colr, file = tmp.file)
  if (dark.theme) {
    par(mar = c(0, 0, 0, 0), bg = "black")
  } else {
    par(mar = c(0, 0, 0, 0))
  }
  plot(cowplot::plot_grid(plotlist = p.list, ncol = ncols))
  dev.off()

  p <- image_read(tmp.file)

  if (method == "viewer") {
    print(p)
    unlink(tmp.file)
  } else if (method == "raster") {
    if (dark.theme) {
      par(mar = c(0, 0, 0, 0), bg = "black")
    } else {
      par(mar = c(0, 0, 0, 0))
    }
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
#' @examples
#'
#' # Load images
#' se <- LoadImages(se)
#'
#' # Overlay the number of unique genes and the number of UMIs per spot on sample 1 HE image on the first two samples
#' MultiFeatureOverlay(se, features = c("nFeature_RNA", "nCount_RNA"), sampleids = 1:2)
#'
#' # Plot selected genes on the first two samples
#' MultiFeatureOverlay(se, features = c("Cck", "Dcn"), sampleids = 1:2)
#'
#' # Plot normalized values on the first two samples
#' se <- SCTransform(se)
#' MultiFeatureOverlay(se, features = c("Cck", "Dcn"), sampleids = 1:2)
#'
#' # Change to scaled data
#' MultiFeatureOverlay(se, features = c("Cck", "Dcn"), sampleids = 1:2, slot = "scale.data", center.zero = TRUE)
#'
#' # Mask images and plot plot the slected genes on the masked images for samples 1 and 2
#' se <- MaskImages(se)
#' MultiFeatureOverlay(se, features = c("Cck", "Dcn"), sampleids = 1:2, type = "masked")
#'
#' @export
#'
#' @seealso \code{\link{ST.FeaturePlot}} and \code{\link{ST.DimPlot}} for how to plot features
#' without the HE image, \code{\link{DimOverlay}} for how to overlay dimensionality reduction output on the
#' HE images and \code{\link{MultiDimOverlay}} for how to overlay dimensionality reduction plots on the HE
#' images in multiple samples.
#'

MultiFeatureOverlay <- function (
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
  pt.alpha = 1,
  shape.by = NULL,
  palette = NULL,
  cols = NULL,
  rev.cols = FALSE,
  center.zero = FALSE,
  channels.use = NULL,
  verbose = FALSE,
  dark.theme = FALSE,
  ...
) {

  # Check to see if Staffli object is present
  if (!"Staffli" %in% names(object@tools)) stop("Staffli object is missing from Seurat object. Cannot plot without coordinates", call. = FALSE)
  st.object <- object@tools$Staffli

  ncols <- ncols %||% ceiling(sqrt(length(x = sampleids)))
  nrows <- ceiling(length(x = sampleids)/ncols)

  p.list <- lapply(sampleids, function(s) {
    FeatureOverlay(object, features = features, sample.index = s, type = type,
                   min.cutoff = min.cutoff, max.cutoff = max.cutoff, slot = slot,
                   blend = blend, pt.size = pt.size, pt.alpha, shape.by = shape.by,
                   palette = palette, cols = cols, rev.cols = rev.cols,
                   grid.ncol = NULL, center.zero = center.zero,
                   channels.use = channels.use, verbose = verbose, dark.theme = dark.theme,... = ...)
  })

  tmp.file <- tempfile(pattern = "", fileext = ".png")

  colf <- ceiling(sqrt(length(x = features)))
  colr <- round(length(x = features)/colf)

  png(width = st.object@xdim*ncols*colf, height = st.object@xdim*nrows*colr, file = tmp.file)
  if (dark.theme) {
    par(mar = c(0, 0, 0, 0), bg = "black")
  } else {
    par(mar = c(0, 0, 0, 0))
  }
  plot(cowplot::plot_grid(plotlist = p.list, ncol = ncols))
  dev.off()

  p <- image_read(tmp.file)

  if (method == "viewer") {
    print(p)
    unlink(tmp.file)
  } else if (method == "raster") {
    if (dark.theme) {
      par(mar = c(0, 0, 0, 0), bg = "black")
    } else {
      par(mar = c(0, 0, 0, 0))
    }
    plot(p)
    unlink(tmp.file)
  } else {
    stop(paste0("Invalid method ", method), call. = FALSE)
  }
}

