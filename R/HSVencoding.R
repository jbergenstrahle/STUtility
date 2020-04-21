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
#' This visualization method should be used only on carefully selected features and you should be
#' aware that color representation of quantitative data can be very misleading. It should only be
#' usde to assess qualitative aspects of the data, for example if you wish to know where 5 "non-overlapping"
#' features are expressed spatially. You should therefore investigate beforehand if the features of interest
#' overlap or, otherwise the results can become very confusing.
#'
#' @section scaling of features:
#' All features are by default scaled independently to a 0 to 1 range which means that the relative
#' differencies between the feature expression levels is not preserved. This is because some features
#' can still be very distinct for a region of interest even though their magnitude of expression is low.
#' If you want to preserve the relative differencies you can set `rescale = FALSE`.
#'
#' @param object Seurat object
#' @param features
#' \itemize{
#'     \item An \code{Assay} feature (e.g. a gene name - "MS4A1")
#'     \item A column name from meta.data (e.g. mitochondrial percentage - "percent.mito")
#' }
#' @param plot.type Select one of 'spots' or 'smooth' [default: 'spots']
#' @param split.hsv Should the HSV colored features be split into separate plots? [default: FALSE]
#' @param rescale Rescale each feature column separately from 0 to 1 range. If set to FALSE, all feature columns
#' will be scaled together from 0 to 1 and preserve the relative differencies
#' @param indices Numeric vector specifying sample indices to include in plot. Default is to show all samples.
#' @param spots Vector of spots to plot (default is all spots)
#' @param min.cutoff,max.cutoff Vector of minimum and maximum cutoff values for each feature,
#'  may specify quantile in the form of 'q##' where '##' is the quantile (eg, 'q1', 'q10')
#' @param slot Which slot to pull expression data from?
#' @param pt.size Adjust point size for plotting
#' @param pt.alpha Adjust opacity of spots.
#' @param pt.border Should a border be drawn around the spots? [default: TRUE]
#' @param add.alpha Adds opacity to spots scaled by feature values. This will disable the pt.alpha parameter
#' @param shape.by If NULL, all points are circles (default). You can specify any spot attribute available in the meta.data slot
#' @param sigma Smoothing bandwidth; only active if \code{plot.type = 'smooth'}. A single positive number, a numeric vector of length 2, or a function that selects the bandwidth automatically [default: 2].
#' See \code{\link{density.ppp}} function from the \code{\link{spatstat}} package for more details.
#' @param highlight.edges Highlights the edges of the tissue. Only active if \code{plot.type = 'smooth'} and if the images have been masked.
#' @param grid.ncol Number of columns for display when combining plots
#' @param dark.theme Use a dark theme for plotting
#' @param theme Add a custom theme to the output ggplot object
#' @param scale.res Integer value setting the resolution of the output raster image. E.g. scale.res = 2 will double the
#' resolution of the output but will also take longer to render. Only active if plot.type is set to 'smooth'.
#' @param verbose Print messages
#' @param ... Extra parameters passed on to \code{\link{STPlot}}
#'
#' @inheritParams STPlot
#' @importFrom cowplot plot_grid
#' @importFrom scales rescale
#' @importFrom ggplot2 ggplot theme theme_void
#' @importFrom zeallot %<-%
#' @importFrom grDevices hsv
#' @importFrom spatstat ppp owin Smooth
#' @importFrom imager imgradient enorm as.cimg
#' @importFrom magick image_crop image_info image_read image_composite image_border image_scale
#'
#' @return A ggplot object
#' @export

HSVPlot <- function (
  object,
  features,
  ncol = NULL,
  plot.type = 'spots',
  split.hsv = FALSE,
  rescale = TRUE,
  indices = NULL,
  spots = NULL,
  min.cutoff = NA,
  max.cutoff = NA,
  slot = "data",
  pt.size = 1,
  pt.alpha = 1,
  pt.border = FALSE,
  add.alpha = FALSE,
  shape.by = NULL,
  sigma = 2,
  highlight.edges = FALSE,
  cols = NULL,
  dark.theme = TRUE,
  grid.ncol = NULL,
  theme = theme_void(),
  scale.res = 1,
  custom.theme = NULL,
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

  # Stop if feature classes are not numeric/integer
  if (!all(data.type %in% c("numeric", "integer"))) {
    stop("Only features of class 'integer' or 'numeric' are allowed ... ")
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
  if (ncol(x = data) < 3) {
    stop("None of the requested features were found: ",
         paste(features, collapse = ", "),
         " in slot ",
         slot,
         call. = FALSE)
  }

  data <- feature.scaler(data, features, min.cutoff, max.cutoff)

  # Subset by index
  if (!is.null(indices)) {
    if (!all(as.character(indices) %in% data[, "sample"])) stop(paste0("Index out of range. "), call. = FALSE)
    data <- data[data[, "sample"] %in% as.character(indices), ]
  } else {
    indices <- unique(data[, "sample"]) %>% as.numeric()
  }

  if (is.null(cols)) {
    # Generate HSV encoded colors
    if (verbose) cat(paste0("Defining Hue for ", length(x = features), " features ... \n"))
    hue_breaks <- seq(0, 1, length.out = length(x = features) + 1)[1:length(x = features)]
    hsv.matrix <- t(matrix(c(hue_breaks, rep(1, length(hue_breaks )), rep(1, length(hue_breaks))), ncol = 3))
    rownames(hsv.matrix) <- c("h", "s", "v")
    ann.cols <- apply(hsv.matrix, 2, function(x) hsv(x[1], x[2], x[3]))
  } else {
    if (length(x = features) != length(x = cols)) stop("Length of features and cols must match ...", call. = FALSE)
    warning("Using user defined colors with opacity. HSV scale will not be used ...", call. = FALSE)
    ann.cols <- cols
    names(cols) <- features
  }

  # Rescale data 0 to 1
  if (rescale) {
    data[, features] <- apply(data[, features], 2, scales::rescale)
  } else {
    data[, features] <- setNames(data.frame(scales::rescale(data[, features] %>% as.matrix() %>% as.numeric()) %>% matrix(ncol = length(x = features))), nm = features)
  }

  # Disable pt.alpha if add.alpha is provided
  if (add.alpha) pt.alpha <- NA

  # Plot HSV encoded feature data
  if (plot.type == 'spots') {
    # Select highest V
    # Select highest V
    d <- create.array.from.feature.vals(data, features, hue_breaks, cols, dark.theme, verbose)

    #red.cols <- data.frame()
    if (verbose) cat("Selecting HSV colors for each spot ... \n")
    data <- create.cols.from.array(data, d, features, cols, split.hsv, dark.theme, add.alpha)

    if (verbose) cat("Plotting features:",
                     ifelse(length(features) == 1, features,  paste0(paste(features[1:(length(features) - 1)], collapse = ", "), " and ", features[length(features)])))
    # Normal visualization -------------------------------------------------------------------------------------

    if (image.type != "empty") {
      dims <- lapply(st.object@dims, function(x) {x[2:3] %>% as.numeric()})
    } else {
      dims <- st.object@limits
    }

    if (!is.null(indices)) dims <- dims[indices]

    # Plot combined HSV
    if (!split.hsv) {
      plot <- STPlot(data, data.type, shape.by, NULL, pt.size, pt.alpha, pt.border = pt.border,
                     palette = "Reds", cols = NULL, ncol = ncol, spot.colors = data$cols,
                     center.zero = F, center.tissue = F, plot.title = "",
                     dims = dims, split.labels = FALSE, dark.theme = dark.theme,
                     pxum = NULL, sb.size = 2.5, custom.theme = custom.theme, ...)
      if (dark.theme) {
        plot <- plot + dark_theme()
      }
      plot <- plot +
        geom_point(data = data.frame(x = rep(-1, length(features)), y = rep(-1, length(features)), features), aes(x, y, colour = features)) +
        scale_color_manual(values = setNames(ann.cols, features))
      return(plot)
    } else {
      plots <- lapply(seq_along(data), function (i) {
        data <- data[[i]]
        plot <- STPlot(data, data.type, shape.by, NULL, pt.size, pt.alpha, pt.border = pt.border,
               palette = NULL, cols = NULL, ncol = ncol, spot.colors = data$cols,
               center.zero = F, center.tissue = F, plot.title = features[i],
               dims = dims, split.labels = FALSE, dark.theme = dark.theme,
               pxum = NULL, sb.size = 2.5, custom.theme = custom.theme, ...)
        if (dark.theme) {
          plot <- plot + dark_theme()
        }
        return(plot)
      })
      ncols <- grid.ncol %||% ceiling(sqrt(length(x = features)))
      nrows <- ceiling(length(x = features)/ncols)
      plot <- cowplot::plot_grid(plotlist = plots, ncol = ncols, nrow = nrows)
      if (dark.theme) plot <- plot + dark_theme()
      return(plot)
    }
  } else if (plot.type == 'smooth') {
    feature.list <- list()
    edges.list <- list()
    for (ftr in features) {
      val.limits <- range(data[, ftr])
      p.list <- list()
      for (i in 1:length(unique(data$sample))) {
        data_subset <- subset(data, sample == i)
        dims <- st.object@rasterlists$processed.masks[[i]] %>% dim()
        if (image.type %in% c('raw', 'masked', 'processed')) {
          extents <- st.object@dims[[i]][2:3] %>% as.numeric()
          data_subset[, c("x", "y")] <- data_subset[, c("x", "y")]/((extents[1]/scale.res)/dims[2])
        } else {
          extents <- st.object@limits[[i]]
          data_subset[, c("x", "y")] <- data_subset[, c("x", "y")]/((extents[1]/scale.res)*scale.res/dims[2])
        }

        ow <- owin(xrange = c(0, dims[2]*scale.res), yrange = c(0, dims[1]*scale.res))
        p <- ppp(x = data_subset[, "x"], y = data_subset[, "y"], window = ow, marks = data_subset[, ftr])
        suppressWarnings({s <- Smooth(p, sigma*scale.res, dimyx = dims*scale.res)})
        m <-  as.matrix(s)
        m[m < 0] <- 0
        m <- m/max(m)

        if (image.type %in% c('processed', 'masked')) {
          msk.type <- paste0(image.type, ".masks")
          msk <- st.object['processed.masks'][[i]]
          if (scale.res != 1) {
            msk <- image_read(msk) %>% image_scale(paste0(st.object@xdim*scale.res)) %>% magick2cimg()
          } else {
            msk <- msk %>% as.cimg()
          }
          if (highlight.edges) {
            edges.list[[i]] <- imgradient(msk, "xy") %>% enorm()
          }
          msk <- msk[, , , 1] %>% as.cimg() %>% threshold()
          m <- t(m) %>% as.cimg()
          masked.m <- m*msk
          p.list[[i]] <- masked.m
        } else {
          p.list[[i]] <- m %>% as.cimg()
        }
      }
      feature.list[[ftr]] <- p.list
    }

    # HSV plot
    hue_breaks <- seq(0, 1, length.out = length(x = features) + 1)[1:length(x = features)]
    rsts <- list()

    if (!split.hsv) {
      for (j in 1:length(unique(data[, "sample"]))) {
        ar <- array(dim = c(rev(dims*scale.res), length(features)))
        n <- 1
        for (i in features) {
          ar[, , n] <- feature.list[[i]][[j]]
          n <- n + 1
        }
        ftr.rst <- apply(ar[, , ], c(1, 2), function(x) {
          if (is.null(cols)) {
            hsvc <- hsv(h = hue_breaks[which.max(x)], s = ifelse(dark.theme, 1, max(x)), v = ifelse(dark.theme, max(x), 1))
          } else {
            hsvc <- cols[which.max(x)]
          }
          if (add.alpha) hsvc <- scales::alpha(hsvc, max(x))
          return(hsvc)
        }) %>% t() %>% as.raster() #%>% as.cimg()

        if (length(edges.list) > 0) {
          ftr.rst[t((edges.list[[j]] > 0)[, , , 1])] <- "#FFFFFF"
        }
        rsts[[j]] <- ftr.rst %>% as.raster()
      }
      ncols <- length(unique(data[, "sample"]))
      nrows <- ceiling(length(unique(data[, "sample"]))/ncols)
    } else {
      for (j in 1:length(unique(data[, "sample"]))) {
        feature.rsts <- list()
        for (i in seq_along(features)) {
          ftr.rst <- sapply(feature.list[[i]][[j]], function(x) {
            if (is.null(cols)) {
              hsvc <- hsv(h = hue_breaks[i], s = ifelse(dark.theme, 1, x), v = ifelse(dark.theme, x, 1))
            } else {
              hsvc <- cols[i]
            }
            if (add.alpha) hsvc <- scales::alpha(hsvc, x)
            return(hsvc)
          }) %>% matrix(nrow = dims[2]*scale.res, ncol = dims[1]*scale.res) %>% t() %>% as.raster() #%>% as.cimg()
          if (length(edges.list) > 0) {
            ftr.rst[t((edges.list[[j]] > 0)[, , , 1])] <- "#FFFFFF"
          }
          feature.rsts[[i]] <-  ftr.rst %>% as.raster()
        }
        rsts[[j]] <- feature.rsts
      }
      rsts <- Reduce(c, rsts)
      # rearrange results
      reord <- rep(seq_along(features), each = 2)
      reord[seq(2, length(reord), 2)] <- reord[seq(2, length(reord), 2)] + length(x = features)
      rsts <- rsts[reord]
      ncols <- length(unique(data[, "sample"]))
      nrows <- length(x = features)
    }

    rsts <- lapply(seq_along(rsts), function(i) {
      im <- rsts[[i]]
      im <- im %>% image_read()
      im <- image_border(im, ifelse(dark.theme, "#000000", "#FFFFFF"), paste(st.object@xdim*scale.res/10, st.object@xdim*scale.res/10, sep = "x"))
      im <- image_annotate(im, text = i, size = round(st.object@xdim/10), color = ifelse(dark.theme, "#FFFFFF", "#000000"))
    })

    tmp.file <- tempfile(pattern = "", fileext = ".png")
    png(width = st.object@xdim*ncols*scale.res, height = st.object@xdim*nrows*scale.res, file = tmp.file)
    par(mfrow = c(nrows, ncols), mar = c(0, 0, 0, 0), bg = ifelse(dark.theme, "black", "white"))
    for (rst in rsts) {
      plot(rst)
    }
    dev.off()

    im <- image_read(tmp.file)
    if (!split.hsv) {
      im <- image_border(im, ifelse(dark.theme, "#000000", "#FFFFFF"), paste0(st.object@xdim*scale.res/2))
    } else {
      im <- image_border(im, ifelse(dark.theme, "#000000", "#FFFFFF"), paste0(st.object@xdim*scale.res))
    }
    tmp.file <- tempfile(pattern = "", fileext = ".png")
    lg <- g_legend(data.frame(x = 1, y = 1, feature = features), data.type = "character", variable = "feature", center.zero = FALSE, cols = ann.cols, val.limits = NULL, dark.theme = dark.theme)

    grobHeight <- function(x) {
      grid::convertHeight(sum(x$heights), "in", TRUE)
    }

    grobWidth <- function(x) {
      grid::convertWidth(sum(x$widths), "in", TRUE)
    }

    ggsave(plot = lg, width = grobWidth(lg), height = grobHeight(lg), filename = tmp.file)
    iminf <- image_info(im)[2:3] %>% as.numeric()
    if (!split.hsv) {
      lgim <- image_read(tmp.file) %>% image_scale(paste0(iminf[2]/5))
    } else {
      lgim <- image_read(tmp.file) %>% image_scale(paste0(iminf[2]/(nrows*2)))
    }
    iminf.lgm <- image_info(lgim)[2:3] %>% as.numeric()
    lgim <- image_crop(lgim, paste0(iminf.lgm[1] - 2, "x", iminf.lgm[2] - 2, "x", 1, "x", 1))
    if (!split.hsv) {
      im <- image_composite(image = im, composite_image = lgim, offset = paste0("+", iminf[1] - st.object@xdim*scale.res/length(features), "+", (iminf[2])/2 - (iminf.lgm[2])/2))
    } else {
      im <- image_composite(image = im, composite_image = lgim, offset = paste0("+", st.object@xdim*ncols*scale.res*1.5, "+", (iminf[2])/2 - (iminf.lgm[2])/2))
    }
    par(mar = c(0, 0, 0, 0), bg = ifelse(dark.theme, "black", "white"))
    plot(im %>% as.raster())
  }
}




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# HSV plots on HE images
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Overlay HSVplot on one selected HE image
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
#' @param ... Extra parameters passed on to \code{\link{ST.ImagePlot}}
#'
#' @inheritParams ST.ImagePlot
#' @inheritParams ST.FeaturePlot
#' @inheritParams HSVPlot
#' @importFrom cowplot plot_grid
#'
#' @return A ggplot object
#'

spatial_hsv_plot <- function (
  object,
  features,
  split.hsv = FALSE,
  sample.index = 1,
  rescale = TRUE,
  spots = NULL,
  type = NULL,
  min.cutoff = NA,
  max.cutoff = NA,
  slot = "data",
  pt.size = 2,
  pt.alpha = 1,
  pt.border = FALSE,
  add.alpha = FALSE,
  shape.by = NULL,
  palette = NULL,
  cols = NULL,
  grid.ncol = NULL,
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

  if (is.null(cols)) {
    # Generate HSV encoded colors
    if (verbose) cat(paste0("Defining Hue for ", length(x = features), " features ... \n"))
    hue_breaks <- seq(0, 1, length.out = length(x = features) + 1)[1:length(x = features)]
    hsv.matrix <- t(matrix(c(hue_breaks, rep(1, length(hue_breaks )), rep(1, length(hue_breaks))), ncol = 3))
    rownames(hsv.matrix) <- c("h", "s", "v")
    ann.cols <- apply(hsv.matrix, 2, function(x) hsv(x[1], x[2], x[3]))
  } else {
    if (length(x = features) != length(x = cols)) stop("Length of features and cols must match ...", call. = FALSE)
    warning("Using user defined colors with opacity. HSV scale will not be used ...", call. = FALSE)
    ann.cols <- cols
    names(cols) <- features
  }

  # Rescale data 0 to 1
  # Add dummy data
  if (is.list(value.scale)) {
    data <- rbind(data, setNames(data.frame(cbind(do.call(cbind, value.scale), matrix(NA, ncol = sum(!colnames(data) %in% features), nrow = 2))), nm = colnames(data)))
  }
  if (rescale) {
    data[, features] <- apply(data[, features], 2, scales::rescale)
  } else {
    data[, features] <- setNames(data.frame(scales::rescale(data[, features] %>% as.matrix() %>% as.numeric()) %>% matrix(ncol = length(x = features))), nm = features)
  }
  data <- na.omit(data)

  # Disable pt.alpha if add.alpha is provided
  if (add.alpha) pt.alpha <- NA

  if (verbose) cat("Plotting features:",
                   ifelse(length(features) == 1, features,  paste0(paste(features[1:(length(features) - 1)], collapse = ", "), " and ", features[length(features)])))

  # Select highest V
  d <- create.array.from.feature.vals(data, features, hue_breaks, cols, dark.theme, verbose)

  #red.cols <- data.frame()
  if (verbose) cat("Selecting HSV colors for each spot ... \n")
  data <- create.cols.from.array(data, d, features, cols, split.hsv, dark.theme, add.alpha)


  # Plot combined HSV
  if (!split.hsv) {
    plot <- ST.ImagePlot(data, data.type, shape.by, NULL, image, dims = imdims,
                         pt.size, pt.alpha, pt.border = pt.border, FALSE, palette,
                         cols, NULL, spot.colors = data$cols,
                         FALSE, plot.title = "", FALSE, dark.theme,
                         pixels.per.um, NULL, custom.theme = custom.theme, ...)
    plot <- plot +
      geom_point(data = data.frame(x = rep(-1, length(features)), y = rep(-1, length(features)), features), aes(x, y, colour = features)) +
      scale_color_manual(values = setNames(ann.cols, features)) +
      theme_void()
    if (dark.theme) {
      plot <- plot + dark_theme()
    }
    return(plot)
  } else {
    plots <- lapply(seq_along(data), function (i) {
      data <- data[[i]]
      plot <- ST.ImagePlot(data, data.type, shape.by, NULL, image, dims = imdims,
                           pt.size, pt.alpha, pt.border = pt.border, add.alpha = FALSE, palette,
                           cols, NULL, spot.colors = data$cols,
                           FALSE, plot.title = features[i], FALSE, dark.theme,
                           pixels.per.um, NULL, custom.theme = custom.theme, ...)
      if (dark.theme) {
        plot <- plot + dark_theme()
      }
      return(plot)
    })
    ncols <- grid.ncol %||% ceiling(sqrt(length(x = features)))
    nrows <- ceiling(length(x = features)/ncols)
    plot <- cowplot::plot_grid(plotlist = plots, ncol = ncols, nrow = nrows)
    if (dark.theme) plot <- plot + dark_theme()
    return(plot)
  }
}


#' Overlay HSV encoded features on HE images
#'
#' Graphs the selected features as a HSVplot on a 2D grid of spots overlaid on top of an HE images.
#' Only numerical features are accepted, e.g. genes or dimensionality reduction output vectors. If you
#' want to draw dimentionality reduction vectors you need to specify the whole names of the vectors, e.g.
#' `features = c("factor_1", "factor_2")` for the two first NMF factors.
#'
#' NOTE that this function draws sample 1 as default, but can take multiple samples as well using the `sampleids argument`.
#'
#' @details It is typically difficult to explore details in the HE image when diplaying multiple samples side by side,
#' so we recommend to draw the plots for one sample at the time. If you have higher resolution images,
#' it could also take significant time to draw the plots.
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
#' @param show.sb Should a scalebar be drawn? [default: TRUE]
#' @param ... Parameters passed to DimOverlay
#'
#' @inheritParams spatial_hsv_plot
#' @inheritParams HSVPlot
#'
#' @examples
#' # Load images
#' se <- se %>% SCTransfrom() %>% LoadImages() %>% RunNMF()
#'
#' # Overlay first two NMF factors on the first two tissue sections
#' HSVPlot(se, features = c("factor_1", "factor_2"), sampleids = 1:2)
#'
#' @export
#'

HSVOverlay <- function (
  object,
  features,
  sampleids = 1,
  rescale = TRUE,
  spots = NULL,
  ncols.features = NULL,
  ncols.samples = NULL,
  type = NULL,
  min.cutoff = NA,
  max.cutoff = NA,
  slot = "data",
  pt.size = 2,
  pt.alpha = 1,
  add.alpha = FALSE,
  shape.by = NULL,
  palette = NULL,
  cols = NULL,
  split.hsv = FALSE,
  dark.theme = FALSE,
  sample.label = TRUE,
  show.sb = TRUE,
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

  ncols.features <- ncols.features %||% length(x = features)
  ncols.samples <- ncols.samples %||% 1

  data <- FetchData(object = object, vars = c(features), cells = spots, slot = slot)
  value.type <- sapply(data, class)
  if (any(!value.type %in% "numeric")) stop("Only numeric features can be plotted with HSVOverlay. \n", call. = FALSE)
  value.scale.list <- lapply(data, range)

  p.list <- lapply(remaining_samples, function(s) {
    spatial_hsv_plot(object = object, features = features, split.hsv = split.hsv,
                     sample.index = s, rescale = rescale, spots = spots, type = type,
                     min.cutoff = min.cutoff, max.cutoff = max.cutoff, slot = slot,
                     pt.size = pt.size, pt.alpha, pt.border = FALSE, add.alpha = add.alpha, shape.by = shape.by,
                     palette = palette, cols = cols, grid.ncol = ncols.features,
                     dark.theme = dark.theme, sample.label = sample.label, show.sb = show.sb,
                     value.scale = value.scale.list, custom.theme = custom.theme, verbose = verbose)#, ... = ...)
  })

  p <- cowplot::plot_grid(plotlist = p.list, ncol = ncols.samples)
  if (dark.theme) p <- p + dark_theme()
  return(p)
}


#' Creates an array of dimensions number_of_spots*3*number_of_features
#'
#' For each feature, a matrix is stored with nSpots number of rows and
#' with the HSV color channels as columns. If dark.theme is set to TRUE,
#' the V channel will be reserved for feature values and the S channel will
#' be set to 1, otherwise the S channel will be resevred for feature values
#' and the V channel will be set to 1.
#'
#' @param data data.frame with feature values
#' @param features feature names
#' @param hue_breaks Hue values (same length as features)
#' @param cols Custom colors
#' @param dark.theme Used to select what channel the feature values should be encoded in
#' @param verbose Print messages

create.array.from.feature.vals <- function (
  data,
  features,
  hue_breaks,
  cols,
  dark.theme,
  verbose
) {
  if (is.null(cols)) {
    d <- array(dim = c(nrow(data), 3, length(x = features)))
    if (verbose) cat("Converting values to HSV colors ... \n")
    for (i in 1:length(features)) {
      ftr <- features[i]
      if (dark.theme) {
        s <- data.frame(h = hue_breaks[i],
                        s = 1,
                        v = data[, ftr, drop = T] %>% as.numeric()) %>% as.matrix()
      } else {
        s <- data.frame(h = hue_breaks[i],
                        s = data[, ftr, drop = T] %>% as.numeric(),
                        v = 1) %>% as.matrix()
      }
      d[, , i] <- s
    }
  } else {
    d <- array(dim = c(nrow(data), 1, length(x = features)))
    if (verbose) cat("Using provided colors ... \n")
    for (i in 1:length(features)) {
      ftr <- features[i]
      s <- data.frame(v = data[, ftr, drop = T] %>% as.numeric()) %>% as.matrix()
      d[, , i] <- s
    }
  }
  return(d)
}


#' Creates HSV colors from an array
#'
#' If split.hsv = FALSE, the feature with the highest value in a spot will define the
#' color for that spot. The intensity of the color will depend on if dark.theme is active and
#' the magnitude of the feature value in that spot.
#'
#' @param data data.frame with feature values
#' @param d array created with \code{create.array.from.vals} function
#' @param features Feature names
#' @param cols Custom colors
#' @param split.hsv Should the features be plotted separately?
#' @param dark.theme Used to select what channel the feature values should be encoded in
#' @param add.alpha Adds opacity to the output colors, defined by the scaled feature values

create.cols.from.array <- function (
  data,
  d,
  features,
  cols,
  split.hsv,
  dark.theme,
  add.alpha
) {
  # If split.hsv is deactivated, get return one data.frame
  if (!split.hsv) {
    if (is.null(cols)) {
      red.cols <- apply(d, 1, function (x) {
        ind <- ifelse(dark.theme, 3, 2)
        max.val <- which.max(x[ind, ])
        hsvc <- hsv(h = x[1, ][which.max(x[ind, ])], s = ifelse(dark.theme, 1, max(x[ind, ])), v = ifelse(dark.theme, max(x[ind, ]), 1))
        if (add.alpha) hsvc <- scales::alpha(hsvc, max(x[ind, ]))
        return(hsvc)
      })
    } else {
      red.cols <- unlist(apply(d, 1, function (x) {
        alpha_col <- cols[which.max(x[1, ])]
        if (add.alpha) {
          alpha_col <- scales::alpha(colour = alpha_col, alpha = max(x[1, ]))
        }
        return(alpha_col)
      }))
    }
    data$cols <- red.cols
    return(data)
  } else {
    # If split.hsv is activated, get return one data.frame
    if (is.null(cols)) {
      full.data <- matrix(ncol = ncol(data), nrow = 0)
      for (i in 1:dim(d)[3]) {
        full.data <- rbind(full.data, cbind(data, setNames(data.frame(d[, , i]), nm = c("h", "s", "v")), variable = features[i]))
      }
      red.cols <- apply(full.data, 1, function (x) {
        hsvc <- hsv(h = x["h"], s = x["s"], v = x["v"])
        if (add.alpha) hsvc <- scales::alpha(hsvc, ifelse(dark.theme, x["v"], x["s"]))
        return(hsvc)
      })
    } else {
      full.data <- matrix(ncol = ncol(data), nrow = 0)
      for (i in 1:dim(d)[3]) {
        full.data <- rbind(full.data, cbind(data, setNames(data.frame(d[, , i]), nm = c("v")), variable = features[i]))
      }
      red.cols <- apply(full.data, 1, function (x) {
        hsvc <- cols[x["variable"]]
        if (add.alpha) hsvc <- scales::alpha(hsvc, ifelse(dark.theme, x["v"], x["s"]))
        return(hsvc)
      })
    }
    full.data$cols <- red.cols
    full.data.split <- split(full.data, full.data$variable)
    return(full.data.split)
  }
}

