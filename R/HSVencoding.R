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
  add.alpha = FALSE,
  shape.by = NULL,
  sigma = 2,
  highlight.edges = FALSE,
  cols = NULL,
  dark.theme = TRUE,
  grid.ncol = NULL,
  verbose = FALSE,
  theme = theme_void(),
  scale.res = 1,
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
    if (is.null(cols)) {
      d <- array(dim = c(nrow(data), 3, length(x = features)))
      if (verbose) cat("Converting values to HSV ... \n")
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

    #red.cols <- data.frame()
    if (verbose) cat("Selecting HSV colors for each spot ... \n")
    if (!split.hsv) {
      if (is.null(cols)) {
        red.cols <- apply(d, 1, function (x) {
          ind <- ifelse(dark.theme, 3, 2)
          max.val <- which.max(x[ind, ])
          hsvc <- hsv(h = x[1, ][which.max(x[ind, ])], s = ifelse(dark.theme, 1, max(x[ind, ])), v = ifelse(dark.theme,  max(x[ind, ]), 1))
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
    } else {
      if (is.null(cols)) {
        full.data <- matrix(ncol = ncol(data), nrow = 0)
        for (i in 1:dim(d)[3]) {
          full.data <- rbind(full.data, cbind(data, setNames(data.frame(d[, , i]), nm = c("h", "s", "v")), variable = features[i]))
        }
        red.cols <- apply(full.data, 1, function (x) {
          hsvc <- hsv(h = x["h"], s = x["s"], v = x["v"])
          if (add.alpha) hsvc <- scales::alpha(hsvc, x["v"])
          return(hsvc)
        })
      } else {
        full.data <- matrix(ncol = ncol(data), nrow = 0)
        for (i in 1:dim(d)[3]) {
          full.data <- rbind(full.data, cbind(data, setNames(data.frame(d[, , i]), nm = c("v")), variable = features[i]))
        }
        red.cols <- apply(full.data, 1, function (x) {
          hsvc <- cols[x["variable"]]
          if (add.alpha) hsvc <- scales::alpha(hsvc, x["v"])
          return(hsvc)
        })
      }

      full.data$cols <- red.cols
      full.data.split <- split(full.data, full.data$variable)
    }

    if (verbose) cat("Plotting features:",
                     ifelse(length(features) == 1, features,  paste0(paste(features[1:(length(features) - 1)], collapse = ", "), " and ", features[length(features)])))
    # Normal visualization -------------------------------------------------------------------------------------

    if (image.type != "empty") {
      dims <- lapply(iminfo(st.object), function(x) {x[2:3] %>% as.numeric()})
    } else {
      dims <- st.object@limits
    }

    if (!is.null(indices)) dims <- dims[indices]

    # Plot combined HSV
    if (!split.hsv) {
      plot <- STPlot(data, data.type, shape.by, NULL, pt.size, pt.alpha,
                     palette = NULL, cols = NULL, ncol, spot.colors = data$cols,
                     center.zero = F, center.tissue = F, plot.title = "",
                     dims, FALSE, theme = theme, dark.theme = dark.theme,
                     pxum = NULL, sb.size = 2.5, ...)
      if (dark.theme) {
        plot <- plot + dark_theme()
      }
      plot <- plot +
        geom_point(data = data.frame(x = rep(-1, length(features)), y = rep(-1, length(features)), features), aes(x, y, colour = features)) +
        scale_color_manual(values = setNames(ann.cols, features))
      suppressWarnings({print(plot)})
    } else {
      plots <- lapply(full.data.split, function (data) {
        plot <- STPlot(data, data.type, shape.by, NULL, pt.size, pt.alpha,
               palette = NULL, cols = NULL, ncol, spot.colors = data$cols,
               center.zero = F, center.tissue = F, plot.title = "",
               dims, FALSE, theme = theme, dark.theme = dark.theme,
               pxum = NULL, sb.size = 2.5, ...)
        if (dark.theme) {
          plot <- plot + dark_theme()
        }
        plot <- plot +
          geom_point(data = data.frame(x = rep(-1, length(features)), y = rep(-1, length(features)), features), aes(x, y, colour = features)) +
          scale_color_manual(values = setNames(ann.cols, features))
      })
      ncols <- grid.ncol %||% ceiling(sqrt(length(x = features)))
      nrows <- ceiling(length(x = features)/ncols)
      plot <- cowplot::plot_grid(plotlist = plots, ncol = ncols, nrow = nrows)
      if (dark.theme) plot <- plot + dark_theme()
      suppressWarnings({print(plot)})
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

