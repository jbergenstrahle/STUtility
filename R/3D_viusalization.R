#' Create 3D stack
#'
#' Creates a stack of point patterns with aligned data and stores it in a
#' `Staffli` object.
#'
#' @param object Seurat object
#' @param limit Cut-off threshold for segmentation of points. Has to be a value between 0-1.
#' @param maxnum Maximum number of points to store for each section [default: 5e5]. If you are analyzing many samples
#' you will have to decrease tsis number, oterwise tere will too many points to draw.
#' @param verbose Print messages
#'
#' @inheritParams rasterize_scatter
#'
#' @return Seurat object
#'
#' @importFrom akima interp
#' @importFrom raster raster
#'
#' @export

Create3DStack <- function (
  object,
  limit = 0.4,
  nx = 200,
  verbose = TRUE
) {

  # Check limit
  if (!0 < limit & limit < 1) stop("limit has to been in the 0-1 range ... \n", call. = FALSE)

  # Check to see if Staffli object is present
  if (!"Staffli" %in% names(object@tools)) stop("Staffli object is missing from Seurat object... \n", call. = FALSE)
  st.object <- GetStaffli(object)

  # Check if processed images are present
  if (!"processed" %in% rasterlists(st.object)) stop("No processed images available in Seurat object ... \n", call. = FALSE)

  # Switch resolution
  if (st.object@xdim < 1500) {
    object <- SwitchResolution(object, xdim = 1500, verbose = verbose)
    st.object <- GetStaffli(object)
  }

  # Obtain number of cells to interpolate over
  coords <- do.call(rbind, lapply(seq_along(st.object@samplenames), function(i) {
    s <- st.object@samplenames[i]
    coords <- subset(st.object[[]], sample == s)[, c("warped_x", "warped_y")]
    dims.raw <- iminfo(st.object)[[i]][2:3] %>% as.numeric()
    dims.scaled <- scaled.imdims(st.object)[[i]]
    sf.xy <- dims.raw[2]/dims.scaled[1]
    coords <- coords/sf.xy
    coords$z <- i
    return(coords)
  }))

  if (verbose) cat("Running approximative segmentation of nuclei ... \n")
  scatters <- do.call(rbind, lapply(seq_along(st.object@samplenames), function(i) {
    s <- st.object@samplenames[i]
    scatter <- scatter_HE(st.object, type = "processed", sample.index = s, maxnum = maxnum, limit = limit, edges = FALSE)
    scatter$z <- i
    return(scatter)
  }))

  max.x = max(coords[, 1])
  min.x = min(coords[, 1])
  max.y = max(coords[, 2])
  min.y = min(coords[, 2])
  res = (max.x - min.x) / nx
  r = raster(xmn = min.x, ymn = min.y, xmx = max.x, ymx = max.y, res = res)
  r[] = 0

  # Obtain grid cells
  if (verbose) cat("Assigning nuclei positions to grid cells ... \n")
  section.input <- rasterize_scatter(scatters, r, nx)

  st.object@scatter.data <- section.input
  object@tools$Staffli <- st.object
  return(object)
}

# TODO: fix cutoff thresholds

#' Plots the values of a feature in 3D
#'
#' @param object Seurat object
#' @param e.anterior <- spots Subset spots to plot
#' @param dims Dimensions to plot
#' @param reduction Reduction object to pull data from
#' @param slot Which slot to pull the data from? [default: 'data']
#' @param zcoords Vector of z coordinates with the same length as the number of sections in the dataset [default: 1:#sections]
#' @param min.cutoff,max.cutoff Vector of minimum and maximum cutoff values for each feature, may specify quantile in the form of 'q##' where '##'
#' is the quantile (eg, 'q1', 'q10'). This can be useful if you have outlier values that skew the colorscale in the plot. For example, if you specify
#' 'q1', you will trim of values below the 1st percentile. [default: no cuttoffs]
#' @param blend Scale and blend expression values to visualize coexpression of two dims (this options will override other coloring parameters).
#' See 'Blending values' below for a more thourough description.
#' @param pt.size Sets the size of points in the 3D plot
#' @param pt.alpha Sets the opacity of the points
#' @param cols Colors used to create a colorscale
#' @param add.margins Add margins along z axis to push sections closer to each other
#' @param channels.use Color channels to use for blending. Has to be a character vector of length 2 or 3 with "red", "green" and "blue"
#' color names specified [default: c("red", "green", "blue)]
#' @param scene Give the scene a name to allow for multiple subplots
#' @param return.data return the data.frame with x,y coordinates and interpolated values
#' instead of plotting
#' @param dark.theme Draws the plot with a dark theme
#'
#' @importFrom plotly plot_ly add_markers layout
#'
#' @export

DimPlot3D <- function (
  object,
  spots = NULL,
  dims = 1,
  reduction = NULL,
  mode = c("cloud", "spots"),
  zcoords = NULL,
  slot = 'data',
  min.cutoff = NA,
  max.cutoff = NA,
  blend = FALSE,
  pt.size = NULL,
  pt.alpha = 1,
  cols = NULL,
  add.margins = 0,
  channels.use = NULL,
  scene = "scene1",
  return.data = FALSE,
  dark.theme = FALSE,
  verbose = FALSE
) {

  # Check to see if Staffli object is present
  if (!"Staffli" %in% names(object@tools)) stop("Staffli object is missing from Seurat object ... \n", call. = FALSE)
  st.object <- GetStaffli(object)

  # Set mode
  mode <- match.arg(mode, c("cloud", "spots"))

  # Set point size if not provides
  pt.size <- pt.size %||% ifelse(mode == "cloud", 0.8, 5)

  # Check to see if Staffli object is present
  if (mode == "cloud") {
    if (!length(x = st.object@scatter.data) > 0) stop("3D stack is missing. Run Create3DStack() first ... \n", call. = FALSE)
    scatter.data <- st.object@scatter.data
  }

  reduction <- reduction %||% {
    default.reductions <- c('umap', 'tsne', 'pca', 'ica')
    object.reductions <- FilterObjects(object = object, classes.keep = 'DimReduc')
    reduc.use <- min(which(x = default.reductions %in% object.reductions))
    default.reductions[reduc.use]
  }

  if (!is.null(spots) & mode == "cloud") stop("Subsetting of spots only works if mode is set to 'spots' ... \n", call. = FALSE)
  spots <- spots %||% colnames(object)

  # prepare data
  signs <- sign(dims); dims <- abs(dims)
  #spots <- spots %||% colnames(x = object)
  data <- Embeddings(object = object[[reduction]])[spots, dims, drop = FALSE]
  if (verbose) cat(paste0("Selected ", length(spots), " spots"))
  data <- as.data.frame(x = t(t(data)*signs))
  dims <- paste0(Key(object = object[[reduction]]), dims)

  coords <- do.call(rbind, lapply(seq_along(st.object@samplenames), function(i) {
    s <- st.object@samplenames[i]
    dims.raw <- as.numeric(st.object@dims[[i]][2:3])
    dims.scaled <- dim(st.object["raw"][[i]])
    sf.xy <- dims.raw[1]/dims.scaled[2]
    coords <- subset(st.object[[]], sample == s)[, c("warped_x", "warped_y")]/sf.xy
    coords$z <- i
    return(coords)
  }))
  coords <- coords[spots, ]

  # Use zcoords to convert z axis
  def_z <- unique(coords$z)
  if (!is.null(zcoords)) {
    if (length(x = def_z) != length(x = zcoords)) stop(paste0("zcoords (length = ", length(zcoords), ") has to be the same length as the number of samples ", length(def_z), " in the data"), call. = FALSE)
    coords$z <- zcoords[coords$z]
    if (mode == "cloud") scatter.data$z <- zcoords[scatter.data$z]
  }

  # Add values
  data <- cbind(coords, data)

  # Scale values
  data <- feature.scaler(data = data, features = dims, min.cutoff = min.cutoff, max.cutoff = max.cutoff)

  # Split data
  if (mode == "cloud") {
    data.list <- split(data, data$z)
    section.input.list <- split(scatter.data, scatter.data$z)
    nxy <- colnames(scatter.data)[1:2] %>% as.numeric()
  }

  xmax <- lapply(st.object@dims, function(d) {d[2] %>% as.numeric()}) %>% unlist() %>% max()

  # Check if blend option is set
  if (blend) {
    if (!length(dims) %in% c(2, 3)) {
      stop("Color blending only works for 2 or 3 dimensions ... \n", call. = FALSE)
    } else {
      channels.use <- channels.use %||% c("red", "green", "blue")[1:length(dims)]
      if (verbose) cat(paste0("Blending colors for dimensions ",
                              paste0(ifelse(length(dims) == 2, paste0(dims[1], " and ", dims[2]), paste0(dims[1], dims[2], " and ", dims[2]))),
                              ": \n", paste(paste(dims, channels.use, sep = ":"), collapse = "\n")))
      if (mode == "cloud") {
        interpolated.data <- do.call(rbind, lapply(seq_along(data.list), function(i) {
          data <- data.list[[i]]
          section.input <- section.input.list[[i]]
          interpolated_data.list <- list()
          for (i in 1:length(dims)) {
            df <- interpolate_2D_data(setNames(data[, c(1:3, i + 3)], nm = c("warped_x", "warped_y", "z", "value")), section.input, nx = nxy[1], ny = nxy[2], xy.range = apply(coords[, 1:2], 2, range))
            colnames(df) <- c("warped_x", "warped_y", "z", dims[i])
            interpolated_data.list[[i]] <- df
          }
          A <- interpolated_data.list[[1]]
          for (i in 2:length(dims)) {
            A <- setNames(cbind(A, interpolated_data.list[[i]][, 4]), nm = c(colnames(A), dims[i]))
          }
          colored.data <- apply(A[, dims], 2, scales::rescale)
          colored.data[is.na(colored.data)] <- 0
          spot.colors <- ColorBlender(colored.data, channels.use)
          return(setNames(data.frame(A[, 1:3], spot.colors), nm = c("x", "y", "z", "spot.colors")))
        }))

        p <- plot_ly(interpolated.data,
                     scene = scene,
                     x = ~xmax - x, y = ~y, z = ~z,
                     marker = list(color = interpolated.data$spot.colors,
                                   showscale = FALSE,
                                   size = pt.size,
                                   opacity = pt.alpha)) %>%
          add_markers() %>%
          layout(title = paste(paste(channels.use, dims, sep = ": "), collapse = "; "), paper_bgcolor = ifelse(dark.theme, 'rgb(0, 0, 0)', 'rgb(255, 255, 255)'),
                 scene = list(zaxis = list(title = 'z', range = c(-add.margins, max(interpolated.data$z) + add.margins), showticks = FALSE)))
        names(p$x$layoutAttrs[[1]]) <- c("title", "paper_bgcolor", scene)
        return(p)
      } else {
        data <- setNames(data, nm = c("x", "y", "z", dims))
        colored.data <- apply(data[, dims], 2, scales::rescale)
        colored.data[is.na(colored.data)] <- 0
        spot.colors <- ColorBlender(colored.data, channels.use)
        data$spot.colors <- spot.colors

        p <- plot_ly(data,
                     scene = scene,
                     x = ~xmax - x, y = ~y, z = ~z,
                     marker = list(color = data$spot.colors,
                                   showscale = FALSE,
                                   size = pt.size,
                                   opacity = pt.alpha)) %>%
          add_markers() %>%
          layout(title = paste(paste(channels.use, dims, sep = ": "), collapse = "; "), paper_bgcolor = ifelse(dark.theme, 'rgb(0, 0, 0)', 'rgb(255, 255, 255)'),
                 scene = list(zaxis = list(title = 'z', range = c(-add.margins, max(data$z) + add.margins), showticks = FALSE)))
        names(p$x$layoutAttrs[[1]]) <- c("title", "paper_bgcolor", scene)
        return(p)
      }
    }
  } else {
    if (length(dims) > 1) stop("Only one feature can be plotted at the time when blend = FALSE ... \n", call. = FALSE)

    # Set range for colors and create colorscale
    if (!is.null(cols)) {
      cs <- cscale(cols)
    } else {
      cs <- "Jet"
    }

    if (mode == "cloud") {
      # Run interpolation
      interpolated.data <- do.call(rbind, lapply(seq_along(data.list), function(i) {
        data <- setNames(data.list[[i]], c("warped_x", "warped_y", "z", "value"))
        section.input <- setNames(section.input.list[[i]], c("warped_x", "warped_y", "z", "value"))
        interpolated_data <- interpolate_2D_data(data, section.input, nx = nxy[1], ny = nxy[2], xy.range = apply(coords[, 1:2], 2, range))
      }))

      if (return.data) return(interpolated.data)

      p <- plot_ly(na.omit(interpolated.data),
                   scene = scene,
                   x = ~xmax - x, y = ~y, z = ~z,
                   marker = list(color = ~val,
                                 showscale = TRUE,
                                 colorscale = cs,
                                 size = pt.size,
                                 opacity = pt.alpha)) %>%
        add_markers() %>%
        layout(title = dims, paper_bgcolor = ifelse(dark.theme, 'rgb(0, 0, 0)', 'rgb(255, 255, 255)'),
               scene = list(zaxis = list(title = 'z', range = c(-add.margins, max(interpolated.data$z) + add.margins), showticks = FALSE)))
      names(p$x$layoutAttrs[[1]]) <- c("title", "paper_bgcolor", scene)
      return(p)
    } else {
      data <- setNames(data, c("x", "y", "z", "value"))

      if (return.data) return(interpolated.data)

      p <- plot_ly(data,
                   scene = scene,
                   x = ~xmax - x, y = ~y, z = ~z,
                   marker = list(color = ~value,
                                 showscale = TRUE,
                                 colorscale = cs,
                                 size = pt.size,
                                 opacity = pt.alpha)) %>%
        add_markers() %>%
        layout(title = dims, paper_bgcolor = ifelse(dark.theme, 'rgb(0, 0, 0)', 'rgb(255, 255, 255)'),
               scene = list(zaxis = list(title = 'z', range = c(-add.margins, max(data$z) + add.margins), showticks = FALSE)))
      names(p$x$layoutAttrs[[1]]) <- c("title", "paper_bgcolor", scene)
      return(p)
    }
  }
}


#' Plots the values of a feature in 3D
#'
#' @param object Seurat object
#' @param spots Subset spots to plot
#' @param features Features to plot
#' @param slot Which slot to pull the data from? [default: 'data']
#' @param zcoords Vector of z coordinates with the same length as the number of sections in the dataset [default: 1:#sections]
#' @param min.cutoff,max.cutoff Vector of minimum and maximum cutoff values for each feature, may specify quantile in the form of 'q##' where '##'
#' is the quantile (eg, 'q1', 'q10'). This can be useful if you have outlier values that skew the colorscale in the plot. For example, if you specify
#' 'q1', you will trim of values below the 1st percentile. [default: no cuttoffs]
#' @param blend Scale and blend expression values to visualize coexpression of two features (this options will override other coloring parameters).
#' See 'Blending values' below for a more thourough description.
#' @param pt.size Sets the size of points in the 3D plot
#' @param pt.alpha Sets the opacity of the points
#' @param cols Colors used to create a colorscale
#' @param add.margins Add margins along z axis to push sections closer to each other
#' @param channels.use Color channels to use for blending. Has to be a character vector of length 2 or 3 with "red", "green" and "blue"
#' color names specified [default: c("red", "green", "blue)]
#' @param scene Give the scene a name to allow for multiple subplots
#' @param return.data return the data.frame with x,y coordinates and interpolated values
#' instead of plotting
#' @param dark.theme Draws the plot with a dark theme
#'
#' @importFrom plotly plot_ly add_markers layout
#'
#' @export

FeaturePlot3D <- function (
  object,
  spots = NULL,
  features,
  mode = c("cloud", "spots"),
  zcoords = NULL,
  slot = 'data',
  min.cutoff = NA,
  max.cutoff = NA,
  blend = FALSE,
  pt.size = NULL,
  pt.alpha = 1,
  cols = NULL,
  add.margins = 0,
  channels.use = NULL,
  scene = "scene1",
  return.data = FALSE,
  dark.theme = FALSE,
  verbose = FALSE
) {

  # Check to see if Staffli object is present
  if (!"Staffli" %in% names(object@tools)) stop("Staffli object is missing from Seurat object ... \n", call. = FALSE)
  st.object <- GetStaffli(object)

  # Set mode
  mode <- match.arg(mode, c("cloud", "spots"))

  # Set point size if not provides
  pt.size <- pt.size %||% ifelse(mode == "cloud", 0.8, 5)

  # Check to see if Staffli object is present
  if (mode == "cloud") {
    if (!length(x = st.object@scatter.data) > 0) stop("3D stack is missing. Run Create3DStack() first ... \n", call. = FALSE)
    scatter.data <- st.object@scatter.data
  }

  # Get value
  if (!is.null(spots) & mode == "cloud") stop("Subsetting of spots only works if mode is set to 'spots' ... \n", call. = FALSE)
  spots <- spots %||% colnames(object)
  values <- FetchData(object = object, vars = features, slot = slot)[spots, , drop = T]

  coords <- do.call(rbind, lapply(seq_along(st.object@samplenames), function(i) {
    s <- st.object@samplenames[i]
    dims.raw <- as.numeric(st.object@dims[[i]][2:3])
    dims.scaled <- dim(st.object["raw"][[i]])
    sf.xy <- dims.raw[1]/dims.scaled[2]
    coords <- subset(st.object[[]], sample == s)[, c("warped_x", "warped_y")]/sf.xy
    coords$z <- i
    return(coords)
  }))
  coords <- coords[spots, ]

  # Use zcoords to convert z axis
  def_z <- unique(coords$z)
  if (!is.null(zcoords)) {
    if (length(x = def_z) != length(x = zcoords)) stop(paste0("zcoords (length = ", length(zcoords), ") has to be the same length as the number of samples ", length(def_z), " in the data"), call. = FALSE)
    coords$z <- zcoords[coords$z]
    if (mode == "cloud") scatter.data$z <- zcoords[scatter.data$z]
  }

  # Add values
  data <- setNames(data.frame(coords, values), c(colnames(coords), features))

  # Scale values
  data <- feature.scaler(data = data, features = features, min.cutoff = min.cutoff, max.cutoff = max.cutoff)

  # Split data
  if (mode == "cloud") {
    data.list <- split(data, data$z)
    section.input.list <- split(scatter.data, scatter.data$z)
    nxy <- colnames(scatter.data)[1:2] %>% as.numeric()
  }

  xmax <- lapply(st.object@dims, function(d) {d[2] %>% as.numeric()}) %>% unlist() %>% max()

  # Check if blend option is set
  if (blend) {
    if (!length(features) %in% c(2, 3)) {
      stop("Color blending only works for 2 or 3 features ... \n", call. = FALSE)
    } else {
      channels.use <- channels.use %||% c("red", "green", "blue")[1:length(features)]
      if (verbose) cat(paste0("Blending colors for features ",
                              paste0(ifelse(length(features) == 2, paste0(features[1], " and ", features[2]), paste0(features[1], features[2], " and ", features[2]))),
                              ": \n", paste(paste(features, channels.use, sep = ":"), collapse = "\n")))
      if (mode == "cloud") {
        interpolated.data <- do.call(rbind, lapply(seq_along(data.list), function(i) {
          data <- data.list[[i]]
          section.input <- section.input.list[[i]]
          interpolated_data.list <- list()
          for (i in 1:length(features)) {
            df <- interpolate_2D_data(setNames(data[, c(1:3, i + 3)], nm = c("warped_x", "warped_y", "z", "value")), section.input, nx = nxy[1], ny = nxy[2], xy.range = apply(coords[, 1:2], 2, range))
            colnames(df) <- c("warped_x", "warped_y", "z", features[i])
            interpolated_data.list[[i]] <- df
          }
          A <- interpolated_data.list[[1]]
          for (i in 2:length(features)) {
            A <- setNames(cbind(A, interpolated_data.list[[i]][, 4]), nm = c(colnames(A), features[i]))
          }
          colored.data <- apply(A[, features], 2, scales::rescale)
          colored.data[is.na(colored.data)] <- 0
          spot.colors <- ColorBlender(colored.data, channels.use)
          return(setNames(data.frame(A[, 1:3], spot.colors), nm = c("x", "y", "z", "spot.colors")))
        }))

        p <- plot_ly(interpolated.data,
                     scene = scene,
                     x = ~xmax - x, y = ~y, z = ~z,
                     marker = list(color = interpolated.data$spot.colors,
                                   showscale = FALSE,
                                   size = pt.size,
                                   opacity = pt.alpha)) %>%
          add_markers() %>%
          layout(title = paste(paste(channels.use, features, sep = ": "), collapse = "; "), paper_bgcolor = ifelse(dark.theme, 'rgb(0, 0, 0)', 'rgb(255, 255, 255)'),
                 scene = list(zaxis = list(title = 'z', range = c(-add.margins, max(interpolated.data$z) + add.margins), showticks = FALSE)))
        names(p$x$layoutAttrs[[1]]) <- c("title", "paper_bgcolor", scene)
        return(p)
      } else {
        data <- setNames(data, nm = c("x", "y", "z", features))
        colored.data <- apply(data[, features], 2, scales::rescale)
        colored.data[is.na(colored.data)] <- 0
        spot.colors <- ColorBlender(colored.data, channels.use)
        data$spot.colors <- spot.colors

        p <- plot_ly(data,
               scene = scene,
               x = ~xmax - x, y = ~y, z = ~z,
               marker = list(color = data$spot.colors,
                             showscale = FALSE,
                             size = pt.size,
                             opacity = pt.alpha)) %>%
          add_markers() %>%
          layout(title = paste(paste(channels.use, features, sep = ": "), collapse = "; "), paper_bgcolor = ifelse(dark.theme, 'rgb(0, 0, 0)', 'rgb(255, 255, 255)'),
                 scene = list(zaxis = list(title = 'z', range = c(-add.margins, max(data$z) + add.margins), showticks = FALSE)))
        names(p$x$layoutAttrs[[1]]) <- c("title", "paper_bgcolor", scene)
        return(p)
      }
    }
  } else {
    if (length(features) > 1) stop("Only one feature can be plotted at the time when blend = FALSE ... \n", call. = FALSE)

    # Set range for colors and create colorscale
    if (!is.null(cols)) {
      cs <- cscale(cols)
    } else {
      cs <- "Jet"
    }

    if (mode == "cloud") {
      # Run interpolation
      interpolated.data <- do.call(rbind, lapply(seq_along(data.list), function(i) {
        data <- data.list[[i]]
        data <- setNames(data, c("warped_x", "warped_y", "z", "value"))
        section.input <- section.input.list[[i]]
        interpolated_data <- interpolate_2D_data(data, section.input, nx = nxy[1], ny = nxy[2], xy.range = apply(coords[, 1:2], 2, range))
      }))
      interpolated.data <- setNames(interpolated.data, c("x", "y", "z", "val"))

      if (return.data) return(interpolated.data)

      p <- plot_ly(na.omit(interpolated.data),
                   scene = scene,
                   x = ~xmax - x, y = ~y, z = ~z,
                   marker = list(color = ~val,
                                 showscale = TRUE,
                                 colorscale = cs,
                                 size = pt.size,
                                 opacity = pt.alpha)) %>%
        add_markers() %>%
        layout(title = features, paper_bgcolor = ifelse(dark.theme, 'rgb(0, 0, 0)', 'rgb(255, 255, 255)'),
               scene = list(zaxis = list(title = 'z', range = c(-add.margins, max(interpolated.data$z) + add.margins), showticks = FALSE)))
      names(p$x$layoutAttrs[[1]]) <- c("title", "paper_bgcolor", scene)
      return(p)
    } else {
      colnames(data) <- c("x", "y", "z", "value")

      if (return.data) return(data)

      p <- plot_ly(data,
                   scene = scene,
                   x = ~xmax - x, y = ~y, z = ~z,
                   marker = list(color = ~value,
                                 showscale = TRUE,
                                 colorscale = cs,
                                 size = pt.size,
                                 opacity = pt.alpha)) %>%
        add_markers() %>%
        layout(title = features, paper_bgcolor = ifelse(dark.theme, 'rgb(0, 0, 0)', 'rgb(255, 255, 255)'),
               scene = list(zaxis = list(title = 'z', range = c(-add.margins, max(data$z) + add.margins), showticks = FALSE)))
      names(p$x$layoutAttrs[[1]]) <- c("title", "paper_bgcolor", scene)
     return(p)
    }
  }
}


#' Plots the values of a feature in 3D
#'
#' @param object Seurat object
#' @param features
#' \itemize{
#'     \item An \code{Assay} feature (e.g. a gene name - "MS4A1")
#'     \item A column name from meta.data (e.g. mitochondrial percentage - "percent.mito")
#'     \item A column name from dimensionality reduction output (e.g. principal component 1 - "PC_1")
#' }
#' @param slot Which slot to pull the data from? [default: 'data']
#' @param zcoords Vector of z coordinates with the same length as the number of sections in the dataset [default: 1:#sections]
#' @param spots Vector of spots to plot (default is all spots)
#' @param min.cutoff,max.cutoff Vector of minimum and maximum cutoff values for each feature, may specify quantile in the form of 'q##' where '##'
#' is the quantile (eg, 'q1', 'q10'). This can be useful if you have outlier values that skew the colorscale in the plot. For example, if you specify
#' 'q1', you will trim of values below the 1st percentile. [default: no cuttoffs]
#' @param blend Scale and blend expression values to visualize coexpression of two features (this options will override other coloring parameters).
#' See 'Blending values' below for a more thourough description.
#' @param pt.size Sets the size of points in the 3D plot
#' @param pt.alpha Sets the opacity of the points
#' @param cols Character vector of colors with equal length to the number of features. These colors will override
#' the selection of HSV colors and can therefore not be encoded in the same way. Instead of tuning the saturation/value
#' parameters we can add an alpha channel, making spots with values close to zero completely transparent. This has to
#' be activated by setting `add.alpha = TRUE`
#' @param add.alpha Adds opacity to colors.
#' @param add.margins Add margins along z axis to push sections closer to each other
#' @param channels.use Color channels to use for blending. Has to be a character vector of length 2 or 3 with "red", "green" and "blue"
#' color names specified [default: c("red", "green", "blue)]
#' @param scene Give the scene a name to allow for multiple subplots
#' @param return.data return the data.frame with x,y coordinates and interpolated values
#' instead of plotting
#' @param dark.theme Draws the plot with a dark theme
#' @param verbose Print messages
#'
#' @importFrom plotly plot_ly add_markers layout
#'
#' @export

HSVPlot3D <- function (
  object,
  features,
  mode = c("cloud", "spots"),
  zcoords = NULL,
  spots = NULL,
  min.cutoff = NA,
  max.cutoff = NA,
  slot = "data",
  blend = FALSE,
  pt.size = NULL,
  pt.alpha = 1,
  cols = NULL,
  add.alpha = FALSE,
  add.margins = 0,
  channels.use = NULL,
  scene = "scene1",
  return.data = FALSE,
  dark.theme = FALSE,
  verbose = FALSE
) {

  # Check to see if Staffli object is present
  if (!"Staffli" %in% names(object@tools)) stop("Staffli object is missing from Seurat object ... \n", call. = FALSE)
  st.object <- GetStaffli(object)

  # Set mode
  mode <- match.arg(mode, c("cloud", "spots"))

  # Set point size if not provides
  pt.size <- pt.size %||% ifelse(mode == "cloud", 0.8, 5)

  # Collect data
  spots <- spots %||% colnames(x = object)
  data <- FetchData(object = object, vars = features, cells = spots, slot = slot)
  data.type <- unique(sapply(data, class))

  # Scale values
  data <- feature.scaler(data = data, features = features, min.cutoff = min.cutoff, max.cutoff = max.cutoff)

  # Stop if feature classes are not numeric/integer
  if (!all(data.type %in% c("numeric", "integer"))) {
    stop("Only features of class 'integer' or 'numeric' are allowed ... ")
  }

  # Check to see if Staffli object is present
  if (mode == "cloud") {
    if (verbose) cat("Plotting features in 'cloud' mode ... \n")
    if (!length(x = st.object@scatter.data) > 0) stop("3D stack is missing. Run Create3DStack() first ... \n", call. = FALSE)
    scatter.data <- st.object@scatter.data
  }

  coords <- do.call(rbind, lapply(seq_along(st.object@samplenames), function(i) {
    s <- st.object@samplenames[i]
    dims.raw <- as.numeric(st.object@dims[[i]][2:3])
    dims.scaled <- dim(st.object["raw"][[i]])
    sf.xy <- dims.raw[1]/dims.scaled[2]
    coords <- subset(st.object[[]], sample == s)[, c("warped_x", "warped_y")]/sf.xy
    coords$z <- i
    return(coords)
  }))

  # Use zcoords to convert z axis
  def_z <- unique(coords$z)
  if (!is.null(zcoords)) {
    if (verbose) cat("Converting z values using provided 'zcoords' ... \n")
    if (length(x = def_z) != length(x = zcoords)) stop(paste0("zcoords (length = ", length(zcoords), ") has to be the same length as the number of samples ", length(def_z), " in the data"), call. = FALSE)
    coords$z <- zcoords[coords$z]
    if (mode == "cloud") scatter.data$z <- zcoords[scatter.data$z]
  }

  # Add values
  data <- cbind(coords, data)

  # Split data
  if (mode == "cloud") {
    data.list <- split(data, data$z)
    section.input.list <- split(scatter.data, scatter.data$z)
    nxy <- colnames(scatter.data)[1:2] %>% as.numeric()
  }

  xmax <- lapply(st.object@dims, function(d) {d[2] %>% as.numeric()}) %>% unlist() %>% max()

  # Generate HSV encoded colors
  if (is.null(cols)) {
    if (verbose) cat(paste0("Defining hue colors for ", length(x = features), " features ... \n"))
    hue_breaks <- seq(0, 1, length.out = length(x = features) + 1)[1:length(x = features)]
    hsv.matrix <- t(matrix(c(hue_breaks, rep(1, length(hue_breaks )), rep(1, length(hue_breaks))), ncol = 3))
    rownames(hsv.matrix) <- c("h", "s", "v")
    ann.cols <- apply(hsv.matrix, 2, function(x) hsv(x[1], x[2], x[3]))
    # Define HSV conversion function
    hsv.func <- ifelse(dark.theme,
                       function(x) hsv(h = x[1, ][which.max(x[3, ])], s = 1, v = max(x[3, ])),
                       function(x) hsv(h = x[1, ][which.max(x[3, ])], v = 1, s = max(x[3, ])))
  } else {
    if (length(x = features) != length(x = cols)) stop("Length of features and cols must match ...", call. = FALSE)
    warning("Using user defined colors with opacity. HSV scale will not be used ...", call. = FALSE)
    ann.cols <- cols
  }

  if (mode == "cloud") {
    if (verbose) cat("Interpolating data across 2D point patterns ... \n")
    if (!is.null(cols) & !add.alpha) warning("add.alpha should be set to TRUE when using custom colors ... \n", call. = FALSE)
    interpolated.data <- do.call(rbind, lapply(seq_along(data.list), function(i) {
      data <- data.list[[i]]
      section.input <- section.input.list[[i]]
      interpolated_data.list <- list()
      for (i in 1:length(features)) {
        df <- interpolate_2D_data(setNames(data[, c(1:3, i + 3)], nm = c("warped_x", "warped_y", "z", "value")), section.input, nx = nxy[1], ny = nxy[2], xy.range = apply(coords[, 1:2], 2, range))
        colnames(df) <- c("warped_x", "warped_y", "z", features[i])
        interpolated_data.list[[i]] <- df
      }
      A <- interpolated_data.list[[1]]
      for (i in 2:length(features)) {
        A <- setNames(cbind(A, interpolated_data.list[[i]][, 4]), nm = c(colnames(A), features[i]))
      }
      A[is.na(A)] <- 0

      # Create hsv matrix and select highest v
      if (is.null(cols)) {
        if (verbose) cat("Converting values to HSV ... \n")
        d <- array(dim = c(nrow(A), 3, length(x = features)))
        for (i in 1:length(features)) {
          ftr <- features[i]
          s <- data.frame(h = hue_breaks[i],
                          s = 1,
                          v = scales::rescale(A[, ftr, drop = T] %>% as.numeric())) %>% as.matrix()
          d[, , i] <- s
        }
        red.cols <- unlist(apply(d, 1, function (x) {
          hsvc <- hsv.func(x)
          if (add.alpha)  hsvc <- scales::alpha(colour = hsvc, alpha = max(x[3, ]))
          return(hsvc)
        }))
      } else {
        if (verbose) cat("Using provided colors ... \n")
        d <- array(dim = c(nrow(A), 1, length(x = features)))
        for (i in 1:length(features)) {
          ftr <- features[i]
          s <- data.frame(v = scales::rescale(A[, ftr, drop = T] %>% as.numeric())) %>% as.matrix()
          d[, , i] <- s
        }
        red.cols <- unlist(apply(d, 1, function (x) {
          alpha_col <- cols[which.max(x[1, ])]
          if (add.alpha) {
            alpha_col <- scales::alpha(colour = alpha_col, alpha = max(x[1, ]))
          }
          return(alpha_col)
        }))
      }

      return(setNames(data.frame(A[, 1:3], red.cols, stringsAsFactors = F), nm = c("x", "y", "z", "spot.colors")))
    }))

    p <- plot_ly(interpolated.data,
                 scene = scene,
                 x = ~xmax - x,
                 y = ~y,
                 z = ~z,
                 marker = list(color = interpolated.data$spot.colors,
                               showscale = FALSE,
                               size = pt.size,
                               opacity = pt.alpha)) %>%
      add_markers() %>%
      layout(title = paste(features, collapse = "; "), paper_bgcolor = ifelse(dark.theme, "rgb(0, 0, 0)", "#FFF"),
             scene = list(zaxis = list(title = 'z', range = c(-add.margins, max(interpolated.data$z) + add.margins), showticks = FALSE)))
    names(p$x$layoutAttrs[[1]]) <- c("title", "paper_bgcolor", scene)
    return(p)
  } else {
    data <- setNames(data, nm = c("x", "y", "z", features))

    if (is.null(cols)) {
      d <- array(dim = c(nrow(data), 3, length(x = features)))
      if (verbose) cat("Converting values to HSV ... \n")
      for (i in 1:length(features)) {
        ftr <- features[i]
        s <- data.frame(h = hue_breaks[i],
                        s = 1,
                        v = scales::rescale(data[, ftr, drop = T] %>% as.numeric())) %>% as.matrix()
        d[, , i] <- s
      }

      red.cols <- unlist(apply(d, 1, function (x) {
        hsvc <- hsv.func(x)
        if (add.alpha)  hsvc <- scales::alpha(colour = hsvc, alpha = max(x[3, ]))
        return(hsvc)
      }))
    } else {
      if (verbose) cat("Using provided colors ... \n")
      d <- array(dim = c(nrow(data), 1, length(x = features)))
      for (i in 1:length(features)) {
        ftr <- features[i]
        s <- data.frame(v = scales::rescale(data[, ftr, drop = T] %>% as.numeric())) %>% as.matrix()
        d[, , i] <- s
      }
      red.cols <- unlist(apply(d, 1, function (x) {
        alpha_col <- cols[which.max(x[1, ])]
        if (add.alpha) {
          alpha_col <- scales::alpha(colour = alpha_col, alpha = max(x[1, ]))
        }
        return(alpha_col)
      }))
    }

    data$spot.colors <- red.cols
    p <- plot_ly(data,
                 scene = scene,
                 x = ~xmax - x, y = ~y, z = ~z,
                 marker = list(color = data$spot.colors,
                               showscale = FALSE,
                               size = pt.size,
                               opacity = pt.alpha)) %>%
      layout(title = paste(features, collapse = "; "), paper_bgcolor = ifelse(dark.theme, "rgb(0, 0, 0)", "#FFF"),
             scene = list(zaxis = list(title = 'z', range = c(-add.margins, max(data$z) + add.margins), showticks = FALSE))) %>%
      add_markers()

    names(p$x$layoutAttrs[[1]]) <- c("title", "paper_bgcolor", scene)
    return(p)
  }
}



#' Create a plotly compatible colorscale
#'
#' @param cols Vector of colors
cscale <- function (cols) {
  rgbs <- do.call(cbind, lapply(cols, col2rgb))
  breaks <- seq(0, 1, length.out = length(cols))

  colorscale <- lapply(seq_along(breaks), function(i) {
    c(breaks[i], paste0('rgb(', paste(rgbs[, i], collapse = ', '), ')'))
  })
  return(colorscale)
}


#' Assign points to a grid
#'
#' takes a point pattern and a raster object as input and
#' assigns each point to a cell of the raster object.
#'
#' @param scatter data.frame with x,y coordinates for a point pattern
#' @param r A 'raster' object
#' @param nx Number of cells to interpolate over across the x axis
#'
#' @importFrom raster coordinates
#'
#' @return data.frame object with point pattern and associated grid.cell number

rasterize_scatter <- function (
  scatter,
  r,
  nx
) {
  pixel.centers = coordinates(r)
  set1 = scatter[, 1:2]
  set2 = pixel.centers[, 1:2]
  new.set2 <- apply(set1, 1, function(x) which.min(colSums((t(set2) - x)^2)))
  section.input = as.data.frame(cbind(scatter[, 1], scatter[, 2], scatter$z, new.set2))
  y = r@nrows
  colnames(section.input) <- c(nx, y, "z", "grid.cell")
  return(section.input)
}


#' Interpolate value across point pattern
#'
#' Takes a data.frame with spot coordinates and a 'value' column as input
#' and interpolates the 'value' across a grid defined by the dimensions nx, ny
#' and finally assigns the interpolated values to points associated with each grid cell.
#'
#' @param data data.frame object with spot coordinates and a 'value' column with numeric data
#' @param section.input data.frame object with x,y coordinates for a point pattern to interpolate values over
#' @param nx,ny Dimensions of grid
#' @param xy.range Range of x,y coordinates to interpolate over
#'
#' @return data.frame object with point pattern and interpolated values

interpolate_2D_data <- function (
  data,
  section.input,
  nx,
  ny,
  xy.range
){
  x1 = as.numeric(data[, "warped_x", drop = T])
  y1 = as.numeric(data[, "warped_y", drop = T])
  w1 = as.numeric(data[, "value"])

  # Run the interpolation
  s1 =  akima::interp(x = x1, y = y1, z = w1, nx = nx, ny = ny, xo = seq(xy.range[1, 1], xy.range[2, 1], length = nx), yo = seq(xy.range[1, 2], xy.range[2, 2], length = ny))
  mat.1 = s1$z
  mat.1 <- mat.1[, ncol(mat.1):1]
  col = as.numeric(mat.1)

  set = col[section.input[, 4]]
  section.xyz.value = cbind(section.input[, 1:3], set)
  df <- as.data.frame(section.xyz.value)
  colnames(df) <- c("x", "y", "z", "val")
  return(df)
}

