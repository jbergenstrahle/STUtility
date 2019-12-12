#' Create 3D stack
#'
#' Creates a stack of point patterns with aligned data and stores it in a
#' `Staffli` object.
#'
#' @param object Seurat object
#' @param limit Cut-off threshold for segmentation of points. Has to be a value between 0-1.
#' @param verbose Print messages
#'
#' @inheritParams rasterize_scatter
#'
#' @return Seurat object
#'
#' @importFrom akima interp
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
    object <- SwitchResolution(object, xdim = 2e3, verbose = verbose)
  }
  st.object <- GetStaffli(object)

  # Obtain number of cells to interpolate over
  coords <- do.call(rbind, lapply(seq_along(st.object@samplenames), function(i) {
    s <- st.object@samplenames[i]
    coords <- subset(st.object[[]], sample == s)[, c("warped_x", "warped_y")]
    dims.raw <- iminfo(st.object)[[s]][2:3] %>% as.numeric()
    dims.scaled <- scaled.imdims(st.object)[[s]]
    sf.xy <- dims.raw[2]/dims.scaled[1]
    coords <- coords/sf.xy
    coords$z <- i
    return(coords)
  }))

  if (verbose) cat("Running approximative segmentation of nuclei ... \n")
  scatters <- do.call(rbind, lapply(seq_along(st.object@samplenames), function(i) {
    s <- st.object@samplenames[i]
    scatter <- scatter_HE(st.object, type = "processed", sample.index = s, maxnum = 5e4, limit = limit, edges = FALSE)
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
  if (verbose) cat("Assigning nuclei position to grid cells ... \n")
  section.input <- rasterize_scatter(scatters, r, nx)

  st.object@scatter.data <- section.input
  object@tools$Staffli <- st.object
  return(object)
}

#' Plots the values of a feature in 3D
#'
#' @param object Seurat object
#' @param feature Feature to plot
#' @param slot Which slot to pull the data from? [default: 'data']
#'
#' @importFrom plotly plot_ly add_markers layout
#'
#' @export

FeaturePlot3D <- function (
  object,
  features,
  slot = 'data',
  cols = NULL
) {

  # Check to see if Staffli object is present
  if (!"Staffli" %in% names(object@tools)) stop("Staffli object is missing from Seurat object ... \n", call. = FALSE)
  st.object <- GetStaffli(object)

  # Check to see if Staffli object is present
  if (!length(x = st.object@scatter.data) > 0) stop("3D stack is missing. Run Create3DStack() first ... \n", call. = FALSE)
  scatter.data <- st.object@scatter.data

  # Get value
  values <- FetchData(object = object, vars = features, slot = slot)[, , drop = T]

  coords <- do.call(rbind, lapply(seq_along(st.object@samplenames), function(i) {
    s <- st.object@samplenames[i]
    dims.raw <- as.numeric(st.object@dims[[s]][2:3])
    dims.scaled <- dim(st.object["raw"][[i]])
    sf.xy <- dims.raw[1]/dims.scaled[2]
    coords <- subset(st.object[[]], sample == s)[, c("warped_x", "warped_y")]/sf.xy
    coords$z <- i
    return(coords)
  }))

  # Add values
  data <- cbind(coords, value = values)

  # Split data
  data.list <- split(data, data$z)
  section.input.list <- split(scatter.data, scatter.data$z)
  nxy <- colnames(scatter.data)[1:2] %>% as.numeric()

  interpolated.data <- do.call(rbind, lapply(seq_along(data.list), function(i) {
    data <- data.list[[i]]
    section.input <- section.input.list[[i]]
    interpolated_data <- interpolate_2D_data(data, section.input, nx = nxy[1], ny = nxy[2])
  }))

  range.vals <- range(interpolated.data$val, na.rm = TRUE)
  if (!is.null(cols)) {
    cs <- cscale(cols)
  } else {
    cs <- "Jet"
  }

  plot_ly(na.omit(interpolated.data),
          x = ~2e3-x, y = ~y, z = ~z,
          marker = list(color = ~val,
                        showscale = T,
                        colorscale = cs,
                        size = 1.2,
                        opacity = 0.6)) %>%
    add_markers() %>%
    layout(paper_bgcolor = 'rgb(0, 0, 0)', scene = list(zaxis = list(title = 'z', range = c(-2, max(interpolated.data$z) + 2), showticks = FALSE)))
}

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
#' @return data.frame object with point pattern and associated grid.cell number

rasterize_scatter <- function (
  scatter,
  r,
  nx
) {
  tab = table(cellFromXY(r, scatter[, 1:2]))
  r[as.numeric(names(tab))] = tab
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
#'
#' @return data.frame object with point pattern and interpolated values

interpolate_2D_data <- function (
  data,
  section.input,
  nx,
  ny
){
  x1 = as.numeric(data[, "warped_x", drop = T])
  y1 = as.numeric(data[, "warped_y", drop = T])
  w1 = as.numeric(data[, "value"])

  # Run the interpolation
  s1 =  akima::interp(x = x1, y = y1, z = w1, nx = nx, ny = ny)
  mat.1 = s1$z
  mat.1 <- mat.1[, ncol(mat.1):1]
  col = as.numeric(mat.1)

  set = col[section.input[, 4]]
  section.xyz.value = cbind(section.input[, 1:3], set)
  df <- as.data.frame(section.xyz.value)
  colnames(df) <- c("x", "y", "z", "val")
  return(df)
}

