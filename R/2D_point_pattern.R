
#' Function used to generate a 2D set of points from an HE image
#'
#' The defualt segmentation uses the HE image as input and defines any pixel with an intensity value
#' below a threshold to be a point. The number of points can be downsampeld to limit the maximum number.
#'
#' @param object Staffli object
#' @param type Sets the image type to run segmentation on
#' @param sample.index Integer value specifying the index of the sample to be analyzed
#' @param limit Sets the intensity threshold in the interval [0, 1]
#' @param maxnum Integer value specifying the maximum number of points
#' @param edges Extracts the coordinates of the edges instead
#' @param custom.edge.detector Custom function used to detect edges

scatter_HE <- function (
  object,
  type = "masked.masks",
  sample.index = NULL,
  limit = 0.5,
  maxnum = 5e4,
  edges = FALSE,
  custom.edge.detector = NULL
){
  sample.index <- sample.index %||% 1

  if (!type %in% rasterlists(object)) stop(paste0(type, " images not found in Staffli object"), call. = FALSE)

  if (edges) {
    bw.image <- get.edges(object, sample.index, type = type)
    xyset = which(bw.image > limit, arr.ind = TRUE)
  } else {
    if (is.null(custom.edge.detector)) {
      im <- as.cimg(object[type][[sample.index]])
      bw.image <- im[, , 1, 1] %>% as.cimg()
      #f <- ecdf(bw.image)
      #bw.image.eq <- f(bw.image) %>% as.cimg(dim = dim(bw.image))
      #bw.image = grayscale(im)
      if (type %in% c("processed.masks", "masked.masks")) {
        xyset = which(bw.image > limit*255, arr.ind = TRUE)
      } else {
        seg <- imagerExtra::ThresholdAdaptive(bw.image, k = 0.1)
        xyset = which(!seg, arr.ind = TRUE)
      }
    } else {
      im <- as.cimg(object[type][[sample.index]])
      if (class(custom.edge.detector) != "function") stop("Provided custom edge detector is not a function")
      seg <- custom.edge.detector(im)
      if (!"pixset" %in% class(seg)) stop("Provided custom edge detector did not return a pixset object")
      xyset = which(seg, arr.ind = TRUE)
    }
  }
  set.seed(1)
  if (maxnum < nrow(xyset)) {
    xyset <- xyset[sample(1:nrow(xyset), size = maxnum, replace = FALSE), ]
  }
  xyset <- xyset[, 1:2] %>% as.data.frame() %>% setNames(nm = c("x", "y"))
  return(xyset)
}


#' Creates 2D point patterns for a set of images
#'
#' @param object Staffli object
#'
#' @inheritParams scatter_HE

grid.from.staffli <- function (
  object,
  type,
  limit = 0.3,
  maxnum = 1e3,
  edges = FALSE,
  custom.edge.detector = NULL
) {
  l <- setNames(lapply(1:length(names(object)), function(sample.index) {
    scatter <- scatter_HE(object = object, sample.index = sample.index, maxnum = maxnum, limit = limit, type = type, edges = edges, custom.edge.detector = custom.edge.detector)
    if (type %in% c("processed", "processed.masks")) {
      xy.names <- c("warped_x", "warped_y")
    } else {
      xy.names <- c("pixel_x", "pixel_y")
    }
    coords <- subset(object[[]], sample == sample.index)[, xy.names]
    dims.raw <- object@dims[[sample.index]][2:3] %>% as.numeric()
    dims.scaled <- scaled.imdims(object)[[sample.index]]
    sf.xy <- dims.raw[2]/dims.scaled[1]
    coords <- coords/sf.xy
    return(list(scatter = scatter, coords = coords))
  }), nm = names(object))
}

