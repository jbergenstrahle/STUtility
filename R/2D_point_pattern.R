
#' Function used to generate a 2D set of points from an HE image
#'
#' The defualt segmentation uses the HE image as input and defines any pixel with an intensity value
#' below a threshold to be a point. The number of points can be downsampeld to limit the maximum number.
#'
#' @param object Seurat object
#' @param type Sets the image type to run segmentation on
#' @param sammple.index Integer value specifying the index of the sample to be analyzed
#' @param limit Sets the intensity threshold in the interval [0, 1]
#' @param maxnum Integer value specifying the maximum number of points
#' @param edges Extracts the coordinates of the edges instead

scatter_HE <- function (
  object,
  type = "masked",
  sample.index = NULL,
  limit = 0.5,
  maxnum = 5e4,
  edges = FALSE
){
  sample.index <- sample.index %||% 1

  if (!type %in% names(object@tools)) stop(paste0(type, " images not fount in Seurat obejct"), call. = FALSE)

  if (edges) {
    bw.image <- get.edges(object, sample.index, type = type)
    xyset = which(bw.image > limit, arr.ind = TRUE)
  } else {
    bw.image = grayscale(as.cimg(object@tools[[type]][[sample.index]]))
    xyset = which(bw.image < limit*255, arr.ind = TRUE)
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
#' @param obejct Seurat object
#'
#' @inheritParams scatter_HE

grid.from.seu <- function (
  object,
  type,
  limit = 0.3,
  maxnum = 5e4,
  edges = FALSE
) {

  if (!type %in% names(object@tools)) stop(paste0(type, " images not fount in Seurat obejct"), call. = FALSE)

  setNames(lapply(1:length(object@tools[[type]]), function(sample.index) {
    scatter <- scatter_HE(object = se, sample.index = sample.index, maxnum = 1e3, limit = limit, type = type, edges = edges)
    if (type == "processed") {
      xy.names <- c("warped_x", "warped_y")
    } else {
      xy.names <- c("pixel_x", "pixel_y")
    }
    coords <- subset(se[[]], sample == sample.index)[, xy.names]
    dims.raw <- as.numeric(se@tools$dims[[sample.index]][2:3])
    dims.scaled <- dim(se@tools$raw[[sample.index]])
    sf.xy <- dims.raw/rev(dims.scaled)
    coords <- coords/sf.xy
    return(list(scatter = scatter, coords = coords))
  }), nm = names(object@tools$masked))
}