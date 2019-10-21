#' The Staffli Class
#'
#' The Staffli object is a representation of images for Spatial Transcriptomics experiments
#' that contains lists of scaled images in raster format and coordinates for the corresponding
#' spots that can be used to map expression data onto the images.
#'
#' @slot imgs A character vector of paths to the raw HE images.
#' @slot rasterlists A list of lists containing images in 'raster' format
#' @slot transformations A list of 3x3 transformation matrices used to transform (x,y)-coordinates
#' of an aligned image to a reference image.
#' @slot meta.data Contains meta-information about each spot, starting with regular (x,y)-coordinates
#' found in the headers of gene-count matrices obtained with the ST method. The meta.data
#' can also contain (adj_x, adj_y)-coordinates which are defined in the same grid but adjusted using the
#' ST spot detector and (pixel_x, pixel_y)-coordinates which are defined in the coordinate system of the
#' original image. Finally, a "sample" column is used to group the spots by image (capture area).
#' @slot xdim The width of the scaled images in pixels.
#' @slot limits Specifies the limits of the array (e.g. limits = c(100, 100) means that the array
#' is a 100 spots wide and a 100 spots high)
#' @slot dims List of numerical vectors specifying the dimensions of the original images.
#' @slot platform Specify the platform used to generate the ST data [options: 'Visium', '2k', '1k']
#' @slot samplenames Character specifying the samplenames.
#' @slot version Package version.
#'
#' @name Staffli-class
#' @rdname Staffli-class
#' @exportClass Seurat
#'
Staffli <- setClass (
  Class = 'Staffli',
  slots = c(
    imgs = 'character',
    rasterlists = 'list',
    scatter.data = 'data.frame',
    transformations = 'list',
    meta.data = 'data.frame',
    xdim = 'numeric',
    limits = 'numeric',
    dims = 'list',
    platform = 'character',
    samplenames = 'character',
    version = 'package_version'
  )
)


#' Create a Staffli object
#'
#' Create a Staffli object from set of images and associated spot coordinates.
#'
#' @param imgs Character vector specifying paths to images in jpg or png format.
#' @param meta.data Spot-level metadata to add to the Staffli object. Should be a data.frame with
#' required fields 'x' and 'y' which specifies the ST array coorindates and a column calles 'sample'
#' which maps spots to a certain image.
#' @param xdim Specifies the width of the scaled images in pixels [default: 400 pixels]
#' @param platform Specify the platform used to generate the ST data [options: 'Visium', '2k', '1k']
#'
#' @importFrom utils packageVersion
#' @importFrom Matrix colSums
#' @export
#'
CreateStaffliObject <- function (
  imgs = NULL,
  meta.data,
  xdim = 400,
  platform = 'Visium'
) {
  if (!all(c("x", "y", "sample") %in% colnames(meta.data))) stop(paste0("Invalid meta.data; one of 'x', 'y' or 'sample' is missing"), call. = FALSE)

  if (platform == 'Visium') {
    limits <- c(71, 71)
  } else if (platform == '1k') {
    limits <- c(33, 35)
  } else if (platform == '1k') {
    limits <- c(67, 64)
  } else {
    stop(paste0(platform, " is not a valid option ... \n"), call. = FALSE)
  }

  object <- new (
    Class = 'Staffli',
    imgs = imgs,
    meta.data = meta.data,
    xdim = xdim,
    limits = limits,
    platform = platform,
    samplenames = unique(meta.data[, "sample"]),
    version = packageVersion(pkg = 'STUtility')
  )

  return(object)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Staffli methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setMethod (
  f = "names",
  signature = "Staffli",
  definition = function(x) {
    x@samplenames
  }
)


setMethod (
  f = "show",
  signature = "Staffli",
  definition = function(object) {
    cat("An object of class", class(x = object), "\n")
    cat(
      nrow(object@meta.data),
      'features across',
      length(unique(object@meta.data[, "sample"])),
      'samples. \n'
    )
    if (length(object@rasterlists) > 0) {
      cat(
        'Available image representations: \n\t',
        paste(names(object@rasterlists), collapse = ", "),
        '\n'
      )
    }
  }
)


#' @importFrom imager as.cimg
#'
setMethod (
  f = "plot",
  signature = "Staffli",
  definition = function(x, type = NULL, ...) {
    object <- x
    ncols <- ceiling(length(x = object@imgs)); nrows <- ceiling(length(x = object@imgs)/ncols)
    if (length(x@rasterlists) == 0) {
      par(mfrow = c(nrows, ncols))
      for (s in object@samplenames) {
        d <- subset(object@meta.data, sample == s)
        plot(d[, c("x", "y")], ...)
      }
    } else {
      type <- type %||% "raw"
      if (type == "raw") {
        xy <- c("pixel_x", "pixel_y")
      } else if (type %in% c("masked", "masked.masks")) {
        if (!"masked" %in% names(object@rasterlists)) stop("Masked images not available in Staffli object", call. = FALSE)
        xy <- c("pixel_x", "pixel_y")
      } else if (type %in% c("processed", "processed.masks")) {
        if (!"processed" %in% names(object@rasterlists)) stop("Masked images not available in Staffli object", call. = FALSE)
        xy <- c("warped_x", "warped_y")
      }
      par(mfrow = c(nrows, ncols))
      for (s in object@samplenames) {
        d <- subset(object@meta.data, sample == s)
        im <- object@rasterlists[[type]][[s]] %>% as.cimg()
        plot(im, axes = FALSE)
        points(d[, xy]/(object@dims[[s]]$width/object@xdim), ...)
      }
    }
  }
)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Staffli functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @param object Object of class Staffli
#' @param image.paths Paths to HE images. This is only required if image paths are missing in the Seurat object.
#' @param xdim Sets the pixel width for scaling, e.g. 400 (maximum allowed width is 1000 pixels)
#' @param verbose Print messages
#'
#' @importFrom magick image_read
#'
#' @export

LoadImages.Staffli <- function (
  object,
  image.paths = NULL,
  xdim = 400,
  verbose = FALSE,
  time.resolve = TRUE
) {

  # Check that image paths are present
  image.paths <- image.paths %||% object@imgs
  if (length(x = image.paths) == 0) stop("No images provided. Provide images using image.paths \n", call. = FALSE)
  if (length(x = image.paths) != length(unique(object@meta.data[, "sample"]))) stop(paste0("Number of images (", length(x = image.paths), ") must match the number of samples (", length(unique(object@meta.data[, "sample"])), ")\n"), call. = FALSE)

  # Check that the image with is no more than 2000 pixels
  if (xdim > 2000) stop("xdim cannot be larger than 2000")

  if (verbose) cat(paste0("Loading images for ", length(x = object@samplenames), " samples: \n"))

  # Read images using the 'magick' library
  imgs <- c()
  dims <- list()
  for (i in seq_along(object@imgs)) {
    path <- object@imgs[i]
    if (verbose) cat("  Reading ", path , " for sample ", object@meta.data[, "sample"][i], " ... \n", sep = "")
    im <- image_read(path)
    dims <- c(dims, list(image_info(im)))
    if (verbose) {
      info <- dims[[i]]
      width <- as.numeric(info[2]); height <- as.numeric(info[3])
      ydim <- round(height/(width/xdim))
      cat("  Scaling down sample ", unique(object@meta.data[, "sample"])[i], " image from ", paste(width, height, sep = "x"), " pixels to ", paste(xdim, ydim, sep = "x"), " pixels \n", sep = "")
    }
    im <- image_scale(im, paste0(xdim))
    #tmpf <- tempfile()
    #image_write(im, path = tmpf)
    #im <- image_read(tmpf)
    imgs <- c(imgs, list(as.raster(im)))
    if (time.resolve) {
      gc()
      sleepy(5)
    }
  }

  object@dims <- setNames(dims, nm = object@samplenames)
  object@rasterlists$raw <- setNames(imgs, nm = object@samplenames)
  object@xdim <- xdim

  return(object)
}




