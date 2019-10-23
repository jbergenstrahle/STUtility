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
    imgs = 'ANY',
    rasterlists = 'list',
    scatter.data = 'data.frame',
    transformations = 'list',
    meta.data = 'data.frame',
    xdim = 'numeric',
    limits = 'list',
    dims = 'list',
    platforms = 'ANY',
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
  platforms = NULL
) {
  if (!all(c("x", "y", "sample") %in% colnames(meta.data))) stop(paste0("Invalid meta.data; one of 'x', 'y' or 'sample' is missing"), call. = FALSE)
  samples <- unique(meta.data[, "sample"])

  # Define platforms if NULL
  platforms <- platforms %||% rep("Visium", length(x = samples))

  # Check that platforms match samples
  if (length(x = samples) != length(x = platforms)) stop("Length of platforms does not match the number of samples", call. = FALSE)

  limits <- list()
  for (i in seq_along(platforms)) {
    platform <- platforms[i]
    if (platform == 'Visium') {
      limits[[i]] <- c(128, 78)
    } else if (platform == '1k') {
      limits[[i]] <- c(33, 35)
    } else if (platform == '2k') {
      limits[[i]] <- c(67, 64)
    } else {
      stop(paste0(platform, " is not a valid option ... \n"), call. = FALSE)
    }
  }

  object <- new (
    Class = 'Staffli',
    imgs = imgs,
    meta.data = meta.data,
    xdim = xdim,
    limits = setNames(limits, samples),
    platforms = platforms,
    samplenames = samples,
    version = packageVersion(pkg = 'STUtility')
  )

  return(object)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Staffli methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setGeneric("iminfo", function(object) {
  standardGeneric("iminfo")
})


setMethod (
  f = "iminfo",
  signature = "Staffli",
  definition = function(object) {
    object@dims
  }
)


setGeneric("scaled.imdims", function(object, type = "raw") {
  standardGeneric("scaled.imdims")
})


setMethod (
  f = "scaled.imdims",
  signature = "Staffli",
  definition = function(object, type = "raw") {
    lapply(object[type], dim)
  }
)


setMethod (
  f = "names",
  signature = "Staffli",
  definition = function(x) {
    x@samplenames
  }
)


setGeneric("rasterlists", function(object) {
  standardGeneric("rasterlists")
})


setMethod (
  f = "rasterlists",
  signature = "Staffli",
  definition = function(object) {
    names(object@rasterlists)
  }
)


setMethod (
  f = "rasterlists",
  signature = "Seurat",
  definition = function(object) {
    if (!"Staffli" %in% names(object@tools)) stop("Staffli not present in Seurat object ... \n", call. = FALSE)
    st.object <- object@tools$Staffli
    names(st.object@rasterlists)
  }
)


setGeneric("samplenames", function(object) {
  standardGeneric("samplenames")
})


setMethod (
  f = "samplenames",
  signature = "Staffli",
  definition = function(object) {
    object@samplenames
  }
)


setMethod (
  f = "samplenames",
  signature = "Seurat",
  definition = function(object) {
    if (!"Staffli" %in% names(object@tools)) stop("Staffli not present in Seurat object ... \n", call. = FALSE)
    st.object <- object@tools$Staffli
    st.object@samplenames
  }
)


setGeneric("GetStaffli", function(object) {
  standardGeneric("GetStaffli")
})


setMethod (
  f = "GetStaffli",
  signature = "Seurat",
  definition = function(object) {
    if (!"Staffli" %in% names(object@tools)) stop("Staffli not present in Seurat object ... \n", call. = FALSE)
    object@tools$Staffli
  }
)


setMethod (
  f = "[[",
  signature = "Staffli",
  definition = function(x, i, j, drop = F) {
    x@meta.data[i, j, drop]
  }
)


setMethod (
  f = "[[<-",
  signature = "Staffli",
  definition = function(x, i, j, ..., value) {
    x@meta.data[i, j] <- value
    return(x)
  }
)


setMethod (
  f = "[",
  signature = "Staffli",
  definition = function(x, i) {
    x@rasterlists[[i]]
  }
)


setMethod (
  f = "[<-",
  signature = "Staffli",
  definition = function(x, i, ..., value) {
    if (length(x@rasterlists[[i]]) != length(value) | class(x@rasterlists[[i]]) != class(value)) {
      stop("Invalid class or the lists are of different lengths", call. = FALSE)
    }
    x@rasterlists[[i]] <- value
    return(x)
  }
)


setMethod (
  f = "show",
  signature = "Staffli",
  definition = function(object) {
    cat("An object of class", class(x = object), "\n")
    cat(
      nrow(object@meta.data),
      'spots across',
      length(unique(object@meta.data[, "sample"])),
      'samples. \n'
    )
    if (length(object@rasterlists) > 0) {
      cat(
        '\nAvailable image representations: \n  ',
        paste(names(object@rasterlists), collapse = ", "),
        '\n'
      )
    }
  }
)


setMethod (
  f = "plot",
  signature = "Staffli",
  definition = function(x, type = NULL, ...) {
    object <- x
    ncols <- ceiling(sqrt(length(x = names(object)))); nrows <- ceiling(length(x = names(object))/ncols)
    if (length(object@rasterlists) == 0) {
      par(mfrow = c(nrows, ncols))
      for (s in names(object)) {
        d <- subset(object[[]], sample == s)
        plot(d[, "x"], object@limits[[s]][2] - d[, "y"], xlim = c(0, object@limits[[s]][1]), ylim = c(0, object@limits[[s]][2]), ann = FALSE)
      }
    } else {
      # Check that type is OK
      choices <- c("processed", "masked", "raw", "processed.masks", "masked.masks")
      if (!is.null(type)) {
        if (!type %in% names(st.object@rasterlists) | !type %in% choices) stop(paste0("type '", type, "' not present in Seurat object"), call. = FALSE)
      }

      type <- type %||% {
        choices <- c("processed", "masked", "raw", "processed.masks", "masked.masks")
        match.arg(arg = choices, choices = names(object@rasterlists), several.ok = TRUE)[1]
      }

      if (type == "raw") {
        xy <- c("pixel_x", "pixel_y")
      } else if (type %in% c("masked", "masked.masks")) {
        if (!"masked" %in% names(object@rasterlists)) stop("Masked images not available in Staffli object", call. = FALSE)
        xy <- c("pixel_x", "pixel_y")
      } else if (type %in% c("processed", "processed.masks")) {
        if (!"processed" %in% names(object@rasterlists)) stop("Masked images not available in Staffli object", call. = FALSE)
        xy <- c("warped_x", "warped_y")
      }
      par(mfrow = c(nrows, ncols), mar = c(0, 0, 0, 0))
      for (s in names(object)) {
        d <- subset(object[[]], sample == s)
        im <- object[type][[s]] %>% as.cimg()
        plot(im, axes = FALSE)
        points(d[, xy]/(iminfo(object)[[s]]$width/object@xdim), ...)
      }
    }
  }
)




