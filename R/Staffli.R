#' @include generics.R
NULL

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
    pixels.per.um = 'numeric',
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
    version = packageVersion(pkg = 'STutility')
  )

  return(object)
}


#' Subset a Seurat object containing Staffli image data
#'
#' Subsets a Seurat object containing Spatial Transcriptomics data while
#' making sure that the images and the spot coordinates are subsetted correctly.
#' If you use the default \code{\link{subset}} function there is a risk that images
#' are kept in the output Seurat object which will make the STUtility functions
#' crash.
#'
#' @param spots A vector of spots to keep
#' @param features A vector of features to keep
#' @param expression Logical expression indicating features/variables to keep
#' @param idents A vector of identity classes to keep
#' @param ... Extra parameters passed to WhichCells, such as slot, invert, or downsample
#'
#' @rdname SubsetSTData
#' @export
#'
SubsetSTData <- function (
  object,
  expression,
  spots = NULL,
  features = NULL,
  idents = NULL,
  ...
) {

  # Check that a Staffli object is present
  if (!"Staffli" %in% names(object@tools)) {
    stop("This Seurat object does not contain any Staffli image object", call. = FALSE)
  }
  # Obtain Staffli object
  st.object <- GetStaffli(object)

  if (!missing(x = expression)) {
    expression <- deparse(expr = substitute(expr = expression))
    object <- subset(object, subset = expression, features = features, cells = spots, idents = idents, ...)
  } else {
    object <- subset(object, features = features, cells = spots, idents = idents, ...)
  }

  # Check spots of new object
  kept.spots <- colnames(object)

  # Subset Staffli object
  st.meta_data <- st.object[[kept.spots, ]]
  st.object@meta.data <- st.meta_data
  samples <- unique(st.meta_data[, "sample"]) %>% as.numeric()

  convert_s <- 1:length(unique(samples))
  names(convert_s) <- samples
  st.object@meta.data$sample <- paste0(convert_s[st.object@meta.data$sample])
  new_samples <- 1:length(unique(samples))

  # Subset each slot in Staffli object
  if (length(st.object@imgs) > 0) {
    st.object@imgs <- st.object@imgs[samples]
  }
  if (length(st.object@rasterlists) > 0) {
    rl <- st.object@rasterlists
    rl <- lapply(rl, function(ls) {
      setNames(ls[samples], nm = new_samples)
    })
    st.object@rasterlists <- rl
  }
  if (length(st.object@scatter.data) > 0) {
    st.object@scatter.data <- subset(st.object@scatter.data, z %in% samples)
    st.object@scatter.data$z <- convert_s[st.object@scatter.data$z]
  }
  if (length(st.object@transformations) > 0) {
    st.object@transformations <- setNames(st.object@transformations[samples], nm = new_samples)
  }
  if (length(st.object@limits) > 0) {
    st.object@limits <- setNames(st.object@limits[samples], nm = new_samples)
  }
  if (length(st.object@dims) > 0) {
    st.object@dims <- setNames(st.object@dims[samples], nm = new_samples)
  }

  st.object@samplenames <- paste0(new_samples)
  object@tools$Staffli <- st.object
  return(object)
}

# TODO: fix sample names when merging

#' Merge two or more Seurat objects containing Staffli image data
#'
#' Merges Seurat objects containing Spatial Transcriptomics data while
#' making sure that the images and spot coordinates are correctly structures.
#' If you use the default \code{\link{merge}} function you will not be able
#' to use any of the STUtility visualization methods on the output object.
#'
#' @param x Seurat object
#' @param y Seurat object (or list of multiple Seurat obejctsa)
#' @param add.spot.ids A character vector of length(x = c(x, y)). Appends the corresponding
#' values to the start of each objects' spot names.
#' @param merge.data 	Merge the data slots instead of just merging the counts (which requires renormalization).
#' This is recommended if the same normalization approach was applied to all objects.
#' @param project Sets the project name for the Seurat object.
#' @param ... Arguments passed to other methods
#'
#' @rdname MergeSTData
#' @export
#'
MergeSTData <- function (
  x = NULL,
  y = NULL,
  add.spot.ids = NULL,
  merge.data = TRUE,
  idents = NULL,
  project = "SeuratProject",
  ...
) {

  # Check that a Staffli object is present
  if (!"Staffli" %in% names(x@tools)) {
    stop("The first Seurat object does not contain any Staffli image object", call. = FALSE)
  }
  i <- 1
  for (obj in y) {
    if (!"Staffli" %in% names(obj@tools)) {
      stop(paste0("Seurat object ", i, " does not contain any Staffli image object"), call. = FALSE)
    }
    i <- i + 1
  }

  # Obtain Staffli objects
  st.x <- GetStaffli(x)
  st.y <- lapply(y, GetStaffli)

  # Merge seurat data
  object <- merge(x = x, y = y, add.cell.ids = add.spot.ids, merge.data = merge.data, project = project)#, ...)

  # Combine Staffli objects into a list
  st.objects <- c(list(st.x), st.y)

  # Check that version matches
  versions <- unlist(lapply(st.objects, function(x) paste0(x@version)))
  if (length(unique(versions)) > 1) {
    warning(paste0("Different versions; ", paste(versions, collapse = ", "), " have been used to process the data"), call. = FALSE)
  }

  # and check that the same xdim has been used
  xdims.check <- unlist(lapply(st.objects, function(x) x@xdim))
  if (length(unique(xdims.check)) > 1) {
    warning(paste0("Different xdims have been used for the different objects; ", paste(xdims.check, collapse = ", "), ". \nAny loaded images will be removed and a defualt value of 400 pixels in width will be set."), call. = FALSE)
    xdim <- 400
  } else {
    xdim <- unique(xdims.check)
  }

  # Merge meta data
  unique.cols <- Reduce(intersect, lapply(st.objects, function(x) colnames(x[[]])))
  st.meta_data <- dplyr::bind_rows(lapply(st.object, `[[`), )
  st.meta_data <- st.meta_data[, unique.cols]
  rownames(st.meta_data) <- colnames(object)

  if (length(unique(xdims.check)) > 1 & all(c("warped_x", "warepd_y") %in% colnames(st.meta_data))) {
    st.meta_data[, c("warped_x", "warepd_y")] <- NULL
  }

  # Check that all objects have images
  imgs.class <- unlist(lapply(st.objects, function(x) class(x@imgs)))
  if (all(imgs.class == "character")) {
    imgs <- unlist(lapply(st.objects, function(x) x@imgs))
  } else {
    imgs <- NULL
  }

  # Define new sample ids
  samples <- c()
  n <- 0
  for (i in seq_along(st.objects)) {
    obj <- st.objects[[i]]
    unique.samples <- unique(obj[[, "sample", drop = TRUE]])
    convert.samples <-seq_along(unique.samples) + n
    names(convert.samples) <- unique.samples
    samples <- c(samples, convert.samples[obj[[, "sample", drop = TRUE]]])
    n <- convert.samples[length(convert.samples)]
  }
  samples <- samples %>% as.numeric() %>% as.character()

  # Replace new sample column of meta data
  st.meta_data[, "sample"] <- samples
  samplenames <- unique(samples)

  if (!length(unique(xdims.check)) > 1) {
    # Check for transformations
    transf.check <- unlist(lapply(st.objects, function(x) length(x@transformations)))
    if (any(transf.check == 0)) {
      transformations <- list()
    } else {
      transformations <- setNames(Reduce(c, lapply(st.objects, function(x) x@transformations)), nm  = samplenames)
    }

    # Check and combine rasterlists
    rasterlists.check <- Reduce(intersect, lapply(st.objects, function(x) names(x@rasterlists)))
    rasterlists <- list()
    if (!is.null(rasterlists.check)) {
      for (rst in rasterlists.check) {
        rasterlists[[rst]] <- setNames(Reduce(c, lapply(st.objects, function(x) x[rst])), nm  = samplenames)
      }
    }
  } else {
    transformations <- list()
    rasterlists <- list()
  }

  # Merge limits
  limits <- setNames(Reduce(c, lapply(st.objects, function(x) x@limits)), nm = samplenames)

  # Merge dims
  check.dims <- unlist(lapply(st.objects, function(x) length(x@dims)))
  if (any(check.dims == 0)) {
    dims <- list()
  } else {
    dims <- setNames(Reduce(c, lapply(st.objects, function(x) x@dims)), nm = samplenames)
  }

  # Merge platforms
  platforms <- unlist(lapply(st.objects, function(x) x@platforms))

  # Create Staffli object
  m <- CreateStaffliObject(imgs = imgs, meta.data = st.meta_data, xdim = xdim, platforms = platforms)
  m@rasterlists <- rasterlists
  m@transformations <- transformations
  m@limits <- limits
  m@dims <- dims

  object@tools$Staffli <- m
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
        if (!type %in% names(object@rasterlists) | !type %in% choices) stop(paste0("type '", type, "' not present in Seurat object"), call. = FALSE)
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




