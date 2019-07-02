#' Function used to read HE images in jpeg format
#'
#' @param object Seurat object
#' @param image.paths Paths to HE images. This is only required if image paths are missing in the Seurat object.
#' @param xdim Sets the pixel width for scaling, e.g. 400 (maximum allowed width is 1000 pixels)
#'
#' @importFrom magick image_read
#'
#' @export

LoadImages <- function(
  object,
  image.paths = NULL,
  xdim = 400
) {
  # Check that image paths are present
  if (!"imgs" %in% names(object@tools)) {
    stop(paste0("Image paths are not present in Seurat object. Provide image.paths manually."))
  }

  # Check that the image with is no more than 2000 pixels
  if (xdim > 2000) stop("xdim cannot be larger than 2000")

  # Check that a sample column exists in meta data
  if (!"sample" %in% colnames(object[[]])) {
    if (length(object@tools$imgs) > 1 & is.null(image.paths)) {
      stop(paste0("Column 'sample' is missing from meta.data and is required when more than one image path (#", length(object@tools$imgs), ") has been provided."), call. = F)
    } else {
      if (!is.null(image.paths)) {
        if (length(image.paths) == 1) {
          object@tools <- list(imgs = image.paths)
        } else {
          stop("Provided image.paths but there are more images than the number of samples (1). \nIf you have more than 1 sample, make sure to provide a 'sample' column in the meta.data slot.")
        }
      }
    }
    samplenames <- "1"
  } else {
    if (!is.null(image.paths)) {
      if (!length(x = unique(object[["sample"]]) == length(x = image.paths))) {
        stop(paste0("Number of images (",
                    length(image.paths),
                    ") does not match the number of samples (",
                    length(x = unique(object[["sample"]]), ")")), call. = F)
      } else {
        object@tools <- list(imgs = image.paths)
      }
    }

    check_group.var <- length(unique(object[["sample", drop = T]])) == length(object@tools$imgs)

    if (!check_group.var) stop(paste0("Number of elements in 'sample' column does not match the number of images provided \n\n",
                                      "elements: \n", paste(unique(object[["sample", drop = T]]), collapse = ", "), "\n\n",
                                      "images: \n", paste(object@tools$imgs, collapse = ", \n")), call. = F)
    samplenames <- unique(object[["sample", drop = T]])
  }

  imgs <- c()
  for (path in object@tools$imgs) {
    imgs <- c(imgs, image_read(path))
  }

  object@tools$dims <- setNames(lapply(imgs, image_info), nm = samplenames)

  imgs <- setNames(lapply(seq_along(imgs), function(i) {
    image_annotate(image_scale(imgs[[i]], paste0(xdim)), text = samplenames[i], size = round(xdim/10))
  }), nm = samplenames)

  object@tools$pointers <- imgs
  object@tools$rasters <- lapply(imgs, as.raster)
  object@tools$xdim <- xdim

  return(object)
}



#' Function used to plot HE images obtained with \code{\link{LoadImages}}
#'
#' @param object Seurat object
#' @param index Image index
#' @param method Specify display method (raster or viewer)
#' @param ncols Number of columns in output grid of images
#' @inheritParams LoadImages
#'
#' @importFrom magick image_append image_annotate image_scale
#'
#' @export

ImagePlot <- function(
  object,
  index = NULL,
  method = "viewer",
  ncols = NULL
) {

  if (!"pointers" %in% names(object@tools)) {
    stop(paste0("Image paths are not present in Seurat object. Run LoadImages() before plotting."))
  }

  images <- object@tools$pointers
  check_pointers <- any(sapply(lapply(images, function(images) {
    try(image_info(images))
  }), class) == "try-error")

  if (check_pointers) {
    warning(paste0("Image pointer is dead. You cannot save or cache image objects between R sessions. \n",
                   "Rerun ImageRead if you want to set the image sizes manually. \n",
                   "Setting image size to '", object@tools$xdim, "' pixels "), call. = F)
    images <- c()
    for (path in object@tools$imgs) {
      images <- c(images, image_read(path))
    }
    images <- setNames(lapply(seq_along(images), function(i) {
      image_annotate(image_scale(images[[i]], paste0(xdim)), text = labels[i], size = round(xdim/10))
    }), nm = names(images))
    object@tools$pointers <- images
  }

  if (is.null(index)) {
    ncols <- ncols %||% round(sqrt(length(x = images)))
    nrows <- round(length(x = images)/ncols)
    stack <- c()
    for (i in 1:nrows) {
      i <- i - 1
      stack <- c(stack, image_append(Reduce(c, images[(i*ncols + 1):(i*ncols + ncols)])))
    }

    final_img <- image_append(Reduce(c, stack), stack = T)
  } else {
    if (!index %in% 1:length(images)) stop("Image index out of bounds", call. = F)
    final_img <-images[[index]]
  }

  if (method == "raster") {
    plot(as.raster(final_img))
  } else if (method == "viewer") {
    print(final_img)
  } else {
    stop("Invalid display method", call. = F)
  }

  return(object)
}
