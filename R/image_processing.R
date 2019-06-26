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
  group.var = "sample",
  xdim = 200
) {

  if (xdim > 1000) stop("xdim cannot be larger than 1000")

  if (!group.var %in% colnames(object[[]])) stop(paste0("group.var", group.var, " not found in meta.data slot"), call. = F)

  if (!is.null(image.paths)) {
    if (!length(x = unique(object[["sample"]]) == length(x = image.paths))) {
      stop(paste0("Number of images (",
                  length(image.paths),
                  ") does not match the number of samples (",
                  length(x = unique(object[["sample"]]), ")")))
    } else {
      object@tools <- list(imgs = image.paths)
    }
  }

  check_group.var <- length(unique(object[[group.var, drop = T]])) == length(object@tools$imgs)

  if (!check_group.var) stop(paste0("number of elements in group.var does not match the number of images provided \n\n",
                                    "elements: \n", paste(unique(object[[group.var, drop = T]]), collapse = ", "), "\n\n",
                                    "images: \n", paste(object@tools$imgs, collapse = ", \n")))


  if (!"imgs" %in% names(object@tools)) {
    stop(paste0("Image paths are not present in Seurat object. Provide image.paths manually."))
  }

  imgs <- c()
  for (path in object@tools$imgs) {
    imgs <- c(imgs, image_read(path))
  }

  object@tools$dims <- setNames(lapply(imgs, image_info), nm = unique(object[[group.var, drop = T]]))

  imgs <- setNames(lapply(seq_along(imgs), function(i) {
    image_annotate(image_scale(imgs[[i]], paste0(xdim)), text = unique(object[[group.var, drop = T]])[i], size = round(xdim/10))
  }), nm = unique(object[[group.var, drop = T]]))

  object@tools$pointers <- imgs
  object@tools$rasters <- lapply(imgs, as.raster)

  return(object)
}



#' Function used to plot HE images obtained with \code{\link{LoadImages}}
#'
#' @param object Seurat object
#' @param index Image index
#' @inheritParams LoadImages
#'
#' @importFrom magick image_append image_annotate image_scale
#'
#' @export

ImagePlot <- function(
  object,
  index = NULL
) {

  if (!"imgs" %in% names(object@tools)) {
    stop(paste0("Image paths are not present in Seurat object. Run LoadImages before plotting."))
  }

  images <- object@tools$pointers
  check_pointers <- any(sapply(lapply(images, function(images) {
    try(image_info(images))
  }), class) == "try-error")

  if (check_pointers) {
    warning(paste0("Image pointer is dead. You cannot save or cache image objects between R sessions. \n",
                   "Rerun ImageRead if you want to set the image sizes manually. \n",
                   "Setting image size to '200' pixels "), call. = F)
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
    ncols <- round(sqrt(length(x = images)))
    nrows <- round(length(x = images)/ncols)
    stack <- c()
    for (i in 1:nrows) {
      i <- i - 1
      stack <- c(stack, image_append(Reduce(c, images[(i*ncols + 1):(i*ncols + ncols)])))
    }
    print(image_append(Reduce(c, stack), stack = T))
  } else {
    if (!index %in% 1:length(images)) stop("Image index out of bounds", call. = F)
    print(images[[index]])
  }

  return(object)
}
