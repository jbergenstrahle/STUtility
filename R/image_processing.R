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
    image_scale(imgs[[i]], paste0(xdim))
  }), nm = samplenames)

  object@tools$pointers <- imgs
  object@tools$rasters <- setNames(lapply(lapply(seq_along(imgs), function(i) {image_annotate(imgs[[i]], text = samplenames[i], size = round(xdim/10))}), as.raster), nm = samplenames)
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
      image_scale(images[[i]], paste0(xdim))
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


#' Masks the background of HE images stored in a Seurat object
#'
#' Algorithm:
#' 1. Blur image with median filter
#' 2. Equalize image colors using an empirical cumulative distribution function
#' 3. Convert image into SLIC superpixels using the \code{\link{slic}} function
#' 4. Convert image to CIELAB colorspace and perform kmeans clustering (k = 2; inside/outside tissue) to segment image
#' 5. Split segmented image into objects and filter out objects with a small area
#' 6. Keep objects which overlaps with adjusted pixel coordinates
#'
#' @param obejct Seurat object
#' @param median.blur Pixel size of median filter used to blur HE images prior to image segmentation
#' @param object.size.th Threshold used to filter out objects in the background which are not part of the
#' tissue, e.g. bubbles or other debris
#' @param verbose Print messages
#'
#' @inheritParams slic
#'
#' @importFrom imager magick2cimg medianblur sRGBtoLab as.cimg split_connected
#' @importFrom dplyr select
#' @importFrom magrittr %>%
#' @importFrom stats kmeans
#'
#' @return A Seurat object with masked HE images
#'
#' @export

MaskImages <- function(
  object,
  compactness = 1,
  median.blur = 10,
  object.size.th = 0.01,
  verbose = FALSE
) {

  object@tools$rasters <- setNames(lapply(seq_along(object@tools$pointers), function(i) {
    im <- magick2cimg(object@tools$pointers[[i]]) %>% medianblur(median.blur)
    f <- ecdf(im)
    im <- f(im) %>% as.cimg(dim = dim(im))
    if (verbose) {
      cat(paste0("Loaded image ", i, "\n"))
      cat(paste0("Running SLIC algorithm \n"))
    }
    out <- slic(im, nS = object@tools$xdim*1.5, compactness)

    d <- sRGBtoLab(out$sp) %>% as.data.frame(wide = "c") %>%
      select(-x,-y)

    km <- kmeans(d, 2)
    seg <- as.cimg(km$cluster - 1, dim = c(dim(im)[1:2], 1, 1)) %>% medianblur(median.blur)

    # Extract pixel sets
    px <- seg > 0.5
    sp <- split_connected(px)

    # Check object sizes
    size.check <- max(unlist(lapply(sp, function(x) sum(x > 0))))/length(seg) > object.size.th

    if (length(sp) == 0 | !size.check) {
      seg <- 1 - seg
      px <- seg > 0.5
      sp <- split_connected(px)
    }

    # Colect pixel coordinates for spots in meta.data slot
    dims.raw <- as.numeric(object@tools$dims[[i]][2:3])
    dims.scaled <- dim(object@tools$rasters[[i]])
    sf.xy <- dims.raw/dims.scaled
    pixel_xy <- sapply(subset(object[[]], sample == paste0(i))[, c("pixel_x", "pixel_y")]/sf.xy, round)
    pixel_coords <- paste(pixel_xy[, 1], pixel_xy[, 2], sep = "x")

    # Select pixelsets overlapping with pixel coordinates in meta.data slot
    inds_inside_tissue <- data.frame()

    inds_inside_tissue <- do.call(rbind, lapply(sp, function(px) {
      m <- px[, ]; colnames(m) <- 1:ncol(m); rownames(m) <- 1:nrow(m)

      inds <- data.frame(label = as.numeric(m),
                         y.id = rep(1:nrow(m), ncol(m)),
                         x.id = rep(1:ncol(m), each = nrow(m)))
      inds$idx <- 1:nrow(inds)
      inds.under.tissue <- subset(inds, label == 1)

      inds_coords <- paste(inds.under.tissue[, "x.id"], inds.under.tissue[, "y.id"], sep = "x")
      if (length(intersect(pixel_coords, inds_coords)) > 0) {
        return(subset(inds, label == 0))
      } else {
        return(NULL)
      }
    }))

    if (verbose) {
      cat(paste0("Masking background in image ", i , "\n"))
    }

    rst <- t(object@tools$rasters[[i]])
    rst[inds_inside_tissue$idx] <- "00000000"
    rst <- t(rst)
    return(rst)
  }), nm = names(object@tools$rasters))

  return(object)
}


#' SLIC algorithm
#'
#' Converts image into SLIC superpixels
#'
#' @param im Image of class "cimg"
#' @param nS Number of superpixels to return
#' @param compactness Controls scaling ratio for pixel values
#' @param ... Parameters passed to kmeans
#'
#' @importFrom purrr map_dbl
#' @importFrom imager imsplit LabtosRGB

slic <- function(
  im,
  nS,
  compactness = 1,
  ...
) {
  #If image is in colour, convert to CIELAB
  if (spectrum(im) == 3) im <- sRGBtoLab(im)

  #The pixel coordinates vary over 1...width(im) and 1...height(im)
  #Pixel values can be over a widely different range
  #We need our features to have similar scales, so
  #we compute relative scales of spatial dimensions to colour dimensions
  sc.spat <- (dim(im)[1:2]*.28) %>% max #Scale of spatial dimensions
  sc.col <- imsplit(im, "c") %>% map_dbl(sd) %>% max

  #Scaling ratio for pixel values
  rat <- (sc.spat/sc.col)/(compactness*10)


  X <- as.data.frame(im*rat, wide="c") %>% as.matrix
  #Generate initial centers from a grid
  ind <- round(seq(1, nPix(im)/spectrum(im),l=nS))
  #Run k-means
  km <- kmeans(X, X[ind, ], ...)

  #Return segmentation as image (pixel values index cluster)
  seg <- as.cimg(km$cluster,dim=c(dim(im)[1:2], 1, 1))
  #Superpixel image: each pixel is given the colour of the superpixel it belongs to
  sp <- map(1:spectrum(im),~ km$centers[km$cluster, 2+.]) %>% do.call(c, .) %>% as.cimg(dim = dim(im))
  #Correct for ratio
  sp <- sp/rat
  if (spectrum(im) == 3)
  {
    #Convert back to RGB
    sp <- LabtosRGB(sp)
  }
  list(km = km, seg = seg, sp = sp)
}
