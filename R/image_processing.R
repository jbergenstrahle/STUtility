#' Function used to read HE images in jpeg format
#'
#' @param object Seurat object
#' @param image.paths Paths to HE images. This is only required if image paths are missing in the Seurat object.
#' @param xdim Sets the pixel width for scaling, e.g. 400 (maximum allowed width is 1000 pixels)
#' @param verbose Print messages
#'
#' @importFrom magick image_read
#'
#' @export

LoadImages <- function (
  object,
  image.paths = NULL,
  xdim = 400,
  verbose = FALSE,
  time.resolve = TRUE
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
  dims <- list()
  for (i in seq_along(object@tools$imgs)) {
    path <- object@tools$imgs[i]
    if (verbose) cat("Reading ", path , " for sample ", unique(object[["sample", drop = T]])[i], " ... \n", sep = "")
    im <- image_read(path)
    dims <- c(dims, list(image_info(im)))
    if (verbose) {
      info <- dims[[i]]
      width <- as.numeric(info[2]); height <- as.numeric(info[3])
      ydim <- round(height/(width/xdim))
      cat("Scaling down sample ", unique(object[["sample", drop = T]])[i], " image from ", paste(width, height, sep = "x"), " pixels to ", paste(xdim, ydim, sep = "x"), " pixels \n", sep = "")
    }
    im <- image_scale(im, paste0(xdim))
    #tmpf <- tempfile()
    #image_write(im, path = tmpf)
    #im <- image_read(tmpf)
    imgs <- c(imgs, list(as.raster(im)))
    if(time.resolve == TRUE){
      gc()
      sleepy(5)
    }
  }

  object@tools$dims <- setNames(dims, nm = samplenames)
  object@tools$raw <- setNames(imgs, nm = samplenames)
  object@tools$xdim <- xdim

  return(object)
}


# TODO: label raster images, why are images black when put into viewer? Compare with warpimages return ...

#' Function used to plot HE images obtained with \code{\link{LoadImages}}
#'
#' @param object Seurat object
#' @param indices Image indices to select
#' @param type Specify which image to display [options: "raw", "masked" or "processed"]
#' @param method Specify display method (raster or viewer).
#' @param ncols Number of columns in output grid of images
#' @param annotate Will put a unique id in the top left corner
#' @inheritParams LoadImages
#'
#' @importFrom magick image_append image_annotate image_scale
#'
#' @export

ImagePlot <- function (
  object,
  indices = NULL,
  type = NULL,
  method = "viewer",
  ncols = NULL,
  annotate = TRUE
) {

  if (!is.null(type)) {
    if (!type %in% names(object@tools)) stop(paste0("Invalid type ", type), call. = FALSE)
  }

  type <- type %||% {
    choices <- c("processed", "masked", "raw", "processed.masks", "masked.masks")
    choices[min(as.integer(na.omit(match(names(object@tools), choices))))]
  }

  images <- object@tools[[type]]
  indices <- indices %||% {
    seq_along(images)
  }

  if (any(!indices %in% 1:length(images))) stop("Image indices out of bounds: ", paste(setdiff(indices, 1:length(images)), collapse = ", "), call. = F)
  images <- images[indices]

  images <- lapply(images, image_read)
  if (annotate) {
    images <- setNames(lapply(seq_along(images ), function(i) {image_annotate(images[[i]], text = names(images)[i], size = round(object@tools$xdim/10))}), nm = names(images))
  }
  ncols <- ncols %||% round(sqrt(length(x = images)))
  nrows <- ceiling(length(x = images)/ncols)

  if (method == "viewer") {
    stack <- c()
    for (i in 1:nrows) {
      i <- i - 1
      stack <- c(stack, image_append(Reduce(c, images[(i*ncols + 1):(i*ncols + ncols)])))
    }

    final_img <- image_append(Reduce(c, stack), stack = T)
    print(final_img)
  } else if (method == "raster"){
    #par(mar = c(0, 0.2, 0, 0.2), mfrow = c(nrows, ncols))
    layout.matrix <- t(matrix(c(1:length(images), rep(0, nrows*ncols - length(images))), nrow = ncols, ncol = nrows))
    layout(mat = layout.matrix)

    for (rst in lapply(images, as.raster)) {
      par(mar = c(0, 0.2, 0, 0.2))
      plot(rst)
    }
    par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)
  } else {
    stop(paste0("Invalid display method: ", method), call. = F)
  }
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
#' @param iso.blur Sigma value (pixels) for isoblurring of HE images prior to image segmentation
#' @param verbose Print messages
#'
#' @inheritParams slic
#'
#' @importFrom imager magick2cimg medianblur sRGBtoLab as.cimg split_connected add imsplit imappend RGBtoHSV blur_anisotropic HSVtoRGB add
#' @importFrom magick image_read
#' @importFrom dplyr select summarize
#' @importFrom magrittr %>%
#' @importFrom stats kmeans
#' @importFrom purrr modify_at
#'
#' @return A Seurat object with masked HE images
#'
#' @export

MaskImages <- function (
  object,
  compactness = 1,
  verbose = FALSE
) {

  rasters <- list()
  masks <- list()
  centers <- list()

  for (i in seq_along(object@tools$raw)) {
    imr <- image_read(object@tools$raw[[i]])

    # segmentation tests
    im <- magick2cimg(imr)
    im <- threshold(im)
    im[, , , 2] <- TRUE; im[, , , 3] <- TRUE
    im <- isoblur(im, 2)

    if (verbose) {
        cat(paste0("Loaded image ", i, "\n"))
        cat(paste0("Running SLIC algorithm \n"))
    }
    out <- slic(im, nS = object@tools$xdim*1.5, compactness)
    out <- out^4
    d <- sRGBtoLab(out) %>% as.data.frame(wide = "c") %>%
      select(-x, -y)

    km <- kmeans(d, 2)
    seg <- as.cimg(km$cluster - 1, dim = c(dim(im)[1:2], 1, 1)) %>% medianblur(20) %>% threshold()

    # Chech that at least one masked region is found
    sp <- split_connected(seg)
    if (length(sp) == 0) sp <- imlist(seg)
    if (length(sp) == 0) stop(paste0("Masking failed for image ", i), call. = FALSE)

    # Colect pixel coordinates for spots in meta.data slot
    dims.raw <- as.numeric(object@tools$dims[[i]][2:3])
    dims.scaled <- dim(object@tools$raw[[i]])
    sf.xy <- dims.raw/rev(dims.scaled)
    pixel_xy <- sapply(subset(object[[]], sample == paste0(i))[, c("pixel_x", "pixel_y")]/sf.xy, round)
    pixel_coords <- paste(pixel_xy[, 1], pixel_xy[, 2], sep = "x")

    # Check object sizes
    #size.check <- ifelse(length(sp) > 0, max(unlist(lapply(sp, function(x) sum(x > 0))))/length(seg) > object.size.th, FALSE)
    pxs.T <- setNames(data.frame(which(seg == TRUE, arr.ind = TRUE)[, 1:2]), nm = c("x", "y"))
    pxs.F <- setNames(data.frame(which(seg == FALSE, arr.ind = TRUE)[, 1:2]), nm = c("x", "y"))
    seg.pxs.T <- paste(pxs.T$x, pxs.T$y, sep = "x")
    seg.pxs.F <- paste(pxs.F$x, pxs.F$y, sep = "x")
    selected.coords.check <- as.logical(which.max(c(length(intersect(seg.pxs.F, pixel_coords)), length(intersect(seg.pxs.T, pixel_coords)))) - 1)

    if (!selected.coords.check) {
      seg <- !seg
      sp <- split_connected(seg)
      if (length(sp) == 0) sp <- imlist(seg)
    }

    # Select pixelsets overlapping with pixel coordinates in meta.data slot
    inds_inside_tissue <- data.frame()
    keep.sp <- list()

    for (j in seq_along(sp)) {
      px <- sp[[j]]
      m <- t(px[, ]); colnames(m) <- 1:ncol(m); rownames(m) <- 1:nrow(m)

      inds <- data.frame(label = as.numeric(m),
                         y.id = rep(1:nrow(m), ncol(m)),
                         x.id = rep(1:ncol(m), each = nrow(m)))
      inds$idx <- 1:nrow(inds)
      inds.under.tissue <- subset(inds, label == 1)
      inds_coords <- paste(inds.under.tissue[, "x.id"], inds.under.tissue[, "y.id"], sep = "x")
      if (length(intersect(pixel_coords, inds_coords)) > 0) {
        inds_inside_tissue <- rbind(inds_inside_tissue, subset(inds, label == 1))
        keep.sp <- c(keep.sp, list(px))
      } else {
        next
      }
    }

    masks[[i]] <- as.raster(add(keep.sp))

    if (verbose) {
      cat(paste0("Masking background in image ", i , "\n"))
    }

    rst <- object@tools$raw[[i]]
    rst[(1:length(rst))[-inds_inside_tissue$idx]] <- "#FFFFFF"

    #rst <- t(rst)
    rasters[[i]] <- rst

    # Find tissue center
    center <- which(seg[, , 1, 1] > 0, arr.ind = T) %>% as.data.frame() %>% summarize(center.x = mean(row), center.y = mean(col))
    centers[[i]] <- center
  }

  object@tools$masked <- setNames(rasters, nm = names(object@tools$raw))
  object@tools$masked.masks <- setNames(masks, nm = names(object@tools$raw))
  object@tools$centers <- setNames(centers, nm = names(object@tools$raw))

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
#' @importFrom purrr map_dbl map
#' @importFrom imager imsplit LabtosRGB sRGBtoLab spectrum nPix as.cimg
#' @importFrom stats kmeans

slic <- function (
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


  X <- as.data.frame(im*rat, wide = "c") %>% as.matrix
  #Generate initial centers from a grid
  ind <- round(seq(1, nPix(im)/spectrum(im), l = nS))
  #Run k-means
  km <- suppressWarnings({kmeans(X, X[ind, ], ...)})

  #Superpixel image: each pixel is given the colour of the superpixel it belongs to
  sp <- map(1:spectrum(im), ~ km$centers[km$cluster, 2+.]) %>% do.call(c, .) %>% as.cimg(dim = dim(im))
  #Correct for ratio
  sp <- sp/rat
  if (spectrum(im) == 3)
  {
    #Convert back to RGB
    sp <- LabtosRGB(sp)
  }

  return(sp)
}





#' Create transformation function
#'
#' Creates a function that takes x,y values as input and applies a transformation
#'
#' @param tr Forward transformation matrix
#' @param forward Logical: sets algorithm to 'forwar', otherwise 'backward'
#'

generate.map.rot <- function (
  tr,
  forward = FALSE
) {

  if (forward) {
    map.rot <- function (
      x,
      y
    ) {
      xy <- t(cbind(x, y, 1))
      #tr <- combine.tr(center.cur, center.new, alpha = angle, mirror.x, mirror.y)
      xytr <- t(tr%*%xy)
      list(x = xytr[, 1], y = xytr[, 2])
    }
  } else {
    map.rot <- function (
      x,
      y
    ) {
      xy <- t(cbind(x, y, 1))
      #tr <- combine.tr(center.cur, center.new, alpha = angle, mirror.x, mirror.y)
      xytr <- t(solve(tr)%*%xy)
      list(x = xytr[, 1], y = xytr[, 2])
    }
  }

  return(map.rot)
}

#' Function used to add whitespace to image
#'
#' This function adds whitespace to a raster image and forces the height/width
#' to be equal to the image diagonal. This way you can ensure that no part of the tissue
#' is cropped after a rotation transformation
#'
#' @param rst Raster image

add.margins <- function(
  rst
) {
  d <- round(sqrt(ncol(rst)^2 + nrow(rst)^2))
  mar.x <- round((d - ncol(rst))/2); mar.y <- round((d - nrow(rst))/2)
  m <- as.raster(matrix(data = "#FFFFFF", nrow = d, ncol = d))
  m[(mar.y + 1):(mar.y + nrow(rst)), (mar.x + 1):(mar.x + ncol(rst))] <- rst
  m <- t(m)
}


#' Apply warping of x, y coordinates using a affine transformation function
#'
#' @param im Raster image
#' @param map.rot Affine transformation function, see \code{\link{generate.map.rot}}
#'
#' @importFrom imager as.cimg imwarp
#' @importFrom grDevices as.raster

Warp <- function (
  im,
  map.rot,
  mask = FALSE
) {
  im <- imwarp(as.cimg(im), map = map.rot, direction = "backward", interpolation = "cubic")
  if (!mask) {
    copy.im <- as.raster(matrix(data = "#FFFFFF", nrow = ncol(im), ncol = nrow(im)))
    copy.im <- imwarp(as.cimg(copy.im), map = map.rot, direction = "backward", interpolation = "cubic")
    inds <- which(copy.im != 255)
    im[inds] <- 255
    imrst <- as.raster(im)
    tab.im <- table(imrst)
    if (length(tab.im) > 2) {
      imrst[imrst == names(tab.im)[which.max(tab.im)]] <- "#FFFFFF"
    }
  } else {
    imrst <- as.raster(t(ifelse(im[, , 1, 1] > 100, 1, 0)))
  }
  return(imrst)
}

# TODO: give examples, fix transformation of pixel coordinates

#' Warps images using various transformations
#'
#' @param object Seurat object
#' @param transforms List of arguments passed to warp function
#' @param verbose Print messages
#'
#' @return Seurat object with processed imaged
#'
#' @export

WarpImages <- function (
  object,
  transforms,
  verbose
) {
  if (!"masked" %in% names(object@tools)) stop(paste0("Masked images are not present in Seurat object"), call. = FALSE)
  if (!any(c("pixel_x", "pixel_y") %in% colnames(object[[]]))) stop(paste0("Pixel coordinates are missing in Seurat object"), call. = FALSE)

  if (!all(names(transforms) %in% names(object@tools$masked))) stop(paste0("transforms does not match the image labels"))

  processed.images <- setNames(ifelse(rep("processed" %in% names(object@tools), length(object@tools$imgs)), object@tools$processed, object@tools$masked), nm = names(object@tools$masked))
  masks <- object@tools$masked.masks
  processed.masks <- object@tools$masked.masks
  if (all(c("warped_x", "warped_y") %in% colnames(object[[]]))) {
    warped_coords <- object[[c("warped_x", "warped_y", "sample")]]
    spots <- rownames(subset(warped_coords, sample %in% names(transforms)))
    pxset <- object[[c("pixel_x", "pixel_y")]]
    warped_coords[spots, 1:2] <- pxset[spots, ]
    warped_coords <- warped_coords[, 1:2]
  } else {
    warped_coords <- object[[c("pixel_x", "pixel_y")]]
  }


  for (i in names(transforms)) {

    if (verbose) cat(paste0("Loading masked image for sample ", i, " ... \n"))
    m <- object@tools$masked[[i]]

    args <- transforms[[i]]
    center.x <- args[["center.x"]] %||% FALSE
    center.y <- args[["center.y"]] %||% FALSE
    angle <- args[["angle"]] %||% 0
    mirror.x <- args[["mirror.x"]] %||% FALSE
    mirror.y <- args[["mirror.y"]] %||% FALSE

    tr <- combine.tr(as.numeric(object@tools$centers[[i]]), rev(dim(m)/2), alpha = angle, mirror.x = mirror.x, mirror.y = mirror.y)

    map.rot.backward <- generate.map.rot(tr)
    map.rot.forward <- generate.map.rot(tr, forward = TRUE)

    # Obtain scale factors
    dims.raw <- as.numeric(object@tools$dims[[i]][2:3])
    dims.scaled <- dim(object@tools$raw[[i]])
    sf.xy <- dims.raw/rev(dims.scaled)
    pixel_xy <- subset(object[[]], sample == paste0(i))[, c("pixel_x", "pixel_y")]/sf.xy

    # Warp pixels
    if (verbose) cat(paste0("Warping pixel coordinates for ", i, " ... \n"))
    warped_xy <- sapply(setNames(as.data.frame(do.call(cbind, map.rot.forward(pixel_xy$pixel_x, pixel_xy$pixel_y))), nm = c("warped_x", "warped_y"))*sf.xy, round, digits = 1)
    warped_coords[rownames(pixel_xy), 1:2] <- warped_xy

    if (verbose) cat(paste0("Warping image for ", i, " ... \n"))
    processed.images[[i]] <- Warp(m, map.rot.backward)
    msk <- masks[[i]]
    if (verbose) cat(paste0("Warping image mask for ", i, " ... \n"))
    processed.masks[[i]] <- Warp(msk, map.rot.backward, mask = T)
    if (verbose) cat(paste0("Finished alignment for sample", i, " \n\n"))
  }

  object@tools$processed <- processed.images
  object@tools$processed.masks <- processed.masks
  object[[c("warped_x", "warped_y")]] <- warped_coords
  return(object)
}

#' Function to perform system sleep if user having memory issues during image loading
#'
#' @param x time to sleep
#' @keywords internal


sleepy <- function(x)
{
  p1 <- proc.time()
  Sys.sleep(x)
  proc.time() - p1 # The cpu usage should be negligible
}


