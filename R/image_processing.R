#' @include generics.R Staffli.R image_processing_utilities.R
NULL

#' @rdname LoadImages
#' @method LoadImages Staffli
#'
#' @export
#' @return  A Staffli object
#' @examples
#' # Create a new Staffli object and plot images
#' st.obj <- CreateStaffliObject(imgs, meta.data)
#' st.obj <- LoadImages(st.obj, verbose = TRUE)
#' plot(st.obj)
#'

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

  # Add dummy pixel_x, pixel_y columns if not provided
  if (!all(c("pixel_x", "pixel_y") %in% colnames(object[[]]))) {
    object[[, c("pixel_x", "pixel_y")]] <- object[[, c("x", "y")]]
    convert.pixel.coords <- TRUE
  } else {
    convert.pixel.coords <- FALSE
  }

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

    imgs <- c(imgs, list(as.raster(im)))
    if (time.resolve) {
      gc()
      sleepy(5)
    }

    # Convert pixel coords if specified
    if (convert.pixel.coords) {
      if (verbose) cat("converting coords")
      limits <- object@limits[[i]]
      im.limits <- dims[[i]][2:3] %>% as.numeric()
      xy <- object[[object[[, "sample", drop = T]] == paste0(i), c(c("pixel_x", "pixel_y"))]]
      spot_intx <- im.limits[1]/(limits[1] - 1)
      spot_inty <- im.limits[2]/(limits[2] - 1)
      xy <- t(t(xy - 1)*c(spot_intx, spot_inty))
      object[[object[[, "sample", drop = T]] == paste0(i), c(c("pixel_x", "pixel_y"))]] <- xy
    }
  }

  # Compute minimum spot distance
  pixels.per.um <- c()
  for (i in seq_along(object@imgs)) {
    xy <- subset(object[[]], sample == i)[, c("pixel_x", "pixel_y")]
    d <- dist(xy) %>% as.matrix()
    diag(d) <- Inf
    min.distance <- apply(d, 2, min) %>% median() %>% round(digits = 1)
    if (object@platforms[i] == "1k") {
      min.spot.distance <- 200
    } else if (object@platforms[i] == "2k") {
      min.spot.distance <- 141
    } else if (object@platforms[i] == "Visium") {
      min.spot.distance <- 100
    }
   pixels.per.um[i] <- min.distance/min.spot.distance
  }
  names(pixels.per.um) <- object@samplenames

  object@dims <- setNames(dims, nm = object@samplenames)
  object@rasterlists$raw <- setNames(imgs, nm = object@samplenames)
  object@xdim <- xdim
  object@pixels.per.um <- pixels.per.um

  return(object)
}

#' @rdname LoadImages
#' @method LoadImages Seurat
#'
#' @export
#' @return A Seurat object
#' @examples
#'
#' # Load images into a Seurat object and plot images
#' se <- LoadImages(se, verbose = TRUE)
#' ImagePlot(se)

LoadImages.Seurat <- function (
  object,
  image.paths = NULL,
  xdim = 400,
  verbose = FALSE,
  time.resolve = TRUE
) {
  if (!"Staffli" %in% names(object@tools)) stop("Staffli not present in Seurat object ... \n", call. = FALSE)
  st.object <- object@tools$Staffli
  object@tools$Staffli <- LoadImages(object = st.object, image.paths, xdim, verbose, time.resolve)
  return(object)
}



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Plot images
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Function used to plot HE images obtained with \code{\link{LoadImages}}
#'
#' @param object Seurat object
#' @param indices Image indices to select
#' @param type Specify which image to display [options: "raw", "masked" or "processed"]
#' @param method Specify display method (raster or viewer).
#' @param ncols Number of columns in output grid of images
#' @param annotate Will put a unique id in the top left corner
#' @param darken Switches the background to black
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
  annotate = TRUE,
  darken = FALSE
) {
  # obtain Staffli object
  if (!"Staffli" %in% names(object@tools)) stop("Staffli not present in Seurat object ... \n", call. = FALSE)
  st.object <- object@tools$Staffli

  # Check if images are available
  if (length(rasterlists(st.object)) == 0) stop(paste0("No images are available in this Seurat object. Please Run LoadImages() first."), call. = FALSE)

  # Check that type is OK
  choices <- c("processed", "masked", "raw", "processed.masks", "masked.masks")
  if (!is.null(type)) {
    if (!type %in% names(st.object@rasterlists) | !type %in% choices) stop(paste0("type '", type, "' not present in Seurat object"), call. = FALSE)
  }

  type <- type %||% {
    choices <- c("processed", "masked", "raw", "processed.masks", "masked.masks")
    match.arg(arg = choices, choices = names(st.object@rasterlists), several.ok = TRUE)[1]
  }

  images <- st.object@rasterlists[[type]]

  # Use all images if indices are not specified
  indices <- indices %||% {
    seq_along(images)
  }

  # Check if indices are OK
  if (any(!indices %in% 1:length(images))) stop("Image indices out of bounds: ", paste(setdiff(indices, 1:length(images)), collapse = ", "), call. = F)
  images <- images[indices]

  # Switch color
  if (darken & (type %in% c("masked", "processed"))) {
    images <- lapply(images, function(im) {
      im[im == "#FFFFFF"] <- "#000000"
      return(im)
    })
  }

  # Read images
  images <- lapply(images, image_read)

  # Add sample ID
  if (annotate) {
    images <- setNames(lapply(seq_along(images ), function(i) {image_annotate(images[[i]], text = names(images)[i], color = ifelse(darken, "#FFFFFF", "#000000"), size = round(st.object@xdim/10))}), nm = names(images))
  }
  ncols <- ncols %||% ceiling(sqrt(length(x = images)))
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
    graphics::layout(mat = layout.matrix)

    for (rst in lapply(images, as.raster)) {
      if (darken) {
        par(mar = c(0, 0.2, 0, 0.2), bg = "black")
      } else {
        par(mar = c(0, 0.2, 0, 0.2))
      }
      plot(rst)
    }
    par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)
  } else {
    stop(paste0("Invalid display method: ", method), call. = F)
  }
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Mask Images
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @rdname MaskImages
#' @method MaskImages Staffli
#'
#' @export
#' @return  A Staffli object
#' @examples
#' # Create a new Staffli object, mask and plot images
#' st.obj <- CreateStaffliObject(imgs, meta.data)
#' st.obj <- LoadImages(st.obj, verbose = TRUE) %>% MaskImages()
#' plot(st.obj)
#'

MaskImages.Staffli <- function (
  object,
  thresholding = TRUE,
  iso.blur = 2,
  channels.use = NULL,
  compactness = 1,
  add.contrast = NULL,
  verbose = FALSE,
  custom.msk.fkn = NULL
) {

  # obtain Staffli object
  if (!"raw" %in% names(object@rasterlists)) stop("Raw images not present in Staffli object, pelase run LoadImages() first ... \n", call. = FALSE)

  rasters <- list()
  masks <- list()

  for (i in seq_along(object@rasterlists$raw)) {
    imr <- image_read(object@rasterlists$raw[[i]])

    # segmentation tests
    im <- magick2cimg(imr)

    if (!is.null(custom.msk.fkn)) {
      if (class(custom.msk.fkn) != "function") stop(paste0("custom.msk.fkn is not a function"), call. = FALSE)
      seg <- custom.msk.fkn(im)
      if (dim(seg)[4] == 3) seg <- seg[, , , 1] %>% as.pixset()
      if (!"pixset" %in% class(seg)) stop(paste0("output from custom.msk.fkn is not a 'pxset' object"), call. = FALSE)
    } else {
      if (thresholding) im <- threshold(im)

      # Select channels to use for masking if not specified and depending on platform
      if (object@platforms[i] == "Visium") {
        channels.use <- channels.use %||% 1:3
        add.contrast = add.contrast %||% FALSE
      } else if (object@platforms[i] %in% c("1k", "2k")) {
        channels.use <- channels.use %||% 1
        add.contrast = add.contrast %||% TRUE
      }

      rm.channels <- (1:3)[-channels.use]
      for (ind in rm.channels) {
        im[, , , ind] <- TRUE
      }

      im <- isoblur(im, iso.blur)

      if (verbose) {
        cat(paste0("Loaded image ", i, "\n"))
        cat(paste0("Running SLIC algorithm \n"))
      }
      out <- slic(im, nS = object@xdim*1.5, compactness)
      if (add.contrast) out <- out^4
      d <- sRGBtoLab(out) %>% as.data.frame(wide = "c") %>%
        dplyr::select(-x, -y)

      km <- kmeans(d, 2)
      seg <- as.cimg(km$cluster - 1, dim = c(dim(im)[1:2], 1, 1)) %>% medianblur(20) %>% threshold()
    }

    # Chech that at least one masked region is found
    sp <- split_connected(seg)
    if (length(sp) == 0) sp <- imlist(seg)
    if (length(sp) == 0) stop(paste0("Masking failed for image ", i), call. = FALSE)

    # Colect pixel coordinates for spots in meta.data slot
    dims.raw <- as.numeric(object@dims[[i]][2:3])
    dims.scaled <- dim(object@rasterlists$raw[[i]])
    sf.xy <- dims.raw/rev(dims.scaled)
    pixel_xy <- sapply(subset(object@meta.data, sample == paste0(i))[, c("pixel_x", "pixel_y")]/sf.xy, round)
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

    if (length(keep.sp) == 0) stop(paste0("Masking failed for sample ", i, " with zero spots were found under the tissue. Check that the pixel coordinates match the HE images."))
    masks[[i]] <- as.raster(add(keep.sp))

    if (verbose) {
      cat(paste0("Masking background in image ", i , "\n"))
    }

    rst <- object@rasterlists$raw[[i]]
    rst[(1:length(rst))[-inds_inside_tissue$idx]] <- "#FFFFFF"

    rasters[[i]] <- rst
  }

  object@rasterlists$masked <- setNames(rasters, nm = object@samplenames)
  object@rasterlists$masked.masks <- setNames(masks, nm = object@samplenames)

  return(object)
}

#' @rdname MaskImages
#' @method MaskImages Seurat
#'
#' @export
#' @return A Seurat object
#' @examples
#'
#' # Load images into a Seurat objectm, mask and plot images
#' se <- LoadImages(se, verbose = TRUE) %>% MaskImages()
#' ImagePlot(se)

MaskImages.Seurat <- function (
  object,
  thresholding = TRUE,
  iso.blur = 2,
  channels.use = NULL,
  compactness = 1,
  add.contrast = NULL,
  verbose = FALSE,
  custom.msk.fkn = NULL
) {
  if (!"Staffli" %in% names(object@tools)) stop("Staffli not present in Seurat object ... \n", call. = FALSE)
  object@tools$Staffli <- MaskImages(object@tools$Staffli, thresholding, iso.blur, channels.use, compactness, add.contrast, verbose, custom.msk.fkn)
  return(object)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Warp Images
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @rdname WarpImages
#' @method WarpImages Staffli
#'
#' @export
#' @return  A Staffli object
#' @examples
#' # Create a new Staffli object, mask, warp and plot images
#' st.obj <- CreateStaffliObject(imgs, meta.data)
#' transforms <- list("2" = list("mirror.y" = TRUE))
#' st.obj <- LoadImages(st.obj, verbose = TRUE) %>% MaskImages() %>% WarpImages(transforms)
#' plot(st.obj)

WarpImages.Staffli <- function (
  object,
  transforms,
  verbose = FALSE
) {

  # Check if masked images are available
  if (!"masked" %in% rasterlists(object)) stop(paste0("Masked images are not present in Seurat object"), call. = FALSE)

  # Check if the transform list is OK
  if (!all(names(transforms) %in% samplenames(object))) stop(paste0("transforms does not match the sample labels"))

  if ("processed" %in% rasterlists(object)) {
    processed.images <- object["processed"]
    processed.masks <- object["processed.masks"]
  } else {
    processed.images <- object["masked"]
    processed.masks <- object["masked.masks"]
  }

  masked.images <- object["masked"]; masked.masks <- object["masked.masks"]

  if (all(c("warped_x", "warped_y") %in% colnames(object[[]]))) {
    warped_coords <- object[[, c("warped_x", "warped_y", "sample")]]
    spots <- rownames(subset(warped_coords, sample %in% names(transforms)))
    pxset <- object[[, c("pixel_x", "pixel_y")]]
    warped_coords[spots, 1:2] <- pxset[spots, ]
    warped_coords <- warped_coords[, 1:2]
  } else {
    warped_coords <- object[[, c("pixel_x", "pixel_y")]]
  }

  # Create a list of 3x3 identity matrices
  transformations <- setNames(lapply(1:length(object@samplenames), function(i) {diag(c(1, 1, 1))}), nm = names(object))

  for (i in names(transforms)) {

    if (verbose) cat(paste0("Loading masked image for sample ", i, " ... \n"))
    m <- masked.images[[i]]

    args <- transforms[[i]]
    center.x <- args[["center.x"]] %||% FALSE
    center.y <- args[["center.y"]] %||% FALSE
    angle <- args[["angle"]] %||% 0
    mirror.x <- args[["mirror.x"]] %||% FALSE
    mirror.y <- args[["mirror.y"]] %||% FALSE

    center <- apply(which(m == "#FFFFFF", arr.ind = T), 2, mean)
    center.new <- rev(dim(m)/2)
    center.new <- c(ifelse(center.x, center.new[1], center[1]), ifelse(center.y, center.new[2], center[2]))
    tr <- combine.tr(center, center.new, alpha = angle, mirror.x = mirror.x, mirror.y = mirror.y)
    transformations[[i]] <- tr

    map.rot.backward <- generate.map.rot(tr)
    map.rot.forward <- generate.map.rot(tr, forward = TRUE)

    # Obtain scale factors
    dims.raw <- as.numeric(object@dims[[i]][2:3])
    dims.scaled <- dim(object["raw"][[i]])
    sf.xy <- dims.raw[2]/dims.scaled[1]
    pixel_xy <- subset(object[[]], sample == paste0(i))[, c("pixel_x", "pixel_y")]/sf.xy

    # Warp pixels
    if (verbose) cat(paste0("Warping pixel coordinates for ", i, " ... \n"))
    warped_xy <- sapply(setNames(as.data.frame(do.call(cbind, map.rot.forward(pixel_xy$pixel_x, pixel_xy$pixel_y))), nm = c("warped_x", "warped_y"))*sf.xy, round, digits = 1)
    warped_coords[rownames(pixel_xy), 1:2] <- warped_xy

    if (verbose) cat(paste0("Warping image for ", i, " ... \n"))
    processed.images[[i]] <- Warp(m, map.rot.backward)
    msk <- masked.masks[[i]]
    if (verbose) cat(paste0("Warping image mask for ", i, " ... \n"))
    processed.masks[[i]] <- Warp(msk, map.rot.backward, mask = T)
    if (verbose) cat(paste0("Finished alignment for sample", i, " \n\n"))
  }

  object@transformations <- transformations
  object@rasterlists$processed <- processed.images
  object@rasterlists$processed.masks <- processed.masks
  object[[, c("warped_x", "warped_y")]] <- warped_coords

  return(object)
}

#' @rdname WarpImages
#' @method WarpImages Seurat
#'
#' @export
#' @return A Seurat object
#' @examples
#' # Load, mask, warp and plot images in a Seurat object
#' # Mirror y axis in sample '2' and rotate sample '3' 10 degrees
#' transforms <- list("2" = list("mirror.y" = TRUE), "3" = list("angle" = 10))
#' se <- LoadImages(se, verbose = TRUE) %>% MaskImages() %>% WarpImages(transforms)
#' ImagePlot(se)

WarpImages.Seurat <- function (
  object,
  transforms,
  verbose = FALSE
) {

  # obtain Staffli object
  if (!"Staffli" %in% names(object@tools)) stop("Staffli not present in Seurat object ... \n", call. = FALSE)
  st.object <- WarpImages(GetStaffli(object), transforms, verbose)
  object@tools$Staffli <- st.object
  return(object)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Align Images (automatic)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @rdname AlignImages
#' @method AlignImages Staffli
#'
#' @export
#' @return  A Staffli object
#' @examples
#' # Create a new Staffli object, mask, align and plot images
#' st.obj <- CreateStaffliObject(imgs, meta.data)
#' st.obj <- LoadImages(st.obj, verbose = TRUE) %>% MaskImages() %>% AlignImages()
#' plot(st.obj)

AlignImages.Staffli <- function (
  object,
  indices = NULL,
  reference.index = NULL,
  use.masked = FALSE,
  verbose = FALSE
) {

  # Check if masked images are available
  if (!"masked" %in% rasterlists(object)) stop(paste0("Masked images are not present in Seurat object"), call. = FALSE)

  # Check if pixel coordinates are available
  if (!any(c("pixel_x", "pixel_y") %in% colnames(object[[]]))) stop(paste0("Pixel coordinates are missing in Staffli object"), call. = FALSE)

  reference.index <- reference.index %||% 1
  if (verbose) cat(paste0("Selecting image ", reference.index, " as reference for alignment. \n"))

  # Obtain edges of tissue sections
  reference.edge <- get.edges(object, index = reference.index)
  indices <- indices %||% (1:length(object@samplenames))[-reference.index]
  edge.list <- lapply(indices, function(i) {
    get.edges(object, index = i, verbose = verbose)
  })

  xyset.ref <- which(reference.edge > 0, arr.ind = T)
  colnames(xyset.ref) <- c("x", "y")
  xyset <- setNames(lapply(edge.list, function(edge) {
    xy <- which(edge > 0, arr.ind = T)
    colnames(xy) <- c("x", "y")
    return(xy)
  }), nm = indices)

  # Obtain reference image
  im.ref <- as.cimg(object["raw"][[reference.index]])

  # Create a list of 3x3 identity matrices
  transformations <- setNames(lapply(1:length(object@samplenames), function(i) {diag(c(1, 1, 1))}), nm = names(object@samplenames))

  # Collect processed/masked images
  if ("processed" %in% rasterlists(object) & !use.masked) {
    processed.images <- object["processed"]
    processed.masks <- object["processed.masks"]
  } else {
    processed.images <- object["masked"]
    processed.masks <- object["masked.masks"]
  }
  warped_coords <- object[[, c("pixel_x", "pixel_y")]]

  for (i in indices) {
    ima <- as.cimg(object["raw"][[i]])
    ima.msk <- as.cimg(object["masked.masks"][[i]])

    if (verbose) cat(paste0("Processing image ", i, " \n Estimating transformation function ... \n"))
    c(xdim, ydim) %<-% scaled.imdims(object)[[i]]

    # Obtain optimal transform and create map functions
    icps <- find.optimal.transform(xyset.ref, xyset[[paste0(i)]], xdim, ydim)
    tr <- icps$icp$map

    # Collect rotation matrix and convert to 3x3 matrix (column/row 3 will not be used)
    tr <- tr[-3, -3]
    reflect.x <- icps$os[1] == xdim; reflect.y <- icps$os[2] == ydim
    center.new <- apply(xyset[[paste0(i)]], 2, mean)
    if (reflect.x) {
      center.new[1] <- xdim - center.new[1]
    }
    if (reflect.y) {
      center.new[2] <- ydim - center.new[2]
    }
    tr.refl <- combine.tr(center.cur = apply(xyset[[paste0(i)]], 2, mean), center.new = center.new, alpha = 0, mirror.x = reflect.x, mirror.y = reflect.y)
    tr <- tr%*%tr.refl
    transformations[[i]] <- tr
    map.affine.backward <- generate.map.affine(tr)
    map.affine.forward <- generate.map.affine(tr, forward = T)

    # Warp images
    if (verbose) cat(paste0(" Applying rigid transformation ... \n"))
    imat = imwarp(ima, map = map.affine.backward, dir = "backward", interpolation = "cubic")
    imat.msk = imwarp(ima.msk, map = map.affine.backward, dir = "backward", interpolation = "linear")
    inds <- which(imat.msk != 255)

    # Obtain scale factors
    dims.raw <- object@dims[[i]][, c("width", "height")] %>% as.numeric()
    dims.scaled <- scaled.imdims(object)[[i]]
    sf.xy <- dims.raw[2]/dims.scaled[1]
    pixel_xy <- subset(object[[]], sample == paste0(i))[, c("pixel_x", "pixel_y")]/sf.xy

    # Warp coordinates
    warped_xy <- map.affine.forward(x = pixel_xy[, 1], y = pixel_xy[, 2])
    warped_coords[rownames(pixel_xy), ] <- sapply(setNames(as.data.frame(do.call(cbind, warped_xy)*sf.xy), nm = c("x", "y")), round, 2)

    if (verbose) cat(paste0(" Cleaning up background ... \n"))
    imrst <- as.raster(imat)
    imat[inds] <- 255
    imrst <- as.raster(imat)
    tab.im <- table(imrst)
    if (length(tab.im) > 2) {
      imrst[imrst == names(which.max(tab.im))] <- "#FFFFFF"
    }

    if (verbose) cat(paste0(" Image ", i, " alignment complete. \n\n"))
    processed.images[[i]] <- imrst
    processed.masks[[i]] <- as.raster(imat.msk)
  }

  object@rasterlists$processed <- processed.images
  object@rasterlists$processed.masks <- processed.masks
  object@transformations <- transformations
  object[[, c("warped_x", "warped_y")]] <- warped_coords

  return(object)
}

#' @rdname AlignImages
#' @method AlignImages Seurat
#'
#' @export
#' @return  A Seurat object
#' @examples
#' # Load, mask, align and plot images
#' se <- LoadImages(se, verbose = TRUE) %>% MaskImages() %>% AlignImages()
#' ImagePlot(se)

AlignImages.Seurat <- function (
  object,
  indices = NULL,
  reference.index = NULL,
  use.masked = FALSE,
  verbose = FALSE
) {
  # Check if masked images are available
  if (!"masked" %in% rasterlists(object)) stop(paste0("Masked images are not present in Seurat object"), call. = FALSE)
  object@tools$Staffli <- AlignImages(GetStaffli(object), indices, reference.index, use.masked, verbose)
  return(object)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Align Images (manual)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# TODO: make it possible to use other input than "masked.masks"

#' @rdname ManualAlignImages
#' @method ManualAlignImages Staffli
#'
#' @export
#' @return  A Staffli object
#' @examples
#' # Create a new Staffli object, mask, align and plot images (will start an interactive shiny session)
#' st.obj <- CreateStaffliObject(imgs, meta.data)
#' st.obj <- LoadImages(st.obj, verbose = TRUE) %>% MaskImages() %>% ManualAlignImages()
#' plot(st.obj)

ManualAlignImages.Staffli <- function (
  object,
  type = "masked.masks",
  reference.index = 1,
  edges = TRUE,
  verbose = FALSE,
  maxnum = 1e3,
  custom.edge.detector = NULL
) {

  # Check if ultiple samples are available
  if (length(x = object@samplenames) == 1) stop("Only one sample present in the Staffli object. At least 2 samples are required for alignment ... \n", call. = FALSE)

  # Check if images have been masked
  if (!"masked" %in% rasterlists(object)) warning(paste0("It is recommended to mask the images before alignment ... \n"), call. = FALSE)

  # Check that type is OK
  choices <- c("processed", "masked", "raw", "processed.masks", "masked.masks")
  if (!is.null(type)) {
    if (!type %in% rasterlists(object) | !type %in% choices) stop(paste0("type '", type, "' not present in Staffli object"), call. = FALSE)
  }

  if (verbose) cat(paste0("Using '", type, "' images as input for alignment ... \n"))

  # Obtain point sets from each image
  if (class(custom.edge.detector) == "function") edges <- FALSE
  scatters <- grid.from.staffli(object, type = type, edges = edges, maxnum = maxnum, custom.edge.detector = custom.edge.detector)
  fixed.scatter <- scatters[[reference.index]]$scatter
  counter <- NULL
  coords.ls <- NULL
  transformations <-  ifelse(rep(type %in% c("processed", "prossesed.masks"), length(names(object))), object@transformations, lapply(seq_along(names(object)), function(i) diag(c(1, 1, 1))))
  tr.matrices <- lapply(transformations, function(x) diag(c(1, 1, 1)))
  image.dims <- lapply(object[type], dim)


  ui <- fluidPage(
    useShinyjs(),
    fluidRow(
      column(3,
             shiny::hr(),
             actionButton(inputId = "info", label = "Instructions"),
             shiny::hr(),
             sliderInput(
               inputId = "angle",
               label = "Rotation angle",
               value = 0, min = -120, max = 120, step = 0.1
             ),
             sliderInput(
               inputId = "shift_x",
               label = "Move along x axis",
               value = 0, min = -200, max = 200, step = 1
             ),
             sliderInput(
               inputId = "shift_y",
               label = "Move along y axis",
               value = 0, min = -200, max = 200, step = 1
             ),
             sliderInput(
               inputId = "size_spot",
               label = "Change spot size",
               value = 0.5, min = 0, max = 5, step = 0.1
             ),
             sliderInput(
               inputId = "size_ref",
               label = "Change reference point size",
               value = 0.3, min = 0, max = 5, step = 0.05
             ),
             sliderInput(
               inputId = "size_target",
               label = "Change target sample point size",
               value = 0.3, min = 0, max = 5, step = 0.05
             ),
             checkboxInput(inputId = "flip_x",
                           label = "Mirror along x axis",
                           value = FALSE),
             checkboxInput(inputId = "flip_y",
                           label = "Mirror along y axis",
                           value = FALSE),
             selectInput(inputId = "sample", choices = (1:length(scatters))[-reference.index], label = "Select sample", selected = reference.index),
             actionButton("myBtn", "Return aligned data")
      ),

      column(8, plotOutput("scatter")
      )
    )
  )

  server <- function(input, output) {

    rotation_angle <- reactive({
      input$angle
    })

    translation_xy <- reactive({
      trxy <- c(input$shift_x, input$shift_y)
      return(trxy)
    })

    mirror_xy <- reactive({
      mirrxy <- c(input$flip_x, input$flip_y)
      return(mirrxy)
    })

    pt_size <- reactive({
      input$size_spot
    })

    pt_size_ref <- reactive({
      input$size_ref
    })

    pt_size_target <- reactive({
      input$size_target
    })

    coords_list <- reactive({

      # Obtain point set and spot pixel coordinates
      ls <- scatter.coords()
      scatter.t <- ls[[1]]; coords.t <- ls[[2]]

      # Set transformation parameters
      xt.yt <- translation_xy()
      xy.alpha <- rotation_angle()
      mirrxy <-  mirror_xy()

      # Apply reflections
      center <- apply(scatter.t, 2, mean)
      tr.mirror <- mirror(mirror.x = mirrxy[1], mirror.y = mirrxy[2], center.cur = center)

      # Apply rotation
      tr.rotate <- rotate(angle = -xy.alpha, center.cur = center)

      # Apply translation
      tr.translate <- translate(translate.x = xt.yt[1], translate.y = -xt.yt[2])

      # Combine transformations
      tr <- tr.translate%*%tr.rotate%*%tr.mirror


      # Apply transformations
      scatter.t <- t(tr%*%rbind(t(scatter.t), 1))[, 1:2]
      coords.t <- t(tr%*%rbind(t(coords.t), 1))[, 1:2]

      return(list(scatter = scatter.t, coords = coords.t, tr = tr, xylimits = image.dims[[input$sample]]))
    })

    output$scatter <- renderPlot({

      coords.ls <<- coords_list()
      c(scatter.t, coords.t, tr, xylimit) %<-% coords.ls

      d <- round((sqrt(xylimit[1]^2 + xylimit[2]^2) - xylimit[2])/2)

      plot(fixed.scatter[, 1], fixed.scatter[, 2], xlim = c(-d, xylimit[1] + d), ylim = c(d + xylimit[2], -d), xaxt = 'n', yaxt = 'n', ann = FALSE, cex = pt_size_ref())
      points(scatter.t[, 1], scatter.t[, 2], col = "gray", cex = pt_size_target())
      points(coords.t[, 1], coords.t[, 2], col = "red", cex = pt_size())

    }, height = 800, width = 800)

    scatter.coords <- eventReactive(input$sample, {
      reset("angle"); reset("shift_x"); reset("shift_y"); reset("flip_x"); reset("flip_y")
      if (!is.null(counter)) {
        scatters[[counter]] <<- coords.ls[c(1, 2)]
        if (!is.null(tr.matrices[[counter]])) {
          tr.matrices[[counter]] <<- coords.ls[[3]]%*%tr.matrices[[counter]]
          #cat("Sample:", counter, "\n",  tr.matrices[[counter]][1, ], "\n", tr.matrices[[counter]][2, ], "\n", tr.matrices[[counter]][3, ], "\n\n")
        } else {
          tr.matrices[[counter]] <<- coords.ls[[3]]
        }
      }
      scatter <- scatters[[as.numeric(input$sample)]]$scatter
      coords <- scatters[[as.numeric(input$sample)]]$coords
      counter <<- as.numeric(input$sample)
      return(list(scatter, coords))
    })

    observe({
      if(input$myBtn > 0){
        if (!is.null(counter)) {
          scatters[[counter]] <<- coords.ls[c(1, 2)]
          if (!is.null(tr.matrices[[counter]])) {
            tr.matrices[[counter]] <<- coords.ls[[3]]%*%tr.matrices[[counter]]
            cat("Sample:", counter, "\n",  tr.matrices[[counter]][1, ], "\n", tr.matrices[[counter]][2, ], "\n", tr.matrices[[counter]][3, ], "\n\n")
          } else {
            tr.matrices[[counter]] <<- coords.ls[[3]]
          }
        }
        stopApp(tr.matrices)
      }
    })

    observeEvent(input$info, {
      showModal(modalDialog(
        title = "Instructions",
        HTML("The selected sample is highlighted by its coordinates under the tissue <br>",
             "highlighted in red. only rigid transformations are allowed, meaning <br>",
             "rotation, shifts along x/y-axes and reflections.<br><br>",
             "1. Select sample that you want to align to a reference [default: 2]<br>",
             "2. Adjust transformation parameters to fit the sample image to the reference<br>",
             "3. Repeat 1-4 until all samples are aligned<br>",
             "4. Press the 'return aligned data' button to return results"),
        easyClose = TRUE,
        footer = NULL
      ))
    })
  }

  # Returned transformation matrices
  alignment.matrices <- runApp(list(ui = ui, server = server))
  alignment.matrices <- lapply(alignment.matrices, function(tr) {
    tr <- solve(tr)
    return(tr)
  })
  if (verbose) cat(paste("Finished image alignment. \n\n"))
  processed.ids <- which(unlist(lapply(alignment.matrices, function(tr) {!all(tr == diag(c(1, 1, 1)))})))

  # Raise error if none of the samples were processed
  if (length(processed.ids) == 0) stop("None of the samples were processed", call. = FALSE)

  im.type <- ifelse(type %in% c("processed.masks", "masked.masks"), gsub(pattern = ".masks", replacement = "", x = type), type)
  processed.images <- object[im.type]

  msk.type <- ifelse(type %in% c("processed", "masked"), paste0(type, ".masks"), type)
  processed.masks <- masks <- object[msk.type]

  if (type %in% c("processed", "processed.masks")) {
    xy.names <- c("warped_x", "warped_y")
  } else {
    xy.names <- c("pixel_x", "pixel_y")
  }
  warped_coords <- object[[, xy.names]]

  # Select input image
  if (type %in% c("masked", "masked.masks")) {
    selected.input.image <- "masked"
  } else if (type %in% c("processed", "processed.masks")) {
    selected.input.image <- "processed"
  } else {
    selected.input.image <- "raw"
  }

  for (i in processed.ids) {

    if (verbose) cat(paste0("Loading masked image for sample ", i, " ... \n"))
    m <- object[selected.input.image][[i]]

    # Obtain alignment matrix
    tr <- alignment.matrices[[i]]
    transformations[[i]] <- tr%*%transformations[[i]]

    map.rot.backward <- generate.map.affine(tr)
    map.rot.forward <- generate.map.affine(tr, forward = TRUE)

    # Obtain scale factors
    dims.raw <- as.numeric(object@dims[[i]][2:3])
    dims.scaled <- scaled.imdims(object)[[i]]
    sf.xy <- dims.raw[1]/dims.scaled[2]
    pixel_xy <- subset(object[[]], sample == paste0(i))[, c("pixel_x", "pixel_y")]/sf.xy

    # Warp pixels
    if (verbose) cat(paste0("Warping pixel coordinates for ", i, " ... \n"))
    warped_xy <- sapply(setNames(as.data.frame(do.call(cbind, map.rot.forward(pixel_xy$pixel_x, pixel_xy$pixel_y))), nm = c("warped_x", "warped_y"))*sf.xy, round, digits = 1)
    #warped_xy <- sapply(setNames(as.data.frame(do.call(cbind, map.rot.forward(pixel_xy$pixel_x, pixel_xy$pixel_y))), nm = c("warped_x", "warped_y")), round, digits = 1)
    warped_coords[rownames(pixel_xy), 1:2] <- warped_xy

    if (verbose) cat(paste0("Warping image for ", i, " ... \n"))
    processed.images[[i]] <- Warp(m, map.rot.backward)

    msk <- masks[[i]]
    if (verbose) cat(paste0("Warping image mask for ", i, " ... \n"))
    processed.masks[[i]] <- Warp(msk, map.rot.backward, mask = T)
    if (verbose) cat(paste0("Finished alignment for sample ", i, " \n\n"))
  }

  object@transformations <- transformations
  object@rasterlists$processed <- processed.images
  object@rasterlists$processed.masks <- processed.masks
  object[[, c("warped_x", "warped_y")]] <- warped_coords

  return(object)
}


#' @rdname ManualAlignImages
#' @method ManualAlignImages Seurat
#'
#' @export
#' @return  A Seurat object
#' @examples
#' # Load, mask, align and plot images (will start an interactive shiny session)
#' se <- LoadImages(se, verbose = TRUE) %>% MaskImages() %>% ManualAlignImages()
#' ImagePlot(se)

ManualAlignImages.Seurat <- function (
  object,
  type = "masked.masks",
  reference.index = 1,
  edges = TRUE,
  verbose = FALSE,
  maxnum = 1e3,
  custom.edge.detector = NULL
) {
  # Check if masked images are available
  if (!"masked" %in% rasterlists(object)) warning(paste0("Masked images are not present in Seurat object"), call. = FALSE)
  object@tools$Staffli <- ManualAlignImages(GetStaffli(object), type, reference.index, edges, verbose, maxnum, custom.edge.detector)
  return(object)
}
