#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load Images
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
      par(mar = c(0, 0.2, 0, 0.2))
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
#' @param object Seurat object
#' @param iso.blur Sigma value (pixels) for isoblurring of HE images prior to image segmentation
#' @param channels.use Select channel to use for masking (default: 1)
#' @param verbose Print messages
#'
#' @inheritParams slic
#'
#' @importFrom imager magick2cimg medianblur sRGBtoLab as.cimg split_connected add imsplit imappend RGBtoHSV blur_anisotropic HSVtoRGB add threshold isoblur imlist
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
  iso.blur = 2,
  channels.use = 1,
  compactness = 1,
  add.contrast = TRUE,
  verbose = FALSE
) {

  rasters <- list()
  masks <- list()
  #centers <- list()

  for (i in seq_along(object@tools$raw)) {
    imr <- image_read(object@tools$raw[[i]])

    # segmentation tests
    im <- magick2cimg(imr)
    im <- threshold(im)

    rm.channels <- (1:3)[-channels.use]
    for (ind in rm.channels) {
      im[, , , ind] <- TRUE
    }

    #im[, , , 2] <- TRUE; im[, , , 3] <- TRUE
    im <- isoblur(im, iso.blur)

    if (verbose) {
        cat(paste0("Loaded image ", i, "\n"))
        cat(paste0("Running SLIC algorithm \n"))
    }
    out <- slic(im, nS = object@tools$xdim*1.5, compactness)
    if (add.contrast) out <- out^4
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
    #center <- which(seg[, , 1, 1] > 0, arr.ind = T) %>% as.data.frame() %>% summarize(center.x = mean(row), center.y = mean(col))
    #centers[[i]] <- center
  }

  object@tools$masked <- setNames(rasters, nm = names(object@tools$raw))
  object@tools$masked.masks <- setNames(masks, nm = names(object@tools$raw))
  #object@tools$centers <- setNames(centers, nm = names(object@tools$raw))

  return(object)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Warp Images
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# TODO: Fix center

#' Warps images using various transformations
#'
#' @param object Seurat object
#' @param transforms List of arguments passed to warp function
#' @param verbose Print messages
#'
#' @return Seurat object with processed imaged

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


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Align Images (automatic)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Automatic alignment of HE stained tissue images
#'
#' Image alignment or image registration consists in finding a rigid tranformation
#' function that remaps pixels between two images so that two images are aligned.
#' AlignImages allows you to align all images to a reference (any image present in the Seurat object)
#' which simplifies the interpretation of spatial heatmaps and can be useful for creating 3D
#' models. The transformation function is learned using the ICP (Iterative Closest Point) on two point
#' sets which defines the edges of the tissues in the HE images. Note that this alignment works best
#' for tissue sections that are intact, i.e. not cropped or folded. Also, because the method only used
#' the edges of the tissue for alignment, you might end up with strange results if the tissue
#' shape is symmetrical.
#'
#' @param object A Seurat object
#' @param indices Integer: sample indices of images to align with reference
#' @param reference.index Integer: sample index of referece image
#' @param verbose Print messages
#'
#' @importFrom imager as.cimg imwarp
#' @importFrom grDevices as.raster
#'
#' @export

AlignImages <- function (
  object,
  indices = NULL,
  reference.index = NULL,
  verbose = FALSE
) {

  if (!"masked" %in% names(object@tools)) stop(paste0("Masked images are not present in Seurat object"), call. = FALSE)
  if (!any(c("pixel_x", "pixel_y") %in% colnames(object[[]]))) stop(paste0("Pixel coordinates are missing in Seurat object"), call. = FALSE)

  reference.index <- reference.index %||% 1
  if (verbose) cat(paste0("Selecting image ", reference.index, " as reference for alignment. \n"))

  reference.edge <- get.edges(object, index = reference.index)
  indices <- indices %||% (1:length(object@tools$imgs))[-reference.index]
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
  im.ref <- as.cimg(object@tools$raw[[reference.index]])

  # Create empty lists
  transformations <- setNames(lapply(1:length(object@tools$imgs), function(i) {diag(c(1, 1, 1))}), nm = names(object@tools$masked))
  processed.images <- setNames(ifelse(rep("processed" %in% names(object@tools), length(object@tools$imgs)), object@tools$processed, object@tools$masked), nm = names(object@tools$masked))
  processed.masks <- object@tools$masked.masks
  warped_coords <- object[[c("pixel_x", "pixel_y")]]

  for (i in indices) {
    ima <- as.cimg(object@tools$raw[[i]])
    ima.msk <- as.cimg(object@tools$masked.masks[[i]])


    if (verbose) cat(paste0("Processing image ", i, " \n Estimating transformation function ... \n"))
    xdim <- object@tools$xdim
    width <- as.numeric(object@tools$dims[[i]][2]); height <- as.numeric(object@tools$dims[[i]][3])
    ydim <- round(height/(width/xdim))

    # Obtain optimal transform and create map functions
    icps <- find.optimal.transform(xyset.ref, xyset[[paste0(i)]], xdim, ydim)
    tr <- icps$icp$map
    # Collect rotation matrix
    tr <- tr[-3, -3]#; tr[1:2, 3] <- 0
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
    width.raw <- as.numeric(object@tools$dims[[i]][2])
    width.scaled <- dim(object@tools$raw[[i]])[i]
    sf.xy <- width.raw/width.scaled
    pixel_xy <- subset(object[[]], sample == paste0(i))[, c("pixel_x", "pixel_y")]/sf.xy

    # Warp coordinates
    warped_xy <- map.affine.forward(pixel_xy[, 1], pixel_xy[, 2])
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

  object@tools$processed <- processed.images
  object@tools$processed.masks <- processed.masks
  object@tools$transformations <- transformations
  object[[c("warped_x", "warped_y")]] <- warped_coords

  return(object)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Align Images (manual)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Manual alignment of images
#'
#' Creates an interactive shiny application to align images manually
#'
#' @param object Seurat object
#' @param Image type used for alignment
#' @param reference.index Specifies reference sample image for alignment(default: 1)
#' @param edges Uses the tissue edges as points set for alignment
#' @param maxnum Maximum grid number
#' @param verbose Print messages
#'
#' @inheritParams grid.from.seu
#'
#' @importFrom shiny runApp fluidPage fluidRow column sliderInput checkboxInput selectInput actionButton plotOutput reactive renderPlot eventReactive observe stopApp
#' @importFrom shinyjs useShinyjs reset
#'
#' @export

ManualAlignImages <- function (
  object,
  type = NULL,
  reference.index = 1,
  edges = TRUE,
  verbose = FALSE,
  maxnum = 1e3
) {

  # use processed images as input if available
  type <- type %||% {
    if ("processed" %in% names(object@tools)) {
      "processed"
    } else if ("masked" %in% names(object@tools)) {
      "masked"
    } else {
      stop(paste0("No masked images are available in the Seurat object"), call. = FALSE)
    }
  }
  if (verbose) cat(paste0("Using ", type, " images as input for alignment ... \n"))

  # Obtain point sets from each image
  scatters <- grid.from.seu(object, type = type, edges = edges, maxnum = maxnum)
  fixed.scatter <- scatters[[reference.index]]$scatter
  counter <- NULL
  coords.ls <- NULL
  tr.matrices <- ifelse(rep(type %in% c("processed", "prossesed.masks"), length(object@tools$imgs)), se@tools$transformations, lapply(seq_along(object@tools[[type]]), function(i) diag(c(1, 1, 1))))
  #tr.matrices <- lapply(seq_along(object@tools[[type]]), function(i) diag(c(1, 1, 1)))


  ui <- fluidPage(
    useShinyjs(),
    fluidRow(
      column(3,
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
               inputId = "size",
               label = "Change point size",
               value = 0.5, min = 0.1, max = 6, step = 0.1
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
      input$size
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

      return(list(scatter = scatter.t, coords = coords.t, tr = tr))
    })

    output$scatter <- renderPlot({

      coords.ls <<- coords_list()
      scatter.t <- coords.ls[[1]]; coords.t <- coords.ls[[2]]

      d <- round((sqrt(400^2 + 400^2) - 400)/2)

      plot(fixed.scatter[, 1], 400 - fixed.scatter[, 2], xlim = c(-d, 400 + d), ylim = c(-d, 400 + d))
      points(scatter.t[, 1], 400 - scatter.t[, 2], col = "gray")
      points(coords.t[, 1], 400 - coords.t[, 2], col = "red", cex = pt_size())

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

  }

  # Returned transformation matrices
  alignment.matrices <- runApp(list(ui = ui, server = server))
  if (verbose) cat(paste("Finished image alignment. \n\n"))
  processed.ids <- which(unlist(lapply(alignment.matrices, function(tr) {!all(tr == diag(c(1, 1, 1)))})))

  if (length(processed.ids) == 0) stop("None of the samples were processed", call. = FALSE)

  # Create lists for transformation
  transformations <- setNames(ifelse(rep("transformations" %in% names(object@tools) & type != "masked", length(object@tools$imgs)), object@tools$transformations, lapply(1:length(object@tools$imgs), function(i) {diag(c(1, 1, 1))})), nm = names(object@tools$masked))

  type <- type %||% match.arg(several.ok = T, names(object@tools), choices = c("processd", "masked"))[1]
  if (!type %in% names(object@tools)) stop(paste0(type, " images not present in Seurat object"), call. = F)

  im.type <- ifelse(type %in% c("processed.masks", "masked.masks"), gsub(pattern = ".masks", replacement = "", x = type), type)
  processed.images <- object@tools[[im.type]]

  msk.type <- ifelse(type %in% c("processed", "masked"), paste0(type, ".masks"), type)
  processed.masks <- masks <- object@tools[[msk.type]]

  if (type %in% c("processed", "processed.masks")) {
    xy.names <- c("warped_x", "warped_y")
  } else {
    xy.names <- c("pixel_x", "pixel_y")
  }
  warped_coords <- object[[xy.names]]

  for (i in processed.ids) {

    if (verbose) cat(paste0("Loading masked image for sample ", i, " ... \n"))
    m <- object@tools$masked[[i]]

    # Obtain alignment matrix
    tr <- alignment.matrices[[i]]
    transformations[[i]] <- tr%*%transformations[[i]]

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
    if (verbose) cat(paste0("Finished alignment for sample ", i, " \n\n"))
  }

  object@tools$transformations <- transformations
  object@tools$processed <- processed.images
  object@tools$processed.masks <- processed.masks
  object[[c("warped_x", "warped_y")]] <- warped_coords
  return(object)

}
