#' Obtain edges of binaryh mask stores in masked.masks list
#'
#' @param object Seurat object
#' @param index Sample index
#' @param verbose Print messages
#'
#' @importFrom imager imgradient add map_il

get.edges <- function (
  object,
  index,
  verbose = FALSE
) {
  if (verbose) cat(paste0(" Detecting edges of sample ", index, "\n"))
  im <- object@tools$masked.masks[[index]]

  grad <- imgradient(as.cimg(im))
  grad.sq <- grad %>% map_il(~ .^2)

  grad.sq <- imager::add(grad.sq)
  grad.sq <- apply(grad.sq, c(1, 2), max)
  return(grad.sq/max(grad.sq))
}



#' Match two point sets using iterative closest point search
#'
#' @param set1 Point set from image to be aligned with reference
#' @param set2 Point set from reference image
#' @param iterations Number of iterations
#' @param mindist Minimum distance that definces valid points
#' @param type Select type of transform to be applied
#' @param threads number of trheads to use
#' @return list with transformed x, y coordinates and list of transformation matrices
#'
#' @param importFrom Rvcg vcgCreateKDtree vcgSearchKDtree
#' @param importFrom Morpho computeTransform applyTransform

icpmat <- function (
  set1,
  set2,
  iterations,
  mindist = 1e+15,
  type = c("rigid", "similarity", "affine"),
  threads = 1
) {
  set1 <- cbind(set1, 0)
  set2 <- cbind(set2, 0)
  xtmp <- set1
  yKD <- vcgCreateKDtree(set2)
  transformations <- list()

  for (i in 1:iterations) {
    clost <- vcgSearchKDtree(yKD, xtmp, 1, threads = threads)
    good <- which(clost$distance < mindist)
    trafo <- computeTransform(set2[clost$index[good], ], xtmp[good, ], type = type, weights = NULL, centerweight = FALSE)
    xtmp <- applyTransform(xtmp[, ], trafo)
    transformations[[i]] <- trafo
  }
  xtmp <- xtmp[, 1:2]
  return(list(xy = xtmp, map =  Reduce(`%*%`, transformations)))
}



#' Finds optimal transform based on RMSE
#'
#' Tests different types of reflection settings and return the optimal solution
#' based on RMSE between the transformed points and thre reference set
#'
#' @param set1 Point set from image to be aligned with reference
#' @param set2 Point set from reference image
#' @param xdim Width of image
#' @param ydim Height of image
#' @return list with the list of tranformation matrices, reflection coordinates and rmse score
#' for the optimal transformation
#'
#' @importFrom Rvcg vcgKDtree

find.optimal.transform <- function (
  set1,
  set2,
  xdim,
  ydim
) {
  os <- matrix(c(c(0, 1, 0, 1)*xdim, c(0, 0, 1, 1)*ydim), ncol = 2)
  trf <- lapply(1:nrow(os), function(i) {
    p <- set1
    px_dims <- os[i, ]
    p <- t(abs(t(p) - px_dims))
    icpr <- icpmat(p, set2, iterations = 10)
    RMSE <- sqrt(sum(vcgKDtree(set2, icpr$xy, k = 1)$distance^2))
    return(list(icp = icpr, os = px_dims, rmse = RMSE))
  })

  rmses <- unlist(lapply(trf, function(x) x$rmse))
  return(trf[[which.min(rmses)]])
}

#' Apply rigid transformation to a set of points
#'
#' Takes a list of obtained with \link{\code{FindOptimalTransform}} and
#' a matrix of x, y coordinates and returns the transformed x, y coordinates
#'
#' @param icp List containing transformation matrices
#' @param set Matrix of x, y coordinates to be transformed
#' @return Matrix of transformed x, y coordinates

apply.transform <- function (
  map,
  set
) {
  set.new <- cbind(as.matrix(set), 0)
  set.new <- applyTransform(x = set.new, map)
  return(set.new[, 1:2])
}


#' Generate a map function
#'
#' Given two sets of points (2D) this function will create a rigid transformation
#' function that minimized the rmse of the two points sets. Because we use the backward
#' transformation function for image warping, we need to compute the forward transformation
#' function for pixel coordinate warping. If set1 is (x, y) and set2 is the target (x', y')
#' point set, we want to find the backward transformation function M(x', y') -> (x, y).
#' The forward transformation will be the inverse of M, i.e. M^-1(x, y) -> (x', y').
#' The iterative closest point algorithm used to find the best rigid transformation function
#' may not find the best solution if the image is reflected. For this reason we calculate
#' the transofmration function for all different combinations of reflections and select the
#' function with the lowest rmse between the aligned set and reference set.
#'
#' @param set1 Matrix of x, y coordinates for points in image to be aligned
#' @param set2 Matrix of x, y coordinates for points in reference
#' @param forward Logical specifying if the forward map function should be returned
#' @param xdim,ydim Integer values specifying the image dimensions
#' @return A transformation function that takes x and y coordinates as input and outputs a
#' list of warped x, y coordinates

generate.map.affine <- function (
  set1,
  set2,
  forward = FALSE,
  xdim,
  ydim
) {
  icps <- find.optimal.transform(set2, set1, xdim, ydim)
  if (forward) {
    map.affine <- function (x, y) {
      p <- cbind(x, y)
      os <- icps$os
      xy <- apply.transform(map = solve(icps$icp$map), p)
      xy <- t(abs(t(xy) - os))
      list(x = xy[, 1], y = xy[, 2])
    }
  } else {
    map.affine <- function (x, y) {
      p <- cbind(x, y)
      p <- t(abs(t(p) - icps$os))
      xy <- apply.transform(map = icps$icp$map, p)
      list(x = xy[, 1], y = xy[, 2])
    }
  }
  return(map.affine)
}


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
  #plot(as.raster(reference.edge))
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
  warp.functions.forward <- lapply(1:length(object@tools$masked), function(i) NULL)
  warp.functions.backward <- lapply(1:length(object@tools$masked), function(i) NULL)
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
    map.affine.backward <- generate.map.affine(xyset[[paste0(i)]], xyset.ref, xdim = xdim, ydim = ydim)
    map.affine.forward <- generate.map.affine(xyset[[paste0(i)]], xyset.ref, xdim = xdim, ydim = ydim, forward = T)
    warp.functions.backward[[i]] <- map.affine.backward
    warp.functions.forward[[i]] <- map.affine.forward

    # Warp images
    if (verbose) cat(paste0(" Applying rigid transformation ... \n"))
    imat = imwarp(ima, map = map.affine.backward, dir = "backward", interpolation = "cubic")
    imat.msk = imwarp(ima.msk, map = map.affine.backward, dir = "backward", interpolation = "linear")
    inds <- which(imat.msk != 255)

    # Obtain scale factors
    dims.raw <- as.numeric(object@tools$dims[[i]][2:3])
    dims.scaled <- dim(object@tools$raw[[i]])
    sf.xy <- dims.raw/dims.scaled
    pixel_xy <- subset(object[[]], sample == paste0(i))[, c("pixel_x", "pixel_y")]/sf.xy

    # Warp coordinates
    warped_xy <- map.affine.forward(pixel_xy[, 1], pixel_xy[, 2])
    warped_coords[rownames(pixel_xy), ] <- sapply(setNames(as.data.frame(t(t(do.call(cbind, warped_xy))*sf.xy)), nm = c("x", "y")), round, 2)

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
  object@tools$warp.functions.forward <- warp.functions.forward
  object@tools$warp.functions.backward <- warp.functions.backward
  object[[c("warped_x", "warped_y")]] <- warped_coords

  return(object)
}
