#' @include generics.R Staffli.R
NULL

#' Finds edges of a mask and returns a set of x,y coordinates
#' aligned along the edges
#'
#' @method get.edges Staffli
#'
#' @return  A list of data.frames with edge coordinates
#' @examples
#' # Create a new Staffli object, mask, align and plot edges for image 2
#' st.obj <- CreateStaffliObject(imgs, meta.data)
#' edges <- LoadImages(st.obj, verbose = TRUE) %>% MaskImages() %>% get.edges(index = 2)
#' plot(as.raster(edges) %>% as.cimg())

get.edges.Staffli <- function (
  object,
  index = 1,
  verbose = FALSE,
  type = "masked.masks"
) {
  # Check that type is present
  if (!type %in% rasterlists(object)) stop(paste0(type, " invalid or not present in Staffli object ... \n"), call. = FALSE)
  if (verbose) cat(paste0(" Detecting edges of sample ", index, "\n"))

  # Select image by index
  im <- object[type][[index]]

  # Detect gradients
  grad <- imgradient(as.cimg(im))
  grad.sq <- grad %>% map_il(~ .^2)

  # Combine edges, select maximum value and normalize to 0-1 range
  grad.sq <- add(grad.sq)
  grad.sq <- apply(grad.sq, c(1, 2), max)
  grad.sq <- grad.sq/max(grad.sq)

  return(grad.sq)
}

#' Finds edges of a mask and returns a set of x,y coordinates
#' aligned along the edges
#'
#' @method get.edges Seurat
#'
#' @return  A list of data.frames with edge coordinates
#' @examples
#' # Mask, align and plot edges for image 2 from a Seurat object
#' edges <- LoadImages(se, verbose = TRUE) %>% MaskImages() %>% get.edges(index = 2)
#' plot(as.raster(edges) %>% as.cimg())

get.edges.Seurat <- function (
  object,
  index = 1,
  verbose = FALSE,
  type = "masked.masks"
) {
  st.object <- GetStaffli(object)
  edges <- get.edges(st.object, index = index, verbose = verbose, type = type)
  return(edges)
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
#' @importFrom Rvcg vcgCreateKDtree vcgSearchKDtree
#' @importFrom Morpho computeTransform applyTransform

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
#' Takes a list of obtained with \code{\link{FindOptimalTransform}} and
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
#' @param tr results obtained with \code{\link{find.optimal.transform}}
#' @return A transformation function that takes x and y coordinates as input and outputs a
#' list of warped x, y coordinates

generate.map.affine <- function (
  tr, #icps,
  forward = FALSE
) {
  #icps <- find.optimal.transform(set2, set1, xdim, ydim)
  if (forward) {
    map.affine <- function (x, y) {
      p <- cbind(x, y)
      #os <- icps$os
      #xy <- apply.transform(map = solve(tr), p)
      #xy <- t(abs(t(xy) - os))
      xy <- t(solve(tr)%*%t(cbind(p, 1)))
      list(x = xy[, 1], y = xy[, 2])
    }
  } else {
    map.affine <- function (x, y) {
      p <- cbind(x, y)
      #p <- t(abs(t(p) - icps$os))
      #xy <- apply.transform(map = tr, p)
      xy <- t(tr%*%t(cbind(p, 1)))
      list(x = xy[, 1], y = xy[, 2])
    }
  }
  return(map.affine)
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
