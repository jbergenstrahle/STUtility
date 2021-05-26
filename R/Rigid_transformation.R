#' Match two point sets using iterative closest point search
#'
#' @param set1,set2 Point set from image to be aligned with reference and point set from reference image
#' @param iterations Number of iterations
#' @param mindist Minimum distance that definces valid points
#' @param type Select type of transform to be applied
#' @param threads number of trheads to use
#'
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
#' @param set1,set2 Point set from image to be aligned with reference and point set from reference image
#' @param xdim,ydim Width and height of image
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
#' Takes a list obtained with \code{\link{find.optimal.transform}} and
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


#' Creates a transformation matrix that rotates an object
#' in 2D around it's center of mass
#'
#' @param angle Angle given in degrees used for rotation
#' @param center.cur Coordinates of the current center of mass

rotate <- function (
  angle,
  center.cur
) {
  alpha <- 2*pi*(angle/360)
  #center.cur <- c(center.cur, 0)
  #points(center.cur[1], center.cur[2], col = "red")
  tr <- rigid.transl(-center.cur[1], -center.cur[2])
  tr <- rigid.transf(center.cur[1], center.cur[2], alpha)%*%tr
  return(tr)
}


#' Creates a transformation matrix that translates an object
#' in 2D
#'
#' @param translate.x,translate.y translation of x, y coordinates

translate <- function (
  translate.x,
  translate.y
) {
  tr <- rigid.transl(translate.x, translate.y)
  return(tr)
}


#' Creates a transformation matrix that mirrors an object
#' in 2D along either the x axis or y axis around its
#' center of mass
#'
#' @param mirror.x,mirror.y Logical specifying whether or not an
#' object should be reflected
#' @param center.cur Coordinates of the current center of mass
#'

mirror <- function (
  mirror.x = FALSE,
  mirror.y = FALSE,
  center.cur
) {
  center.cur <- c(center.cur, 0)
  tr <- rigid.transl(-center.cur[1], -center.cur[2])
  tr <- rigid.refl(mirror.x, mirror.y)%*%tr
  tr <- rigid.transl(center.cur[1], center.cur[2])%*%tr
  return(tr)
}


#' Stretch along angle
#' 
#' Creates a transformation matrix that stretches an object
#' along a specific axis
#' 
#' @param r stretching factor
#' @param alpha angle
#' @param center.cur Coordinates of the current center of mass
#' 

stretch <- function(r, alpha, center.cur) {
  center.cur <- c(center.cur, 0)
  tr <- rigid.transl(-center.cur[1], -center.cur[2])
  tr <- rigid.rot(alpha, forward = TRUE)%*%tr
  tr <- rigid.stretch(r)%*%tr
  tr <- rigid.rot(alpha, forward = FALSE)%*%tr
  tr <- rigid.transl(center.cur[1], center.cur[2])%*%tr
  return(tr)
}


#' Creates a transformation matrix for rotation
#' 
#' Creates a transformation matrix for clockwise rotation by 'alpha' degrees
#' 
#' @param alpha rotation angle
#' 

rigid.rot <- function (
  alpha = 0,
  forward = TRUE
) {
  alpha <- 2*pi*(alpha/360)
  tr <- matrix(c(cos(alpha), ifelse(forward, -sin(alpha), sin(alpha)), 0, ifelse(forward, sin(alpha), -sin(alpha)), cos(alpha), 0, 0, 0, 1), nrow = 3)
  return(tr)
}


#' Creates a transformation matrix for rotation and translation
#'
#' Creates a transformation matrix for clockwise rotation by 'alpha' degrees
#' followed by a translation with an offset of (h, k). Points are assumed to be
#' centered at (0, 0).
#'
#' @param h Numeric: offset along x axis
#' @param k Numeric: offset along y axis
#' @param alpha rotation angle
#'

rigid.transf <- function (
  h = 0,
  k = 0,
  alpha = 0
) {
  tr <- matrix(c(cos(alpha), -sin(alpha), 0, sin(alpha), cos(alpha), 0, h, k, 1), nrow = 3)
  return(tr)
}

#' Creates a transformation matrix for translation with an offset of (h, k)
#'
#' @param h Numeric: offset along x axis
#' @param k Numeric: offset along y axis
#'

rigid.transl <- function (
  h = 0,
  k = 0
) {
  tr <-  matrix(c(1, 0, 0, 0, 1, 0, h, k, 1), nrow = 3)
  return(tr)
}

#' Creates a transformation matrix for reflection
#'
#' Creates a transformation matrix for reflection where mirror.x will reflect the
#' points along the x axis and mirror.y will reflect thepoints along the y axis.
#' Points are assumed to be centered at (0, 0)
#'
#' @param mirror.x,mirror.y Logical: mirrors x or y axis if set to TRUE

rigid.refl <- function (
  mirror.x,
  mirror.y
) {
  tr <- diag(c(1, 1, 1))
  if (mirror.x) {
    tr[1, 1] <- - tr[1, 1]
  }
  if (mirror.y) {
    tr[2, 2] <- - tr[2, 2]
  }
  return(tr)
}

#' Creates a transformation matrix for stretching
#' 
#' Creates a transformation matrix for stretching by a factor of r 
#' along the x axis.
#' 
#' @param r stretching factor

rigid.stretch <- function (
  r
) {
  tr <- matrix(c(r, 0, 0, 0, 1, 0, 0, 0, 1), ncol = 3)
}


#' Combines rigid tranformation matrices
#'
#' Combines rigid tranformation matrices in the following order:
#' translation of points to origin (0, 0) -> reflection of points
#' -> rotation by alpha degrees and translation of points to new center
#'
#' @param center.cur (x, y) image pixel coordinates specifying the current center of the tissue (stored in slot "tools" as "centers")
#' @param center.new (x, y) image pixel coordinates specifying the new center (image center)
#' @param alpha Rotation angle
#'
#' @inheritParams rigid.transf
#' @inheritParams rigid.transl
#' @inheritParams rigid.refl
#'
#' @examples
#' library(imager)
#' library(tidyverse)
#' im <- load.image("https://upload.wikimedia.org/wikipedia/commons/thumb/f/fd/Aster_Tataricus.JPG/1024px-Aster_Tataricus.JPG")
#' d <- sRGBtoLab(im) %>% as.data.frame(wide="c")%>%
#'   dplyr::select(-x,-y)
#'
#' km <- kmeans(d, 2)
#'
#' # Run a segmentation to extract flower
#' seg <- as.cimg(abs(km$cluster - 2), dim = c(dim(im)[1:2], 1, 1))
#' plot(seg); highlight(seg == 1)
#'
#' # Detect edges
#' dx <- imgradient(seg, "x")
#' dy <- imgradient(seg, "y")
#' grad.mag <- sqrt(dx^2 + dy^2)
#' plot(grad.mag)
#'
#' # Extract points at edges
#' edges.px <- which(grad.mag > max(grad.mag[, , 1, 1])/2, arr.ind = TRUE)
#' points(edges.px, col = "green", cex = 0.1)
#'
#' # Apply transformations to point set
#' tr1 <- combine.tr(center.cur = apply(edges.px[, 1:2], 2, mean),
#'                   center.new = c(1200, 1200), alpha = 90)
#' tr2 <- combine.tr(center.cur = apply(edges.px[, 1:2], 2, mean),
#'                   center.new = c(500, 1200), mirror.x = T, alpha = 30)
#' tr3 <- combine.tr(center.cur = apply(edges.px[, 1:2], 2, mean),
#'                   center.new = c(1200, 500), mirror.y = T, alpha = 270)
#' plot(edges.px, xlim = c(0, 1700), ylim = c(0, 1700), cex = 0.1)
#' points(t(tr1%*%t(edges.px[, 1:3])), cex = 0.1, col = "red")
#' points(t(tr2%*%t(edges.px[, 1:3])), cex = 0.1, col = "yellow")
#' points(t(tr3%*%t(edges.px[, 1:3])), cex = 0.1, col = "blue")
#'
#' @export

combine.tr <- function (
  center.cur,
  center.new,
  alpha,
  mirror.x = FALSE,
  mirror.y = FALSE
) {

  alpha <- 2*pi*(alpha/360)
  center.cur <- c(center.cur, 0)
  tr <- rigid.transl(-center.cur[1], -center.cur[2])

  # reflect
  tr <- rigid.refl(mirror.x, mirror.y)%*%tr

  # rotate and translate to new center
  tr <- rigid.transf(center.new[1], center.new[2], alpha)%*%tr
}
