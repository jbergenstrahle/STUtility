#' Creates a transformation matrix for rotation and translation
#'
#' Creates a transformation matrix for clockwise rotation by 'alpha' degrees
#' followed by a translation with an offset of (h, k). Points are assumed to be
#' centered at (0, 0).
#'
#' @param h Numeric: offset along x axis
#' @param k Numeric: offset along y axis
#'

rigid.transf <- function (
  h = 0,
  k = 0,
  alpha = 0
) {
  tr <- matrix(c(cos(alpha), sin(alpha), 0, -sin(alpha), cos(alpha), 0, h, k, 1), nrow = 3)
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

#' Combines rigid tranformation matrices
#'
#' Combines rigid tranformation matrices in the following order:
#' translation of points to origin (0, 0) -> reflection of points
#' -> rotation by alpha degrees and translation of points to new center
#'
#' @param center.cur (x, y) image pixel coordinates specifying the current center of the tissue (stored in slot "tools" as "centers")
#' @param center.new (x, y) image pixel coordinates specifying the new center (image center)
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
#' tr1 <- combine.tr(center.cur = apply(edges.px[, 1:2], 2, mean), center.new = c(1200, 1200), alpha = 90)
#' tr2 <- combine.tr(center.cur = apply(edges.px[, 1:2], 2, mean), center.new = c(500, 1200), mirror.x = T, alpha = 30)
#' tr3 <- combine.tr(center.cur = apply(edges.px[, 1:2], 2, mean), center.new = c(1200, 500), mirror.y = T, alpha = 270)
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
