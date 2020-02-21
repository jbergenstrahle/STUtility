#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Mask Images
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Masking function t
#'
#' This function takes an image of class \code{cimg} as input and returns
#' an image of class \code{pixset}. The method is typically suitable for images
#' with little background debris. Not suitable for 'Visium' images with fiducials.
#'
#' @param im An image of class \code{cimg}
#'
#' @importFrom imager imgradient enorm isoblur threshold as.cimg as.pixset
#'
#' @export
#'

EdgeMask <- function (
  im
) {
  im <- im %>% imgradient("xy") %>% enorm() %>% isoblur(2)
  im <- im^0.1
  im <- threshold(im)
  im <- im[, , , 1] %>% as.cimg() %>% as.pixset()
}


