#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Switch resolution
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Function used to reload images in different quality and apply image transformations
#'
#' If you want to use histological images in your analyses the recommendation is to
#' load them in low resolution (default 400 pixels width) and apply any transformations
#' that you want, e.g. masking and alignment. Then, whenever you need to increase the
#' resolution of your images you can use this function to reload higher resolution
#' versions of the same images and apply the same transformations.
#'
#' Masking and alignment has been optimized for small images and will be significanlty slower
#' for images of higher resolution. For this reson we recommend users to use this function
#' instead.
#'
#' @param object Seurat object
#' @param xdim Sets the pixel width for scaling, e.g. 400 (maximum allowed width is 2000 pixels)
#' @param verbose Print messages
#'
#' @importFrom magick image_read image_scale image_info
#' @importFrom imager as.cimg magick2cimg
#' @importFrom grDevices as.raster
#'
#' @export
#'
#' @examples
#'
#' # Create Seurat object, load, mask and align images at low resolution
#' se <- InputFromTable(infoTable) %>%
#'    LoadImages() %>%
#'    MaskImages() %>%
#'    AlignImages()
#'
#' # Reload images and apply transformations to higher resolution images, e.g. 2000 pixels width
#' se <- SwitchResolution(se, 2000)
#'

SwitchResolution <- function (
  object,
  xdim,
  verbose = FALSE
) {
  # Run checks
  if (!all(c("masked", "processed") %in% names(object@tools))) stop("No image transformations have been executed. Aborting ...", call. = FALSE)

  # Load images
  if (xdim > 2000) stop("Maximum size allowed is 2000px width. Set a lower number for xdim.", call. = FALSE)
  if (class(xdim) != "numeric") stop("xdim must be a numeric value.", call. = FALSE)
  old.xdim <- object@tools$xdim
  object@tools$xdim <- xdim
  high.res.images <- setNames(lapply(seq_along(object@tools$imgs), function(i) {
    x <- object@tools$imgs[i]
    im <- image_read(path = x) %>% image_scale(paste0(xdim))
    if (verbose) cat(paste0("Loaded image ", i, " at ", paste0(image_info(im)$width), "x", image_info(im)$height), "resolution ... \n")
    as.raster(im)
  }), nm = names(object@tools$raw))

  # Blow up masks
  masked.masks <- object@tools$masked.masks
  masked.masks <- lapply(masked.masks, function(msk) {
    msk <- image_read(msk) %>% image_scale(paste0(xdim))
    if (verbose) cat(paste0("Scaled mask ", i, " to ", paste0(image_info(msk)$width), "x", image_info(msk)$height), "resolution ... \n")
    msk <- msk %>% as.raster()
    return(msk)
  })

  # Apply masks to higher resoltuion images
  masked <- setNames(lapply(seq_along(high.res.images), function(i) {
    im <- high.res.images[[i]] %>% as.cimg()
    msk <- masked.masks[[i]] %>% as.cimg()
    msked.im <- im*msk
    msked.im <- as.raster(msked.im)
    msked.im[msked.im == "#000000"] <- "#FFFFFF"
    if (verbose) cat(paste0("Finished masking image ", i, "... \n"))
    return(msked.im)
  }), nm = names(masked.masks))

  # Transform images if processed images are available
  if ("processed" %in% names(object@tools)) {
    transforms <- object@tools$transformations
    warped_coords <- object[[c("pixel_x", "pixel_y")]]
    processed.images <- masked
    processed.masks <- object@tools$processed.masks
    for (i in seq_along(transforms)) {
      # Obtain scale factors
      dims.raw <- as.numeric(object@tools$dims[[i]][2:3])
      dims.scaled <- dim(high.res.images[[i]])
      sf.xy <- dims.raw/rev(dims.scaled)
      pixel_xy <- subset(object[[]], sample == paste0(i))[, c("pixel_x", "pixel_y")]/sf.xy

      # Define warp functions
      m <- masked[[i]] %>% as.cimg()
      tr <- transforms[[i]]
      s <- xdim/old.xdim
      tr <- matrix(c(1, 1, 1, 1, 1, 1, s, s, 1), ncol = 3)*tr
      map.affine.backward <- generate.map.affine(tr)
      map.affine.forward <- generate.map.affine(tr, forward = TRUE)

      # Warp pixels
      if (verbose) cat(paste0("Warping pixel coordinates for ", i, " ... \n"))
      warped_xy <- sapply(setNames(as.data.frame(do.call(cbind, map.affine.forward(pixel_xy$pixel_x, pixel_xy$pixel_y))), nm = c("warped_x", "warped_y"))*sf.xy, round, digits = 1)
      warped_coords[rownames(pixel_xy), 1:2] <- warped_xy

      if (verbose) cat(paste0("Warping image for ", i, " ... \n"))
      imat <- imwarp(m, map = map.affine.backward, dir = "backward", interpolation = "cubic")
      if (verbose) cat(paste0("Scaling processed image mask for ", i, " ... \n"))
      imat.msk <- image_read(processed.masks[[i]]) %>% image_scale(paste0(xdim)) %>% magick2cimg()*255
      inds <- which(imat.msk != 255)

      if (verbose) cat(paste0(" Cleaning up background ... \n"))
      imrst <- as.raster(imat)
      imat[inds] <- 255
      imrst <- as.raster(imat)
      tab.im <- table(imrst)
      if (length(tab.im) > 2) {
        imrst[imrst == names(which.max(tab.im))] <- "#FFFFFF"
      }

      if (verbose) cat(paste0(" Image ", i, " processing complete. \n\n"))
      processed.images[[i]] <- imrst
      processed.masks[[i]] <- imat.msk
    }
    object@tools$processed <- processed.images
    object@tools$processed.masks <- processed.masks
  }

  # Save values
  object@tools$raw <- high.res.images
  object@tools$masked <- masked
  object@tools$masked.masks <- masked.masks
  object[[c("warped_x", "warped_y")]] <- warped_coords

  return(object)
}

