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
#' @param reference.index Set reference image used for image alignment. This should be the same the reference.index
#' used for \code{\link{AlignImages}} or \code{\link{ManualAlignImages}}.
#' @param verbose Print messages
#'
#' @importFrom magick image_read image_scale image_info
#' @importFrom imager as.cimg magick2cimg imsub
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
  reference.index = 1,
  verbose = FALSE
) {

  # Check if masked images are available
  if (!"masked" %in% rasterlists(object)) stop(paste0("Masked images are not present in Seurat object"), call. = FALSE)
  st.object <- GetStaffli(object)

  if (all(!c("masked", "processed") %in% rasterlists(st.object))) stop("No image transformations have been executed. Aborting ...", call. = FALSE)

  # Load images
  if (xdim > 2000) stop("Maximum size allowed is 2000px width. Set a lower number for xdim.", call. = FALSE)
  if (class(xdim) != "numeric") stop("xdim must be a numeric value.", call. = FALSE)
  old.xdim <- st.object@xdim
  st.object@xdim <- xdim
  high.res.images <- setNames(lapply(seq_along(names(st.object)), function(i) {
    x <- st.object@imgs[i]
    im <- image_read(path = x)
    
    # Check if image has been cropped
    ds <- st.object@dims[[i]]
    if (ds$cropped) {
      crop.geom <- paste0(ds$max_x - ds$min_x, "x", ds$max_y - ds$min_y, "+", ds$min_x, "+", ds$min_y)
      im <- image_crop(im, geometry = crop.geom)
    }
    
    im <- im %>% image_scale(paste0(xdim))
    if (verbose) cat(paste0("Loaded section ", i, " image at ", paste0(image_info(im)$width), "x", image_info(im)$height), "resolution ... \n")
    
    as.raster(im)
  }), nm = names(st.object))

  # Blow up masks
  masked.masks <- st.object["masked.masks"]
  masked.masks <- setNames(lapply(seq_along(masked.masks), function(n) {
    msk <- masked.masks[[n]]
    msk <- image_read(msk) %>% image_scale(paste0(xdim))
    msk <- msk %>% as.raster()

    # Check if dimensions match
    im <- high.res.images[[n]]
    diff <- nrow(msk) - nrow(im)

    if (diff > 0) {
      for (i in 1:diff) {
        msk <- msk[-nrow(msk), ]
      }
    } else if (diff < 0) {
      for (j in 1:abs(diff)) {
        msk <- rbind(as.matrix(msk), rep("#000000ff", ncol(msk))) %>% as.raster()
      }
    }

    if (verbose) cat(paste0("Scaled section ", n, " image mask to ", paste0(ncol(msk)), "x", nrow(msk)), "resolution ... \n")
    return(msk)
  }), nm = names(st.object))

  # Apply masks to higher resolution images
  masked <- setNames(lapply(seq_along(high.res.images), function(i) {
    im <- high.res.images[[i]] %>% as.cimg()
    msk <- masked.masks[[i]] %>% as.cimg()
    msked.im <- im*msk
    msked.im <- as.raster(msked.im)
    msked.im[msked.im == "#000000"] <- "#FFFFFF"
    if (verbose) cat(paste0("Finished masking section ", i, " image ... \n"))
    return(msked.im)
  }), nm = names(st.object))

  # Transform images if processed images are available
  if ("processed" %in% rasterlists(st.object)) {
    transforms <- st.object@transformations
    processed.images <- masked
    processed.masks <- st.object["processed.masks"]
    mref <- processed.images[[reference.index]]
    for (i in seq_along(transforms)) {
      s <- unique(st.object@meta.data$sample)[i]
      # Obtain scale factors
      dims.raw <- as.numeric(st.object@dims[[i]][2:3])
      dims.scaled <- dim(high.res.images[[i]])
      sf.xy <- dims.raw[1]/dims.scaled[2]

      # Define warp functions
      m <- masked[[i]] %>% as.cimg()
      tr <- transforms[[i]]
      s <- xdim/old.xdim
      tr <- matrix(c(1, 1, 1, 1, 1, 1, s, s, 1), ncol = 3)*tr
      map.affine.backward <- generate.map.affine(tr)
      map.affine.forward <- generate.map.affine(tr, forward = TRUE)

      if (verbose) cat(paste0("Warping image for section ", i, " ... \n"))
      imat <- Warp(m, map.affine.backward)
      select_rows <- min(nrow(imat), nrow(mref)); select_cols <- min(ncol(imat), ncol(mref))
      imat <- imat[1:select_rows, 1:select_cols] %>% as.cimg()

      if (verbose) cat(paste0("Scaling processed mask for section ", i, " image ... \n"))
      imat.msk <- image_read(processed.masks[[i]]) %>% image_scale(paste0(xdim)) %>% magick2cimg()*255

      diff <- ncol(imat.msk) - ncol(imat)

      if (diff > 0) {
        imat.msk <- imsub(imat.msk, y <= ncol(imat))
      } else if (diff < 0) {
        imat.msk <- as.array(imat.msk)
        empty_msk <- array(dim = c(dim(imat.msk)[1], ncol(imat.msk) + abs(diff), dim(imat.msk)[3:4]))
        for (j in 1:dim(imat.msk)[4]) {
          empty_msk[, , 1, j] <- cbind(imat.msk[, , 1, j], matrix(rep(0, nrow(imat.msk)*abs(diff)), ncol = abs(diff)))
        }
        imat.msk <- as.cimg(empty_msk)
      }

      if (verbose) cat(paste0("Cleaning up background ... \n"))
      imat.masked <- imat*threshold(imat.msk, thr = 1)
      imat.masked <- imat.masked + (!threshold(imat.msk, thr = 1))*255
      imat.masked[imat.masked > 255] <- 255

      if (verbose) cat(paste0("Sample ", i, " image processing complete. \n\n"))
      processed.images[[i]] <- imat.masked %>% as.raster()
      processed.masks[[i]] <- imat.msk %>% as.raster()
    }
    st.object@rasterlists$processed <- processed.images
    st.object@rasterlists$processed.masks <- processed.masks
  }

  # Save values
  st.object@rasterlists$raw <- high.res.images
  st.object@rasterlists$masked <- masked
  st.object@rasterlists$masked.masks <- masked.masks

  object@tools$Staffli <- st.object

  return(object)
}

