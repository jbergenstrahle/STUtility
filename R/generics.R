# ----------------------------------
# -------- Image processing --------
# ----------------------------------

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load Images
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Function used to read HE images in jpeg or png format
#'
#' @param object Seurat or Staffli object
#' @param image.paths Paths to HE images. This is only required if image paths are missing in the Seurat object.
#' @param xdim Sets the pixel width for scaling, e.g. 400 (maximum allowed width is 2000 pixels)
#' @param crop.to.fiducials Set to TRUE if you want to crop out background from the images outside the fiducials.
#' @param time.resolve Activate to stop R from loading raw images into memory
#' @param verbose Print messages
#' @param crop.scale.factors Numeric vector of length 4 providing a scaling factor for each side of the image.
#' The scaling factors define the number of spot widths away from the capture area to include in the cropped
#' image if  `crop.to.fiducials` is set to TRUE. The default values are c(9, 10, 10, 8) which corresponds to
#' left, top, right and bottom.
#' @importFrom magick image_read geometry_area image_crop
#' @importFrom imager magick2cimg
#'
#' @export
#'

LoadImages <- function (
  object,
  image.paths = NULL,
  xdim = 400,
  crop.to.fiducials = FALSE,
  crop.scale.factors = c(9, 10, 10, 8),
  verbose = TRUE,
  time.resolve = TRUE
) {
  UseMethod(generic = 'LoadImages', object = object)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Mask Images
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Image masking
#'
#' Masks the background of a tissue section in HE images
#'
#' @details
#' The masking method provided uses the SLIC method; Simple Linear Iterative Clustering (for more info, see
#' \url{https://www.r-bloggers.com/superpixels-in-imager/}). This algorithm generates superpixels by
#' clustering pixels based on their color similarity and proximity in the image plane.
#' Before running the SLIC method on the HE images we have noticed that some pre-processing can improve the
#' masking significantly.
#' First of all, a thresholding is applied to the image using a variant of Otsu's method (see \code{\link[imager]{threshold}}.
#' Thereafter, you can select which color channels to use for the masking. it can be beneficial to remove 1 or two color channels from the image
#' using \code{channels.use}. Specifying \code{channels.use = 1} will keep only the first channel before running SLIC.
#' The next step is to apply some blurring efect on the image to "smooth" out speckles in the image. This
#' "smoothing" effect can be adjusted with \code{iso.blur}, where a higher \code{iso.blur} leads to more smoothing.
#' The compactness will adjust the number of superpixels to compute. If you double this number you will get twice
#' as many superpixels. This can sometimes be helpful to get a finer sensitivity of the masking.
#'
#' @section Custom mask function:
#' The `custom.msk.function` gives you the option to use your own masking function specifically designed for your
#' images. The custom function has to take a "cimg" class image as input and return an object of class "pxset".
#'
#' @section Troubleshooting:
#' Masking HE images is a non trivial problem as the tissue morphology comes in many different shapes and
#' colors. The default masking algorithm can fail or perform poorly for a number of different reasons and
#' below is a couple of examples of common problems.
#' \itemize{
#'    \item{
#'      Bubbles and dirt - If you have air bubbles or other no tissue residue in you images, there is a risk that
#'      this will be picked up as tissue by the `MaskImages` function. The function uses the spotfile coordinates
#'      to define what region is outside or inside tissue, but if these are not provided, any part of the image
#'      picked up by the masking algorithm will be interpreted as relevant. If there is no way around this, you
#'      can mask the images manually using an image editing software before loading them into your Seurat object.
#'    }
#'    \item{
#'      Tissue gets masked instead of background - If this happens, it probably means that your pixel coordinates provided
#'      in the spotfiles do not match the HE images provided. For example, if you've run the ST spot detector on a set of
#'      HE images, you have to provide these exact HE images when running `InputFromTable` or otherwise the pixel coordinates
#'      will be incorrect.
#'    }
#'    \item{
#'      Parts of tissue gets masked - The masking algorithm relies on differencies in color intensity to segment out the tissue
#'      area. A common problem is that parts of your tissue is more similar to the background than the rest of the tissue, whcih
#'      typically happens for cell sparse tissue regions such as adipose tisse or connective tissue. If this happens you can try
#'      out different settings for the `thresholding`, `ìso.blur`, `channels.use`, `compactness` or `add.contrast` options or you might have to
#'      write your own masking function and pass it using the `custom.msk.fkn` option.
#'    }
#' }
#'
#' @param object Seurat or Staffli object o pre-procesing step [default: TRUE]
#' @param thresholding Applies thresholding step
#' @param iso.blur Sets the level of smoothing in the pre-procesing step [default: 2]
#' @param channels.use Select channel to use for masking [default: 1 for '1k' and '2k' platforms and 1:3 for 'Visium' platform]
#' @param compactness Scales the number of super-pixels [default: 1]
#' @param custom.msk.fkn Custom masking function that takes an image of class "cimg" as input and returns a mask
#' of class "pixset" outlining the tissue area.
#' @param add.contrast Add contrast to pre-procesing step [default: TRUE platforms and 1:3 for FALSE for 'Visium' platform]
#' @param verbose Print messages
#'
#' @importFrom imager magick2cimg medianblur sRGBtoLab as.cimg split_connected add imsplit imappend RGBtoHSV blur_anisotropic HSVtoRGB add threshold isoblur imlist
#' @importFrom magick image_read
#' @importFrom dplyr select summarize
#' @importFrom magrittr %>%
#' @importFrom stats kmeans
#' @importFrom purrr modify_at
#'
#' @return A Seurat or Staffli object with masked HE images
#'
#' @export

MaskImages <- function (
  object,
  thresholding = FALSE,
  iso.blur = 2,
  channels.use = NULL,
  compactness = 1,
  add.contrast = NULL,
  verbose = FALSE,
  custom.msk.fkn = NULL
) {
  UseMethod(generic = 'MaskImages', object = object)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Warp Images
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Warps images using various transformations
#'
#' This function can be used to apply rigid transformations on masked images including;
#' rotations, reflections and translations. See details below for more information.
#'
#' @details
#' The transformations are controlled by specifying verious rigid transformations
#' for each sample using the transforms variable. The transforms variable should be a
#' named list of lists with one element for each sample that you want to apply a transformation
#' on. For example, if you want to rotate sample 1 90 degrees clockwise and reflect sample 2
#' along its x axis you can pass `transforms <- list("1" = list("angle" = 90), "2" = list("mirror.x" = TRUE))`.
#' When multiple transformations are passed, they are applied in the following order;
#' rotation -> reflection -> translation.
#'
#'
#'
#' @param object Seurat object
#' @param transforms List of arguments passed to warp function. The available options are
#' "angle" [numeric], "mirror.x" [logical], "mirror.y" [logical], "shift.x" [numeric]and "shift.y" [numeric].
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
  UseMethod(generic = 'WarpImages', object = object)
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
#' @param use.masked Setting this to TRUE will ignore any "processed" images available in the object and
#' apply the alignment from scratch.
#' @param verbose Print messages
#'
#' @importFrom imager as.cimg imwarp
#' @importFrom grDevices as.raster
#' @importFrom zeallot %<-%
#'
#' @export

AlignImages <- function (
  object,
  indices = NULL,
  reference.index = NULL,
  use.masked = FALSE,
  verbose = FALSE
) {
  UseMethod(generic = "AlignImages", object = object)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Align Images (manual)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Manual alignment of images
#'
#' Creates an interactive shiny application to align images manually
#'
#' @param object Seurat object
#' @param type Image type to use as input for alignment [default: 'masked.masks']
#' @param reference.index Specifies reference sample image for alignment(default: 1)
#' @param edges Uses the tissue edges as points set for alignment
#' @param maxnum Maximum grid number
#' @param custom.edge.detector Custom function used to detect edges in tissue image. If a function is provided, the
#' edges option will be overridden.
#' @param verbose Print messages
#'
#' @importFrom shiny runApp fluidPage fluidRow column sliderInput checkboxInput selectInput actionButton plotOutput reactive renderPlot eventReactive observe stopApp
#' @importFrom shinyjs useShinyjs reset
#' @importFrom zeallot %<-%
#'
#' @export

ManualAlignImages <- function (
  object,
  type = NULL,
  reference.index = 1,
  edges = TRUE,
  verbose = FALSE,
  maxnum = 1e3,
  fix.axes = FALSE,
  custom.edge.detector = NULL
) {
  UseMethod(generic = "ManualAlignImages", object = object)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Crop Images
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Manual cropping of images
#'
#' Takes a predefined list of crop windows and cuts each image
#'
#' @details
#' Let's say that you want to crop out a smaller piece of secttion 1 that should be 500x500 pixels
#' and should be offset by 400 pixels along the x axis and 400 pixels along the y axis. The y axis
#' starts from the top of the image so the offset will essentially determine the top left corner
#' of the output cropped image. To do this, you can set `crop.geometry.list = list("1" = "500x500+400+400")`.
#' If `xdim` is not specified, the function will collect the predefined `xdim` that was set when
#' running `LoadImages`. The `xdim` value determines the maximum allowed width of the cropped image,
#' so if you for example set a large crop window from the original image, e.g. "2000x2000+1000+1000"
#' the output image will still be downscaled to a width of `xdim` to keep the output image in low
#' resolution.
#'
#' @param object Seurat object
#' @param crop.geometry.list List of crop windows. Each element of the list should be named by a section number that
#' the cropping should be applied to, e.g. "1", "2", "3", etc. The crop window is defined by a string speficying the
#' size of the output window and the offset in pixels; "_width_x_height_+_offsetx_+_offsety_". Note that the cropping
#' will be applied to the original HE images, i.e. the images that were loaded when running `LoadImages`. See details
#' below for more information.
#' @param xdim Maximum width of cropped window
#' @param return.spots.vec Returns a list with `object` as the first element and a list of spots used to convert
#' between old and new ids.
#' @param time.resolve Activate to stop R from loading raw images into memory
#' @param verbose Print messages
#'
#' @importFrom magick image_read image_crop image_info image_scale
#' @importFrom zeallot %<-%
#'
#' @export

CropImages <- function (
  object,
  crop.geometry.list,
  xdim = NULL,
  return.spots.vec = FALSE,
  time.resolve = FALSE,
  verbose = FALSE
) {
  UseMethod(generic = "CropImages", object = object)
}


# ----------------------------------
# --- Image processing utilities ---
# ----------------------------------

#' Obtain edges of binary mask stores in masked.masks list
#'
#' @param object Seurat object
#' @param index Sample index
#' @param verbose Print messages
#' @param Input image type to run edge detection on [default: 'masked.masks']
#'
#' @importFrom imager imgradient add map_il

get.edges <- function (
  object,
  index = 1,
  verbose = FALSE,
  type = "masked.masks"
) {
  UseMethod("get.edges", object = object)
}


# ----------------------------------
# --------- Staffli methods --------
# ----------------------------------

#' Extracts a list of image info from a Staffli object
#'
#' @param object A Staffli object
#' @export

iminfo <- function (object) {
  UseMethod("iminfo", object)
}



#' Extracts scaled dimensions of images in Staffli object
#'
#' @param object A Staffli object
#' @export

scaled.imdims <- function (object, ...) {
  UseMethod("scaled.imdims", object, ...)
}


#' Extracts available image types from Seurat or Staffli objects
#'
#' @param object A Seurat or Staffli object
#' @export

rasterlists <- function (object) {
  UseMethod("rasterlists", object)
}


#' Extracts available samplenames from Seurat or Staffli objects
#'
#' @param object A Seurat or Staffli object
#' @export

samplenames <- function (object) {
  UseMethod("samplenames", object)
}


#' Extracts Staffli from Seurat object
#'
#' @param object A Seurat object
#' @export

GetStaffli <- function (object) {
  UseMethod("GetStaffli", object)
}
