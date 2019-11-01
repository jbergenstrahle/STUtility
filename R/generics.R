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
#' @param verbose Print messages
#'
#' @importFrom magick image_read
#'
#' @export
#'

LoadImages <- function (
  object,
  image.paths = NULL,
  xdim = 400,
  verbose = FALSE,
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
#' First of all, it can be beneficial to remove 1 or two color channels from the image
#' using \code{channels.use}. Specifying \code{channels.use = 1} will keep only the first channel before running
#' SLIC.
#' The next step is to apply some blurring efect
#'
#'
#' @section Masking options:
#'
#'
#' @section Custom mask function:
#' Masking can sometimes very difficult to accomplish without cutting off parts of the tissue or
#' including unwanted parts of the tissue. If the masking fails you also have the option to
#' write your own function that takes a "cimg" class image as input and returns a masked "pxset".
#'
#' @param object Seurat or Staffli object
#' @param iso.blur Sigma value (pixels) for isoblurring of HE images prior to image segmentation
#' @param channels.use Select channel to use for masking (default: 1)
#' @param compactness Scales the number of super-pixels
#' @param custom.msk.fkn Custom masking function that takes an image of class "cimg" as input and returns a mask
#' of class "pixset" outlining the tissue area.
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
  UseMethod(generic = 'MaskImages', object = object)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Warp Images
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Warps images using various transformations
#'
#' @param object Seurat object
#' @param transforms List of arguments passed to warp function
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
  UseMethod(generic = "ManualAlignImages", object = object)
}


# ----------------------------------
# --- Image processing utilities ---
# ----------------------------------

#' Obtain edges of binaryh mask stores in masked.masks list
#'
#' @param object Seurat object
#' @param index Sample index
#' @param verbose Print messages
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

scaled.imdims <- function (object) {
  UseMethod("scaled.imdims", object)
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
