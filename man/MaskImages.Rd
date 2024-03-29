% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/image_processing.R
\name{MaskImages}
\alias{MaskImages}
\alias{MaskImages.Staffli}
\alias{MaskImages.Seurat}
\title{Image masking}
\usage{
MaskImages(
  object,
  thresholding = FALSE,
  iso.blur = 2,
  channels.use = NULL,
  compactness = 1,
  add.contrast = NULL,
  verbose = FALSE,
  custom.msk.fkn = NULL
)

\method{MaskImages}{Staffli}(
  object,
  thresholding = TRUE,
  iso.blur = 2,
  channels.use = NULL,
  compactness = 1,
  add.contrast = NULL,
  verbose = FALSE,
  custom.msk.fkn = NULL
)

\method{MaskImages}{Seurat}(
  object,
  thresholding = TRUE,
  iso.blur = 2,
  channels.use = NULL,
  compactness = 1,
  add.contrast = NULL,
  verbose = FALSE,
  custom.msk.fkn = NULL
)
}
\arguments{
\item{object}{Seurat or Staffli object o pre-procesing step [default: TRUE]}

\item{thresholding}{Applies thresholding step}

\item{iso.blur}{Sets the level of smoothing in the pre-procesing step [default: 2]}

\item{channels.use}{Select channel to use for masking [default: 1 for '1k' and '2k' platforms and 1:3 for 'Visium' platform]}

\item{compactness}{Scales the number of super-pixels [default: 1]}

\item{add.contrast}{Add contrast to pre-procesing step [default: TRUE platforms and 1:3 for FALSE for 'Visium' platform]}

\item{verbose}{Print messages}

\item{custom.msk.fkn}{Custom masking function that takes an image of class "cimg" as input and returns a mask
of class "pixset" outlining the tissue area.}
}
\value{
A Seurat or Staffli object with masked HE images

A Staffli object

A Seurat object
}
\description{
Masks the background of a tissue section in HE images
}
\details{
The masking method provided uses the SLIC method; Simple Linear Iterative Clustering (for more info, see
\url{https://www.r-bloggers.com/superpixels-in-imager/}). This algorithm generates superpixels by
clustering pixels based on their color similarity and proximity in the image plane.
Before running the SLIC method on the HE images we have noticed that some pre-processing can improve the
masking significantly.
First of all, a thresholding is applied to the image using a variant of Otsu's method (see \code{\link[imager]{threshold}}.
Thereafter, you can select which color channels to use for the masking. it can be beneficial to remove 1 or two color channels from the image
using \code{channels.use}. Specifying \code{channels.use = 1} will keep only the first channel before running SLIC.
The next step is to apply some blurring efect on the image to "smooth" out speckles in the image. This
"smoothing" effect can be adjusted with \code{iso.blur}, where a higher \code{iso.blur} leads to more smoothing.
The compactness will adjust the number of superpixels to compute. If you double this number you will get twice
as many superpixels. This can sometimes be helpful to get a finer sensitivity of the masking.
}
\section{Custom mask function}{

The `custom.msk.function` gives you the option to use your own masking function specifically designed for your
images. The custom function has to take a "cimg" class image as input and return an object of class "pxset".
}

\section{Troubleshooting}{

Masking HE images is a non trivial problem as the tissue morphology comes in many different shapes and
colors. The default masking algorithm can fail or perform poorly for a number of different reasons and
below is a couple of examples of common problems.
\itemize{
   \item{
     Bubbles and dirt - If you have air bubbles or other no tissue residue in you images, there is a risk that
     this will be picked up as tissue by the `MaskImages` function. The function uses the spotfile coordinates
     to define what region is outside or inside tissue, but if these are not provided, any part of the image
     picked up by the masking algorithm will be interpreted as relevant. If there is no way around this, you
     can mask the images manually using an image editing software before loading them into your Seurat object.
   }
   \item{
     Tissue gets masked instead of background - If this happens, it probably means that your pixel coordinates provided
     in the spotfiles do not match the HE images provided. For example, if you've run the ST spot detector on a set of
     HE images, you have to provide these exact HE images when running `InputFromTable` or otherwise the pixel coordinates
     will be incorrect.
   }
   \item{
     Parts of tissue gets masked - The masking algorithm relies on differencies in color intensity to segment out the tissue
     area. A common problem is that parts of your tissue is more similar to the background than the rest of the tissue, whcih
     typically happens for cell sparse tissue regions such as adipose tisse or connective tissue. If this happens you can try
     out different settings for the `thresholding`, `ìso.blur`, `channels.use`, `compactness` or `add.contrast` options or you might have to
     write your own masking function and pass it using the `custom.msk.fkn` option.
   }
}
}

\examples{
\dontrun{
# Create a new Staffli object, mask and plot images
st.obj <- CreateStaffliObject(imgs, meta.data)
st.obj <- LoadImages(st.obj, verbose = TRUE) \%>\% MaskImages()
plot(st.obj)
}

\dontrun{
# Load images into a Seurat objectm, mask and plot images
se <- LoadImages(se, verbose = TRUE) \%>\% MaskImages()
ImagePlot(se)
}
}
