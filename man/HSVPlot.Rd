% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HSVencoding.R
\name{HSVPlot}
\alias{HSVPlot}
\title{HSV encoded plots}
\usage{
HSVPlot(
  object,
  features,
  ncol = NULL,
  plot.type = "spots",
  split.hsv = FALSE,
  rescale = TRUE,
  indices = NULL,
  spots = NULL,
  min.cutoff = NA,
  max.cutoff = NA,
  slot = "data",
  pt.size = 1,
  pt.alpha = 1,
  pt.border = FALSE,
  add.alpha = FALSE,
  sigma = 2,
  highlight.edges = FALSE,
  cols = NULL,
  dark.theme = TRUE,
  grid.ncol = NULL,
  theme = theme_void(),
  scale.res = 1,
  verbose = FALSE,
  ...
)
}
\arguments{
\item{object}{Seurat object}

\item{features}{\itemize{
    \item An \code{Assay} feature (e.g. a gene name - "MS4A1")
    \item A column name from meta.data (e.g. mitochondrial percentage - "percent.mito")
}}

\item{ncol}{Number of columns to arrange the samples into. This can for example be useful to adjust if you want to visualize the samples
in just in one row or one column.}

\item{plot.type}{Select one of 'spots' or 'smooth' [default: 'spots']}

\item{split.hsv}{Should the HSV colored features be split into separate plots? [default: FALSE]}

\item{rescale}{Rescale each feature column separately from 0 to 1 range. If set to FALSE, all feature columns
will be scaled together from 0 to 1 and preserve the relative differencies}

\item{indices}{Numeric vector specifying sample indices to include in plot. Default is to show all samples.}

\item{spots}{Vector of spots to plot (default is all spots)}

\item{min.cutoff, max.cutoff}{Vector of minimum and maximum cutoff values for each feature,
may specify quantile in the form of 'q##' where '##' is the quantile (eg, 'q1', 'q10')}

\item{slot}{Which slot to pull expression data from?}

\item{pt.size}{Adjust point size for plotting}

\item{pt.alpha}{Adjust opacity of spots.}

\item{pt.border}{Should a border be drawn around the spots? [default: TRUE]}

\item{add.alpha}{Adds opacity to spots scaled by feature values. This will disable the pt.alpha parameter}

\item{sigma}{Smoothing bandwidth; only active if \code{plot.type = 'smooth'}. A single positive number, a numeric vector of length 2, or a function that selects the bandwidth automatically [default: 2].
See \code{\link{density.ppp}} function from the \code{\link{spatstat.core}} package for more details.}

\item{highlight.edges}{Highlights the edges of the tissue. Only active if \code{plot.type = 'smooth'} and if the images have been masked.}

\item{cols}{A vector of colors to use for colorscale, e.g. \code{cols = c("blue", "white", "red")} will
create a gradient color scale going from blue to white to red. This options will deactivate the \code{palette}
option.}

\item{dark.theme}{Use a dark theme for plotting}

\item{grid.ncol}{Number of columns for display when combining plots}

\item{theme}{Add a custom theme to the output ggplot object}

\item{scale.res}{Integer value setting the resolution of the output raster image. E.g. scale.res = 2 will double the
resolution of the output but will also take longer to render. Only active if plot.type is set to 'smooth'.}

\item{verbose}{Print messages}

\item{...}{Extra parameters passed on to \code{\link{STPlot}}}
}
\value{
A ggplot object
}
\description{
Using an HSV encoding of feature values, this functions can be used to color
code expression profiles of multiple features and visualize spatially.
}
\details{
Using RGB encoding, we can show up to 3 features at the same time in the
"red", "green" and "blue" color channels. Whenever two or three features overlap,
the color will be a mix of the three channels, e.g. 50% green and 50% red will give a yellow color.
This strategy is very effective when looking at features values with significant
overlap but is limited to show maximum three features.

If we want to show more than three features in the same plot, this
function provides a strategy to do this as long as the overlap between features
is relatively low. First, a color is assigned to each of N features by cutting
the hue (H) into N values with an even interval. The feature values (e.g. gene expression)
are then scaled to a 0-1 range which is encoded in the Value channel (V).
For each spot, the color with the highest V is selected meaning that only the
feature with the highest value will be shown in the plot. This strategy works well
for features with no or very little overlap but gets cluttered when to many
features are included.

This visualization method should be used only on carefully selected features and you should be
aware that color representation of quantitative data can be very misleading. It should only be
usde to assess qualitative aspects of the data, for example if you wish to know where 5 "non-overlapping"
features are expressed spatially. You should therefore investigate beforehand if the features of interest
overlap or, otherwise the results can become very confusing.
}
\section{scaling of features}{

All features are by default scaled independently to a 0 to 1 range which means that the relative
differencies between the feature expression levels is not preserved. This is because some features
can still be very distinct for a region of interest even though their magnitude of expression is low.
If you want to preserve the relative differencies you can set `rescale = FALSE`.
}

