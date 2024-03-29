% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Rigid_transformation.R
\name{find.optimal.transform}
\alias{find.optimal.transform}
\title{Finds optimal transform based on RMSE}
\usage{
find.optimal.transform(set1, set2, xdim, ydim)
}
\arguments{
\item{set1, set2}{Point set from image to be aligned with reference and point set from reference image}

\item{xdim, ydim}{Width and height of image}
}
\value{
list with the list of tranformation matrices, reflection coordinates and rmse score
for the optimal transformation
}
\description{
Tests different types of reflection settings and return the optimal solution
based on RMSE between the transformed points and thre reference set
}
