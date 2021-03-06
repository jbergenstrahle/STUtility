% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/3D_viusalization.R
\name{interpolate_2D_data}
\alias{interpolate_2D_data}
\title{Interpolate value across point pattern}
\usage{
interpolate_2D_data(data, section.input, nx, ny, xy.range)
}
\arguments{
\item{data}{data.frame object with spot coordinates and a 'value' column with numeric data}

\item{section.input}{data.frame object with x,y coordinates for a point pattern to interpolate values over}

\item{nx, ny}{Dimensions of grid}

\item{xy.range}{Range of x,y coordinates to interpolate over}
}
\value{
data.frame object with point pattern and interpolated values
}
\description{
Takes a data.frame with spot coordinates and a 'value' column as input
and interpolates the 'value' across a grid defined by the dimensions nx, ny
and finally assigns the interpolated values to points associated with each grid cell.
}
