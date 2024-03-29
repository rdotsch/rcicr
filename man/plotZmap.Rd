% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotZmap.R
\name{plotZmap}
\alias{plotZmap}
\title{Plots a Z-map}
\usage{
plotZmap(
  zmap,
  bgimage = "",
  sigma,
  threshold = 3,
  mask = NULL,
  decoration = T,
  targetpath = "zmaps",
  filename = "zmap",
  size = 512,
  ...
)
}
\arguments{
\item{zmap}{A matrix containing z-scores that map onto a given base image. zmap and baseimage must have the same dimensions.}

\item{bgimage}{A matrix containing the grayscale image to use as a background. This should be either the base image or the final CI. If not this argument is not given, only the Z-map will be drawn.}

\item{sigma}{The sigma of the smoothing that was applied to the CI to create the Z-map.}

\item{threshold}{Integer specifying the threshold z-score (default: 3). Z-scores below the threshold will not be plotted on the z-map.}

\item{mask}{Optional. A boolean matrix with the same dimensions as zmap. If a cell evaluates to TRUE, the corresponding zmap pixel will be masked. Can also be the filename (as a string) of a black and white PNG image to be converted to a matrix (black = masked).}

\item{decoration}{Optional boolean specifying whether the Z-map should be plotted with margins, text (sigma, threshold) and a scale (default: TRUE).}

\item{targetpath}{String specifying path to save the Z-map PNG to.}

\item{filename}{Optional string to specify a file name for the Z-map PNG.}

\item{size}{Integer specifying the width and height of the PNG image (default: 512).}

\item{...}{Additional arguments to be passed to raster::plot. Only applied when decoration is TRUE.}
}
\value{
Nothing. It writes a Z-map image.
}
\description{
Plots a Z-map given a matrix of z-scores that maps onto a specified base image.
}
\details{
This function takes in a matrix of z-scores (as returned by generateCI) and an Rdata file containing a base image. It returns a Z-map image in PNG format.
Unlisted additional arguments will be passed to raster::plot. For example, a different color palette can be specified using the \code{col} argument. See raster::plot for details.
}
