% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/batchGenerateCI2IFC.R
\name{batchGenerateCI2IFC}
\alias{batchGenerateCI2IFC}
\title{Generates multiple 2IFC classification images by participant or condition}
\usage{
batchGenerateCI2IFC(
  data,
  by,
  stimuli,
  responses,
  baseimage,
  rdata,
  save_as_png = TRUE,
  targetpath = "./cis",
  antiCI = FALSE,
  scaling = "autoscale",
  constant = 0.1,
  label = ""
)
}
\arguments{
\item{data}{Data frame}

\item{by}{String specifying column name that specifies the smallest unit (participant, condition) to subset the data on and calculate CIs for.}

\item{stimuli}{String specifying column name in data frame that contains the stimulus numbers of the presented stimuli.}

\item{responses}{String specifying column name in data frame that contains the responses coded 1 for original stimulus selected and -1 for inverted stimulus selected.}

\item{baseimage}{String specifying which base image was used. Not the file name, but the key used in the list of base images at time of generating the stimuli.}

\item{rdata}{String pointing to .RData file that was created when stimuli were generated. This file contains the contrast parameters of all generated stimuli.}

\item{save_as_png}{Boolean stating whether to additionally save the CI as PNG image.}

\item{targetpath}{Optional string specifying path to save PNGGs to (default: ./cis).}

\item{antiCI}{Optional boolean specifying whether antiCI instead of CI should be computed.}

\item{scaling}{Optional string specifying scaling method: \code{none}, \code{constant}, \code{matched}, \code{independent}, or \code{autoscale} (default).}

\item{constant}{Optional number specifying the value used as constant scaling factor for the noise (only works for \code{scaling='constant'}).}

\item{label}{Optional string to insert in file names of PNGGs to make them easier to identify.}
}
\value{
List of classification image data structures (which are themselves lists of pixel matrix of classification noise only, scaled classification noise only, base image only and combined).
}
\description{
Generate classification image for 2 images forced choice reverse correlation task.
}
\details{
This funcions saves the classification images by participant or condition as PNG to a folder and returns the CIs.
}
