#' Generates multiple classification images by participant or condition
#'
#' Generate classification image for any reverse correlation task that displays independently generated alternatives.
#'
#' This function saves the classification images by participant or condition as PNG to a folder and returns the CIs.
#'
#' @export
#' @import dplyr
#' @param data Data frame
#' @param by String specifying column name that specifies the smallest unit (participant, condition) to subset the data on and calculate CIs for.
#' @param stimuli String specifying column name in data frame that contains the stimulus numbers of the presented stimuli.
#' @param responses String specifying column name in data frame that contains the responses coded 1 for original stimulus selected and -1 for inverted stimulus selected.
#' @param baseimage String specifying which base image was used. Not the file name, but the key used in the list of base images at time of generating the stimuli.
#' @param rdata String pointing to .RData file that was created when stimuli were generated. This file contains the contrast parameters of all generated stimuli.
#' @param save_as_png Boolean stating whether to additionally save the CI as PNG image.
#' @param targetpath Optional string specifying path to save PNGs to (default: ./cis).
#' @param label Optional string to insert in file names of PNGs to make them easier to identify.
#' @param antiCI Optional boolean specifying whether antiCI instead of CI should be computed.
#' @param scaling Optional string specifying scaling method: \code{none}, \code{constant},  \code{independent} or \code{autoscale} (default).
#' @param constant Optional number specifying the value used as constant scaling factor for the noise (only works for \code{scaling='constant'}).
#' @return List of classification image data structures (which are themselves lists of pixel matrix of classification noise only, scaled classification noise only, base image only and combined).
batchGenerateCI <- function(data, by, stimuli, responses, baseimage, rdata, save_as_png=TRUE, targetpath='./cis', label='', antiCI=FALSE, scaling='autoscale', constant=0.1) {

  if (scaling == 'autoscale') {
    doAutoscale <- TRUE
    scaling <- 'none'
  } else {
    doAutoscale <- FALSE
  }

  pb <- dplyr::progress_estimated(length(unique(data[,by])))
  cis <- list()

  for (unit in unique(data[,by])) {

    # Update progress bar
    pb$tick()$print()

    # Get subset of data
    unitdata <- data[data[,by] == unit, ]

    # Specify filename for CI PNG
    if (label == '') {
      filename <- paste0(baseimage, '_', by, '_', unitdata[1,by])
    } else {
      filename <- paste0(baseimage, '_', label, '_', by, '_', unitdata[1,by])
    }

    # Compute CI with appropriate settings for this subset (Optimize later so rdata file is loaded only once)
    cis[[filename]] <- generateCI(
      stimuli=unitdata[,stimuli],
      responses=unitdata[,responses],
      baseimage=baseimage,
      rdata=rdata,
      save_as_png=save_as_png,
      filename=paste0(filename),
      targetpath=targetpath,
      antiCI=antiCI,
      scaling=scaling,
      scaling_constant=constant,
      participants=NA)
  }

  if (doAutoscale) {
    cis <- autoscale(cis, save_as_pngs=save_as_png, targetpath=targetpath)
  }

  pb$stop()
  return(cis)

}
