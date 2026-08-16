#' Generates multiple 2IFC classification images by participant or condition
#'
#' Generate classification image for 2 images forced choice reverse correlation task.
#'
#' This function saves the classification images by participant or condition as PNG to a folder and returns the CIs.
#'
#' @export
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @param data Data frame
#' @param by String specifying column name that specifies the smallest unit (participant, condition) to subset the data on and calculate CIs for.
#' @param stimuli String specifying column name in data frame that contains the stimulus numbers of the presented stimuli.
#' @param responses String specifying column name in data frame that contains the responses coded 1 for original stimulus selected and -1 for inverted stimulus selected.
#' @param baseimage String specifying which base image was used. Not the file name, but the key used in the list of base images at time of generating the stimuli.
#' @param rdata String pointing to .RData file that was created when stimuli were generated. This file contains the contrast parameters of all generated stimuli.
#' @param save_as_png Boolean stating whether to additionally save the CI as PNG image.
#' @param targetpath String specifying the directory to save PNGs to. Required when \code{save_as_png = TRUE}; there is no default path. It is created if it does not exist. Use \code{tempdir()} if you only want to try the function out.
#' @param label Optional string to insert in file names of PNGs to make them easier to identify.
#' @param antiCI Optional boolean specifying whether antiCI instead of CI should be computed.
#' @param scaling Optional string specifying scaling method: \code{none}, \code{constant}, \code{matched}, \code{independent}, or \code{autoscale} (default).
#' @param constant Optional number specifying the value used as constant scaling factor for the noise (only works for \code{scaling='constant'}).
#' @return List of classification image data structures (which are themselves lists of pixel matrix of classification noise only, scaled classification noise only, base image only and combined).
#' @examples
#' # a synthetic square grayscale image stands in for a real base face photo
#' base_face <- tempfile(fileext = ".png")
#' png::writePNG(matrix(runif(32 * 32), 32, 32), base_face)
#'
#' stimulus_path <- tempdir()
#' generateStimuli2IFC(
#'   base_face_files = list(face = base_face),
#'   n_trials = 6,
#'   img_size = 32,
#'   stimulus_path = stimulus_path,
#'   seed = 1,
#'   ncores = 1,
#'   nscales = 1,
#'   save_as_png = FALSE
#' )
#' rdata_file <- list.files(stimulus_path, pattern = "\\.Rdata$", full.names = TRUE)[1]
#'
#' # two "participants", three trials each
#' data <- data.frame(
#'   participant = rep(c("p1", "p2"), each = 3),
#'   stimulus = 1:6,
#'   response = sample(c(1, -1), 6, replace = TRUE)
#' )
#'
#' cis <- suppressWarnings(batchGenerateCI2IFC(
#'   data = data, by = "participant", stimuli = "stimulus", responses = "response",
#'   baseimage = "face", rdata = rdata_file, save_as_png = FALSE
#' ))
batchGenerateCI2IFC <- function(data, by, stimuli, responses, baseimage, rdata, save_as_png = TRUE, targetpath, antiCI = FALSE, scaling = 'autoscale', constant = 0.1, label = '') {

  # targetpath is required, not defaulted: a default path writes to the user's
  # filespace uninvited, which CRAN policy does not allow.
  if (save_as_png && missing(targetpath)) {
    stop(paste0('save_as_png is TRUE but no targetpath was given. Supply ',
                'targetpath = <a directory> to say where the PNGs should go, ',
                'or set save_as_png = FALSE to compute the classification ',
                'images without writing them. Use tempdir() if you only want ',
                'to try the function out.'))
  }

  if (scaling == 'autoscale') {
    doAutoscale <- TRUE
    scaling <- 'none'
  } else {
    doAutoscale <- FALSE
  }

  cis <- list()

  # Remove by = NA's from data
  data <- data[!is.na(data[, by]), ]

  by.levels <- unique(data[, by])
  # dplyr::progress_estimated() is deprecated; use the base R progress bar
  pb <- txtProgressBar(min = 0, max = length(by.levels), style = 3)
  pb_i <- 0

  for (unit in by.levels) {

    # Update progress bar
    pb_i <- pb_i + 1
    setTxtProgressBar(pb, pb_i)

    # Get subset of data
    unitdata <- data[data[, by] == unit, ]

    # Specify filename for CI PNG
    if (label == '') {
      filename <- paste0(baseimage, '_', by, '_', unitdata[1, by])
    } else {
      filename <- paste0(baseimage, '_', label, '_', by, '_', unitdata[1, by])
    }

    # Compute CI with appropriate settings for this subset (Optimize later so rdata file is loaded only once)
    cis[[filename]] <- generateCI2IFC(
      stimuli = unitdata[, stimuli],
      responses = unitdata[, responses],
      baseimage = baseimage,
      rdata = rdata,
      save_as_png = save_as_png,
      filename = filename,
      targetpath = targetpath,
      antiCI = antiCI,
      scaling = scaling,
      constant = constant
    )
  }

  if (doAutoscale) {
    cis <- autoscale(cis, save_as_pngs = save_as_png, targetpath = targetpath)
  }

  close(pb)
  return(cis)

}
