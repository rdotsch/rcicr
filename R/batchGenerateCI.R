#' Generates multiple classification images by participant or condition
#'
#' Generate one classification image per participant or condition, for any reverse correlation task.
#'
#' Splits \code{data} by the \code{by} column, calls \code{\link{generateCI}} for each part, and returns the CIs. By default each CI is also saved as a PNG.
#'
#' @export
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @param data Data frame with one row per trial.
#' @param by Name of the column that splits the data into units, such as participants or conditions. One CI is computed per unit.
#' @param stimuli Name of the column holding the stimulus numbers of the presented stimuli.
#' @param responses Name of the column holding the responses: 1 where the original stimulus was chosen, -1 where the inverted one was.
#' @param baseimage String naming the base image: not its file name, but its key in the \code{base_face_files} list passed to \code{\link{generateStimuli2IFC}}.
#' @param rdata Path to the \code{.Rdata} file written when the stimuli were generated. It holds the contrast parameters of every stimulus.
#' @param save_as_png Boolean: also save the CI as a PNG image.
#' @param targetpath Directory to save PNGs to. Required when \code{save_as_png = TRUE}; there is no default. The directory is created if it does not exist; to just try the function out, use \code{tempdir()}.
#' @param label Optional string added to the PNG file names to make them easier to identify.
#' @param antiCI Boolean: compute the anti-CI, the classification image with its sign flipped, instead of the CI.
#' @param scaling Scaling method: \code{none}, \code{constant}, \code{matched}, \code{independent} or \code{autoscale} (default). \code{autoscale} computes the CIs unscaled and then puts them on one scale with \code{\link{autoscale}}.
#' @param constant Scaling constant for the noise. Used only when \code{scaling = 'constant'}.
#' @return Named list with one classification image per unit. Each is itself a list of pixel matrices: the raw noise (\code{ci}), the scaled noise (\code{scaled}), the base image (\code{base}) and the two combined (\code{combined}).
#' @examples
#' # a synthetic square grayscale image stands in for a real base face photo
#' base_face <- tempfile(fileext = ".png")
#' png::writePNG(matrix(runif(32 * 32), 32, 32), base_face)
#'
#' stimulus_path <- tempfile("stimuli")
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
#' cis <- suppressWarnings(batchGenerateCI(
#'   data = data, by = "participant", stimuli = "stimulus", responses = "response",
#'   baseimage = "face", rdata = rdata_file, save_as_png = FALSE
#' ))
batchGenerateCI <- function(data, by, stimuli, responses, baseimage, rdata, save_as_png = TRUE, targetpath, label = '', antiCI = FALSE, scaling = 'autoscale', constant = 0.1) {

  # targetpath is required, not defaulted: a default path writes to the user's
  # filespace uninvited, which CRAN policy does not allow.
  if (save_as_png && missing(targetpath)) {
    stop(paste0('save_as_png is TRUE but no targetpath was given. Supply ',
      'targetpath = <a directory> to say where the PNGs should go, ',
      'or set save_as_png = FALSE to compute the classification ',
      'images without writing them. Use tempdir() if you only want ',
      'to try the function out.'
    ))
  }

  if (scaling == 'autoscale') {
    doAutoscale <- TRUE
    scaling <- 'none'
  } else {
    doAutoscale <- FALSE
  }

  # Match batchGenerateCI2IFC(): rows without a grouping value cannot name
  # an output CI and should not be turned into a spurious NA group.
  # Columns are read with [[ ]], never [, col]: on a tibble the latter stays a
  # one-column tibble, which the loop and the == below treat as one unit,
  # mixing every group into a single CI (#336).
  data <- data[!is.na(data[[by]]), , drop = FALSE]

  # dplyr::progress_estimated() is deprecated; use the base R progress bar
  pb <- txtProgressBar(min = 0, max = length(unique(data[[by]])), style = 3)
  cis <- list()
  pb_i <- 0

  for (unit in unique(data[[by]])) {

    # Update progress bar
    pb_i <- pb_i + 1
    setTxtProgressBar(pb, pb_i)

    # Get subset of data
    unitdata <- data[data[[by]] == unit, , drop = FALSE]

    # Specify filename for CI PNG
    if (label == '') {
      filename <- paste0(baseimage, '_', by, '_', unitdata[[by]][1])
    } else {
      filename <- paste0(baseimage, '_', label, '_', by, '_', unitdata[[by]][1])
    }

    # Compute CI with appropriate settings for this subset (Optimize later so rdata file is loaded only once)
    cis[[filename]] <- generateCI(
      stimuli = unitdata[[stimuli]],
      responses = unitdata[[responses]],
      baseimage = baseimage,
      rdata = rdata,
      save_as_png = save_as_png,
      filename = paste0(filename),
      targetpath = targetpath,
      antiCI = antiCI,
      scaling = scaling,
      scaling_constant = constant,
      participants = NA
    )
  }

  if (doAutoscale) {
    cis <- autoscale(cis, save_as_pngs = save_as_png, targetpath = targetpath)
  }

  close(pb)
  return(cis)

}
