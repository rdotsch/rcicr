#' Computes cumulative trial CIs correlations with final/target CI
#'
#' Computes cumulative trial CIs correlations with final/target CI.
#'
#' Use for instance for plotting curves of trial-final/target CI correlations to estimate how many trials are necessary in your task
#'
#' @export
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom stats cor
#' @param stimuli Vector with stimulus numbers (should be numeric) that were presented in the order of the response vector. Stimulus numbers must match those in file name of the generated stimuli.
#' @param responses Vector specifying the responses in the same order of the stimuli vector, coded 1 for original stimulus selected and -1 for inverted stimulus selected.
#' @param baseimage String specifying which base image was used. Not the file name, but the key used in the list of base images at time of generating the stimuli.
#' @param rdata String pointing to .RData file that was created when stimuli were generated. This file contains the contrast parameters of all generated stimuli.
#' @param targetci List Target CI object generated with rcicr functions to correlate cumulative CIs with.
#' @param step Step size in sequence of trials to compute correlations with.
#' @return Vector containing correlation between cumulative CI and final/target CI.
#' @examples
#' \donttest{
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
#' responses <- sample(c(1, -1), 6, replace = TRUE)
#' correlations <- suppressWarnings(computeCumulativeCICorrelation(
#'   stimuli = 1:6, responses = responses, baseimage = "face", rdata = rdata_file
#' ))
#' }
computeCumulativeCICorrelation <- function(stimuli, responses, baseimage, rdata, targetci=list(), step=1) {

  # Coerce to plain vectors: tibble columns stay one-column tibbles rather than
  # dropping to vectors (see the same handling in generateCI()).
  stimuli <- unlist(stimuli, use.names = FALSE)
  responses <- unlist(responses, use.names = FALSE)

  # Load parameter file (created when generating stimuli)
  load(rdata)

  # Check whether critical variables have been loaded
  if (!exists('s', envir=environment(), inherits=FALSE) & !exists('p', envir=environment(), inherits=FALSE) ) {
    stop('File specified in rdata argument did not contain s or p variable.')
  }

  if (!exists('base_faces', envir=environment(), inherits=FALSE)) {
    stop('File specified in rdata argument did not contain base_faces variable.')
  }

  if (!exists('stimuli_params', envir=environment(), inherits=FALSE)) {
    stop('File specified in rdata argument did not contain stimuli_params variable.')
  }

  # Convert s to p (if rdata file originates from pre-0.3.3)
  if (exists('s', envir=environment(), inherits=FALSE)) {
    p <- list(patches=s$sinusoids, patchIdx=s$sinIdx, noise_type='sinusoid')
    rm(s)
  }


  # Get base image
  base <- base_faces[[baseimage]]
  if (is.null(base)) {
    stop(paste0('File specified in rdata argument did not contain any reference to base image label: ', baseimage, ' (NOTE: file contains references to the following base image label(s): ', paste(names(base_faces), collapse=', '), ')'))
  }


  # Retrieve parameters of actually presented stimuli (this will work with non-consecutive stims as well)
  params <- stimuli_params[[baseimage]][stimuli,]

  # Check whether parameters were found in this .rdata file
  if (length(params) == 0) {
    stop(paste0('No parameters found for base image: ', base))
  }

  # Compute final classification image if necessary
  if (length(targetci) == 0) {
    finalCI <- generateCINoise(params, responses, p)
  } else {
    finalCI <- targetci$ci
  }

  # Compute correlations with final CI with cumulative CI
  # dplyr::progress_estimated() is deprecated; use the base R progress bar
  pb <- txtProgressBar(min = 0, max = length(responses), style = 3)

  correlations <- vector()
  corcounter <- 1
  for (trial in seq(1,length(responses), step)) {
    setTxtProgressBar(pb, trial)

    cumCI <- generateCINoise(params[1:trial,], responses[1:trial], p)
    correlations[corcounter] <- cor(as.vector(cumCI), as.vector(finalCI))
    corcounter <- corcounter + 1
  }
  close(pb)

  # Return correlations
  return(correlations)
}
