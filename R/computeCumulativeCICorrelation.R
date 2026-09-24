#' Computes cumulative trial CIs correlations with final/target CI
#'
#' Correlates the CI built from the first trials with the final or target CI, adding trials one step at a time.
#'
#' Plot the resulting curve to estimate how many trials your task needs.
#'
#' @section Repeated presentations of the same stimulus:
#' This function takes the trials in the order they were presented and does not average repeated
#' presentations of a stimulus. \code{\link{generateCI}} does average them, per unique stimulus,
#' before building its classification image. Averaging here would discard the presentation order
#' that a cumulative curve is about.
#'
#' Without a \code{targetci}, the final CI is built here from the same trials as the curve. When
#' the evaluated trials reach the last one, as they always do at the default \code{step = 1}, the
#' curve's last point compares that CI with itself and is exactly 1. That shows self-consistency,
#' not convergence. A larger \code{step} can stop short, because trials are taken at
#' \code{seq(1, length(responses), step)}: with six responses and \code{step = 2}, the last trial
#' evaluated is the fifth, and the curve ends at whatever that partial CI correlates to (0.97 in
#' one such set, not 1).
#'
#' This assumes the CI compared against varies at all. Responses that cancel exactly, with every
#' presentation of a stimulus answered both ways, average to a CI that is zero everywhere. A
#' correlation with a constant is undefined, so then \strong{every} point on the curve is
#' \code{NA}. Such a curve means the responses carry no net signal, not that the call failed.
#'
#' A \code{targetci} with masked pixels (\code{\link{generateCI}} stores \code{NA} in every pixel
#' a \code{mask} excludes) is correlated over the unmasked pixels only. If the mask covers
#' \emph{every} pixel, no pairs remain and the whole curve is \code{NA}, as above.
#'
#' If every stimulus was shown equally often, the final CI computed here is identical to the one
#' \code{generateCI} returns. If not, the two weight the data differently (each trial equally
#' here, each unique stimulus equally there) and they diverge: on an 8-trial set with counts
#' 4/2/1/1 they correlate at 0.77.
#'
#' So to see how the CI approaches the one you will report, pass that CI as
#' \code{targetci = generateCI(...)} instead of relying on the default, which is always built
#' without a \code{mask}.
#'
#' @export
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom stats cor
#' @param stimuli Numeric vector of stimulus numbers, one per response and in the same order. Each must be a positive whole number no larger than the number of trials saved for the selected base image. Numbers may repeat and need not be consecutive. Factors, characters and logicals are rejected: if your data hold stimulus labels, check them against the generated stimulus filenames before converting them to numbers.
#' @param responses Vector of responses in the same order as \code{stimuli}: 1 where the original stimulus was chosen, -1 where the inverted one was.
#' @param baseimage String naming the base image: not its file name, but its key in the \code{base_face_files} list passed to \code{\link{generateStimuli2IFC}}.
#' @param rdata Path to the \code{.Rdata} file written when the stimuli were generated. It holds the contrast parameters of every stimulus.
#' @param targetci Optional target CI to correlate the cumulative CIs with, as returned by \code{\link{generateCI}}. Without it, the final CI of these trials is used.
#' @param step Number of trials added between successive correlations.
#' @return Vector of correlations between each cumulative CI and the final or target CI.
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
#' responses <- sample(c(1, -1), 6, replace = TRUE)
#' correlations <- suppressWarnings(computeCumulativeCICorrelation(
#'   stimuli = 1:6, responses = responses, baseimage = "face", rdata = rdata_file
#' ))
computeCumulativeCICorrelation <- function(stimuli, responses, baseimage, rdata, targetci = list(), step = 1) {

  # Coerce to plain vectors: tibble columns stay one-column tibbles rather than
  # dropping to vectors (see the same handling in generateCI()).
  trials <- coerceTrialVectors(stimuli, responses, NA)
  stimuli <- trials$stimuli
  responses <- trials$responses

  # load() assigns straight into this function's frame, so any object stored in
  # the .Rdata file silently overwrites an argument of the same name - the same
  # hazard handled in generateCI() and generateReferenceDistribution2IFC(). No
  # field written today collides with these arguments, so this changes nothing
  # now; it is here because the collision that actually bit us (the z-map
  # `sigma`, fixed in #146) was created from the .Rdata side, by adding a field,
  # not by adding an argument.
  #
  # Captured after the unlist() above, so the coerced vectors are what gets
  # restored rather than the original tibble columns.
  .args <- captureArgs(environment())

  # Load parameter file (created when generating stimuli)
  loadRdata(rdata, environment())

  list2env(.args, envir = environment())

  # Check whether critical variables have been loaded
  if (!exists('s', envir = environment(), inherits = FALSE) && !exists('p', envir = environment(), inherits = FALSE)) {
    stop('File specified in rdata argument did not contain s or p variable.', rdataWriterNote(environment()))
  }

  if (!exists('base_faces', envir = environment(), inherits = FALSE)) {
    stop('File specified in rdata argument did not contain base_faces variable.', rdataWriterNote(environment()))
  }

  if (!exists('stimuli_params', envir = environment(), inherits = FALSE)) {
    stop('File specified in rdata argument did not contain stimuli_params variable.', rdataWriterNote(environment()))
  }

  # Convert s to p (if rdata file originates from pre-0.3.3)
  if (exists('s', envir = environment(), inherits = FALSE)) {
    p <- list(patches = s$sinusoids, patchIdx = s$sinIdx, noise_type = 'sinusoid')
    rm(s)
  }


  # Get base image
  base <- base_faces[[baseimage]]
  if (is.null(base)) {
    stop(paste0('File specified in rdata argument did not contain any reference to base image label: ', baseimage, ' (NOTE: file contains references to the following base image label(s): ', paste(names(base_faces), collapse = ', '), ')'))
  }


  # Retrieve parameters of actually presented stimuli (this will work with
  # non-consecutive stims as well). drop = FALSE keeps a single presented stimulus
  # a one-row matrix rather than a vector, so the params[1:trial, ] slice in the
  # cumulative loop below stays valid; without it a length-1 `stimuli` aborted
  # with "incorrect number of dimensions" regardless of parameter count.
  validateStimulusIds(stimuli, nrow(stimuli_params[[baseimage]]))
  params <- stimuli_params[[baseimage]][stimuli, , drop = FALSE]

  # Check whether parameters were found in this .rdata file
  if (length(params) == 0) {
    stop(paste0('No parameters found for base image: ', baseimage))
  }

  # Truncate a pre-0.3.0 parameter set from 4096 to 4092, exactly as generateCI()
  # does, so this function can read the same old files. Without it the extra four
  # unused contrasts reach generateNoiseImage() as a length mismatch and abort.
  # See generateCI() and ChangeLog 0.3.0-29 for why 4096 was over-allocated.
  if (ncol(params) == 4096) {
    params <- params[, 1:4092, drop = FALSE]
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
  for (trial in seq(1, length(responses), step)) {
    setTxtProgressBar(pb, trial)

    cumCI <- generateCINoise(params[1:trial, ], responses[1:trial], p)
    correlations[corcounter] <- cor(as.vector(cumCI), as.vector(finalCI),
                                    use = 'pairwise.complete.obs')
    corcounter <- corcounter + 1
  }
  close(pb)

  # Return correlations
  return(correlations)
}
