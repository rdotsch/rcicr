#' Generate classification image noise pattern based on set of stimuli (matrix: trials, parameters), responses (vector), and sinusoid
#'
#' @export
#' @param stimuli Matrix with one row per trial, each row containing the 4092 parameters for the original stimulus.
#' @param responses Vector containing the response to each trial (1 if participant selected original, -1 if participant selected inverted;
#' this can be changed into a scale).
#' @param p 3D patch matrix (generated using \code{generateNoisePattern()}).
#' @return The classification image as pixel matrix.
#' @examples
#' p <- generateNoisePattern(img_size = 32, nscales = 1)
#' nparams <- max(p$patchIdx)
#'
#' # two trials, one chosen original (1) and one inverted (-1)
#' stimuli <- matrix(runif(2 * nparams, -1, 1), nrow = 2)
#' responses <- c(1, -1)
#'
#' ci <- generateCINoise(stimuli, responses, p)
generateCINoise <- function(stimuli, responses, p) {

  weighted <- stimuli * responses

  # Only aggregate if more than one stimulus/response row
  if (is.null(dim(weighted))) {
    params <- weighted
  } else {
    # Compute mean and return to original variance
    params <- colMeans(weighted) #* sqrt(length(responses))
  }

  return(generateNoiseImage(params, p))
}
