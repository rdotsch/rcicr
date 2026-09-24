#' Generate classification image noise pattern based on set of stimuli (matrix: trials, parameters), responses (vector), and sinusoid
#'
#' @export
#' @param stimuli Matrix with one parameter row per response. A parameter vector is also accepted for a single trial with exactly one response.
#' @param responses Vector with the response to each trial: 1 where the participant chose the original, -1 where they chose the inverted image. Other values, such as ratings on a scale, weight the trials accordingly.
#' @param p Noise basis, as returned by \code{\link{generateNoisePattern}}.
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

  n_trials <- if (is.null(dim(stimuli))) 1L else nrow(stimuli)
  if ((!is.null(dim(stimuli)) && !is.matrix(stimuli)) ||
        length(stimuli) == 0L || length(responses) != n_trials) {
    stop('stimuli must have one parameter row per response, or be a single-trial vector with one response.')
  }

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
