#' Generate single noise image based on parameter vector
#'
#' @export
#' @param params Vector with each value specifying the contrast of each patch in noise.
#' @param p 3D patch matrix (generated using \code{generateNoisePattern()}).
#' @return The noise pattern as pixel matrix.
#' @examples
#' p <- generateNoisePattern(img_size = 32, nscales = 2)
#'
#' # one contrast weight per patch index
#' params <- rnorm(max(p$patchIdx))
#'
#' noise <- generateNoiseImage(params, p)
generateNoiseImage <- function(params, p) {

  # Normalise a pre-0.3.3 noise pattern before anything reads p$patchIdx. This
  # rename used to happen *after* the length check below, so an old-style p
  # hit max(NULL) and always aborted - meaning the backward compatibility this
  # code provides never actually worked.
  if ('sinusoids' %in% names(p)) {
    p <- list(patches = p$sinusoids, patchIdx = p$sinIdx, noise_type = 'sinusoid')
  }

  # Abort stimulus generation if number of params doesn't equal number of patches
  if (length(params) != max(p$patchIdx)) {

    if ((length(params) == max(p$patchIdx) + 1) && (min(p$patchIdx) == 0)) {
      # Some versions of dependencies created patch indices starting with 0, latest dependencies
      # start counting at 1. Fix this.

      warning("Rdata patch indices start at 0, whereas parameters are used from position 1. Due to this mismatch, one sinusoid will not be shown in resulting CI.")

    } else {
      stop("Stimulus generation aborted: number of parameters doesn't equal number of patches!")

    }
  }

  # Average the weighted patches over the third dimension, i.e. one mean per
  # pixel across all patch layers. rowMeans() with dims=2 does exactly that.
  # It replaces apply(..., 1:2, mean), which was ~31x slower at the averaging
  # step (1.47s -> 0.05s at 512px, nscales=5); the whole call is ~6x faster,
  # the rest being the weighted array above, which both forms have to build.
  # Note that dims MUST be given: rowMeans() on a 3-D array defaults to dims=1,
  # which collapses dimensions 2 and 3 and silently returns a vector that
  # array() would then recycle. See PR #122.
  noise <- rowMeans(p$patches * array(params[p$patchIdx], dim(p$patches)), dims = 2)
  return(noise)

}
