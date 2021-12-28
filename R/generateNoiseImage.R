#' Generate single noise image based on parameter vector
#'
#' @export
#' @param params Vector with each value specifying the contrast of each patch in noise.
#' @param p 3D patch matrix (generated using \code{generateNoisePattern()}).
#' @return The noise pattern as pixel matrix.
#' @examples
#' #params <- rnorm(4092) # generates 4092 normally distributed random values
#' #s <- generateNoisePattern(img_size=256)
#' #noise <- generateNoiseImage(params, p)
generateNoiseImage <- function(params, p) {

  # Abort stimulus generation if number of params doesn't equal number of patches
  if (length(params) != max(p$patchIdx)) {

    if ((length(params) == max(p$patchIdx) + 1)& (min(p$patchIdx) == 0)) {
      # Some versions of dependencies created patch indices starting with 0, latest dependencies
      # start counting at 1. Fix this.

      warning("Rdata patch indices start at 0, whereas parameters are used from position 1. Due to this mismatch, one sinusoid will not be shown in resulting CI.")

    } else {
      stop("Stimulus generation aborted: number of parameters doesn't equal number of patches!")

    }
  }

  if ('sinusoids' %in% names(p)) {
    # Pre 0.3.3 noise pattern, rename for appropriate use
    p <- list(patches=p$sinusoids, patchIdx=p$sinIdx, noise_type='sinusoid')
  }

  noise <- apply(p$patches * array(params[p$patchIdx], dim(p$patches)), 1:2, mean)
  return(noise)

}