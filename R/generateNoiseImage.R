#' Generate single noise image based on parameter vector
#'
#' @export
#' @param params Vector of contrast weights, one per patch in the noise basis.
#' @param p Noise basis, as returned by \code{\link{generateNoisePattern}}.
#' @return The noise pattern as pixel matrix.
#' @examples
#' p <- generateNoisePattern(img_size = 32, nscales = 2)
#'
#' # one contrast weight per patch index
#' params <- rnorm(max(p$patchIdx))
#'
#' noise <- generateNoiseImage(params, p)
generateNoiseImage <- function(params, p) {

  # A pre-0.3.3 noise pattern, normalised before anything reads p$patchIdx
  # (NEWS 1.1.0).
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

  # One mean per pixel across the patch layers. dims = 2 is required: on a 3-D
  # array rowMeans() defaults to dims = 1 and returns a vector, which array()
  # would silently recycle. ~31x faster than apply(..., 1:2, mean) (#122).
  noise <- rowMeans(p$patches * array(params[p$patchIdx], dim(p$patches)), dims = 2)
  return(noise)

}

# The basis as the render loops use it, with its patch indices held as integer:
# params[patchIdx] then skips converting every index, which made a 512-pixel
# render about a quarter faster and halves what each parallel worker is sent
# (analyses/worker-memory.md). A working copy only; the .Rdata file keeps the
# type it was saved with.
renderingBasis <- function(p) {
  for (idx in intersect(c('patchIdx', 'sinIdx'), names(p))) {
    if (is.double(p[[idx]]) && all(p[[idx]] == trunc(p[[idx]])) &&
          all(abs(p[[idx]]) <= .Machine$integer.max)) {
      storage.mode(p[[idx]]) <- 'integer'
    }
  }
  p
}
