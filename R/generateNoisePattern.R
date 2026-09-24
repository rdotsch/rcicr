#' Generate sinusoid noise pattern
#'
#' @export
#' @param img_size Width and height of the noise pattern, in pixels.
#' @param nscales Number of spatial scales (default: 5). Each additional scale adds a higher spatial frequency.
#' @param noise_type Noise pattern type: \code{sinusoid} (default) or \code{gabor}.
#' @param sigma Sigma of the Gabor patches when \code{noise_type = 'gabor'} (default: 25).
#' @param pre_0.3.0 Boolean: build the noise pattern the way rcicr did before version 0.3.0. Leave it at \code{FALSE} unless you need to recreate that old behaviour.
#' @return The noise basis: a list holding the 3D array of patches (\code{patches}), an index
#' array of the same size mapping each pixel to its contrast parameter (\code{patchIdx}), and
#' the \code{noise_type} and \code{generator_version}.
#' @examples
#' generateNoisePattern(256)
generateNoisePattern <- function(img_size = 512, nscales = 5, noise_type = 'sinusoid', sigma = 25, pre_0.3.0 = FALSE) { # nolint: object_name_linter.
  # Settings of sinusoids
  orientations <- c(0, 30, 60, 90, 120, 150)
  phases <- c(0, pi / 2)
  scales <- 2^(0:(nscales - 1))

  # Size of patches per scale
  mg <- matlab::meshgrid(1:img_size, 1:img_size, 1:length(scales)) # nolint: seq_linter.
  patchSize <- mg$x / mg$y

  # Number of patch layers needed
  nrPatches <- length(scales) * length(orientations) * length(phases)

  # Preallocate memory
  patches <- matlab::zeros(c(img_size, img_size, nrPatches))
  patchIdx <- matlab::zeros(c(img_size, img_size, nrPatches))

  # Counters
  if (pre_0.3.0) {
    co <- 0 # patch layer counter
    idx <- 0 # contrast index counter
  } else {
    co <- 1 # patch layer counter
    idx <- 1 # contrast index counter
  }

  for (scale in scales) {
    for (orientation in orientations) {
      for (phase in phases) {
        # Generate single patch
        size <- patchSize[scale, img_size]

        if (noise_type == 'gabor') {
          p <- generateGabor(size, 1.5, orientation, phase, sigma, 1)
        } else {
          p <- generateSinusoid(size, 2, orientation, phase, 1)
        }

        # Repeat to fill scale
        patches[, , co] <- matlab::repmat(p, scale)

        # Create index matrix
        for (col in 1:scale) {
          for (row in 1:scale) {

            # Insert absolute index for later contrast weighting
            patchIdx[(size * (row - 1) + 1):(size * row), (size * (col - 1) + 1):(size * col), co] <- idx

            # Update contrast counter
            idx <- idx + 1

          }
        }

        # Update layer counter
        co <- co + 1

      }
    }
  }

  return(list(patches = patches, patchIdx = patchIdx, noise_type = noise_type, generator_version = utils::packageVersion('rcicr')))
}
