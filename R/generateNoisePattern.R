#' Generate sinusoid noise pattern
#'
#' @export
#' @param img_size Width and height of the noise pattern, in pixels. It must be divisible by \code{2^(nscales - 1)}, because the finest scale tiles the image with that many patches per side.
#' @param nscales Number of spatial scales (default: 5). Each additional scale adds a higher spatial frequency. \code{img_size} must be divisible by \code{2^(nscales - 1)}.
#' @param noise_type Noise pattern type: \code{sinusoid} (default) or \code{gabor}.
#' @param sigma Sigma of the Gabor patches when \code{noise_type = 'gabor'} (default: 25).
#' @param pre_0.3.0 Boolean: build the noise pattern the way rcicr did before version 0.3.0. Leave it at \code{FALSE} unless you need to recreate that old behaviour.
#' @return The noise basis: a list holding the 3D array of patches (\code{patches}), an index
#' array of the same size mapping each pixel to its contrast parameter (\code{patchIdx}), and
#' the \code{noise_type} and \code{generator_version}.
#' @examples
#' generateNoisePattern(256)
generateNoisePattern <- function(img_size = 512, nscales = 5, noise_type = 'sinusoid', sigma = 25, pre_0.3.0 = FALSE) { # nolint: object_name_linter.
  validateTiling(img_size, nscales)

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

# The finest scale tiles the image with 2^(nscales - 1) patches per side; a
# fractional patch size is truncated and the tiling then fails with a
# replacement-length error that does not say why (#339).
validateTiling <- function(img_size, nscales) {
  tiles <- 2^(nscales - 1)
  # Only a value that can be checked is; anything else keeps its old behaviour.
  if (!isTRUE(img_size %% tiles != 0)) return(invisible(NULL))
  sizes <- c(floor(img_size / tiles), ceiling(img_size / tiles)) * tiles
  sizes <- sizes[sizes > 0]
  max_scales <- 1
  while (img_size %% 2^max_scales == 0) max_scales <- max_scales + 1
  stop("img_size (", img_size, ") must be divisible by 2^(nscales - 1) = ", tiles,
       " for nscales = ", nscales, ": the finest scale tiles the image with ", tiles,
       " patches per side. Use img_size ", paste(sizes, collapse = " or "),
       ", or nscales = ", max_scales, if (max_scales > 1) " or fewer", ".", call. = FALSE)
}
