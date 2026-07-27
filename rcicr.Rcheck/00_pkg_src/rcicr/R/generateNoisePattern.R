#' Generate sinusoid noise pattern
#'
#' @export
#' @import matlab
#' @param img_size Integer specifying size of the noise pattern in number of pixels.
#' @param nscales Integer specifying the number of incremental spatial scales. Defaults to 5. Higher numbers will add higher spatial frequency scales.
#' @param noise_type String specifying noise pattern type (defaults to \code{sinusoid}; other options: \code{gabor}).
#' @param sigma Number specifying the sigma of the Gabor patch if noise_type is set to \code{gabor} (defaults to 25).
#' @param pre_0.3.0 Boolean specifying whether the noise pattern should be created in a way compatible with older versions of rcicr (< 0.3.0). If you are starting a new project, you should keep this at the default setting (FALSE). There is no reason to set this to TRUE, with the sole exception to recreate behavior of rcicr prior to version 0.3.0.
#' @return List with two elements: the 3D noise matrix with size \code{img_size}, and an indexing
#' matrix with the same size to easily change contrasts.
#' @examples
#' generateNoisePattern(256)
generateNoisePattern <- function(img_size=512, nscales=5, noise_type='sinusoid', sigma=25, pre_0.3.0=FALSE) {
  # Settings of sinusoids
  orientations <- c(0, 30, 60, 90, 120, 150)
  phases <- c(0, pi/2)
  scales <- 2^(0:(nscales-1))

  # Size of patches per scale
  mg <- matlab::meshgrid(1:img_size, 1:img_size,1:length(scales))
  patchSize = mg$x / mg$y

  # Number of patch layers needed
  nrPatches = length(scales) * length(orientations) * length(phases)

  # Preallocate memory
  patches = matlab::zeros(c(img_size, img_size, nrPatches))
  patchIdx = matlab::zeros(c(img_size, img_size, nrPatches))

  # Counters
  if (pre_0.3.0) {
    co = 0 # patch layer counter
    idx = 0 # contrast index counter
  } else {
    co = 1 # patch layer counter
    idx = 1 # contrast index counter
  }

  for (scale in scales) {
    for (orientation in orientations) {
      for (phase in phases) {
        # Generate single patch
        size <- patchSize[scale, img_size]

        if (noise_type=='gabor') {
          p <- generateGabor(size, 1.5, orientation, phase, sigma, 1)
        } else {
          p <- generateSinusoid(size, 2, orientation, phase, 1)
        }

        # Repeat to fill scale
        patches[,,co] <- matlab::repmat(p, scale)

        # Create index matrix
        for (col in 1:scale) {
          for (row in 1:scale) {

            # Insert absolute index for later contrast weighting
            patchIdx[(size * (row-1) + 1) : (size * row), (size * (col-1) + 1) : (size * col), co] = idx

            # Update contrast counter
            idx = idx + 1

          }
        }

        # Update layer counter
        co = co + 1

      }
    }
  }

  return(list(patches=patches, patchIdx=patchIdx, noise_type=noise_type, generator_version=utils::packageVersion('rcicr')))
}
