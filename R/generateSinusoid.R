#' Generate single sinusoid patch
#'
#' @export
#' @import matlab
#' @param img_size Integer specifying size of sinusoid patch in number of pixels.
#' @param cycles Integer specifying number of cycles sinusoid should span.
#' @param angle Value specifying the angle (rotation) of the sinusoid.
#' @param phase Value specifying phase of sinusoid.
#' @param contrast Value between -1.0 and 1.0 specifying contrast of sinusoid.
#' @return The sinusoid image with size \code{img_size}.
#' @examples
#' generateSinusoid(512, 2, 90, pi/2, 1.0)
generateSinusoid <- function(img_size, cycles, angle, phase, contrast) {

  # Generates an image matrix containing a sinusoid, angle (in degrees) of 0 will give vertical, 90 horizontally oriented sinusoid
  angle <- deg2rad(angle)
  sinepatch = matlab::repmat(matlab::linspace(0, cycles, img_size), img_size, 1)
  sinusoid <- (sinepatch * cos(angle) + t(sinepatch) * sin(angle)) * 2 * pi
  sinusoid <- contrast * sin(sinusoid + phase)
  return(sinusoid)
}
