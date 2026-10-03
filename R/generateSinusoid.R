#' Generate single sinusoid patch
#'
#' @export
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
  # A ramp from 0 to `cycles` along each row. As matlab::linspace() and
  # repmat() built it (#208), a single pixel holds `cycles`, not 0.
  ramp <- if (img_size < 2) cycles else seq(0, cycles, length.out = img_size)
  sinepatch <- if (img_size == 1) ramp else matrix(ramp, img_size, img_size, byrow = TRUE)
  sinusoid <- (sinepatch * cos(angle) + t(sinepatch) * sin(angle)) * 2 * pi
  sinusoid <- contrast * sin(sinusoid + phase)
  return(sinusoid)
}
