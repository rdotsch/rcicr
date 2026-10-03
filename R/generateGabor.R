#' Generate single gabor patch
#'
#' @export
#' @param img_size Integer specifying size of gabor patch in number of pixels.
#' @param cycles Integer specifying number of cycles the sinusoid should span.
#' @param angle Value specifying the angle (rotation) of the sinusoid.
#' @param phase Value specifying phase of sinusoid.
#' @param sigma Value specifying the standard deviation, in pixels, of the
#'   Gaussian mask applied on top of the sinusoid.
#' @param contrast Value between -1.0 and 1.0 specifying contrast of sinusoid.
#' @return The gabor patch image with size \code{img_size}.
#' @examples
#' generateGabor(512, 2, 90, pi/2, 25, 1.0)
generateGabor <- function(img_size, cycles, angle, phase, sigma, contrast) {

  s <- generateSinusoid(img_size, cycles, angle, phase, contrast)
  # Coordinates from -0.5 to 0.5, as scales::rescale() gave them (#208): 0 for
  # a single pixel.
  x0 <- if (img_size == 1) 0 else (seq_len(img_size) - 1) / (img_size - 1) - 0.5
  gauss_x <- matrix(x0, img_size, img_size, byrow = TRUE)
  gauss_y <- matrix(x0, img_size, img_size)
  gauss_mask <- exp(-(((gauss_x^2) + (gauss_y^2)) / (2 * (sigma / img_size)^2)))
  return(gauss_mask * s)

}
