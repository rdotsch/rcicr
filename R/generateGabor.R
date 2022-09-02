#' Generate single gabor patch
#'
#' @export
#' @import matlab
#' @param img_size Integer specifying size of gabor patch in number of pixels.
#' @param cycles Integer specifying number of cycles the sinusoid should span.
#' @param angle Value specifying the angle (rotation) of the sinusoid.
#' @param phase Value specifying phase of sinusoid.
#' @param sigma of guassian mask on top of sinusoid.
#' @param contrast Value between -1.0 and 1.0 specifying contrast of sinusoid.
#' @return The sinusoid image with size \code{img_size}.
#' @examples
#' generateSinusoid(512, 2, 90, pi/2, 1.0)
generateGabor <- function(img_size, cycles, angle, phase, sigma, contrast) {

  s <- generateSinusoid(img_size, cycles, angle, phase, contrast)
  x0 <- scales::rescale(1:img_size, to= c(-.5,.5))
  gauss <- matlab::meshgrid(x0, x0)
  gauss_mask = exp( -(((gauss$x^2)+(gauss$y^2)) / (2* (sigma/img_size)^2)) )
  return(gauss_mask * s)

}