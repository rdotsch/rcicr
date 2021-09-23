
#' Convert angle in degrees to radians
#' @export
#' @param deg Angle in degrees
#' @examples
#' deg2rad(180)
deg2rad <- function(deg) {
  (deg * pi) / (180)
}