#' Determines optimal scaling constant for a list of CIs
#'
#' @export
#' @import matlab
#' @import png
#' @param cis List of cis, each of which are a list containing the pixel matrices of at least the noise pattern (\code{$ci}) and if the noise patterns need to be written to PNGs, also the base image (\code{$base}).
#' @param save_as_pngs Boolean, when set to true, the autoscaled noise patterns will be combined with their respective base images and saved as PNGs (using the key of the list as name).
#' @param targetpath Optional string specifying path to save PNGs to (default: ./cis).
#' @return The input \code{cis} list, with a \code{$scaled} pixel matrix added to each element. The
#' scaling constant that was determined is printed to the console, not returned.
#' @examples
#' cis <- list(
#'   participant1 = list(ci = matrix(runif(64, -0.2, 0.2), 8, 8), base = matrix(0.5, 8, 8)),
#'   participant2 = list(ci = matrix(runif(64, -0.3, 0.3), 8, 8), base = matrix(0.5, 8, 8))
#' )
#' scaled_cis <- autoscale(cis, save_as_pngs = FALSE)
autoscale <- function(cis, save_as_pngs=TRUE, targetpath='./cis') {

  # Get range of each ci
  ranges <- matlab::zeros(length(names(cis)), 2)
  for (ciname in names(cis)) {
    ranges[which(ciname==names(cis)), ] <- range(cis[[ciname]]$ci)
  }

  # Determine the lowest possible scaling factor constant
  if (abs(min(ranges[,1])) > max(ranges[,2])) {
    constant <- abs(min(ranges[,1]))
  }  else {
    constant <- max(ranges[,2])
  }

  write(paste0("Using scaling factor constant:", constant), stdout())

  # Scale all noise patterns
  for (ciname in names(cis)) {
    cis[[ciname]]$scaled <-  (cis[[ciname]]$ci + constant) / (2*constant)

    # Combine and save to PNG if necessary
    if (save_as_pngs) {
      ci <- (cis[[ciname]]$scaled + cis[[ciname]]$base) / 2

      dir.create(targetpath, recursive=T, showWarnings = F)

      png::writePNG(ci, paste0(targetpath, '/', ciname, '_autoscaled.png'))
    }

  }

  return(cis)
}
