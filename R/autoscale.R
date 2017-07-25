#' Determines optimal scaling constant for a list of CIs
#'
#' @export
#' @import matlab
#' @import png
#' @param cis List of cis, each of which are a list containing the pixel matrices of at least the noise pattern (\code{$ci}) and if the noise patterns need to be written to PNGs, also the base image (\code{$base}).
#' @param save_as_pngs Boolean, when set to true, the autoscaled noise patterns will be combined with their respective base images and saved as PNGs (using the key of the list as name).
#' @param targetpath Optional string specifying path to save PNGs to (default: ./cis).
#' @return List of scaled noise patterns and determined scaling factor.
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
