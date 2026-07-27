#' Determines optimal scaling constant for a list of CIs
#'
#' @export
#' @import matlab
#' @import png
#' @param cis List of cis, each of which are a list containing the pixel matrices of at least the noise pattern (\code{$ci}) and if the noise patterns need to be written to PNGs, also the base image (\code{$base}).
#' @param save_as_pngs Boolean, when set to true, the autoscaled noise patterns will be combined with their respective base images and saved as PNGs (using the key of the list as name).
#' @param targetpath Optional string specifying path to save PNGs to (default: ./cis).
#' @return The input \code{cis} list, with the \code{$scaled} pixel matrix of each element
#' replaced by its autoscaled version, and \code{$combined} refreshed to match (where a
#' \code{$base} image is present). The scaling constant that was determined is printed to
#' the console, not returned.
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

    # Refresh $combined to match the new $scaled. It used to be left untouched,
    # so a CI that came from batchGenerateCI() (which scales with 'none' before
    # handing over to autoscale) kept a $combined built from the *unscaled*
    # noise -- a nearly invisible overlay -- while its $scaled was correct. The
    # two fields silently disagreed, and anyone plotting the returned $combined
    # got the un-autoscaled image. The PNGs written below were always right,
    # which is why this went unnoticed.
    if (!is.null(cis[[ciname]]$base)) {
      cis[[ciname]]$combined <- (cis[[ciname]]$scaled + cis[[ciname]]$base) / 2
    }

    # Save to PNG if necessary
    if (save_as_pngs) {
      dir.create(targetpath, recursive=T, showWarnings = F)

      png::writePNG(cis[[ciname]]$combined,
                    paste0(targetpath, '/', ciname, '_autoscaled.png'))
    }

  }

  return(cis)
}
