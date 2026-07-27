#' Determines optimal scaling constant for a list of CIs
#'
#' @export
#' @import matlab
#' @import png
#' @param cis List of cis, each of which are a list containing the pixel matrices of at least the noise pattern (\code{$ci}) and if the noise patterns need to be written to PNGs, also the base image (\code{$base}).
#' @param save_as_pngs Boolean, when set to true, the autoscaled noise patterns will be combined with their respective base images and saved as PNGs (using the key of the list as name).
#' @param targetpath Optional string specifying path to save PNGs to (default: ./cis).
#' @return The input \code{cis} list, with the \code{$scaled} pixel matrix of each element
#' replaced by its autoscaled version. The scaling constant that was determined is printed
#' to the console, not returned.
#'
#' \strong{Look at \code{$scaled}, not \code{$combined}.} \code{$combined} is left exactly
#' as it was passed in: autoscaling deliberately does not disturb it, so a combination made
#' before this call survives unchanged. Only \code{$scaled} reflects the autoscaled result.
#' If you want the autoscaled noise shown over the base image, build it yourself with
#' \code{(ci$scaled + ci$base) / 2} — that is exactly what \code{save_as_pngs = TRUE}
#' writes to disk.
#'
#' This matters most after \code{\link{batchGenerateCI}} or
#' \code{\link{batchGenerateCI2IFC}}, which scale with \code{'none'} before handing over to
#' this function: their \code{$combined} is therefore an overlay of the \emph{unscaled}
#' noise and will look almost blank, while \code{$scaled} is the image you want.
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

    # Note that $combined is deliberately NOT updated here. It stays as the
    # caller supplied it, so whatever combination was made before autoscaling
    # survives this call unchanged. Use $scaled to see the autoscaled result;
    # that is the field this function exists to produce, and it is what the
    # PNG below is built from. Do not "fix" this by rewriting $combined -- it
    # would silently change what existing analysis scripts plot.
    if (save_as_pngs) {
      ci <- (cis[[ciname]]$scaled + cis[[ciname]]$base) / 2

      dir.create(targetpath, recursive=T, showWarnings = F)

      png::writePNG(ci, paste0(targetpath, '/', ciname, '_autoscaled.png'))
    }

  }

  return(cis)
}
