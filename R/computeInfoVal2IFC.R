#' Computes Informational Value
#'
#' Computes Informational Value for a single CI in a 2IFC task.
#'
#' The Informational Value metric can be considered as a z-score that quantifies the signal
#' present in a classification image. The higher the Informational Value, the more signal. It is
#' possible to use a cut-off such as z = 1.96 to select classification images with significant
#' signal under alpha = 0.05.
#'
#' Informational Value is computed by simulating random responding under identical task parameters to
#' an empirical dataset. The metric quantifies how unlikely it is to observe these data under the
#' null-hypothesis that there is no signal (i.e., that there is only random responding).
#'
#' For more information see Brinkman, Goffin, Aarts, van Haren, & Dotsch (in prep).
#'
#' @export
#' @importFrom stats mad median
#' @param ci A classification image object (list-type) as returned by generateCI
#' @param rdata String pointing to .RData file that was created when stimuli were generated. This file contains the contrast parameters of all generated stimuli and possibly its corresponding reference distribution generated with generateReferenceDistribution().
#' @return Informational value (z-score)

computeInfoVal2IFC <- function(ci, rdata) {

  # Load parameter file (created when generating stimuli)
  load(rdata)

  # Check whether reference norms are present. If not, re-generate
  if (!exists("reference_norms", envir=environment(), inherits=FALSE)) {

    # Reference norms not present in rdata file, re-generate
    generateReferenceDistribution2IFC(rdata)

    # Re-load rdata file
    load(rdata)

    write("Note that now that this simulated reference distribution has been saved to the .Rdata file, the next time you call computeInfoVal2IFC(), it will not need to be computed again.")

  } else {

    write("Reference distribution found in rdata file.", stdout())

  }

  # Compute informational value metric
  infoVal <- (norm(ci$ci) - median(reference_norms) ) / mad(reference_norms)

  write( paste0("Informational value: z = ", infoVal, " (reference median = ", median(reference_norms), "; MAD = ", mad(reference_norms), "; iterations = ", length(reference_norms),  ")"), stdout() )

  return(infoVal)

}