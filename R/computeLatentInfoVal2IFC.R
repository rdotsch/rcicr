#' Informational value of a latent-space classification image
#'
#' How much signal a latent classification image carries, measured against a
#' null distribution built by permuting the responses that produced it.
#'
#' The direction a participant produces has some length whatever they did, so
#' the picture alone cannot separate a consistent participant from one pressing
#' keys at random, and \code{latent_scaling = 'sd'} removes the difference
#' entirely by construction. This gives the number that does separate them.
#'
#' The null is built by re-drawing the responses at random over the same trials
#' and recomputing the direction, so it holds the stimuli fixed and asks only
#' what the responses contributed. It needs no rendering, which is why it is
#' cheap here where the equivalent for the pixel pipeline is expensive enough
#' that \code{\link{computeInfoVal2IFC}} caches it in the stimulus file.
#'
#' Length is measured in standard deviations of the generator's training faces,
#' summed across dimensions, so a generator whose components have very different
#' scales does not have its widest dimension decide the answer on its own.
#'
#' This is not the formula \code{\link{computeInfoVal2IFC}} uses. That one
#' compares against a simulated reference distribution of noise images and its
#' definition is fixed by the published erratum. The two numbers are not
#' comparable with each other.
#'
#' @export
#' @param latent_ci The list returned by \code{\link{generateLatentCI}}.
#' @param rdata String pointing to the .Rdata file written by
#'   \code{\link{generateStimuliLatent2IFC}}, the same one the classification
#'   image was computed from.
#' @param stimuli Vector of stimulus numbers, in the order they were presented:
#'   the same vector passed to \code{\link{generateLatentCI}}.
#' @param iter Number of permutations in the null distribution.
#' @param response_seed Integer seeding the permutations, for a reproducible
#'   null. Defaults to NULL, which leaves the caller's random number stream
#'   alone.
#' @return A list with the informational value, the observed length, the
#'   proportion of the null at least as long as the observed direction, and the
#'   null's median and median absolute deviation.
#' @examples
#' # This function is part of the experimental latent-space module.
#' options(rcicr.experimental = TRUE)
#'
#' faces <- replicate(6, tempfile(fileext = ".png"))
#' for (i in seq_along(faces)) {
#'   png::writePNG(matrix(runif(16 * 16), 16, 16), faces[i])
#' }
#' generator <- latentGeneratorPCA(faces, n_components = 3, img_size = 16)
#'
#' stimuli <- generateStimuliLatent2IFC(
#'   generator, n_trials = 20, stimulus_path = tempdir(),
#'   seed = 1, save_as_png = FALSE
#' )
#' ci <- generateLatentCI(
#'   stimuli = 1:20, responses = rep(c(1, -1), 10),
#'   rdata = stimuli$rdata, save_as_png = FALSE
#' )
#'
#' computeLatentInfoVal2IFC(ci, stimuli$rdata, 1:20, iter = 200)$infoVal
computeLatentInfoVal2IFC <- function(latent_ci, rdata, stimuli, iter = 10000, response_seed = NULL) {

  requireExperimental('computeLatentInfoVal2IFC')

  stored <- loadLatentStimulusParams(rdata)

  if (!identical(latent_ci$generator_fingerprint, stored$generator_spec$fingerprint)) {
    msg <- paste0(
      'This classification image was not computed from these stimuli: its ',
      'generator fingerprint does not match the one in the stimulus file.'
    )
    stop(msg, call. = FALSE)
  }

  params <- selectLatentParams(stored$latent_params, unlist(stimuli, use.names = FALSE))
  latent_sd <- stored$generator_spec$latent_sd

  observed <- latentNorm(latent_ci$direction, latent_sd)

  # response_seed, not seed: the stimuli were already drawn under the stimulus
  # file's own seed, and reusing that name would suggest this reproduces them.
  if (!is.null(response_seed)) {
    set.seed(response_seed)
  }

  n_trials <- nrow(params)
  reference <- numeric(iter)
  for (i in seq_len(iter)) {
    responses <- ((stats::runif(n_trials) > 0.5) * 2) - 1
    reference[i] <- latentNorm(weightedLatentMean(params, responses), latent_sd)
  }

  reference_median <- stats::median(reference)
  reference_mad <- stats::mad(reference)

  return(list(
    infoVal = (observed - reference_median) / reference_mad,
    observed_norm = observed,
    p = mean(reference >= observed),
    reference_median = reference_median,
    reference_mad = reference_mad,
    iter = iter
  ))
}

# Length in units of the generator's own per-dimension spread, so the answer
# does not depend on how the latent space happens to be scaled.
latentNorm <- function(direction, latent_sd) {
  return(sqrt(sum((direction / latent_sd)^2)))
}
