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
#' The null shuffles the responses the participant actually gave over the same
#' trials and recomputes the direction, so it holds both the stimuli and the
#' responses fixed and asks only what their pairing contributed. Shuffling
#' rather than re-drawing matters whenever a participant pressed one key more
#' often than the other: fresh balanced responses would be a null that an
#' imbalanced participant could never have produced, and their bias would come
#' back as signal. It needs no rendering, which is why it is
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
#' @param responses Vector of responses, the same vector passed to
#'   \code{\link{generateLatentCI}}. These are what the null shuffles, so
#'   passing anything else makes the number meaningless.
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
#'   generator, n_trials = 20,
#'   stimulus_path = file.path(tempdir(), "rcicr_latent_infoval_example"),
#'   seed = 1, save_as_png = FALSE
#' )
#' ci <- generateLatentCI(
#'   stimuli = 1:20, responses = rep(c(1, -1), 10),
#'   rdata = stimuli$rdata, save_as_png = FALSE
#' )
#'
#' computeLatentInfoVal2IFC(ci, stimuli$rdata, 1:20, rep(c(1, -1), 10), iter = 200)$infoVal
computeLatentInfoVal2IFC <- function(latent_ci, rdata, stimuli, responses, iter = 10000, response_seed = NULL) {

  requireExperimental('computeLatentInfoVal2IFC')

  stored <- loadLatentStimulusParams(rdata)

  # The stimulus set, not just the generator. Two experiments run through one
  # generator have identical generator fingerprints by design, so checking only
  # that would accept a classification image from one experiment scored against
  # the other's perturbations, which produces a plausible and meaningless
  # number.
  if (!identical(latent_ci$stimulus_fingerprint, stimulusFingerprint(stored))) {
    msg <- paste0(
      'This classification image was not computed from these stimuli. Pass the ',
      'stimulus file that generateLatentCI() was given, and the same stimuli ',
      'and responses.'
    )
    stop(msg, call. = FALSE)
  }

  # Aggregated the way generateLatentCI() aggregates when pooling, so the
  # observed direction and the null are computed over the same trials.
  trials <- coerceTrialVectors(stimuli, responses, NA)
  aggregated <- aggregateResponses(trials$stimuli, trials$responses)

  params <- selectLatentParams(stored$latent_params, aggregated$stimuli)
  latent_sd <- stored$generator_spec$latent_sd
  observed_responses <- aggregated$responses

  observed <- latentNorm(latent_ci$direction, latent_sd)

  # response_seed, not seed: the stimuli were already drawn under the stimulus
  # file's own seed, and reusing that name would suggest this reproduces them.
  if (!is.null(response_seed)) {
    set.seed(response_seed)
  }

  # A permutation of what the participant actually answered, not fresh coin
  # flips. Drawing balanced signs would compare an imbalanced response vector
  # against a null it could never have produced, so a participant who pressed
  # one key more often than the other would score as carrying signal. The
  # multiset of responses is held fixed and only its pairing with the stimuli
  # is broken, which is the question being asked.
  reference <- numeric(iter)
  for (i in seq_len(iter)) {
    permuted <- sample(observed_responses)
    reference[i] <- latentNorm(weightedLatentMean(params, permuted), latent_sd)
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

# Identifies the perturbations themselves, so a classification image can be tied
# to the stimulus set it came from and not merely to the generator that rendered
# it. Derived rather than stored, so it costs the .Rdata contract nothing.
stimulusFingerprint <- function(stored) {
  return(paste0(stored$generator_spec$fingerprint, '|',
                fingerprintNumeric(stored$latent_params, stored$base_latent)))
}

# Length in units of the generator's own per-dimension spread, so the answer
# depends neither on how the latent space is scaled nor on how many dimensions
# it has. The root mean square rather than the sum, matching what
# applyLatentScaling() means by a displacement of so many standard deviations:
# the two were written with different conventions, and the disagreement made a
# search step overshoot by a factor of sqrt(latent_dim).
#
# Informational value is unaffected by the change, being a ratio in which any
# constant factor cancels.
latentNorm <- function(direction, latent_sd) {
  return(sqrt(mean((direction / latent_sd)^2)))
}
