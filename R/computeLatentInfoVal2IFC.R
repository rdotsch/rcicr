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
#' Each iteration rebuilds the direction through the same computation that
#' produced the observed one, so the cost grows with the number of trials as
#' well as with \code{iter}: a 300-trial study takes a few seconds at the
#' default. That is still far below what the pixel-noise equivalent costs, which
#' has to render.
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
#' @param participants Vector of participant identifiers, the same vector passed
#'   to \code{\link{generateLatentCI}}. Defaults to NA, matching its default.
#'   The null is built through the same estimator the classification image was,
#'   so a per-participant image is compared against a per-participant null.
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
#' computeLatentInfoVal2IFC(
#'   ci, stimuli$rdata, 1:20, rep(c(1, -1), 10), iter = 200
#' )$infoVal
computeLatentInfoVal2IFC <- function(latent_ci, rdata, stimuli, responses, participants = NA, iter = 10000, response_seed = NULL) {

  requireExperimental('computeLatentInfoVal2IFC')

  stored <- loadLatentStimulusParams(rdata)

  checkParticipants(participants, stimuli, 'computeLatentInfoVal2IFC')
  checkStimuliAreTrials(stimuli, 'computeLatentInfoVal2IFC')

  # The trials and the answers, not just the stimulus file. Fingerprinting the
  # file alone still accepts a classification image built from a different
  # subset of its trials, or from different responses to the same ones, and
  # scores it against a null it has nothing to do with.
  if (!identical(latent_ci$analysis_inputs,
                 analysisInputs(stored, stimuli, responses, participants))) {
    msg <- paste0(
      'This classification image was not computed from these inputs. Pass the ',
      'stimulus file, stimuli, responses and participants that ',
      'generateLatentCI() was given.'
    )
    stop(msg, call. = FALSE)
  }

  latent_sd <- stored$generator_spec$latent_sd
  observed <- latentNorm(latent_ci$direction, latent_sd)
  trial_responses <- unlist(responses, use.names = FALSE)

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
  #
  # Permuted at trial level and then put through latentDirection(), the same
  # function that produced the observed direction. Permuting after aggregation
  # would shuffle means rather than answers, and pooling the trials of a
  # per-participant classification image would compare it against a null built
  # by a different estimator.
  reference <- numeric(iter)
  for (i in seq_len(iter)) {
    permuted <- latentDirection(stored$latent_params, stimuli,
                                permuteResponses(trial_responses, participants),
                                participants)
    reference[i] <- latentNorm(permuted$direction, latent_sd)
  }

  reference_median <- stats::median(reference)
  reference_mad <- stats::mad(reference)

  return(list(
    infoVal = standardiseAgainstNull(observed, reference, reference_median, reference_mad),
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

# Shuffled within each participant, not across all of them. A global shuffle
# moves answers between people, so a participant who always pressed one key and
# one who always pressed the other would produce a null full of mixed answers
# that neither of them could have given, and the difference between their key
# biases would be reported as latent signal. Within a participant the multiset
# each person actually supplied is preserved.
permuteResponses <- function(responses, participants) {
  if (all(is.na(participants))) {
    return(sample(responses))
  }

  participants <- unlist(participants, use.names = FALSE)
  permuted <- responses
  for (pid in unique(participants)) {
    rows <- participants == pid
    permuted[rows] <- sample(responses[rows])
  }

  return(permuted)
}

# (observed - median) / mad is 0 / 0 when every permutation gives the same
# norm, which is not a corner case here: it is exactly what a participant who
# pressed one key throughout produces, since their responses have one
# arrangement. The limit is what is reported, so the number stays readable in
# the direction it means, and p is unaffected either way.
standardiseAgainstNull <- function(observed, reference, reference_median, reference_mad) {
  if (reference_mad > 0) {
    return((observed - reference_median) / reference_mad)
  }

  # A zero MAD does not mean a null with no spread. It means more than half the
  # permutations landed on one value, which repeated stimuli make ordinary:
  # stimuli c(1, 1, 2) answered c(1, 1, -1) put two thirds of the arrangements on
  # one weighting. Reporting an infinity there would put no scale at all beside a
  # perfectly finite tail probability.
  #
  # The standard deviation is a scale for that case; the limit is only for a null
  # that genuinely cannot vary, which is the all-one-answer case the permutation
  # null exists to handle.
  spread <- stats::sd(reference)
  if (spread > 0) {
    return((observed - reference_median) / spread)
  }

  if (isTRUE(all.equal(observed, reference_median))) {
    return(0)
  }

  return(sign(observed - reference_median) * Inf)
}

# The stimulus set plus which of its trials were analysed and what was answered.
# Participants are folded in as their factor levels, so a per-participant
# classification image cannot be scored against a pooled null.
analysisInputs <- function(stored, stimuli, responses, participants) {
  trials <- coerceTrialVectors(stimuli, responses, participants)
  pids <- if (all(is.na(trials$participants))) NA_character_ else as.character(trials$participants)

  # The trials themselves, not summary statistics of them. Six statistics over a
  # vector of 1s and -1s carry almost nothing: the number of each is fixed by the
  # sum, so only the index-weighted sum separates the arrangements and collisions
  # are easy to write down -- c(1, -1, -1, 1, -1, -1) and c(-1, 1, 1, -1, -1, -1)
  # agree on all six. The vectors are a few hundred numbers, so keeping them is
  # both exact and cheaper than any digest base R could build.
  #
  # Participants as characters, so a factor and the character vector it came from
  # are the same analysis.
  # As doubles, because the storage type is not part of the analysis: stimuli
  # written as 1:20 are integer and the same values read back from a CSV are
  # double, and identical() separates them. Both select and weight the same
  # trials, so a classification image must not be refused against its own data
  # for having made a round trip through a file.
  return(list(stimulus = stimulusFingerprint(stored),
              stimuli = as.numeric(trials$stimuli),
              responses = as.numeric(trials$responses),
              participants = pids))
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
