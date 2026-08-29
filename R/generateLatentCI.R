#' Computes a classification image in a generative model's latent space
#'
#' Turn the responses from a latent-space 2IFC task into a direction in latent
#' space, and render the face at the end of it.
#'
#' The arithmetic is the one \code{\link{generateCINoise}} performs on sinusoid
#' contrast weights, applied to latent perturbations instead: each trial's
#' perturbation is weighted by the response, and the weighted perturbations are
#' averaged. The result is a direction, which is added to the base latent to
#' give the classification latent, and rendered.
#'
#' \subsection{What the direction estimates}{
#'
#' Perturbations are drawn with a covariance set by the generator's
#' \code{latent_sd}, and a response-weighted mean under Gaussian sampling
#' estimates the participant's preference direction multiplied by that
#' covariance. The recovered direction is therefore pulled towards the
#' dimensions along which faces vary most, which is what keeps the rendered
#' result on the face manifold: dividing the covariance back out would give a
#' direction that fits the responses equally well and renders as a face nobody
#' has. Two participants' directions computed from the same generator are
#' weighted identically, so they remain comparable to each other.
#' }
#'
#' \subsection{Scaling}{
#'
#' A raw direction shrinks as \code{1 / sqrt(n_trials)}, so at any realistic
#' trial count rendering it unscaled gives a face barely distinguishable from
#' the base. \code{latent_scaling} decides what to do about that, and the choice
#' is the same trade-off the pixel pipeline's \code{scaling} argument makes.
#'
#' \code{'sd'}, the default, rescales the direction so its root-mean-square
#' displacement is \code{scaling_constant} standard deviations of the training
#' faces. Every classification image is then equally visible, and none is
#' comparable with another: the magnitude has been normalised away, so a
#' participant who answered at random and one who answered consistently produce
#' equally strong-looking faces. Use \code{\link{computeLatentInfoVal2IFC}} to
#' tell those apart, never the picture.
#'
#' \code{'constant'} divides the direction by \code{scaling_constant}, the same
#' number for every classification image, so magnitudes stay comparable across
#' participants and conditions. This is the option to use when the images are
#' shown side by side as evidence.
#'
#' \code{'none'} renders the direction as computed.
#' }
#'
#' @export
#' @param stimuli Vector of stimulus numbers, in the order they were presented.
#' @param responses Vector of responses, coded 1 when the original image was
#'   chosen and -1 when the inverted image was. The same coding
#'   \code{\link{generateCI}} takes.
#' @param rdata String pointing to the .Rdata file written by
#'   \code{\link{generateStimuliLatent2IFC}}.
#' @param targetpath String specifying the directory to save the classification
#'   image to. Required unless \code{save_as_png} is FALSE; there is no default
#'   path.
#' @param generator The generator that made the stimuli. Defaults to NULL, which
#'   rebuilds it from the .Rdata file. That works for a generator built by
#'   \code{\link{latentGeneratorPCA}}, which travels inside the file; a
#'   generator too large to store has to be passed back in here, and one whose
#'   fingerprint does not match the file is refused.
#' @param participants Vector of participant identifiers, one per trial.
#'   Defaults to NA, which pools every trial into one classification image. When
#'   given, a direction is computed per participant and those are averaged, so
#'   participants who ran more trials do not count for more.
#' @param latent_scaling String specifying how the direction is scaled before
#'   rendering: \code{sd} (default), \code{constant} or \code{none}. See
#'   Scaling above.
#' @param scaling_constant Number the scaling uses. Under \code{sd} it is a
#'   number of standard deviations (default 2); under \code{constant} it is the
#'   divisor.
#' @param antiCI Boolean specifying whether the anti-classification image is
#'   rendered instead: the same direction, negated.
#' @param save_as_png Boolean specifying whether to write the classification
#'   image to disk (default: TRUE).
#' @param filename Name for the written image, without extension. Defaults to
#'   the label recorded in the .Rdata file.
#' @return A list with the classification latent, the direction before and after
#'   scaling, the rendered classification image, the base latent and its render,
#'   the per-participant directions when \code{participants} was given, and
#'   fingerprints of the generator and of the stimulus set it came from.
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
#'   stimulus_path = file.path(tempdir(), "rcicr_latent_ci_example"),
#'   seed = 1, save_as_png = FALSE
#' )
#'
#' ci <- generateLatentCI(
#'   stimuli = 1:20,
#'   responses = rep(c(1, -1), 10),
#'   rdata = stimuli$rdata,
#'   targetpath = file.path(tempdir(), "rcicr_latent_ci_example")
#' )
#' dim(ci$ci_image)
generateLatentCI <- function(stimuli, responses, rdata, targetpath, generator = NULL, participants = NA, latent_scaling = 'sd', scaling_constant = 2, antiCI = FALSE, save_as_png = TRUE, filename = '') {

  requireExperimental('generateLatentCI')

  if (save_as_png && missing(targetpath)) {
    msg <- paste0(
      'No targetpath was given. Supply targetpath = <a directory> to say where ',
      'the classification image should go, or save_as_png = FALSE to compute ',
      'it without writing anything.'
    )
    stop(msg, call. = FALSE)
  }

  trials <- coerceTrialVectors(stimuli, responses, participants)

  stored <- loadLatentStimulusParams(rdata)
  generator <- resolveGenerator(generator, stored$generator_spec, 'generateLatentCI')

  if (all(is.na(trials$participants))) {
    aggregated <- aggregateResponses(trials$stimuli, trials$responses)
    trials$stimuli <- aggregated$stimuli
    trials$responses <- aggregated$responses
  }

  params <- selectLatentParams(stored$latent_params, trials$stimuli)

  pid_directions <- NULL
  if (all(is.na(trials$participants))) {
    direction <- weightedLatentMean(params, trials$responses)
  } else {
    pid_directions <- computeParticipantDirections(params, trials$responses, trials$participants)
    direction <- colMeans(pid_directions)
  }

  if (antiCI) {
    direction <- -direction
  }

  scaled <- applyLatentScaling(direction, latent_scaling, scaling_constant,
                               generator$latent_sd)

  latent_ci <- stored$base_latent + scaled
  ci_image <- renderUnchecked(generator, latent_ci, validate = FALSE)[1, , ]
  base_image <- renderUnchecked(generator, stored$base_latent, validate = FALSE)[1, , ]

  if (save_as_png) {
    label <- if (filename == '') latentLabel(rdata) else filename
    saveToImage(label, ci_image, targetpath, filename, antiCI)
  }

  return(list(
    latent_ci = latent_ci,
    direction = direction,
    scaled_direction = scaled,
    ci_image = ci_image,
    base_latent = stored$base_latent,
    base_image = base_image,
    pid_directions = pid_directions,
    generator_fingerprint = generator$fingerprint,
    stimulus_fingerprint = stimulusFingerprint(stored)
  ))
}

# The response-weighted mean, the same operation generateCINoise() performs on
# sinusoid contrast weights. A single trial arrives as a vector rather than a
# one-row matrix, as it does there.
weightedLatentMean <- function(params, responses) {
  weighted <- params * responses

  if (is.null(dim(weighted))) {
    return(weighted)
  }

  return(colMeans(weighted))
}

# One direction per participant, so a participant who ran twice the trials does
# not count twice in the average.
computeParticipantDirections <- function(params, responses, participants) {
  pids <- sort(unique(participants))
  directions <- matrix(0, nrow = length(pids), ncol = ncol(params),
                       dimnames = list(as.character(pids), NULL))

  for (i in seq_along(pids)) {
    rows <- participants == pids[i]
    directions[i, ] <- weightedLatentMean(params[rows, , drop = FALSE], responses[rows])
  }

  return(directions)
}

# Indexing by stimulus number works for non-consecutive stimuli too, matching
# selectStimulusParams() on the pixel-noise side.
selectLatentParams <- function(latent_params, stimuli) {
  if (any(stimuli < 1) || any(stimuli > nrow(latent_params))) {
    msg <- paste0(
      'stimuli must be numbers between 1 and ', nrow(latent_params),
      ', the number of trials in this stimulus file. The range given is ',
      min(stimuli), ' to ', max(stimuli), '.'
    )
    stop(msg, call. = FALSE)
  }

  return(latent_params[stimuli, , drop = FALSE])
}

applyLatentScaling <- function(direction, latent_scaling, scaling_constant, latent_sd) {
  if (latent_scaling == 'none') {
    return(direction)
  }

  if (latent_scaling == 'constant') {
    return(direction / scaling_constant)
  }

  if (latent_scaling == 'sd') {
    # The displacement measured in standard deviations of the training faces, so
    # scaling_constant reads as "move the face this many SDs".
    in_sd <- sqrt(mean((direction / latent_sd)^2))

    if (in_sd == 0) {
      return(direction)
    }

    return(direction * (scaling_constant / in_sd))
  }

  warning(paste0('Unknown latent_scaling: ', latent_scaling,
                 '. Returning the direction unscaled.'))

  return(direction)
}

# The generator that made the stimuli: the one the caller passed, checked
# against the file, or the one the file carries.
resolveGenerator <- function(generator, spec, caller) {
  if (!is.null(generator)) {
    validateGenerator(generator, probe = FALSE)
    matchGenerator(generator, spec, caller)
    return(generator)
  }

  rebuilt <- generatorFromSpec(spec)

  if (is.null(rebuilt)) {
    msg <- paste0(
      'These stimuli were made with a ', spec$kind, ' generator, which is too ',
      'large to store in the stimulus file, so ', caller, '() cannot rebuild ',
      'it. Pass the same generator as generator = <your generator>. Only a ',
      'generator built by latentGeneratorPCA() travels inside the file.'
    )
    stop(msg, call. = FALSE)
  }

  return(rebuilt)
}

# The label the stimuli were written under, recovered from the file name the
# .Rdata was saved as, so a classification image is named after its stimulus set
# without the label having to be passed in again.
latentLabel <- function(rdata) {
  return(sub('_seed_.*$', '', basename(rdata)))
}
