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
#'   \code{\link{generateCI}} takes, but checked here rather than assumed:
#'   anything else is an error. A 0/1 coding is what several experiment
#'   programs write and it raises no error further down, so it would otherwise
#'   turn the estimator into a quantity with no meaning and still render a
#'   plausible face.
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
#'   the per-participant directions when \code{participants} was given,
#'   fingerprints of the generator and of the stimulus set it came from, and the
#'   trials it was computed over, which \code{\link{computeLatentInfoVal2IFC}}
#'   requires in order to score it against a null built from the same data.
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

  stored <- loadLatentStimulusParams(rdata)
  generator <- resolveGenerator(generator, stored$generator_spec, 'generateLatentCI')

  # Length before coding: a caller who passed the wrong number of responses is
  # better told that than told their values are out of range.
  coerceTrialVectors(stimuli, responses, participants)
  checkParticipants(participants, stimuli, 'generateLatentCI')
  checkStimuliAreTrials(stimuli, 'generateLatentCI')
  checkResponseCoding(responses, 'generateLatentCI')

  estimate <- latentDirection(stored$latent_params, stimuli, responses, participants)
  direction <- estimate$direction
  pid_directions <- estimate$pid_directions

  # Both, so that the group direction stays the average of the participant
  # directions returned beside it. The pixel pipeline negates its parameters
  # before computing anything, which has the same effect there; here the
  # estimator is shared with the permutation null, which antiCI must not touch,
  # so the negation happens on the way out instead.
  if (antiCI) {
    direction <- -direction
    if (!is.null(pid_directions)) {
      pid_directions <- -pid_directions
    }
  }

  scaled <- applyLatentScaling(direction, latent_scaling, scaling_constant,
                               generator$latent_sd)

  latent_ci <- stored$base_latent + scaled
  ci_image <- renderUnchecked(generator, latent_ci, validate = FALSE)[1, , ]
  base_image <- renderUnchecked(generator, stored$base_latent, validate = FALSE)[1, , ]

  if (save_as_png) {
    label <- if (filename == '') latentOutputName(stored, rdata) else filename
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
    stimulus_fingerprint = stimulusFingerprint(stored),
    analysis_inputs = analysisInputs(stored, stimuli, responses, participants)
  ))
}

# Responses are 1 and -1 and nothing else, checked wherever a caller supplies
# them. A 0/1 coding is the one that matters: it is what plenty of experiment
# software writes, it raises no error anywhere downstream, and it turns the
# response-weighted mean into a quantity with no meaning while still rendering a
# perfectly plausible face. NA reaches the render and comes out as a blank
# image, which is at least visible, but not until much later.
#
# Stricter than the pixel-noise pipeline, which documents that responses can be
# a scale. That latitude is not carried over here deliberately: this module is
# new, so nothing depends on it, and a scale can be added as a decision rather
# than inherited as an accident.
# One participant identifier per trial, or the NA sentinel meaning "pool
# everything". coerceTrialVectors() checks stimuli against responses and not this,
# so a short vector reaches `participants == pid` and recycles: two identifiers
# for twenty trials silently deal the trials out alternately and give a group
# classification image over participants nobody ran.
#
# The check lives here rather than in coerceTrialVectors(), which the pixel
# pipeline shares. Tightening it there would turn data that generateCI() has
# always accepted into an error, which is the one thing this package does not do
# to an existing script; #294 tracks that side.
# Whole numbers, because a matrix subscript truncates rather than refuses:
# latent_params[1.5, ] is trial 1, so a fractional trial id gives a
# classification image built from a perturbation the caller did not name, with an
# informational value that agrees with it.
#
# Checked here rather than only at the subscript, because pooling aggregates the
# trials first and a non-finite id does not survive that intact.
checkStimuliAreTrials <- function(stimuli, caller) {
  values <- unlist(stimuli, use.names = FALSE)

  if (!is.numeric(values) || !all(is.finite(values)) || any(values != trunc(values))) {
    offending <- values[!is.finite(values) | values != trunc(values)]
    msg <- paste0(
      caller, '() takes whole trial numbers as stimuli, indexing the ',
      'perturbations in the stimulus file. These are not: ',
      paste(utils::head(unique(offending), 3), collapse = ', '), '.'
    )
    stop(msg, call. = FALSE)
  }

  return(invisible(TRUE))
}

checkParticipants <- function(participants, stimuli, caller) {
  values <- unlist(participants, use.names = FALSE)

  # All-NA is the sentinel meaning "pool every trial". A vector that is only
  # partly NA is not: computeParticipantDirections() takes its levels from
  # sort(unique(participants)), which drops NA, so those trials would be left
  # out of every per-participant direction and out of the group average of them,
  # while the length check above passed and nothing said so.
  if (all(is.na(values))) {
    return(invisible(TRUE))
  }

  # Length first, as with responses: a caller who supplied the wrong number of
  # identifiers is better told that than told some of them are missing.
  n_participants <- length(values)
  n_trials <- length(unlist(stimuli, use.names = FALSE))

  if (n_participants != n_trials) {
    msg <- paste0(
      caller, '() needs one participant identifier per trial. It was given ',
      n_participants, ' for ', n_trials, ' trials. Leave participants at its ',
      'default to pool every trial into one classification image.'
    )
    stop(msg, call. = FALSE)
  }

  if (anyNA(values)) {
    msg <- paste0(
      caller, '() was given participant identifiers for some trials and NA for ',
      sum(is.na(values)), ' of ', length(values), '. Those trials would be ',
      'dropped from every participant\'s direction without appearing anywhere ',
      'in the result. Name a participant for every trial, or leave participants ',
      'at its default to pool them all into one classification image.'
    )
    stop(msg, call. = FALSE)
  }

  return(invisible(TRUE))
}

checkResponseCoding <- function(responses, caller) {
  values <- unlist(responses, use.names = FALSE)

  if (anyNA(values) || !all(values %in% c(-1, 1))) {
    offending <- unique(values[is.na(values) | !(values %in% c(-1, 1))])
    msg <- paste0(
      caller, '() takes responses coded 1 (the original image was chosen) and ',
      '-1 (the inverted image was), one per trial. It was given ',
      paste(utils::head(offending, 3), collapse = ', '), '.'
    )
    stop(msg, call. = FALSE)
  }

  return(invisible(TRUE))
}

# The whole path from trials to a direction, in one place because
# computeLatentInfoVal2IFC() has to build its null through exactly this
# computation. Two implementations of it drifted apart once already: the null
# pooled where the classification image had averaged per participant, and it
# permuted responses that had already been collapsed over repeated stimuli, so
# it was scoring the observed direction against a differently-built one.
latentDirection <- function(latent_params, stimuli, responses, participants) {
  trials <- coerceTrialVectors(stimuli, responses, participants)

  # Pooling collapses repeated presentations of a stimulus to their mean
  # response, so this has to happen after any permutation rather than before it.
  if (all(is.na(trials$participants))) {
    aggregated <- aggregateResponses(trials$stimuli, trials$responses)
    params <- selectLatentParams(latent_params, aggregated$stimuli)

    return(list(direction = weightedLatentMean(params, aggregated$responses),
                pid_directions = NULL))
  }

  params <- selectLatentParams(latent_params, trials$stimuli)
  pid_directions <- computeParticipantDirections(params, trials$responses,
                                                 trials$participants)

  return(list(direction = colMeans(pid_directions),
              pid_directions = pid_directions))
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
  # Whole numbers, because a matrix subscript truncates rather than refuses:
  # latent_params[1.5, ] is trial 1, so a fractional trial id would give a
  # classification image built from a perturbation the caller did not name, and
  # the informational value would agree with it.
  if (!all(is.finite(stimuli)) || any(stimuli != trunc(stimuli))) {
    offending <- stimuli[!is.finite(stimuli) | stimuli != trunc(stimuli)]
    msg <- paste0(
      'stimuli must be whole trial numbers. These are not: ',
      paste(utils::head(unique(offending), 3), collapse = ', '), '.'
    )
    stop(msg, call. = FALSE)
  }

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

# A single positive finite number, checked where it is applied rather than at the
# entry, so latent_scaling = 'none' does not reject a constant it never reads.
#
# Zero divides the direction by zero under 'constant'. A negative turns the
# classification image around under either mode, which is what antiCI is for --
# doing it through the scaling constant leaves the number reading as a size while
# meaning a size and a reversal, and nothing downstream can tell.
checkScalingConstant <- function(scaling_constant, latent_scaling) {
  if (length(scaling_constant) != 1 || !is.finite(scaling_constant) ||
        scaling_constant <= 0) {
    msg <- paste0(
      'scaling_constant must be a single positive number when latent_scaling ',
      'is "', latent_scaling, '". It is ',
      paste(format(scaling_constant), collapse = ', '),
      '. To turn the classification image around, use antiCI = TRUE.'
    )
    stop(msg, call. = FALSE)
  }

  return(invisible(TRUE))
}

applyLatentScaling <- function(direction, latent_scaling, scaling_constant, latent_sd) {
  if (latent_scaling == 'none') {
    return(direction)
  }

  if (latent_scaling == 'constant') {
    checkScalingConstant(scaling_constant, latent_scaling)

    return(direction / scaling_constant)
  }

  if (latent_scaling == 'sd') {
    checkScalingConstant(scaling_constant, latent_scaling)

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

# The label the stimuli were written under, so a classification image is named
# after its stimulus set without the label having to be passed in again.
#
# Taken from the file's own `label` field. Recovering it from the file name
# instead means splitting on '_seed_', which a label may itself contain --
# "my_seed_test" came back as "my" -- so the image was named after part of its
# set, and two sets whose labels shared a first segment wrote the same file.
# Falling back to the file name for a stimulus file written before `label` was
# stored, which is the append-only contract working as intended.
# The default output name for a classification image. The label alone does not
# identify a stimulus set -- two sets made at different seeds share it, and the
# default label means most sets share it -- so analysing both into one targetpath
# wrote ci_rcic_latent.png twice and the second replaced the first.
#
# The seed completes the identity: two sets at the same label and seed have the
# same perturbations, and generateStimuliLatent2IFC() already refuses to write
# them into one directory.
latentOutputName <- function(stored, rdata) {
  label <- latentLabel(stored, rdata)

  if (is.null(stored$seed) || length(stored$seed) != 1 || is.na(stored$seed)) {
    return(label)
  }

  return(paste(label, stored$seed, sep = '_'))
}

latentLabel <- function(stored, rdata) {
  if (!is.null(stored$label) && !is.na(stored$label) && nzchar(stored$label)) {
    return(stored$label)
  }

  return(sub('_seed_.*$', '', basename(rdata)))
}
