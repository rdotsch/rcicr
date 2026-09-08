# Preparing generateCI()'s inputs: the steps between the arguments the caller
# passed and the parameter matrix the CI is actually computed from. All pure --
# no I/O, no reliance on anything load() put in a frame -- so each can be tested
# on values handed to it directly.

# Coerce the trial vectors to plain vectors and check they describe the same
# trials.
#
# Data read with readr or manipulated with dplyr comes back as a tibble, where
# tbl[, "col"] stays a one-column tibble rather than dropping to a vector the
# way df[, "col"] does. That made aggregate() fail downstream with the opaque
# message "arguments must have same length".
#
# Input: stimuli, responses, participants (NA when not used)
# Output: list of the three, coerced
coerceTrialVectors <- function(stimuli, responses, participants) {
  stimuli <- coerceStimulusIds(stimuli)
  responses <- unlist(responses, use.names = FALSE)
  if (!all(is.na(participants))) {
    participants <- unlist(participants, use.names = FALSE)
  }

  if (length(stimuli) != length(responses)) {
    stop(paste0('stimuli and responses must have the same length (stimuli: ',
      length(stimuli), ', responses: ', length(responses), ').'
    ))
  }

  if (!all(is.na(participants)) && length(participants) != length(stimuli)) {
    stop(paste0('participants must have one ID per trial (participants: ',
      length(participants), ', trials: ', length(stimuli), ').'
    ))
  }

  validateStimulusIds(stimuli)

  return(list(stimuli = stimuli, responses = responses,
    participants = participants
  ))
}

# Check before unlist() can turn factor codes or logicals into numeric indices.
coerceStimulusIds <- function(stimuli) {
  if (is.list(stimuli)) {
    return(unlist(lapply(stimuli, coerceStimulusIds), use.names = FALSE))
  }
  if (is.factor(stimuli) || !typeof(stimuli) %in% c('integer', 'double')) {
    stop('stimuli must contain numeric stimulus numbers, not factors, characters or logicals.')
  }
  return(unlist(stimuli, use.names = FALSE))
}

# Used before aggregation and again with the selected base's saved trial count.
validateStimulusIds <- function(stimuli, n_trials = Inf) {
  if (is.factor(stimuli) || !typeof(stimuli) %in% c('integer', 'double') ||
        length(stimuli) == 0L || any(!is.finite(stimuli)) ||
        any(stimuli < 1 | stimuli != floor(stimuli))) {
    stop('stimuli must contain at least one finite, positive, whole-number stimulus ID.')
  }
  if (any(stimuli > n_trials)) {
    stop(paste0('stimuli contains an ID beyond the ', n_trials,
      ' trials saved for the selected base image.'
    ))
  }
  return(invisible(NULL))
}

# Look up the base image the CI is built on.
#
# Input: base_faces list from the .Rdata, base image label
# Output: the base image matrix
selectBaseImage <- function(base_faces, baseimage) {
  base <- base_faces[[baseimage]]
  if (is.null(base)) {
    stop(paste0('File specified in rdata argument did not contain any ',
      'reference to base image label: ', baseimage, ' (NOTE: file ',
      'contains references to the following base image label(s): ',
      paste(names(base_faces), collapse = ', '), ')'
    ))
  }

  return(base)
}

# Collapse repeated presentations of the same stimulus, so each unique stimulus
# gets equal weight regardless of how many times it was presented. Used only
# when `participants` is NA; the per-participant path weights within
# participants instead. generateCI()'s roxygen documents what this changes for
# unbalanced designs.
#
# Input: stimuli, responses
# Output: list of both, one entry per unique stimulus
aggregateResponses <- function(stimuli, responses) {
  aggregated <- aggregate(responses, by = list(stimuli = stimuli), FUN = mean)

  return(list(stimuli = aggregated$stimuli, responses = aggregated$x))
}

# Retrieve the noise parameters of the stimuli that were actually presented.
# Indexing by the stimulus numbers works for non-consecutive stimuli too.
#
# Input: stimuli_params list from the .Rdata, base image label, stimulus numbers
# Output: parameter matrix (or vector, for a single trial)
selectStimulusParams <- function(stimuli_params, baseimage, stimuli) {
  if (length(stimuli_params[[baseimage]]) == 0L) {
    stop(paste0('No parameters found for base image: ', baseimage))
  }
  validateStimulusIds(stimuli, nrow(stimuli_params[[baseimage]]))
  params <- stimuli_params[[baseimage]][stimuli, ]

  if (length(params) == 0) {
    stop(paste0('No parameters found for base image: ', baseimage))
  }

  # Files written before rcicr 0.3.0 hold 4096 parameters per trial where only
  # 4092 patch indices exist: 6 orientations x 2 phases x sum(4^0..4^4) is 4092,
  # while pre-0.3.0 drew a round 4096 and never referred to the last 4. See
  # ChangeLog, 0.3.0-29.
  if (!is.vector(params)) {
    if (ncol(params) == 4096) {
      params <- params[, 1:4092]
    }
  } else {
    # In case we only have a single trial as input. This tested
    # `length(params) == 4092` and then truncated to 4092 -- a no-op that could
    # never fire on the 4096-parameter input it exists for, so a single-trial CI
    # from a pre-0.3.0 file died in generateNoiseImage() with "number of
    # parameters doesn't equal number of patches". The multi-trial branch above
    # was always correct.
    if (length(params) == 4096) {
      params <- params[1:4092]
    }
  }

  return(params)
}
