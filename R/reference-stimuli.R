# References over a subset of the saved stimuli (#349). Brinkman et al. (2019)
# require the reference to use the stimuli the CI was built from; the default
# reference uses every saved stimulus, so a CI from a subset needs its own.
# Matching the reference to the design is the researcher's call, so nothing
# here changes a default: these paths run only when reference_stimuli is given.

# NULL for the full saved set, so an explicit full set takes the unchanged
# default path; otherwise sorted distinct integers. IDs often arrive as doubles,
# and c(1, 2, 3) is not identical() to seq_len(3).
canonicalReferenceStimuli <- function(reference_stimuli, n_trials) {
  if (is.null(reference_stimuli)) return(NULL)
  ids <- tryCatch({
    ids <- coerceStimulusIds(reference_stimuli)
    validateStimulusIds(ids, n_trials)
    ids
  }, error = function(e) stop('reference_stimuli: ', conditionMessage(e), call. = FALSE))
  repeated <- unique(ids[duplicated(ids)])
  if (length(repeated) > 0L) {
    shown <- paste(utils::head(sort(repeated), 5), collapse = ', ')
    if (length(repeated) > 5) shown <- paste0(shown, ', ...')
    stop('reference_stimuli lists stimulus ', shown, ' more than once. A reference ',
         'assumes one response per stimulus from one responder, so it cannot describe ',
         'repeated presentations or several participants.', call. = FALSE)
  }
  ids <- sort(as.integer(ids))
  if (identical(ids, seq_len(n_trials))) NULL else ids
}

validTrialCount <- function(n_trials) {
  is.numeric(n_trials) && length(n_trials) == 1L && is.finite(n_trials) && n_trials >= 1 &&
    n_trials == trunc(n_trials)
}

# The file's contents in an environment of their own, whatever the file's
# base-image layout. selectReferenceBase() keeps one only for independent bases.
referenceSource <- function(rdata, selection) {
  if (selection$independent) return(selection$source)
  source <- new.env(parent = emptyenv())
  loadRdata(rdata, source)
  source
}

# NULL when the default reference applies: no reference_stimuli, or all of them.
subsetReferenceFor <- function(rdata, selection, reference_stimuli) {
  if (is.null(reference_stimuli)) return(NULL)
  source <- referenceSource(rdata, selection)
  if (!validTrialCount(source$n_trials)) stop('The stimulus file must contain a positive integer n_trials.')
  ids <- canonicalReferenceStimuli(reference_stimuli, source$n_trials)
  if (is.null(ids)) return(NULL)
  list(source = source, reference_stimuli = ids)
}

referenceLabel <- function(source, selection) {
  if (selection$independent) selection$baseimage else names(source$stimuli_params)[1]
}

subsetReferenceCache <- function(source) {
  cache <- source$reference_norms_by_stimuli
  if (!is.null(cache) && !is.list(cache)) stop('reference_norms_by_stimuli must be a list.')
  if (is.null(cache)) list() else cache
}

# Shared-parameter files key their entries on a NULL base: every base has the
# same noise, so one reference serves them all.
subsetReferenceIndex <- function(cache, reference_stimuli, base_key) {
  for (i in seq_along(cache)) {
    if (identical(cache[[i]]$reference_stimuli, reference_stimuli) &&
          identical(cache[[i]]$baseimage, base_key)) {
      return(i)
    }
  }
  NA_integer_
}

# The default reference's arithmetic over the selected stimuli. The stream
# still replays the generator's draws for every saved trial, so a subset
# reference is reproducible from the file alone.
generateSubsetReference <- function(source, rdata, selection, reference_stimuli, iter,
                                    ncores, response_seed, save_rdata, reference_method) {
  if (length(iter) != 1L || !is.finite(iter) || iter < 1 || iter != trunc(iter)) {
    stop('iter must be a positive integer.')
  }
  label <- referenceLabel(source, selection)
  base_key <- if (selection$independent) label else NULL
  if (is.null(response_seed)) requireStimulusSeed(source$seed, rdata, base_key, subset = TRUE)
  write('Building the reference from the saved noise, please wait...', stdout())
  if (iter < 10000) warning('You should set iter >= 10000 for InfoVal statistic to be reliable')
  write('Computing reference distribution, please wait...', stdout())
  norms <- referenceNorms(source, label, reference_stimuli, iter, ncores, response_seed,
                          reference_method)
  if (save_rdata) {
    cache <- subsetReferenceCache(source)
    entry <- list(
      reference_stimuli = reference_stimuli, baseimage = base_key, norms = norms,
      response_seed = response_seed, source = 'saved_noise',
      fingerprint = referenceSnapshot(norms), method = reference_method
    )
    i <- subsetReferenceIndex(cache, reference_stimuli, base_key)
    if (is.na(i)) i <- length(cache) + 1L
    cache[[i]] <- entry
    source$reference_norms_by_stimuli <- cache
    # Saved from the file's own environment, never from a caller's frame, so
    # the file keeps every field it had and gains only this one.
    saveRdataSafely(ls(source, all.names = TRUE), rdata, source)
  }
  norms
}

subsetReference <- function(rdata, iter, force_gen_ref_dist, response_seed, source, selection,
                            reference_stimuli, reference_method, report) {
  label <- referenceLabel(source, selection)
  base_key <- if (selection$independent) label else NULL
  report(source$n_trials, reference_stimuli)
  cache <- subsetReferenceCache(source)
  i <- subsetReferenceIndex(cache, reference_stimuli, base_key)
  entry <- if (is.na(i)) NULL else cache[[i]]
  norms <- resolveReferenceNorms(entry, rdata, iter, force_gen_ref_dist, response_seed,
                                 base_key, seedless = is.null(source$seed),
                                 reference_stimuli = reference_stimuli,
                                 reference_method = reference_method)
  if (!is.numeric(norms) || !length(norms) || any(!is.finite(norms))) {
    stop('Invalid cached reference for these reference_stimuli. ',
         'Use force_gen_ref_dist = TRUE to regenerate it.')
  }
  # Few stimuli give few distinct norms: one gives one, so the MAD is 0 and the
  # InfoVal would be infinite or NaN rather than a number.
  if (mad(norms) == 0) {
    stop('The reference over these ', length(reference_stimuli), ' stimuli has a MAD of 0, ',
         'so no InfoVal can be computed from it: too few stimuli for random responses to ',
         'give distinct norms.', call. = FALSE)
  }
  base_note <- if (is.null(base_key)) '' else paste0('baseimage = ', base_key, '; ')
  list(median = median(norms), mad = mad(norms), iter = length(norms),
       note = paste0(base_note, 'reference over ', length(reference_stimuli), ' of ',
                     source$n_trials, ' stimuli; '))
}

# A message, never a warning: whether the design and the reference match is
# the researcher's decision, and the number is returned either way. A CI
# without the attribute (an older version's, or one built by hand) has nothing
# to compare.
reportTrialDesign <- function(target_ci, n_trials, reference_stimuli) {
  issue <- trialDesignIssue(target_ci, n_trials, reference_stimuli)
  if (is.null(issue)) return(invisible(NULL))
  if (issue$averaged) {
    message('This classification image ', issue$built_from, '. ', uncalibratedDesignAdvice())
  } else {
    message('This classification image was built from ', issue$built_from, ' of the ',
            n_trials, ' saved stimuli, but the reference is built over ', issue$reference, '. ',
            mismatchedDesignAdvice('reference_stimuli = attr(<your CI>, "trial_design")$stimuli'))
  }
  invisible(NULL)
}

# NULL when the CI's recorded design matches its reference, or nothing is
# recorded to compare.
trialDesignIssue <- function(target_ci, n_trials, reference_stimuli) {
  design <- attr(target_ci, 'trial_design', exact = TRUE)
  if (is.null(design) || !validTrialCount(n_trials)) return(NULL)
  if (isTRUE(design$repeated) || isTRUE(design$n_participants > 1L)) {
    built_from <- if (isTRUE(design$n_participants > 1L)) {
      paste0('averages ', design$n_participants, ' participants')
    } else {
      'averages repeated presentations of the same stimuli'
    }
    return(list(averaged = TRUE, built_from = built_from))
  }
  used <- if (is.null(reference_stimuli)) seq_len(n_trials) else reference_stimuli
  if (identical(design$stimuli, used)) return(NULL)
  list(averaged = FALSE, built_from = length(design$stimuli), reference = length(used))
}

uncalibratedDesignAdvice <- function() {
  paste0('Every InfoVal reference in rcicr assumes one response per stimulus from one ',
         'responder, and none is defined for this design, so this InfoVal is not calibrated ',
         'for it. Where each participant saw every stimulus once, compute InfoVal for each ',
         'participant\'s own CI instead.')
}

mismatchedDesignAdvice <- function(call) {
  paste0('Brinkman et al. (2019) require the reference to use the stimuli the CI was ',
         'built from. To score it that way, pass ', call, '.')
}
