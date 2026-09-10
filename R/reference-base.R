selectReferenceBase <- function(rdata, baseimage) {
  source <- new.env(parent = emptyenv())
  load(rdata, envir = source)
  params <- source$stimuli_params
  labels <- names(params)
  if (!is.null(baseimage) && (!is.character(baseimage) || length(baseimage) != 1L ||
                                is.na(baseimage) || !baseimage %in% labels)) {
    stop('baseimage must be one saved base label: ', paste(labels, collapse = ', '), '.')
  }
  independent <- length(params) > 1L &&
    !all(vapply(params, function(x) identical(unname(x), unname(params[[1]])), logical(1)))
  if (!independent) return(list(independent = FALSE))
  if (is.null(baseimage)) {
    stop('Saved base images have different noise parameters. Supply baseimage = <label> ',
         'for the CI being scored. Available labels: ', paste(labels, collapse = ', '), '.')
  }
  list(independent = TRUE, baseimage = baseimage, source = source)
}

# The saved parameter matrix for one base image, normalized: selectStimulusParams()
# drops the four columns a pre-0.3.0 file holds and never indexes, so its width is
# also the number of draws the generator spent on each trial.
# Ties a `reference_norms_source` marker to the norms it describes.
#
# The marker alone can outlive its data: an older rcicr re-saves every object it
# loaded, so a release that predates the marker will preserve it while replacing
# `reference_norms` with a rebuilt-basis distribution of its own. The file then
# looks vouched-for and is not. A fingerprint that travels with the marker turns
# that into a mismatch, and the norms are refreshed as if unmarked.
#
# Values drawn from the vector rather than computed from it, and numeric rather
# than text. `format()` honours `OutDec` and `scipen`; `sum()` accumulates in
# long double where the platform has one, so it can round differently for the
# same bytes on a build without it. Either would make a cache fail to match
# itself and be rebuilt on every call. Selection and `identical()` involve no
# arithmetic and no formatting, so the fingerprint travels with the file.
referenceFingerprint <- function(norms) {
  n <- length(norms)
  c(n, norms[c(1L, (n + 1L) %/% 2L, n)])
}

# A cache is trusted only on positive evidence that it describes this file's
# saved noise: the marker, and a fingerprint bound to the values it marks.
vouchedReference <- function(norms, marker, fingerprint) {
  identical(marker, 'saved_noise') &&
    identical(fingerprint, referenceFingerprint(norms))
}

# Said by both storage shapes, so they cannot drift apart.
#
# It reports a migration the caller did not ask for, so it must not be the thing
# that fails their call: under options(warn = 2) a warning is an error, and the
# InfoVal -- already computed, and the corrected one -- would never be returned.
# On a read-only archive there is no cache to write either, so every later call
# would rebuild and abort again and the file could never be scored. Degrading to
# a message under that setting keeps the value reachable and still says what
# happened; ordinary warning behaviour is unchanged.
warnReferenceSuperseded <- function() {
  msg <- paste0('This stimulus file carried a reference distribution ',
                'from before rcicr built references from the saved noise, and ',
                'rebuilding it from the saved noise gave different values. The ',
                'InfoVal this returns supersedes any computed from this file before.')
  tryCatch(warning(msg, call. = FALSE), error = function(e) message(msg))
}

# An automatic refresh is not something the caller asked for, so it must not
# move their random stream. Before the refresh existed, a cached file was a
# plain cache hit that consumed nothing; regenerating seeds and draws, which
# would silently change every later sample() in an old analysis script. An
# explicit regeneration keeps its documented stream behaviour.
preserveRandomStream <- function(expr) {
  had <- exists('.Random.seed', envir = globalenv(), inherits = FALSE)
  before <- if (had) get('.Random.seed', envir = globalenv(), inherits = FALSE)
  on.exit({
    if (had) {
      # nolint next: object_name_linter. R owns this name, not this package.
      assign('.Random.seed', before, envir = globalenv())
    } else if (exists('.Random.seed', envir = globalenv(), inherits = FALSE)) {
      rm('.Random.seed', envir = globalenv())
    }
  }, add = TRUE)
  expr
}

savedReferenceParams <- function(source, baseimage) {
  n_trials <- source$n_trials
  if (length(n_trials) != 1L || !is.finite(n_trials) || n_trials < 1 || n_trials != trunc(n_trials)) {
    stop('The stimulus file must contain a positive integer n_trials.')
  }
  params <- selectStimulusParams(source$stimuli_params, baseimage, seq_len(n_trials))
  matrix(params, nrow = n_trials)
}

# The reference distribution is a function of the saved noise alone, so it is
# built from the stored parameters and basis. Re-generating the stimuli instead
# reopened every base image, which an archived or moved experiment no longer has
# (#301), and rebuilt the basis from fields older files do not carry.
referenceNoise <- function(source, baseimage, ncores) {
  n_trials <- source$n_trials
  params <- savedReferenceParams(source, baseimage)
  p <- if (exists('s', envir = source, inherits = FALSE)) source$s else source$p
  if (is.null(p)) stop('The stimulus file does not contain its saved noise basis (p or s).')
  pb <- txtProgressBar(min = 0, max = n_trials, style = 3)
  on.exit(close(pb), add = TRUE)
  cl <- startBackend(ncores)
  on.exit(stopClusterSafely(cl), add = TRUE)
  noise <- foreach::foreach(trial = seq_len(n_trials), .combine = 'cbind',
    .packages = 'rcicr', .options.snow = progressOption(pb, cl)
  ) %dopar% {
    if (is.null(cl)) setTxtProgressBar(pb, trial)
    as.vector(generateNoiseImage(params[trial, ], p))
  }
  matrix(noise, ncol = n_trials)
}

# Put the random stream where the simulated responses expect to find it.
#
# generateStimuli2IFC() seeds on the stimulus seed and then spends one draw per
# parameter per trial, and the reference's responses have always been drawn from
# whatever that left behind -- which is what makes an InfoVal reproducible from
# the stimulus file under a fixed RNGkind(), and is documented as a guarantee on
# ?generateReferenceDistribution2IFC. Since the stimuli
# are no longer re-generated, that consumption is replayed here instead. The
# width of the saved matrix is the count, so nothing has to be assumed about a
# file that does not record its nscales.
seedResponseStream <- function(source, baseimage, response_seed) {
  if (!is.null(response_seed)) {
    set.seed(response_seed)
    return(invisible(NULL))
  }
  nparams <- ncol(savedReferenceParams(source, baseimage))
  set.seed(source$seed)
  for (trial in seq_len(source$n_trials)) runif(nparams)
  invisible(NULL)
}

generateBaseReference <- function(selection, rdata, iter, ncores, response_seed, save_rdata) {
  source <- selection$source
  if (length(iter) != 1L || !is.finite(iter) || iter < 1 || iter != trunc(iter)) {
    stop('iter must be a positive integer.')
  }
  stimuli <- referenceNoise(source, selection$baseimage, ncores)
  seedResponseStream(source, selection$baseimage, response_seed)
  if (iter < 10000) warning('You should set iter >= 10000 for InfoVal statistic to be reliable')
  pb <- txtProgressBar(min = 0, max = iter, style = 3)
  on.exit(close(pb), add = TRUE)
  norms <- numeric(iter)
  for (i in seq_len(iter)) {
    responses <- ((runif(source$n_trials) > 0.5) * 2) - 1
    ci <- (as.matrix(stimuli) %*% as.matrix(responses)) / ncol(stimuli)
    norms[i] <- norm(ci, 'f')
    setTxtProgressBar(pb, i)
  }
  if (save_rdata) {
    cache <- source$reference_norms_by_base
    if (is.null(cache)) cache <- list()
    if (!is.list(cache)) stop('reference_norms_by_base must be a list.')
    cache[[selection$baseimage]] <- list(
      norms = norms, response_seed = response_seed,
      source = 'saved_noise', fingerprint = referenceFingerprint(norms)
    )
    source$reference_norms_by_base <- cache
    save(list = ls(source, all.names = TRUE), file = rdata, envir = source)
  }
  invisible(norms)
}

computeBaseInfoVal <- function(target_ci, rdata, iter, force_gen_ref_dist, response_seed, selection) {
  cache <- selection$source$reference_norms_by_base
  if (!is.null(cache) && !is.list(cache)) stop('reference_norms_by_base must be a list.')
  entry <- cache[[selection$baseimage]]

  # Same rule as the shared path, for the same reason: an entry written before
  # references came from the saved noise cannot be shown to describe it, so it
  # is refreshed once rather than certified. An entry recording a response_seed
  # is a null someone asked for and is left alone.
  previous <- entry$norms
  forced_by_caller <- force_gen_ref_dist || !is.null(response_seed)
  stale <- !is.null(entry) && is.null(entry$response_seed) &&
    !vouchedReference(previous, entry$source, entry$fingerprint)
  inherited_iter <- FALSE
  if (stale && !forced_by_caller) {
    iter <- length(previous)
    inherited_iter <- iter < 10000
  }
  readonly_refresh <- stale && !forced_by_caller && !writableFile(rdata)

  if (forced_by_caller || stale || is.null(entry)) {
    simulate <- function() {
      withCallingHandlers(
        generateReferenceDistribution2IFC(
          rdata, iter = iter, response_seed = response_seed,
          save_rdata = is.null(response_seed) && !readonly_refresh,
          baseimage = selection$baseimage
        ),
        warning = function(cond) {
          if (inherited_iter && grepl('iter >= 10000', conditionMessage(cond), fixed = TRUE)) {
            invokeRestart('muffleWarning')
          }
        }
      )
    }
    norms <- if (stale && !forced_by_caller) preserveRandomStream(simulate()) else simulate()
    if (readonly_refresh) {
      note <- paste0('Rebuilt this base image\'s reference distribution from its saved noise, ',
                     'but ', rdata, ' is not writable, so the rebuilt values were used without ',
                     'being stored. The next call will rebuild them again.')
      write(note, stdout())
    }
    if (stale && is.null(response_seed) && !identical(norms, previous)) {
      warnReferenceSuperseded()
    }
  } else {
    norms <- previous
    if (!is.numeric(norms) || !length(norms) || any(!is.finite(norms))) {
      stop('Invalid cached reference for baseimage ', selection$baseimage,
           '. Use force_gen_ref_dist = TRUE to regenerate it.')
    }
  }
  cinorm <- norm(matrix(target_ci[['ci']]), 'f')
  info_val <- (cinorm - median(norms)) / mad(norms)
  write(paste0('Informational value: z = ', info_val, ' (baseimage = ', selection$baseimage,
               '; ci norm = ', cinorm, '; reference median = ', median(norms),
               '; MAD = ', mad(norms), '; iterations = ', length(norms), ')'), stdout())
  return(info_val)
}

# A seam, not a convenience: an automatic refresh has to know whether it may
# write before it starts, and a test needs to deny that without depending on
# permission bits the test's own user may outrank.
writableFile <- function(path) {
  unname(file.access(path, mode = 2)) == 0L
}
