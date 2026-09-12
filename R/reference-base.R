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

# Retain every value: identical() needs no hash dependency, formatting, or arithmetic.
referenceSnapshot <- function(norms) {
  list(norms = norms)
}

vouchedReference <- function(norms, marker, fingerprint) {
  identical(marker, 'saved_noise') &&
    identical(fingerprint, referenceSnapshot(norms))
}

# Both cache layouts use this policy; their generators own the storage details.
resolveReferenceNorms <- function(entry, rdata, iter, force_gen_ref_dist,
                                  response_seed, baseimage = NULL) {
  forced <- force_gen_ref_dist || !is.null(response_seed)
  stale <- !is.null(entry) && is.null(entry$response_seed) &&
    !vouchedReference(entry$norms, entry$source, entry$fingerprint)
  if (!forced && !stale && !is.null(entry)) {
    write('Using reference distribution found in rdata file.', stdout())
    return(entry$norms)
  }

  # An unsolicited refresh must preserve the precision and RNG state of a cache hit.
  automatic <- stale && !forced
  if (automatic) iter <- length(entry$norms)
  readonly <- automatic && !writableFile(rdata)
  save_rdata <- is.null(response_seed) && !readonly
  simulate <- function() {
    withCallingHandlers(
      generateReferenceDistribution2IFC(
        rdata, iter = iter, response_seed = response_seed,
        save_rdata = save_rdata, baseimage = baseimage
      ),
      warning = function(cond) {
        # The inherited count is not a new choice the caller can act on.
        if (automatic && iter < 10000 &&
              grepl('iter >= 10000', conditionMessage(cond), fixed = TRUE)) {
          invokeRestart('muffleWarning')
        }
      }
    )
  }
  norms <- if (automatic) preserveRandomStream(simulate()) else simulate()

  if (readonly) {
    label <- if (is.null(baseimage)) 'the reference' else paste0('the reference for baseimage ', baseimage)
    write(paste0('Rebuilt ', label, ' from saved noise, but ', rdata,
                 ' is not writable, so the rebuilt values were used without being stored. ',
                 'The next call will rebuild them again.'), stdout())
  } else if (save_rdata) {
    write('The reference distribution has been saved to the .Rdata file for reuse.', stdout())
  } else {
    write(paste0('Reference distribution simulated with response_seed = ', response_seed,
                 '. This independent draw has deliberately not been saved.'), stdout())
  }
  if (stale && is.null(response_seed) && !identical(norms, entry$norms)) {
    message('This stimulus file carried a reference distribution ',
            'from before rcicr built references from the saved noise, and ',
            'rebuilding it from the saved noise gave different values. The ',
            'InfoVal this returns supersedes any computed from this file before.')
  }
  norms
}

# A cache hit consumed no random numbers before automatic migration existed.
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

# Replay one normalized parameter matrix's draws to preserve the historical response stream.
# selectStimulusParams() removes the four unused columns of pre-0.3.0 files.
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
      source = 'saved_noise', fingerprint = referenceSnapshot(norms)
    )
    source$reference_norms_by_base <- cache
    save(list = ls(source, all.names = TRUE), file = rdata, envir = source)
  }
  invisible(norms)
}

computeBaseInfoVal <- function(target_ci, rdata, iter, force_gen_ref_dist, response_seed, selection) {
  cache <- selection$source$reference_norms_by_base
  if (!is.null(cache) && !is.list(cache)) stop('reference_norms_by_base must be a list.')
  norms <- resolveReferenceNorms(cache[[selection$baseimage]], rdata, iter,
                                 force_gen_ref_dist, response_seed, selection$baseimage)
  if (!is.numeric(norms) || !length(norms) || any(!is.finite(norms))) {
    stop('Invalid cached reference for baseimage ', selection$baseimage,
         '. Use force_gen_ref_dist = TRUE to regenerate it.')
  }
  cinorm <- norm(matrix(target_ci[['ci']]), 'f')
  info_val <- (cinorm - median(norms)) / mad(norms)
  write(paste0('Informational value: z = ', info_val, ' (baseimage = ', selection$baseimage,
               '; ci norm = ', cinorm, '; reference median = ', median(norms),
               '; MAD = ', mad(norms), '; iterations = ', length(norms), ')'), stdout())
  return(info_val)
}

# Named so tests can model read-only archives even when running as root.
writableFile <- function(path) {
  unname(file.access(path, mode = 2)) == 0L
}
