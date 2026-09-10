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
# the stimulus file alone, and is documented as a guarantee. Since the stimuli
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
    cache[[selection$baseimage]] <- list(norms = norms, response_seed = response_seed)
    source$reference_norms_by_base <- cache
    save(list = ls(source, all.names = TRUE), file = rdata, envir = source)
  }
  invisible(norms)
}

computeBaseInfoVal <- function(target_ci, rdata, iter, force_gen_ref_dist, response_seed, selection) {
  cache <- selection$source$reference_norms_by_base
  if (!is.null(cache) && !is.list(cache)) stop('reference_norms_by_base must be a list.')
  entry <- cache[[selection$baseimage]]
  if (!is.null(response_seed) || force_gen_ref_dist || is.null(entry)) {
    norms <- generateReferenceDistribution2IFC(rdata, iter = iter,
                                               response_seed = response_seed, save_rdata = is.null(response_seed),
                                               baseimage = selection$baseimage)
  } else {
    norms <- entry$norms
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
