selectReferenceBase <- function(rdata, baseimage) {
  source <- new.env(parent = emptyenv())
  loadRdata(rdata, source)
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
                                  response_seed, baseimage = NULL, seedless = FALSE,
                                  reference_stimuli = NULL, reference_method = 'gram') {
  forced <- force_gen_ref_dist || !is.null(response_seed)
  stale <- !is.null(entry) && is.null(entry$response_seed) &&
    !vouchedReference(entry$norms, entry$source, entry$fingerprint)
  if (!forced && !stale && !is.null(entry)) {
    write('Using reference distribution found in rdata file.', stdout())
    # A message, not a warning: under warn = 2 the notice must not cost the
    # caller the value (see "A cached reference is trusted only on positive
    # evidence" in DECISIONS.md).
    if (seedless && is.null(entry$response_seed)) {
      message(rdata, ' has no stimulus seed, so its stored reference distribution cannot ',
              'be regenerated from it. It is used as stored. ',
              seededReferenceAdvice(rdata, baseimage, !is.null(reference_stimuli)))
    }
    return(entry$norms)
  }

  # An unsolicited refresh must preserve the precision and RNG state of a cache hit.
  automatic <- stale && !forced
  if (automatic) iter <- length(entry$norms)
  # Only InfoVal scoring reaches here, where the cache is an optimization: the
  # norms are simulated before the save is attempted, so an archive that cannot
  # be written must not cost the caller the number it already has. A direct
  # generateReferenceDistribution2IFC(save_rdata = TRUE) still errors, because
  # there the save is the request. A seeded draw is never stored, so it is not
  # a save that failed and writability does not describe it.
  wanted_save <- is.null(response_seed)
  readonly <- wanted_save && !writableFile(rdata)
  save_rdata <- wanted_save && !readonly
  simulate <- function() {
    withCallingHandlers(
      generateReferenceDistribution2IFC(
        rdata, iter = iter, response_seed = response_seed,
        save_rdata = save_rdata, baseimage = baseimage, reference_stimuli = reference_stimuli,
        reference_method = reference_method
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
    write(paste0('Built ', label, ' from saved noise, but ', rdata,
                 ' is not writable, so the values were used without being stored. ',
                 'The next call will build them again.'), stdout())
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
referenceNoise <- function(source, baseimage, ncores, reference_stimuli = NULL) {
  params <- savedReferenceParams(source, baseimage)
  if (!is.null(reference_stimuli)) params <- params[reference_stimuli, , drop = FALSE]
  n_trials <- nrow(params)
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

# The norms of `iter` CIs built from random responses over the selected saved
# stimuli, the responses continuing the stream seedResponseStream() sets up.
# Nothing before the first response draw consumes random numbers.
#
# "gram" uses the stimulus Gram matrix G: norm(S r / n) = sqrt(t(r) G r) / n
# needs no pixels-sized product per draw. "images" is the calculation of
# rcicr 1.5.0 and earlier, bit for bit. They differ by rounding only
# (analyses/gram-reference-accuracy.md, #354).
referenceNorms <- function(source, label, reference_stimuli, iter, ncores, response_seed,
                           reference_method) {
  reportReferenceMethod(reference_method)
  if (reference_method == 'gram') {
    params <- savedReferenceParams(source, label)
    if (!is.null(reference_stimuli)) params <- params[reference_stimuli, , drop = FALSE]
    gram <- stimulusGram(params, referenceBasis(source))
    rm(params)
    seedResponseStream(source, label, response_seed)
    return(gramNorms(gram, iter))
  }
  noise <- referenceNoise(source, label, ncores, reference_stimuli)
  seedResponseStream(source, label, response_seed)
  renderedNorms(noise, iter)
}

# A message, never a warning: it changes no number. Only computing a reference
# reaches here, so a stored one is reused without it.
reportReferenceMethod <- function(reference_method) {
  if (reference_method == 'gram') {
    message('InfoVal reference computed with reference_method = "gram". References from rcicr ',
            '1.5.0 and earlier used "images"; pass reference_method = "images" to reproduce them ',
            'bit for bit.')
  }
  invisible(NULL)
}

referenceBasis <- function(source) {
  p <- if (exists('s', envir = source, inherits = FALSE)) source$s else source$p
  if (is.null(p)) stop('The stimulus file does not contain its saved noise basis (p or s).')
  if ('sinusoids' %in% names(p)) p <- list(patches = p$sinusoids, patchIdx = p$sinIdx)
  p
}

# The arithmetic of every reference before #354, kept for the route that uses it.
renderedNorms <- function(noise, iter) {
  pb <- txtProgressBar(min = 0, max = iter, style = 3)
  on.exit(close(pb), add = TRUE)
  norms <- numeric(iter)
  for (i in seq_len(iter)) {
    responses <- ((runif(ncol(noise)) > 0.5) * 2) - 1
    ci <- (noise %*% as.matrix(responses)) / ncol(noise)
    norms[i] <- norm(ci, 'f')
    setTxtProgressBar(pb, i)
  }
  norms
}

# G = t(S) S for the noise S = P t(X), with P the basis generateNoiseImage()
# averages: one patch per pixel per layer, over the number of layers. Built
# from whichever is smaller, the rendered S or the parameters' cross-product
# t(P) P.
#
# Both are sums over pixels, so P is built a block of image columns at a time
# and the sum accumulated: no full-size copy of the basis is held, nor, for the
# rendered build, the full S. Each block's t(P) P is sparse (about 1 in 1,000
# entries at 512px), so the cross-product build stays sparse until G.
stimulusGram <- function(params, p) {
  d <- dim(p$patches)
  npix <- d[1] * d[2]
  nparams <- max(p$patchIdx)
  if (ncol(params) != nparams) {
    # generateNoiseImage()'s check, once rather than once per trial.
    if (ncol(params) == nparams + 1 && min(p$patchIdx) == 0) {
      warning('Rdata patch indices start at 0, whereas parameters are used from position 1. Due to this mismatch, one sinusoid will not be shown in resulting CI.')
      params <- params[, seq_len(nparams), drop = FALSE]
    } else {
      stop("Stimulus generation aborted: number of parameters doesn't equal number of patches!")
    }
  }
  from_noise <- gramFromRenderedNoise(npix, nrow(params), nparams)
  # At most 16,384 pixels of basis per block, and for the rendered build about
  # 32 MB of noise. Summing blocks holds G two more times, so when G is about
  # as large as the noise itself, one block costs less.
  block_pixels <- if (from_noise) min(2^14, 2^22 %/% nrow(params)) else 2^14
  if (from_noise && npix <= nrow(params) + block_pixels) block_pixels <- npix
  columns <- max(1L, block_pixels %/% d[1])
  total <- NULL
  for (first in seq(1L, d[2], by = columns)) {
    block <- first:min(d[2], first + columns - 1L)
    index <- p$patchIdx[, block, , drop = FALSE]
    # A 0 index is a cell of the layer a 0-based file never wrote; its patch
    # value is 0, so dropping it drops nothing (DECISIONS.md, "4096 -> 4092").
    keep <- index != 0
    basis <- Matrix::sparseMatrix(i = rep(seq_len(d[1] * length(block)), d[3])[keep],
                                  j = index[keep],
                                  x = p$patches[, block, , drop = FALSE][keep] / d[3],
                                  dims = c(d[1] * length(block), nparams))
    part <- if (from_noise) {
      tcrossprod(as.matrix(Matrix::tcrossprod(params, basis)))
    } else {
      Matrix::crossprod(basis)
    }
    total <- if (is.null(total)) part else total + part
  }
  if (from_noise) total else as.matrix(params %*% total %*% t(params))
}

# The rendered noise is pixels x trials; the parameters' cross-product is
# parameters x parameters. Either gives G; build it from the smaller.
gramFromRenderedNoise <- function(npix, n_trials, nparams) {
  npix * n_trials <= nparams^2
}

# One runif(n * k) consumes the stream exactly as k successive runif(n) calls.
# runif(), not rbinom(), in both routes: see DECISIONS.md, "purrr::rbernoulli()
# was replaced with runif(), not rbinom()".
gramNorms <- function(gram, iter, block = 1000L) {
  n <- nrow(gram)
  pb <- txtProgressBar(min = 0, max = iter, style = 3)
  on.exit(close(pb), add = TRUE)
  norms <- numeric(iter)
  done <- 0L
  while (done < iter) {
    k <- min(block, iter - done)
    responses <- matrix(((runif(n * k) > 0.5) * 2) - 1, nrow = n, ncol = k)
    norms[done + seq_len(k)] <- sqrt(pmax(colSums(responses * (gram %*% responses)), 0)) / n
    done <- done + k
    setTxtProgressBar(pb, done)
  }
  norms
}

# A missing or NULL stimulus seed leaves no stream to replay: set.seed(NULL)
# reseeds from the clock, so the default reference could never be reproduced
# (#334). Checked before referenceNoise(), the slow step.
requireStimulusSeed <- function(seed, rdata, baseimage = NULL, subset = FALSE) {
  if (is.null(seed)) {
    stop(rdata, ' has no stimulus seed to replay, so its default reference distribution ',
         'could not be reproduced. ', seededReferenceAdvice(rdata, baseimage, subset), call. = FALSE)
  }
  invisible(NULL)
}

seededReferenceAdvice <- function(rdata, baseimage = NULL, subset = FALSE) {
  base_arg <- if (is.null(baseimage)) '' else paste0(', baseimage = ', encodeString(baseimage, quote = '"'))
  if (subset) base_arg <- paste0(base_arg, ', reference_stimuli = <the same stimuli>')
  paste0('Store a reproducible one with generateReferenceDistribution2IFC(',
         encodeString(rdata, quote = '"'), ', response_seed = <n>', base_arg,
         '); later computeInfoVal2IFC() calls reuse it. If the file cannot be written, ',
         'pass response_seed = <n> to computeInfoVal2IFC() instead.')
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

generateBaseReference <- function(selection, rdata, iter, ncores, response_seed, save_rdata,
                                  reference_method) {
  source <- selection$source
  if (length(iter) != 1L || !is.finite(iter) || iter < 1 || iter != trunc(iter)) {
    stop('iter must be a positive integer.')
  }
  if (is.null(response_seed)) requireStimulusSeed(source$seed, rdata, selection$baseimage)
  if (iter < 10000) warning('You should set iter >= 10000 for InfoVal statistic to be reliable')
  norms <- referenceNorms(source, selection$baseimage, NULL, iter, ncores, response_seed,
                          reference_method)
  if (save_rdata) {
    cache <- source$reference_norms_by_base
    if (is.null(cache)) cache <- list()
    if (!is.list(cache)) stop('reference_norms_by_base must be a list.')
    cache[[selection$baseimage]] <- list(
      norms = norms, response_seed = response_seed,
      source = 'saved_noise', fingerprint = referenceSnapshot(norms), method = reference_method
    )
    source$reference_norms_by_base <- cache
    saveRdataSafely(ls(source, all.names = TRUE), rdata, source)
  }
  invisible(norms)
}

baseReference <- function(rdata, iter, force_gen_ref_dist, response_seed, selection,
                          reference_method, report) {
  report(selection$source$n_trials, NULL)
  cache <- selection$source$reference_norms_by_base
  if (!is.null(cache) && !is.list(cache)) stop('reference_norms_by_base must be a list.')
  norms <- resolveReferenceNorms(cache[[selection$baseimage]], rdata, iter,
                                 force_gen_ref_dist, response_seed, selection$baseimage,
                                 seedless = is.null(selection$source$seed),
                                 reference_method = reference_method)
  if (!is.numeric(norms) || !length(norms) || any(!is.finite(norms))) {
    stop('Invalid cached reference for baseimage ', selection$baseimage,
         '. Use force_gen_ref_dist = TRUE to regenerate it.')
  }
  list(median = median(norms), mad = mad(norms), iter = length(norms),
       note = paste0('baseimage = ', selection$baseimage, '; '))
}

# Named so tests can model read-only archives even when running as root.
writableFile <- function(path) {
  unname(file.access(path, mode = 2)) == 0L
}
