# The generator contract.
#
# R cannot run a StyleGAN, so the method and the renderer are separated: every
# function in the latent module works on latent vectors and calls out through
# this contract to turn them into pixels. That keeps the arithmetic testable
# with no GPU, no Python and no network, and it lets a lab plug in whatever
# generator it already has.
#
# A generator is a named list carrying the class attribute below. The fields:
#
#   kind        backend identifier, e.g. 'pca'
#   latent_dim  length of one latent vector
#   img_size    rendered images are img_size by img_size
#   space       which space the latents live in: 'pca', 'z', 'w', 'w+'
#   latent_mean centre of the space, length latent_dim
#   latent_sd   per-dimension scale, so a caller's sigma is in SD units
#   render      function(latents) taking a matrix of latents, one row per
#               stimulus, returning an array[nrow(latents), img_size,
#               img_size] of greyscale values in [0, 1]
#   fingerprint string identifying this renderer across sessions
#   state       backend-private
#
# render() is batched rather than one-latent-at-a-time because an external
# backend pays its startup cost per call, where a per-trial loop would dominate
# the runtime of a 770-trial study.

rcicrGenerator <- function(kind, latent_dim, img_size, space, latent_mean, latent_sd, render, fingerprint, state = list()) {
  structure(
    list(
      kind = kind,
      latent_dim = as.integer(latent_dim),
      img_size = as.integer(img_size),
      space = space,
      latent_mean = as.numeric(latent_mean),
      latent_sd = as.numeric(latent_sd),
      render = render,
      fingerprint = fingerprint,
      state = state
    ),
    class = c(paste0('rcicr_generator_', kind), 'rcicr_generator')
  )
}

# A generator is checked at the entry point rather than trusted, so a backend
# that returns the wrong shape fails on the call instead of 30 trials into a
# render loop with a non-conformable-arrays error from inside a worker.
#
# probe = FALSE skips the single test render, for a caller that has already
# rendered through this generator in the same call.
validateGenerator <- function(generator, probe = TRUE) {
  if (!inherits(generator, 'rcicr_generator')) {
    msg <- paste0(
      'Expected a generator built by one of rcicr\'s ',
      'latentGenerator* functions. Got an object of class ',
      paste(class(generator), collapse = '/'), '.'
    )
    stop(msg, call. = FALSE)
  }

  required <- c(
    'kind', 'latent_dim', 'img_size', 'space', 'latent_mean', 'latent_sd',
    'render', 'fingerprint'
  )
  missing_fields <- setdiff(required, names(generator))
  if (length(missing_fields) > 0) {
    msg <- paste0(
      'This generator is missing the field(s) ',
      paste(missing_fields, collapse = ', '), '.'
    )
    stop(msg, call. = FALSE)
  }

  checkPositiveScalar(generator$latent_dim, 'latent_dim')
  checkPositiveScalar(generator$img_size, 'img_size')

  if (!is.function(generator$render)) {
    stop('This generator\'s render field is not a function.', call. = FALSE)
  }

  for (field in c('latent_mean', 'latent_sd')) {
    value <- generator[[field]]
    if (length(value) != generator$latent_dim) {
      msg <- paste0(
        'This generator\'s ', field, ' has length ', length(value),
        ', but its latent_dim is ', generator$latent_dim, '.'
      )
      stop(msg, call. = FALSE)
    }
    if (!all(is.finite(value))) {
      msg <- paste0('This generator\'s ', field, ' contains non-finite values.')
      stop(msg, call. = FALSE)
    }
  }

  if (any(generator$latent_sd < 0)) {
    stop('This generator\'s latent_sd contains negative values.', call. = FALSE)
  }

  if (!is.character(generator$fingerprint) || length(generator$fingerprint) != 1) {
    stop('This generator\'s fingerprint is not a single string.', call. = FALSE)
  }

  if (probe) {
    validateRender(generator)
  }

  return(invisible(generator))
}

# One render of the centre of the space, checked against the declared contract.
validateRender <- function(generator) {
  probe <- renderLatent(generator, matrix(generator$latent_mean, nrow = 1), validate = FALSE)

  expected <- c(1L, generator$img_size, generator$img_size)
  if (!is.array(probe) || length(dim(probe)) != 3 || !identical(as.integer(dim(probe)), expected)) {
    got <- if (is.array(probe)) {
      paste(dim(probe), collapse = ' by ')
    } else {
      paste('an object of class', paste(class(probe), collapse = '/'))
    }
    msg <- paste0(
      'This generator\'s render() must return an array of ',
      'nrow(latents) by img_size by img_size. Rendering one latent returned ', got,
      ', where ', paste(expected, collapse = ' by '), ' was expected.'
    )
    stop(msg, call. = FALSE)
  }

  if (!is.numeric(probe) || !all(is.finite(probe))) {
    stop('This generator\'s render() returned non-finite pixel values.', call. = FALSE)
  }

  if (min(probe) < 0 || max(probe) > 1) {
    msg <- paste0(
      'This generator\'s render() returned pixel values outside ', '[0, 1] (range ',
      format(min(probe)), ' to ', format(max(probe)),
      '). Clamp or rescale inside render().'
    )
    stop(msg, call. = FALSE)
  }

  return(invisible(TRUE))
}

# Render latents through a generator, accepting a single vector as well as a
# matrix so callers holding one latent do not each write the same as.matrix.
renderLatent <- function(generator, latents, validate = TRUE) {
  if (validate) {
    validateGenerator(generator, probe = FALSE)
  }

  if (is.null(dim(latents))) {
    latents <- matrix(latents, nrow = 1)
  }
  latents <- as.matrix(latents)

  if (ncol(latents) != generator$latent_dim) {
    msg <- paste0(
      'Latents have ', ncol(latents), ' dimension(s), but this ',
      'generator\'s latent_dim is ', generator$latent_dim, '.'
    )
    stop(msg, call. = FALSE)
  }

  return(generator$render(latents))
}

checkPositiveScalar <- function(value, name) {
  if (length(value) != 1 || !is.finite(value) || value < 1) {
    msg <- paste0(
      'A generator\'s ', name, ' must be a single positive ', 'number. It is ',
      paste(format(value), collapse = ', '), '.'
    )
    stop(msg, call. = FALSE)
  }

  return(invisible(TRUE))
}

# No hashing function exists in base R and none of the Imports provides one, so
# a fingerprint is a set of summary statistics at full precision. It detects a
# generator that is not the one the stimuli were made with, which is what it is
# for. It is not a cryptographic digest and does not survive being described as
# one: a determined caller can construct a collision.
fingerprintNumeric <- function(...) {
  parts <- lapply(list(...), function(x) {
    x <- as.numeric(x)
    if (length(x) == 0) {
      return('0')
    }
    stats <- c(length(x), sum(x), sum(x^2), sum(x * seq_along(x)), min(x), max(x))
    paste(format(stats, digits = 17, scientific = TRUE), collapse = ',')
  })

  return(paste(unlist(parts), collapse = ':'))
}
