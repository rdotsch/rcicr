# Helpers for reading .Rdata stimulus files: guarding the caller's frame
# against load(), and reporting which rcicr wrote a file.

# Snapshot a function's arguments so they can be restored after load(), which
# assigns straight into the calling frame and silently overwrites an argument
# that shares a name with an object in the .Rdata file.
#
# Only arguments that are *required* and were not supplied are skipped. mget()
# forces the promise, and for those it raises "argument is missing, with no
# default" -- which happens whenever a wrapper forwards its own missing
# argument, as batchGenerateCI() does with targetpath. Skipping them is also
# correct: an argument with no value cannot be overwritten into a wrong one.
#
# Defaulted arguments that were not supplied must NOT be skipped, even though
# missing() reports them missing too. Their default is the value the function
# goes on to use, and it is exactly as vulnerable to being replaced by the
# .Rdata file as one the caller passed explicitly.
captureArgs <- function(env) {
  fmls <- formals(sys.function(sys.parent()))
  nms <- names(fmls)
  required <- vapply(fmls, function(d) identical(d, quote(expr = )), logical(1))
  absent <- vapply(nms,
                   function(nm) eval(bquote(missing(.(as.name(nm)))), env),
                   logical(1))
  mget(nms[!(required & absent)], envir = env)
}

# Which fields an .Rdata file has follows from which rcicr wrote it, so the
# writing version is what turns "this file is missing stimuli_params" into
# "regenerate it, or install the version that made it". Appended to the
# validation errors, which run on a frame `load()` has just written into.
#
# p$generator_version is preferred because the top-level field was a hardcoded
# '0.4.0' from 0.4.0 through 1.1.0 -- see DECISIONS.md, "`generator_version` in
# old `.Rdata` files is not trustworthy". Either may be a character string (old
# files) or a package_version (since 1.2.0).
rdataWriterNote <- function(env) {
  field <- function(name, from = env) {
    if (exists(name, envir = from, inherits = FALSE)) get(name, envir = from, inherits = FALSE)
  }

  p <- field('p')
  version <- if (is.list(p)) p$generator_version
  if (is.null(version)) version <- field('generator_version')

  if (is.null(version) || !length(version)) {
    # Not "so it predates 0.4.0": a truncated file, a hand-rewritten one, or one
    # rcicr never wrote also has no version field, and this runs on a file
    # already known to be broken. State the absence and stop there.
    return(paste0(' The file records no writing version. rcicr has recorded one since 0.4.0,',
                  ' so what wrote this file is unknown.'))
  }

  version <- tryCatch(as.character(numeric_version(version)),
                      error = function(e) as.character(version))

  if (identical(version, '0.4.0')) {
    return(paste0(' The file reports rcicr 0.4.0, which is what every version from',
                  ' 0.4.0 through 1.1.0 recorded, so what wrote it is unknown.'))
  }

  paste0(' The file was written by rcicr ', version, '.')
}

# Read a stimulus .Rdata file into a frame of its own and return the four
# objects the CI pipeline uses.
#
# The point is where it loads. load() assigns into whatever environment it is
# given, so reading the file directly inside generateCI() put every field of it
# in scope alongside that function's arguments, where a shared name silently
# won: the noise `sigma` stored since 1.1.0 replaced the z-map blur `sigma`, and
# the value the caller passed was ignored. That is why generateCI() no longer
# calls captureArgs(); computeInfoVal2IFC() and computeCumulativeCICorrelation()
# still load into their own frames and still need it.
#
# A dedicated environment rather than this function's own frame, because the
# frame is not argument-free either: it holds `rdata`, and an older
# generateReferenceDistribution2IFC() saved its own `rdata` argument into the
# file, so that name occurs in real files. Loading into the frame would be safe
# only for as long as nothing read `rdata` after the load -- a hazard narrowed
# rather than removed.
loadStimulusParams <- function(rdata) {
  env <- new.env(parent = emptyenv())
  load(rdata, envir = env)

  has <- function(name) exists(name, envir = env, inherits = FALSE)
  take <- function(name) get(name, envir = env, inherits = FALSE)

  if (!has('s') && !has('p')) {
    stop('File specified in rdata did not contain s or p variable.', rdataWriterNote(env))
  }

  if (!has('base_faces')) {
    stop('File specified in rdata did not contain base_faces variable.', rdataWriterNote(env))
  }

  if (!has('stimuli_params')) {
    stop('File specified in rdata did not contain stimuli_params variable.', rdataWriterNote(env))
  }

  if (!has('img_size')) {
    stop('File specified in rdata did not contain img_size variable.', rdataWriterNote(env))
  }

  # Convert s to p (if rdata file originates from pre-0.3.3). Checked before p
  # because that is the precedence the in-frame version had: it overwrote a
  # loaded p whenever s was also present.
  p <- if (has('s')) {
    s <- take('s')
    list(patches = s$sinusoids, patchIdx = s$sinIdx, noise_type = 'sinusoid')
  } else {
    take('p')
  }

  return(list(p = p, base_faces = take('base_faces'),
    stimuli_params = take('stimuli_params'), img_size = take('img_size')
  ))
}

# The reader for the latent module's stimulus files.
#
# Like loadStimulusParams() this loads into an isolated environment rather than
# the caller's frame, so none of the saved names -- several of which are also
# argument names, deliberately, because they mean the same thing -- can
# overwrite an argument. Nothing here needs captureArgs() as a result.
loadLatentStimulusParams <- function(rdata) {
  env <- new.env(parent = emptyenv())
  load(rdata, envir = env)

  has <- function(name) exists(name, envir = env, inherits = FALSE)
  take <- function(name) get(name, envir = env, inherits = FALSE)

  # A pixel-noise file loaded here would otherwise fail on a missing field and
  # send the reader looking for a corrupt file rather than the wrong function.
  if (has('stimuli_params') && !has('latent_params')) {
    msg <- paste0(
      'File specified in rdata is a pixel-noise stimulus file, not a latent ',
      'one. Use generateCI() for it, or generateStimuliLatent2IFC() to make ',
      'latent stimuli.'
    )
    stop(msg, call. = FALSE)
  }

  for (field in c('latent_params', 'base_latent', 'generator_spec', 'img_size')) {
    if (!has(field)) {
      msg <- paste0('File specified in rdata did not contain ', field, ' variable.')
      stop(msg, rdataWriterNote(env))
    }
  }

  return(list(
    latent_params = take('latent_params'),
    base_latent = take('base_latent'),
    generator_spec = take('generator_spec'),
    img_size = take('img_size'),
    latent_sigma = if (has('latent_sigma')) take('latent_sigma') else NA_real_,
    # The label the stimuli were written under. Read rather than parsed back out
    # of the file name: a label containing "_seed_" -- which nothing forbids --
    # truncates at the first occurrence, so the classification image is named
    # after part of its set and two sets can produce the same file name.
    label = if (has('label')) take('label') else NA_character_
  ))
}
