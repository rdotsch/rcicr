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
