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
  loadRdata(rdata, env)

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

# load() for a stimulus file, which is only ever damaged by a save that was
# interrupted. Every read of one goes through here, so that such a file points
# to the backup saveRdataSafely() left instead of failing with a bare load
# error before any save could reach its own recovery.
loadRdata <- function(file, envir) {
  tryCatch(load(file, envir = envir), error = function(e) {
    backup <- rdataBackupPath(file)
    if (file.exists(backup)) {
      stop(conditionMessage(e), '\n', restoreAdvice(file, backup), call. = FALSE)
    }
    stop(e)
  })
}

# Write a stimulus file in place, keeping a verified backup until the write is
# complete. The file is never replaced by renaming a new one over it: that
# would give it the R process's owner and group and drop its ACLs, which base R
# cannot restore. Power loss is not covered; base R has no fsync().
#
# The phases matter because an interruption must undo different things in
# each: before the backup exists there is nothing to restore, and once save()
# has returned there is nothing to roll back.
saveRdataSafely <- function(names, file, envir) {
  backup <- rdataBackupPath(file)

  # Phase 0: checks; nothing has been written.
  prefix <- paste0(basename(file), '.rcicr-staging-')
  leftovers <- list.files(dirname(file), all.files = TRUE)
  leftovers <- leftovers[startsWith(leftovers, prefix)]
  if (length(leftovers) > 0) {
    # Not deleted: one could belong to another session saving right now.
    warning('Incomplete copies of ', file, ' left by an interrupted save can be deleted: ',
            paste(file.path(dirname(file), leftovers), collapse = ', '), call. = FALSE)
  }
  if (file.exists(backup)) {
    if (!rdataLoads(file)) {
      stop(file, ' does not load. ', restoreAdvice(file, backup), call. = FALSE)
    }
    # The file is complete, possibly already the newer version, so the
    # backup is stale; restoring it could put older contents back.
    unlink(backup)
    message('Removed ', backup, ', left by an earlier save of a file that loads.')
  }
  staging <- tempfile(pattern = prefix, tmpdir = dirname(file))
  if (nchar(basename(staging), type = 'bytes') > 255) {
    warning(file, ' was saved without a backup: its name is too long to add one beside it.',
            call. = FALSE)
    return(invisible(writeRdata(names, file, envir)))
  }

  # Phase 1: back up. An interruption removes the staging copy; the original is untouched.
  if (!createStaging(staging)) {
    if (.Platform$OS.type == 'unix' && !writableDir(dirname(file))) {
      warning(file, ' was saved without a backup: its directory does not allow new files.',
              call. = FALSE)
      return(invisible(writeRdata(names, file, envir)))
    }
    stop('Could not create a backup beside ', file, ' (the disk may be full), ',
         'so it was not saved and is unchanged.', call. = FALSE)
  }
  phase <- 1
  on.exit(if (phase == 1) unlink(staging), add = TRUE)
  # Owner-only, whatever the original allows: the copy takes this process's
  # group. On Windows this sets only the read-only attribute, so the copy has
  # its folder's permissions.
  makePrivate(staging)
  if (.Platform$OS.type == 'unix' && file.mode(staging) != as.octmode('600')) {
    stop('Could not restrict the backup of ', file, ' to its owner, so it was not saved ',
         'and is unchanged.', call. = FALSE)
  }
  if (!isTRUE(copyInto(file, staging)) || !sameContents(file, staging)) {
    stop('Could not back up ', file, ' (the disk may be full), so it was not saved ',
         'and is unchanged.', call. = FALSE)
  }
  if (!renameFile(staging, backup)) {
    stop('Could not back up ', file, ', so it was not saved and is unchanged.', call. = FALSE)
  }

  # Phase 2: save. An error or interrupt restores the original from the backup.
  phase <- 2
  on.exit(if (phase == 2) restoreFromBackup(file, backup), add = TRUE)
  writeRdata(names, file, envir)

  # Phase 3: committed, so never rolled back.
  phase <- 3
  if (!removeFile(backup)) {
    warning(file, ' was saved, but its backup could not be deleted: ', backup,
            '. Delete it; it holds the previous contents.', call. = FALSE)
  }
  invisible(NULL)
}

rdataBackupPath <- function(file) paste0(file, '.rcicr-backup')

# The command is shown to be pasted into R, so both paths are encoded as R
# string literals: a quote or a Windows backslash would otherwise break it.
restoreAdvice <- function(file, backup) {
  lit <- function(x) encodeString(x, quote = '"')
  paste0('An interrupted save left a backup of it at ', backup, '. Restore it with ',
         'file.copy(', lit(backup), ', ', lit(file), ', overwrite = TRUE, copy.mode = FALSE), ',
         'check that the file loads again, and only then delete the backup. Copy rather ',
         'than rename, which would change the file\'s owner and permissions.')
}

restoreFromBackup <- function(file, backup) {
  if (!sameContents(file, backup)) copyInto(backup, file)
  if (sameContents(file, backup)) {
    unlink(backup)
  } else {
    warning(file, ' could not be restored after the failed save. ', restoreAdvice(file, backup),
            call. = FALSE)
  }
}

rdataLoads <- function(file) {
  tryCatch({
    suppressWarnings(load(file, envir = new.env(parent = emptyenv())))
    TRUE
  }, error = function(e) FALSE)
}

sameContents <- function(a, b) {
  isTRUE(unname(tools::md5sum(a)) == unname(tools::md5sum(b)))
}

# Named so tests can model a full disk, a refused rename, a partial write or a
# held file, none of which can be produced on demand.
writeRdata <- function(names, file, envir) save(list = names, file = file, envir = envir)
createStaging <- function(path) suppressWarnings(file.create(path))
makePrivate <- function(path) Sys.chmod(path, '600', use_umask = FALSE)
copyInto <- function(from, to) file.copy(from, to, overwrite = TRUE, copy.mode = FALSE)
renameFile <- function(from, to) file.rename(from, to)
removeFile <- function(path) suppressWarnings(file.remove(path))
writableDir <- function(path) unname(file.access(path, mode = 2)) == 0L
