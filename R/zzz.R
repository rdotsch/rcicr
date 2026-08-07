# Default number of CPU cores for the parallel loops.
#
# Normally detectCores() - 1, which is what this package has always used. Under
# R CMD check it returns 2: CRAN policy caps packages at two cores in examples,
# tests and vignettes, and R sets _R_CHECK_LIMIT_CORES_ while checking. The
# environment variable is only ever set by the checker, so this does not change
# what any user gets at the console -- it only keeps the package inside CRAN's
# limit while it is being checked.
default_ncores <- function() {
  if (nzchar(Sys.getenv("_R_CHECK_LIMIT_CORES_"))) {
    return(2L)
  }
  max(1L, parallel::detectCores() - 1)
}

# Register a foreach backend and return the cluster to stop afterwards, or NULL
# when running serially.
#
# With ncores == 1 there is nothing to parallelise, but the old code built a
# cluster anyway: makeCluster(1) starts a second R process, makes it
# library(rcicr), ships each iteration to it and waits. That is pure overhead,
# and it dominated R CMD check -- the test suite reported [8s/126s], eight
# seconds of CPU against 126 elapsed, across 22 cluster spawns.
#
# registerDoSEQ() runs the very same %dopar% loop in the current process, so no
# loop body needed to change. Numeric output is unaffected: neither parallel
# loop in this package draws random numbers (stimulus parameters are drawn under
# set.seed() before the loop starts), so there is no worker RNG stream to
# diverge from the sequential one.
startBackend <- function(ncores) {
  if (is.null(ncores) || is.na(ncores) || ncores < 2) {
    foreach::registerDoSEQ()
    return(NULL)
  }

  cl <- parallel::makeCluster(ncores, outfile = "")
  doParallel::registerDoParallel(cl)
  cl
}

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

# CRAN Note avoidance


if (getRversion() >= "2.15.1")
  utils::globalVariables(
    c(
      # Suppress checking notes for variables loaded at runtime from .RData files
      "p", "s", "base_faces", "stimuli_params", "img_size", "base_face_files", "n_trials", "seed", "noise_type", "reference_norms", "reference_norms_seed",

      # Suppress checking notes for variables in foreach loop (parallel runs)
      "obs"
      )
  )
