# Parallel execution helpers: choosing a core count, registering and
# tearing down the foreach backend, and reporting progress from the parent.

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
#
# doSNOW rather than doParallel because only it honours .options.snow, whose
# progress callback runs in the *parent*. Ticking a progress bar from inside the
# loop body updates each worker's private copy, which nobody sees -- issue #178.
startBackend <- function(ncores) {
  if (is.null(ncores) || is.na(ncores) || ncores < 2) {
    foreach::registerDoSEQ()
    return(NULL)
  }

  cl <- parallel::makeCluster(ncores, outfile = "")
  doSNOW::registerDoSNOW(cl)
  cl
}

# Progress callback for a %dopar% loop, as .options.snow expects it.
#
# Empty list when running serially: registerDoSEQ() ignores .options.snow
# entirely, so the serial path keeps ticking the bar from inside the loop body
# and the two mechanisms never both run on the same loop.
progressOption <- function(pb, cl) {
  if (is.null(cl)) return(list())
  list(progress = function(n) setTxtProgressBar(pb, n))
}

# Stop a parallel cluster if it is still running.
# Intended for on.exit(), so that workers are released even when an error
# interrupts a foreach loop. Without that, the socket connections leak and R
# reports "closing unused connections" warnings later (issue #50).
# Callers set their cluster variable to NULL after a normal stopCluster(), which
# turns the registered on.exit() call into a no-op.
# Input: cluster object, or NULL
# Output: nothing
stopClusterSafely <- function(cl) {
  if (!is.null(cl)) {
    try(parallel::stopCluster(cl), silent = TRUE)
  }
  invisible(NULL)
}
