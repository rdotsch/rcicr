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

# CRAN Note avoidance


if (getRversion() >= "2.15.1")
  utils::globalVariables(
    c(
      # Suppress checking notes for variables loaded at runtime from .RData files
      "p", "s", "base_faces", "stimuli_params", "img_size", "base_face_files", "n_trials", "seed", "noise_type", "reference_norms",

      # Suppress checking notes for variables in foreach loop (parallel runs)
      "obs"
      )
  )
