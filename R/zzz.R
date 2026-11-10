# CRAN Note avoidance

if (getRversion() >= "2.15.1")
  utils::globalVariables(
    c(
      # Suppress checking notes for variables loaded at runtime from .RData files
      "p", "s", "base_faces", "stimuli_params", "img_size", "base_face_files", "n_trials", "seed", "noise_type", "reference_norms", "reference_norms_seed",

      # Suppress checking notes for variables in foreach loop (parallel runs)
      "obs", "trial"
    )
  )

# Under Ubuntu's intel-mkl without MKL_THREADING_LAYER=GNU, MKL's OpenMP and
# R's conflict and base R's matrix products come out wrong, mostly doubled.
# Products of this size were wrong in almost every column there; a correct
# BLAS is within about 1e-15 of plain R arithmetic, which no BLAS reaches.
blas_product_ok <- function() {
  n <- 40L
  k <- 20L
  x <- matrix(sin(seq_len(n * n)), n)
  a <- x + t(x)
  b <- matrix(ifelse(cos(seq_len(n * k) * 1.7) > 0, 1, -1), n, k)
  got <- a %*% b
  want <- vapply(seq_len(k), function(j) colSums(a * b[, j]), numeric(n))
  all(abs(got - want) <= 1e-8 * max(abs(want)))
}

.onAttach <- function(libname, pkgname) {
  if (!blas_product_ok()) {
    packageStartupMessage(
      "rcicr: this R's BLAS returned a wrong matrix product in a check run when rcicr was attached. ",
      "Reference norms and InfoVal values computed in this session will be wrong. ",
      "With Intel MKL from Debian or Ubuntu, setting the environment variable ",
      "MKL_THREADING_LAYER=GNU before starting R fixed this in our tests."
    )
  }
}
