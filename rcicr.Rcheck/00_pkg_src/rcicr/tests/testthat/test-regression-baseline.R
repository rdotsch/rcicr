# Golden-master / characterization tests.
#
# These pin the numeric output of the core pipeline under NORMAL (default)
# parameters, as it behaved BEFORE the BACKLOG.md P0 bug fixes were applied.
# Their job is to answer one question for researchers:
#
#   "Does fixing these bugs change the classification images I already have?"
#
# For the default configuration the answer must be NO, and these tests enforce
# that. Unlike test-known-bugs.R, these are expected to PASS at all times. If a
# bug fix turns one of them red, that fix changes results under normal usage
# and MUST be documented in NEWS.md under "Reproducibility impact" before it
# is merged - it is not a test to be casually updated.
#
# Values were captured on the pre-fix code (main @ 1487369, R 4.3.3) with
# nscales = 5 (default), noise_type = "sinusoid" (default), seed = 1.
#
# NOTE the explicit base:: prefixes below. The matlab package (imported by
# rcicr) exports its own sum(), which follows MATLAB semantics and returns
# COLUMN sums for a matrix rather than a single total. Under devtools::test()
# that masks base::sum() here, silently turning a scalar assertion into a
# 64-element one. Package code is currently unaffected (its only sum() calls
# are on vectors, where the two agree), but do not drop these prefixes.

# helper: build a stimulus set at default nscales/noise_type, return .Rdata path
baseline_rdata <- function(dir, img_size = 64, n_trials = 20, seed = 1) {
  input <- file.path(dir, "in")
  output <- file.path(dir, "out")
  dir.create(input)
  dir.create(output)
  bf <- file.path(input, "base.png")
  make_square_png(bf, size = img_size, seed = seed)

  generateStimuli2IFC(
    base_face_files = list(face = bf),
    n_trials = n_trials,
    img_size = img_size,
    stimulus_path = output,
    seed = seed,
    ncores = 1,
    save_as_png = FALSE
    # nscales and noise_type deliberately left at their defaults
  )
  list.files(output, pattern = "\\.Rdata$", full.names = TRUE)[1]
}

test_that("the default pipeline produces bit-stable classification images", {
  tmp <- withr::local_tempdir()
  rdata_path <- baseline_rdata(tmp)

  set.seed(42)
  responses <- sample(c(1, -1), 20, replace = TRUE)
  ci <- generateCI(1:20, responses, "face", rdata_path, save_as_png = FALSE)

  # Whole-image summaries: catch any change in the noise basis, the parameter
  # draw, the weighting, or the scaling.
  expect_equal(base::sum(ci$ci), -0.3441236120, tolerance = 1e-8)
  expect_equal(stats::sd(as.vector(ci$ci)), 0.0131844947, tolerance = 1e-8)

  expect_equal(base::sum(ci$scaled), 2044.1127686365, tolerance = 1e-8)
  expect_equal(base::sum(ci$base), 2024.5058823529, tolerance = 1e-8)
  expect_equal(base::sum(ci$combined), 2034.3093254947, tolerance = 1e-8)

  # Individual pixels: catch changes that preserve the summaries but move
  # values around (e.g. a transposed or reordered noise basis).
  expect_equal(ci$ci[1, 1], -0.002182945530, tolerance = 1e-8)
  expect_equal(ci$ci[32, 32], 0.005652215824, tolerance = 1e-8)
  expect_equal(ci$ci[64, 64], -0.004785719128, tolerance = 1e-8)
})

test_that("the noise basis itself is stable at default nscales", {
  # generateNoisePattern() is the foundation of every stimulus and every CI.
  # If this changes, every previously generated stimulus set is invalidated.
  p <- generateNoisePattern(img_size = 32, nscales = 5)

  expect_equal(dim(p$patches), c(32, 32, 60))
  expect_equal(max(p$patchIdx), 4092)
  expect_equal(base::sum(p$patches), 3859.937507680037, tolerance = 1e-8)
  expect_equal(stats::sd(as.vector(p$patches)), 0.7043160907, tolerance = 1e-8)
})

test_that("InfoVal is stable for a given CI and reference distribution", {
  # Guards the metric researchers actually report. Uses a pre-seeded reference
  # distribution so this pins the infoVal arithmetic itself (Euclidean norm,
  # median, mad) rather than the simulation.
  tmp <- withr::local_tempdir()
  rdata_path <- make_fixture_rdata(tmp, img_size = 32, n_trials = 6, nscales = 1, seed = 1)
  seed_reference_norms(rdata_path, n = 50, seed = 1)

  target_ci <- list(ci = matrix(withr::with_seed(2, rnorm(32 * 32)), 32, 32))

  invisible(utils::capture.output(
    iv <- computeInfoVal2IFC(target_ci, rdata_path)
  ))

  e <- new.env()
  load(rdata_path, envir = e)
  expected <- (norm(matrix(target_ci$ci), "f") - median(e$reference_norms)) /
    mad(e$reference_norms)

  expect_equal(iv, expected, tolerance = 1e-10)
  # Pinned literal: catches a change in the formula itself, which expected
  # (recomputed the same way) would not.
  expect_equal(iv, 13.4008068693, tolerance = 1e-8)
})
