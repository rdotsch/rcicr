# Golden-master / characterization tests.
#
# These pin the numeric output of the core pipeline under NORMAL (default)
# parameters, as it behaved BEFORE the P0 bug fixes were applied.
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
# The explicit base:: prefixes below are kept deliberately. rcicr no longer
# imports matlab wholesale, so nothing masks sum() today - but matlab::sum()
# returns COLUMN sums for a matrix rather than a single total, which would turn
# a scalar assertion here into a 64-element one without failing. Prefixed, that
# cannot happen again. test-namespace-imports.R guards the import itself.

# helper: build a stimulus set at default nscales/noise_type, return .Rdata path

# Skipped on CRAN: the golden master exists to stop *us* changing researchers'
# numbers between commits. It is meaningless as a check on an end user's machine
# and it is expensive. It keeps running on GitHub Actions, which is where it
# actually does its job.
skip_on_cran()

baseline_rdata <- function(dir, img_size = 64, n_trials = 20, seed = 1) {
  input <- file.path(dir, "in")
  output <- file.path(dir, "out")
  dir.create(input)
  dir.create(output)
  bf <- file.path(input, "base.png")
  make_square_png(bf, size = img_size, seed = seed) # nolint: object_usage_linter.

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

test_that("the default pipeline produces stable z-maps", {
  # Unlike its neighbours, this block does NOT pin pre-P0-fix behaviour, and
  # cannot: z-maps were wrong before 1.2.0 fixed the sigma collision (an .Rdata
  # field named `sigma` shadowed generateCI()'s z-map blur argument, so the
  # blur ran at 25 rather than the documented 3). The only sane baseline is
  # current, post-fix output. Captured 2026-07-28 on R 4.3.3.
  #
  # Its second job is cross-platform. A test on one machine cannot check that
  # two platforms agree -- there is no second machine in scope -- but pinned
  # constants executed on every platform turn that question into "does each
  # platform match this stored number?". This file runs on macOS as well as
  # Linux (the R-CMD-check jobs set NOT_CRAN=true, so skip_on_cran() above does
  # not fire there), which is what makes these literals a cross-platform check
  # rather than a Linux-only one. Nothing else covers the z-map that way: the
  # release gate compares versions, not platforms, and runs on ubuntu only.
  #
  # The tolerance is deliberately 1e-8 rather than exact. Cross-platform
  # floating-point noise at the last bits is expected and is not what this is
  # hunting; a real divergence in blur() or the per-pixel t-tests would be far
  # larger.
  tmp <- withr::local_tempdir()
  rdata_path <- baseline_rdata(tmp)

  # Same stimulus set and responses as the classification-image test above, so
  # a failure here against a pass there isolates the z-map step itself.
  set.seed(42)
  responses <- sample(c(1, -1), 20, replace = TRUE)

  zmap_at <- function(method) {
    ci <- generateCI(
      1:20, responses, "face", rdata_path, save_as_png = FALSE,
      zmap = TRUE, zmapmethod = method, sigma = 3, threshold = 0,
      # zmapdecoration = FALSE is required, not stylistic: the default TRUE
      # dies inside base R with "figure margins too large" below 256px, and
      # this fixture is 64px (issue #177).
      zmapdecoration = FALSE, zmaptargetpath = file.path(tmp, "zmaps"),
      n_cores = 1
    )
    ci$zmap
  }

  quick <- zmap_at("quick")
  ttest <- zmap_at("t.test")

  expect_equal(dim(quick), c(64, 64))
  expect_equal(dim(ttest), c(64, 64))
  # threshold = 0 keeps every pixel, so these are whole z-maps rather than
  # sparse ones and the summaries below cover all of them.
  expect_false(anyNA(quick))
  expect_false(anyNA(ttest))

  # 'quick' ends in scale(), so its sum is 0 and its sd is 1 *by construction*.
  # They are asserted as a property -- that the standardisation happened -- and
  # deliberately not counted as value checks: they would hold for any
  # standardised matrix whatsoever, including a completely wrong one. The
  # statistics that actually carry information about this z-map follow.
  expect_equal(base::sum(quick), 0, tolerance = 1e-8)
  expect_equal(stats::sd(as.vector(quick)), 1, tolerance = 1e-8)

  expect_equal(base::sum(abs(quick)), 3312.8075360544, tolerance = 1e-8)
  expect_equal(min(quick), -3.1778878721, tolerance = 1e-8)
  expect_equal(max(quick), 2.6937556237, tolerance = 1e-8)
  expect_equal(quick[1, 1], -0.314076657789, tolerance = 1e-8)
  expect_equal(quick[32, 32], -0.210475043121, tolerance = 1e-8)
  expect_equal(quick[64, 64], -0.420994859868, tolerance = 1e-8)

  # 't.test' is the more valuable of the two here: one t-test per pixel, so it
  # exercises considerably more numerical machinery than the blur and is the
  # likelier of the pair to expose an architecture difference. It is not
  # standardised, so its sum and sd are real values.
  expect_equal(base::sum(ttest), 1.9976499050, tolerance = 1e-8)
  expect_equal(stats::sd(as.vector(ttest)), 1.1110921150, tolerance = 1e-8)
  expect_equal(base::sum(abs(ttest)), 3640.3369637317, tolerance = 1e-8)
  expect_equal(min(ttest), -3.7394405368, tolerance = 1e-8)
  expect_equal(max(ttest), 4.1456488445, tolerance = 1e-8)
  expect_equal(ttest[1, 1], -0.217696650289, tolerance = 1e-8)
  expect_equal(ttest[32, 32], 0.511619338519, tolerance = 1e-8)
  expect_equal(ttest[64, 64], -0.350640218903, tolerance = 1e-8)

  # The two methods must not be silently the same function: a wiring bug that
  # routed both to one implementation would satisfy every literal above only if
  # the literals were also wrong, but this states it directly.
  expect_false(isTRUE(all.equal(quick, ttest, check.attributes = FALSE)))
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
