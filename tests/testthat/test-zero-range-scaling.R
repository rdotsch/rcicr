# Scaling a classification image that has no range.
#
# Responses that cancel exactly give an all-zero CI. It is a result, not a
# failure, so the range-based methods render it neutral rather than dividing by
# a zero range and returning NaN (issue #303).

# Two presentations of each stimulus with opposite responses cancel exactly.
zero_stimuli <- c(1, 1, 2, 2)
zero_responses <- c(1, -1, 1, -1)

test_that("a zero-signal CI scales to a finite image under every method", {
  tmp <- withr::local_tempdir()
  rdata <- make_fixture_rdata(tmp, img_size = 32, n_trials = 2, nscales = 1, seed = 1)

  for (scaling in c("none", "constant", "matched", "independent")) {
    ci <- suppressWarnings(generateCI(zero_stimuli, zero_responses, "base", rdata,
                                      save_as_png = FALSE, n_cores = 1, scaling = scaling))

    expect_true(all(ci$ci == 0), info = scaling)
    expect_true(all(is.finite(ci$scaled)), info = scaling)
    expect_true(all(is.finite(ci$combined)), info = scaling)
  }
})

test_that("independent scaling renders a zero-signal CI at 0.5, and says so", {
  tmp <- withr::local_tempdir()
  rdata <- make_fixture_rdata(tmp, img_size = 32, n_trials = 2, nscales = 1, seed = 1)

  expect_warning(
    ci <- generateCI(zero_stimuli, zero_responses, "base", rdata,
                     save_as_png = FALSE, n_cores = 1, scaling = "independent"),
    regexp = "no range to scale"
  )

  expect_equal(unique(as.vector(ci$scaled)), 0.5)
  # 'constant' scaling already answers 0.5 for this CI; the two now agree.
  constant_ci <- generateCI(zero_stimuli, zero_responses, "base", rdata,
                            save_as_png = FALSE, n_cores = 1, scaling = "constant")
  expect_equal(ci$scaled, constant_ci$scaled)
})

test_that("matched scaling renders a zero-signal CI at the base image's midrange", {
  tmp <- withr::local_tempdir()
  base_png <- file.path(tmp, "dim.png")
  # Range [0.2, 0.6], and contrast maximization off so it survives: the
  # midrange is 0.4, which tells this policy apart from a literal 0.5.
  withr::with_seed(1, {
    png::writePNG(matrix(runif(32 * 32, 0.2, 0.6), 32, 32), base_png)
  })
  suppressWarnings(generateStimuli2IFC(
    base_face_files = list(dim = base_png), n_trials = 2, img_size = 32,
    stimulus_path = tmp, seed = 1, nscales = 1, ncores = 1,
    maximize_baseimage_contrast = FALSE, save_as_png = FALSE
  ))
  rdata <- list.files(tmp, pattern = "\\.Rdata$", full.names = TRUE)

  e <- new.env()
  load(rdata, envir = e)
  base <- e$base_faces$dim
  expected <- min(base) + (max(base) - min(base)) / 2

  expect_warning(
    ci <- generateCI(zero_stimuli, zero_responses, "dim", rdata,
                     save_as_png = FALSE, n_cores = 1, scaling = "matched"),
    regexp = "no range to scale"
  )

  expect_equal(unique(as.vector(ci$scaled)), expected)
  # Not vacuous: this base image's midrange is not the 0.5 the other methods use.
  expect_false(isTRUE(all.equal(expected, 0.5)))
})

test_that("a masked zero-signal CI keeps its NA and gains no NaN", {
  tmp <- withr::local_tempdir()
  rdata <- make_fixture_rdata(tmp, img_size = 32, n_trials = 2, nscales = 1, seed = 1)

  mask <- matrix(1, 32, 32)
  mask[1:16, ] <- 0

  expect_warning(
    ci <- generateCI(zero_stimuli, zero_responses, "base", rdata, save_as_png = FALSE,
                     n_cores = 1, scaling = "independent", mask = mask),
    regexp = "no range to scale"
  )

  expect_equal(sum(is.na(ci$ci)), 512)
  expect_true(all(is.na(ci$scaled[mask == 0])))
  expect_equal(unique(as.vector(ci$scaled[mask == 1])), 0.5)
})

test_that("autoscale renders an all-zero list neutral and leaves a mixed list alone", {
  tmp <- withr::local_tempdir()
  rdata <- make_fixture_rdata(tmp, img_size = 32, n_trials = 2, nscales = 1, seed = 1)

  zero <- generateCI(zero_stimuli, zero_responses, "base", rdata,
                     save_as_png = FALSE, n_cores = 1, scaling = "none")
  signal <- generateCI(c(1, 2), c(1, -1), "base", rdata,
                       save_as_png = FALSE, n_cores = 1, scaling = "none")
  expect_gt(max(abs(signal$ci)), 0)

  expect_warning(
    all_zero <- autoscale(list(a = zero, b = zero), save_as_pngs = FALSE),
    regexp = "exactly zero"
  )
  expect_equal(unique(as.vector(all_zero$a$scaled)), 0.5)
  expect_equal(unique(as.vector(all_zero$b$scaled)), 0.5)
  # Only $scaled is written here, degenerate or not.
  expect_equal(all_zero$a$combined, zero$combined)

  # A list holding signal never reached the degenerate branch and must not now.
  mixed <- autoscale(list(a = zero, b = signal), save_as_pngs = FALSE)
  expect_equal(unique(as.vector(mixed$a$scaled)), 0.5)
  expect_true(all(is.finite(mixed$b$scaled)))
  # autoscale() deliberately does not touch $combined.
  expect_equal(mixed$a$combined, zero$combined)
})

test_that("the guards leave every CI that has a range untouched", {
  base <- matrix(seq(0, 1, length.out = 16), 4, 4)
  ci <- matrix(c(-0.2, 0.1, 0.3, -0.05), 4, 4)

  expect_equal(rcicr:::applyScaling(base, ci, "none", 0.1), ci)
  expect_equal(rcicr:::applyScaling(base, ci, "constant", 0.5), (ci + 0.5) / 1)
  expect_equal(
    rcicr:::applyScaling(base, ci, "independent", 0.1),
    (ci + max(abs(ci))) / (2 * max(abs(ci)))
  )
  expect_equal(
    rcicr:::applyScaling(base, ci, "matched", 0.1),
    min(base) + (max(base) - min(base)) * (ci - min(ci)) / (max(ci) - min(ci))
  )
})

test_that("the two guards fire on different CIs", {
  base <- matrix(seq(0, 1, length.out = 16), 4, 4)
  uniform <- matrix(0.3, 4, 4)

  # A uniform *non-zero* CI has no range for 'matched' to map...
  expect_warning(matched <- rcicr:::applyScaling(base, uniform, "matched", 0.1),
                 regexp = "no range to scale")
  expect_equal(unique(as.vector(matched)), 0.5)

  # ...but 'independent' divides by its magnitude, which is not zero, so it
  # renders at the top of the range exactly as before and must not warn.
  expect_no_warning(independent <- rcicr:::applyScaling(base, uniform, "independent", 0.1))
  expect_equal(unique(as.vector(independent)), 1)
})

test_that("an entirely masked CI is left as it was", {
  base <- matrix(seq(0, 1, length.out = 16), 4, 4)
  all_na <- matrix(NA_real_, 4, 4)

  # No values at all is a masking mistake rather than a zero-signal result, and
  # is not what the degenerate-range guards answer.
  scaled <- suppressWarnings(rcicr:::applyScaling(base, all_na, "independent", 0.1))
  expect_true(all(is.na(scaled)))
})
