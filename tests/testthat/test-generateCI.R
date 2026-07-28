test_that("unknown baseimage key errors with a helpful message", {
  tmp <- withr::local_tempdir()
  rdata_path <- make_fixture_rdata(tmp, img_size = 32, n_trials = 6, nscales = 1, seed = 1)

  set.seed(1)
  responses <- sample(c(1, -1), 6, replace = TRUE)

  expect_error(
    generateCI(
      stimuli = 1:6, responses = responses, baseimage = "doesnotexist",
      rdata = rdata_path, save_as_png = FALSE
    ),
    "did not contain any"
  )
})

test_that("generateCI returns a list of correctly-sized matrices with no NAs", {
  tmp <- withr::local_tempdir()
  rdata_path <- make_fixture_rdata(tmp, img_size = 32, n_trials = 6, nscales = 1, seed = 1)

  set.seed(1)
  responses <- sample(c(1, -1), 6, replace = TRUE)

  ci <- generateCI(
    stimuli = 1:6, responses = responses, baseimage = "base",
    rdata = rdata_path, save_as_png = FALSE
  )

  expect_named(ci, c("ci", "scaled", "base", "combined"))
  for (field in names(ci)) {
    expect_equal(dim(ci[[field]]), c(32, 32))
    expect_false(anyNA(ci[[field]]))
  }
})

test_that("antiCI negates the classification image for a single-response-set case", {
  tmp <- withr::local_tempdir()
  rdata_path <- make_fixture_rdata(tmp, img_size = 32, n_trials = 6, nscales = 1, seed = 1)

  set.seed(1)
  responses <- sample(c(1, -1), 6, replace = TRUE)

  ci <- generateCI(
    stimuli = 1:6, responses = responses, baseimage = "base",
    rdata = rdata_path, save_as_png = FALSE
  )
  anti_ci <- generateCI(
    stimuli = 1:6, responses = responses, baseimage = "base",
    rdata = rdata_path, save_as_png = FALSE, antiCI = TRUE
  )

  expect_equal(anti_ci$ci, -ci$ci)
})

test_that("generateCI masks the intended region, from a matrix or a PNG alike", {
  # The existing matrix-mask test in test-fixed-bugs.R covers the region for one
  # input form. The PNG form goes through the same applyMask() but a different
  # import path, and nothing pinned it: a polarity slip there would invert which
  # half of a face survives, silently and plausibly.
  tmp <- withr::local_tempdir()
  rdata_path <- make_fixture_rdata(tmp, img_size = 32, n_trials = 6, nscales = 1, seed = 1)

  mask <- matrix(1, 32, 32) # 1 = kept
  mask[1:10, 1:10] <- 0     # 0 = masked, as black is in a PNG

  mask_file <- file.path(tmp, "mask.png")
  png::writePNG(mask, mask_file)

  masked_ci <- function(m) {
    set.seed(1)
    generateCI(
      stimuli = 1:6, responses = sample(c(1, -1), 6, replace = TRUE),
      baseimage = "base", rdata = rdata_path, save_as_png = FALSE, mask = m
    )$ci
  }

  from_matrix <- masked_ci(mask)
  from_png <- masked_ci(mask_file)

  # The masked corner is dropped and the rest survives, for both input forms...
  for (ci in list(from_matrix, from_png)) {
    expect_true(all(is.na(ci[1:10, 1:10])))
    expect_false(anyNA(ci[11:32, 11:32]))
  }

  # ...and the two forms agree exactly, so neither can drift in polarity.
  expect_equal(from_matrix, from_png)
})
