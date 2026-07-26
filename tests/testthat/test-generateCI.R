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
