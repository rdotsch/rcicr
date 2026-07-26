test_that("generateNoisePattern returns the expected list structure", {
  p <- generateNoisePattern(img_size = 32, nscales = 2)
  expect_named(p, c("patches", "patchIdx", "noise_type", "generator_version"))
})

test_that("number of patch layers scales with nscales", {
  for (nscales in c(1, 2, 3)) {
    p <- generateNoisePattern(img_size = 16, nscales = nscales)
    expect_equal(dim(p$patches)[3], nscales * 6 * 2)
  }
})

test_that("max patch index matches the documented nparams formula", {
  expected <- list("1" = 12, "2" = 60, "5" = 4092)
  for (nscales in c(1, 2, 5)) {
    p <- generateNoisePattern(img_size = 16, nscales = nscales)
    expect_equal(max(p$patchIdx), expected[[as.character(nscales)]])
    expect_equal(max(p$patchIdx), sum(6 * 2 * (2^(0:(nscales - 1)))^2))
  }
})

test_that("noise_type is passed through", {
  p <- generateNoisePattern(16, nscales = 1, noise_type = "gabor")
  expect_equal(p$noise_type, "gabor")
})

test_that("generator_version matches the installed package version", {
  p <- generateNoisePattern(16, nscales = 1)
  expect_equal(as.character(p$generator_version), as.character(utils::packageVersion("rcicr")))
})

test_that("pre_0.3.0 controls whether patch indices start at 0 or 1", {
  p_legacy <- generateNoisePattern(16, nscales = 1, pre_0.3.0 = TRUE)
  expect_equal(min(p_legacy$patchIdx), 0)

  p_current <- generateNoisePattern(16, nscales = 1, pre_0.3.0 = FALSE)
  expect_equal(min(p_current$patchIdx), 1)
})

test_that("patches and patchIdx contain no NAs", {
  p <- generateNoisePattern(16, nscales = 1)
  expect_false(anyNA(p$patches))
  expect_false(anyNA(p$patchIdx))
})
