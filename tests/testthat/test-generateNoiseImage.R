test_that("generateNoiseImage returns a matrix of the expected size", {
  p <- generateNoisePattern(16, nscales = 1)
  img <- generateNoiseImage(rep(1, max(p$patchIdx)), p)
  expect_equal(dim(img), c(16, 16))
})

test_that("mismatched params length errors", {
  p <- generateNoisePattern(16, nscales = 1)
  expect_error(
    generateNoiseImage(rep(1, max(p$patchIdx) - 5), p),
    "Stimulus generation aborted"
  )
})

test_that("legacy 0-indexed patchIdx warns but still returns a valid image", {
  p0 <- generateNoisePattern(16, nscales = 1, pre_0.3.0 = TRUE)
  expect_equal(min(p0$patchIdx), 0)

  params <- rep(1, max(p0$patchIdx) + 1)
  expect_warning(
    result <- generateNoiseImage(params, p0),
    "patch indices start at 0"
  )
  expect_equal(dim(result), c(16, 16))
  expect_false(anyNA(result))
})

# NOTE: support for the pre-0.3.3 sinusoids/sinIdx list shape is currently
# broken; test-known-bugs.R holds a failing test for the intended behaviour.

test_that("all-zero params gives an all-zero image", {
  p <- generateNoisePattern(16, nscales = 1)
  result <- generateNoiseImage(rep(0, max(p$patchIdx)), p)
  expect_equal(result, matrix(0, 16, 16))
})
