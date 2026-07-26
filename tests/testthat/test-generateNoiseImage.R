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

test_that("legacy sinusoids/sinIdx list shape currently errors (known issue)", {
  # generateNoiseImage()'s validation checks `p$patchIdx` before the block
  # that renames `sinusoids`/`sinIdx` to `patches`/`patchIdx`, so a
  # pre-0.3.3-style `p` list is rejected instead of being supported as the
  # roxygen comment ("Pre 0.3.3 noise pattern, rename for appropriate use")
  # implies. Documenting current (broken) behavior rather than the
  # documented intent - not fixed as part of this task.
  p <- generateNoisePattern(16, nscales = 1)
  p_legacy <- list(sinusoids = p$patches, sinIdx = p$patchIdx)
  params <- rep(1, max(p$patchIdx))
  suppressWarnings(expect_error(
    generateNoiseImage(params, p_legacy),
    "Stimulus generation aborted"
  ))
})

test_that("all-zero params gives an all-zero image", {
  p <- generateNoisePattern(16, nscales = 1)
  result <- generateNoiseImage(rep(0, max(p$patchIdx)), p)
  expect_equal(result, matrix(0, 16, 16))
})
