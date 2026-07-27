test_that("vector stimuli path (single trial) matches generateNoiseImage directly", {
  p <- generateNoisePattern(16, nscales = 1)
  nparams <- max(p$patchIdx)
  stim_vec <- rep(1, nparams)

  expect_equal(
    generateCINoise(stim_vec, 1, p),
    generateNoiseImage(stim_vec * 1, p)
  )
})

test_that("multi-row stimuli are averaged after weighting by responses", {
  p <- generateNoisePattern(16, nscales = 1)
  nparams <- max(p$patchIdx)
  stim <- matrix(c(rep(1, nparams), rep(-1, nparams)), nrow = 2, byrow = TRUE)
  resp <- c(1, -1)

  expect_equal(
    generateCINoise(stim, resp, p),
    generateNoiseImage(colMeans(stim * resp), p)
  )
})

test_that("negating responses negates the result for a single trial", {
  p <- generateNoisePattern(16, nscales = 1)
  nparams <- max(p$patchIdx)
  stim_vec <- rep(1, nparams)

  expect_equal(
    generateCINoise(stim_vec, -1, p),
    -generateCINoise(stim_vec, 1, p)
  )
})
