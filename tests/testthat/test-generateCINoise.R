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


test_that("CI pixels equal the mean signed trial images for scaled responses", {
  p <- list(
    patches = array(c(2, -1, 4, 3, 5, -2, 1, 6, -3, 2, 7, -4), c(2, 3, 2)),
    patchIdx = array(c(1, 2, 3, 1, 2, 3, 3, 1, 2, 2, 3, 1), c(2, 3, 2))
  )
  stimuli <- rbind(c(1, 4, -2), c(3, -1, 5), c(-4, 2, 1), c(2, 6, -3))
  responses <- c(1, -1, 0.5, 0)
  expected <- matrix(0, 2, 3)
  # Expand each trial into pixels before combining responses: no production
  # helper or parameter-averaging expression participates in the oracle.
  for (row in 1:2) {
    for (column in 1:3) {
      for (trial in 1:4) {
        for (layer in 1:2) {
          expected[row, column] <- expected[row, column] +
            responses[trial] * stimuli[trial, p$patchIdx[row, column, layer]] *
              p$patches[row, column, layer] / (4 * 2)
        }
      }
    }
  }
  expect_equal(generateCINoise(stimuli, responses, p), expected, tolerance = 1e-14)
})
