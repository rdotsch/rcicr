test_that("computeCumulativeCICorrelation returns one correlation per step, ending at 1", {
  tmp <- withr::local_tempdir()
  rdata_path <- make_fixture_rdata(tmp, img_size = 32, n_trials = 6, nscales = 1, seed = 1)

  set.seed(1)
  responses <- sample(c(1, -1), 6, replace = TRUE)

  correlations <- suppressWarnings(
    computeCumulativeCICorrelation(
      stimuli = 1:6, responses = responses, baseimage = "base",
      rdata = rdata_path, step = 1
    )
  )

  expect_length(correlations, 6)
  # With no targetci supplied, the final CI is the one built from all trials,
  # so the cumulative correlation at the last trial must be 1.
  expect_equal(correlations[6], 1, tolerance = 1e-8)
})
