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

test_that("the self-computed final CI matches generateCI() only under equal repeat counts", {
  skip_if_not_installed("withr")

  # This function does not aggregate repeated presentations and generateCI() does,
  # which is deliberate -- see DECISIONS.md, "computeCumulativeCICorrelation() does
  # not aggregate repeated stimuli". The difference is invisible where each stimulus
  # was presented the same number of times and real where it was not, so pinning
  # only one of those halves would not distinguish the intended behaviour from an
  # accidental one.
  tmp <- withr::local_tempdir()
  rdata_path <- make_fixture_rdata(tmp, img_size = 32, n_trials = 8, nscales = 2, seed = 1)

  e <- new.env()
  load(rdata_path, envir = e)

  self_ci <- function(stimuli, responses) {
    generateCINoise(e$stimuli_params$base[stimuli, , drop = FALSE], responses, e$p)
  }
  aggregated_ci <- function(stimuli, responses) {
    suppressWarnings(generateCI(
      stimuli = stimuli, responses = responses, baseimage = "base",
      rdata = rdata_path, save_as_png = FALSE, scaling = "none"
    ))$ci
  }

  equal_stim <- c(1, 1, 2, 2)
  equal_resp <- c(1, -1, 1, 1)
  expect_equal(self_ci(equal_stim, equal_resp), aggregated_ci(equal_stim, equal_resp))

  # Unequal counts: the two weight the data differently -- each trial equally here,
  # each unique stimulus equally there -- so they must differ, and by more than a
  # rounding artefact.
  uneq_stim <- c(1, 1, 1, 1, 2, 2, 3, 4)
  uneq_resp <- c(1, -1, 1, 1, -1, 1, 1, -1)
  uneq_self <- self_ci(uneq_stim, uneq_resp)
  uneq_aggregated <- aggregated_ci(uneq_stim, uneq_resp)

  expect_false(isTRUE(all.equal(uneq_self, uneq_aggregated)))
  expect_gt(max(abs(uneq_self - uneq_aggregated)), 1e-3)
  expect_lt(stats::cor(as.vector(uneq_self), as.vector(uneq_aggregated)), 0.99)
})

test_that("the endpoint of 1 holds only when the evaluated trials reach the last one", {
  skip_if_not_installed("withr")

  # Trials are taken at seq(1, length(responses), step), which need not include the
  # final trial. Where it does -- always at the default step = 1 -- the last point
  # compares the final CI with itself and is 1. Where it does not, the curve ends
  # on a partial CI and lands below 1, so documenting the endpoint as unconditional
  # would be wrong.
  tmp <- withr::local_tempdir()
  rdata_path <- make_fixture_rdata(tmp, img_size = 32, n_trials = 6, nscales = 1, seed = 1)

  responses <- withr::with_seed(2, sample(c(1, -1), 6, replace = TRUE))
  curve <- function(step) {
    suppressWarnings(computeCumulativeCICorrelation(
      stimuli = 1:6, responses = responses, baseimage = "base",
      rdata = rdata_path, step = step
    ))
  }

  step_1 <- curve(1)
  expect_length(step_1, 6)
  expect_equal(step_1[length(step_1)], 1, tolerance = 1e-8)

  # step = 2 evaluates trials 1, 3 and 5, so the sixth response never enters a
  # cumulative CI even though it is in the final CI being compared against. Three
  # points rather than six is what shows the sequence stopped short.
  step_2 <- curve(2)
  expect_length(step_2, 3)
  expect_lt(step_2[length(step_2)], 1 - 1e-6)
})
