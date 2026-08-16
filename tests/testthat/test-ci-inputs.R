# Direct tests for the input-preparation helpers generateCI() is built from.
#
# These assert each helper's own contract on values handed to it. The pipeline
# behaviour they implement is pinned end to end elsewhere and stays there --
# test-fixed-bugs.R for the 4096-parameter truncation, test-generateCI-paths.R
# for the aggregation, test-generateCI.R for the base-image lookup.

# --------------------------------------------------------------------------
# coerceTrialVectors
# --------------------------------------------------------------------------

test_that("coerceTrialVectors drops one-column tibbles to plain vectors", {
  skip_if_not_installed("tibble")

  tbl <- tibble::tibble(stim = 1:4, resp = c(1, -1, 1, -1))
  out <- rcicr:::coerceTrialVectors(tbl[, "stim"], tbl[, "resp"], NA)

  # The point of the coercion: a tibble column stays a one-column tibble under
  # [, "col"], which aggregate() later rejects with "arguments must have same
  # length". Assert it is a bare vector, not merely that it has four elements.
  expect_null(dim(out$stimuli))
  expect_null(dim(out$responses))
  expect_type(out$stimuli, "integer")
  expect_equal(out$stimuli, 1:4)
  expect_equal(out$responses, c(1, -1, 1, -1))
})

test_that("coerceTrialVectors leaves plain vectors and NA participants alone", {
  out <- rcicr:::coerceTrialVectors(1:3, c(1, -1, 1), NA)

  expect_equal(out$stimuli, 1:3)
  expect_equal(out$responses, c(1, -1, 1))
  # NA must survive as the sentinel: generateCI() branches on
  # all(is.na(participants)), so unlisting it into something else would send
  # every call down the per-participant path.
  expect_true(all(is.na(out$participants)))
  expect_length(out$participants, 1)
})

test_that("coerceTrialVectors coerces participants only when they are given", {
  skip_if_not_installed("tibble")

  tbl <- tibble::tibble(pid = c("a", "a", "b"))
  out <- rcicr:::coerceTrialVectors(1:3, c(1, -1, 1), tbl[, "pid"])

  expect_null(dim(out$participants))
  expect_equal(out$participants, c("a", "a", "b"))
})

test_that("coerceTrialVectors rejects mismatched stimuli and responses", {
  expect_error(rcicr:::coerceTrialVectors(1:4, c(1, -1), NA),
    "same length \\(stimuli: 4, responses: 2\\)"
  )
  # The counts in the message are the useful part, so check they track the
  # input rather than being a fixed string.
  expect_error(rcicr:::coerceTrialVectors(1:2, c(1, -1, 1), NA),
    "stimuli: 2, responses: 3"
  )
})

# --------------------------------------------------------------------------
# selectBaseImage
# --------------------------------------------------------------------------

test_that("selectBaseImage returns the requested base image", {
  faces <- list(one = matrix(0.1, 2, 2), two = matrix(0.9, 2, 2))

  expect_equal(rcicr:::selectBaseImage(faces, "two"), matrix(0.9, 2, 2))
})

test_that("selectBaseImage names the labels the file does contain", {
  faces <- list(alpha = matrix(0, 2, 2), beta = matrix(0, 2, 2))

  # Listing what is available is what makes this error actionable -- the usual
  # cause is a typo or a file from a different stimulus set.
  err <- expect_error(rcicr:::selectBaseImage(faces, "gamma"),
    "did not contain any reference to base image label: gamma"
  )
  expect_match(conditionMessage(err), "alpha, beta")
})

# --------------------------------------------------------------------------
# aggregateResponses
# --------------------------------------------------------------------------

test_that("aggregateResponses averages repeats and returns one row per stimulus", {
  out <- rcicr:::aggregateResponses(
    stimuli   = c(1, 1, 2, 3, 3, 3),
    responses = c(1, -1, 1, 1, 1, -1)
  )

  expect_equal(out$stimuli, c(1, 2, 3))
  # Stimulus 1 saw 1 and -1 (mean 0); stimulus 3 saw two 1s and a -1 (mean 1/3).
  expect_equal(out$responses, c(0, 1, 1 / 3))
})

test_that("aggregateResponses is a no-op on already-unique stimuli", {
  out <- rcicr:::aggregateResponses(c(3, 1, 2), c(1, -1, 1))

  # aggregate() sorts by the grouping variable, so the responses must be
  # carried along with it rather than staying in input order.
  expect_equal(out$stimuli, c(1, 2, 3))
  expect_equal(out$responses, c(-1, 1, 1))
})

# --------------------------------------------------------------------------
# selectStimulusParams
# --------------------------------------------------------------------------

test_that("selectStimulusParams picks the rows of the stimuli presented", {
  params <- list(base = matrix(1:20, nrow = 4))
  out <- rcicr:::selectStimulusParams(params, "base", c(1, 3))

  expect_equal(out, params$base[c(1, 3), ])
})

test_that("selectStimulusParams handles non-consecutive and repeated stimuli", {
  params <- list(base = matrix(1:20, nrow = 4))

  expect_equal(rcicr:::selectStimulusParams(params, "base", c(4, 1)),
    params$base[c(4, 1), ]
  )
  expect_equal(rcicr:::selectStimulusParams(params, "base", c(2, 2)),
    params$base[c(2, 2), ]
  )
})

test_that("selectStimulusParams truncates 4096-parameter files to 4092", {
  # Pre-0.3.0 rcicr drew 4096 contrasts per trial where only 4092 patch indices
  # exist. Both branches truncate; the single-trial one silently did not until
  # #235's sibling fix, which is pinned end to end in test-fixed-bugs.R.
  wide <- list(base = matrix(seq_len(3 * 4096), nrow = 3))

  multi <- rcicr:::selectStimulusParams(wide, "base", 1:2)
  expect_equal(ncol(multi), 4092)
  expect_equal(multi, wide$base[1:2, 1:4092])

  single <- rcicr:::selectStimulusParams(wide, "base", 2)
  expect_length(single, 4092)
  expect_equal(single, wide$base[2, 1:4092])
})

test_that("selectStimulusParams leaves a 4092-parameter file untouched", {
  # The truncation must fire on 4096 and nothing else: trimming a current file
  # would drop real parameters.
  ok <- list(base = matrix(seq_len(2 * 4092), nrow = 2))

  expect_equal(ncol(rcicr:::selectStimulusParams(ok, "base", 1:2)), 4092)
  expect_length(rcicr:::selectStimulusParams(ok, "base", 1), 4092)
})

test_that("selectStimulusParams errors when the base image has no parameters", {
  empty <- list(base = matrix(numeric(0), nrow = 0, ncol = 0))

  expect_error(rcicr:::selectStimulusParams(empty, "base", integer(0)),
    "No parameters found for base image: base"
  )
})
