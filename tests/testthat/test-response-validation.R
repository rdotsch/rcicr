# Responses that cannot weight a trial stop the call instead of turning the CI
# NA or weighting FALSE trials 0 (#372).

quietly <- function(expr) {
  out <- NULL
  utils::capture.output(out <- suppressMessages(suppressWarnings(expr)))
  out
}

responses8 <- c(1, -1, -1, 1, 1, 1, -1, 1)

test_that("missing, infinite and non-numeric responses stop generateCI()", {
  rdata <- make_fixture_rdata(withr::local_tempdir(), n_trials = 8)
  ci <- function(responses, participants = NA) {
    quietly(generateCI(1:8, responses, "base", rdata, participants = participants,
                       save_as_png = FALSE, n_cores = 1))
  }
  gap <- responses8
  gap[c(3, 6)] <- c(NA, NaN)
  for (participants in list(NA, rep(1:2, each = 4))) {
    expect_error(ci(gap, participants), "no finite value for 2 of 8 trials (trials 3, 6)",
                 fixed = TRUE)
    expect_error(ci(replace(responses8, 8, Inf), participants), "(trial 8)", fixed = TRUE)
    expect_error(ci(factor(responses8), participants), "These are factor")
    expect_error(ci(as.character(responses8), participants), "These are character")
    expect_error(ci(responses8 > 0, participants), "ifelse(responses, 1, -1)", fixed = TRUE)
  }
  # The recoding the logical message gives reproduces the numeric CI exactly.
  expect_identical(ci(ifelse(responses8 > 0, 1, -1))$ci, ci(responses8)$ci)
  # Ratings still weight trials.
  expect_false(anyNA(ci(c(3, -2, 1, 0, 2, -1, 1, 3))$ci))
})

test_that("computeCumulativeCICorrelation() rejects them too", {
  rdata <- make_fixture_rdata(withr::local_tempdir(), n_trials = 8)
  expect_error(quietly(computeCumulativeCICorrelation(1:8, replace(responses8, 2, NA), "base", rdata)),
               "(trial 2)", fixed = TRUE)
})

test_that("the batch functions name the groups with missing responses before computing any CI", {
  rdata <- make_fixture_rdata(withr::local_tempdir(), n_trials = 8)
  data <- data.frame(pp = rep(c("a", "b"), each = 4), stim = 1:8, resp = responses8)
  data$resp[c(6, 7)] <- NA
  local_mocked_bindings(generateCI = function(...) stop("computed"))
  for (batch in list(batchGenerateCI, batchGenerateCI2IFC)) {
    expect_error(quietly(batch(data = data, by = "pp", stimuli = "stim", responses = "resp",
                               baseimage = "base", rdata = rdata, save_as_png = FALSE)),
                 "2 rows have no finite response in column resp (pp b)", fixed = TRUE)
  }
  data$resp <- data$resp > 0
  expect_error(quietly(batchGenerateCI2IFC(data = data, by = "pp", stimuli = "stim",
                                           responses = "resp", baseimage = "base",
                                           rdata = rdata, save_as_png = FALSE)),
               "These are logical")
})
