# --------------------------------------------------------------------------
# The .Rdata validation errors name the version that wrote the file (#181)
#
# Which fields a file carries follows from which rcicr saved it, so the version
# is what turns "this file is missing stimuli_params" into "regenerate it, or
# install the version that made it". The field it reads cannot be trusted at
# face value: the top-level generator_version was a hardcoded '0.4.0' from 0.4.0
# through 1.1.0, so these tests pin the three things it must say apart -- a real
# version, that sentinel, and no version at all.
#
# Reached by removing a field from a good fixture, as the other "did not
# contain" tests in test-error-paths.R do, so the file stays realistic.
# --------------------------------------------------------------------------

# Builds a fixture missing `stimuli_params`, with whatever version fields the
# caller asks for. `p_version = NULL` drops it, so the fallback is exercised.
broken_rdata <- function(dir, top_version = "keep", p_version = "keep") {
  rdata <- make_fixture_rdata(dir) # nolint: object_usage_linter.
  e <- new.env()
  load(rdata, envir = e)

  if (!identical(p_version, "keep")) {
    p <- get("p", envir = e)
    p$generator_version <- p_version
    assign("p", p, envir = e)
  }
  if (!identical(top_version, "keep")) {
    if (is.null(top_version)) rm("generator_version", envir = e) else assign("generator_version", top_version, envir = e)
  }
  rm("stimuli_params", envir = e)

  save(list = ls(e), file = rdata, envir = e)
  rdata
}

failure_message <- function(rdata) {
  tryCatch(
    {
      generateCI(
        stimuli = 1:6, responses = rep(c(1, -1), 3), baseimage = "base",
        rdata = rdata, save_as_png = FALSE
      )
      NA_character_
    },
    error = conditionMessage
  )
}

test_that("the missing-field error names the version that wrote the file", {
  skip_if_not_installed("withr")

  msg <- failure_message(broken_rdata(withr::local_tempdir()))

  expect_match(msg, "did not contain stimuli_params", fixed = TRUE)
  expect_match(msg, paste("written by rcicr", utils::packageVersion("rcicr")), fixed = TRUE)
})

test_that("p$generator_version wins over the top-level 0.4.0 sentinel", {
  skip_if_not_installed("withr")

  # A file as 1.1.0 wrote it: the top-level field says 0.4.0 whatever made it,
  # while p carried the truth all along.
  msg <- failure_message(
    broken_rdata(withr::local_tempdir(), top_version = "0.4.0", p_version = "1.1.0")
  )

  expect_match(msg, "written by rcicr 1.1.0", fixed = TRUE)
  expect_false(grepl("unknown", msg, fixed = TRUE))
})

test_that("a file that only claims 0.4.0 is reported as unknown, not as 0.4.0", {
  skip_if_not_installed("withr")

  msg <- failure_message(
    broken_rdata(withr::local_tempdir(), top_version = "0.4.0", p_version = "0.4.0")
  )

  expect_match(msg, "every version from 0.4.0 through 1.1.0", fixed = TRUE)
  expect_match(msg, "unknown", fixed = TRUE)
  # The claim must not be repeated as if it were the answer.
  expect_false(grepl("written by rcicr 0.4.0", msg, fixed = TRUE))
})

test_that("a file with no version field reports the absence, not an age", {
  skip_if_not_installed("withr")

  # This fixture is a *current* file with its version fields stripped, which is
  # the case that makes the point: an absent field is equally a truncated file, a
  # hand-rewritten one, or one rcicr never wrote, so the message may not conclude
  # the file is old.
  msg <- failure_message(
    broken_rdata(withr::local_tempdir(), top_version = NULL, p_version = NULL)
  )

  expect_match(msg, "records no writing version", fixed = TRUE)
  expect_match(msg, "what wrote this file is unknown", fixed = TRUE)
  expect_false(grepl("predates", msg, fixed = TRUE))
})

test_that("an unparseable version is reported rather than raising from the guard", {
  skip_if_not_installed("withr")

  # numeric_version() rejects this. The guard exists to explain a bad file, so
  # it must not fail on one.
  msg <- failure_message(
    broken_rdata(withr::local_tempdir(), top_version = "not-a-version", p_version = NULL)
  )

  expect_match(msg, "did not contain stimuli_params", fixed = TRUE)
  expect_match(msg, "not-a-version", fixed = TRUE)
})

test_that("computeCumulativeCICorrelation names the writing version too", {
  skip_if_not_installed("withr")

  rdata <- broken_rdata(withr::local_tempdir())

  expect_error(
    computeCumulativeCICorrelation(
      stimuli = 1:6, responses = rep(c(1, -1), 3), baseimage = "base", rdata = rdata
    ),
    paste("written by rcicr", utils::packageVersion("rcicr")),
    fixed = TRUE
  )
})
