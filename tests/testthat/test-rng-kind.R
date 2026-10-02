# The stimulus file records the RNG kind its seed belongs to, and references
# replay under it (#315); generateStimuli2IFC() gives the caller's stream back
# (#189).

local_rng_kind <- function(kind, envir = parent.frame()) {
  previous <- RNGkind()
  withr::defer(suppressWarnings(do.call(RNGkind, as.list(previous))), envir = envir)
  RNGkind(kind)
}

saved_field <- function(rdata, name) {
  e <- new.env()
  load(rdata, envir = e)
  get0(name, envir = e, inherits = FALSE)
}

reference <- function(rdata, ...) {
  norms <- NULL
  utils::capture.output(norms <- suppressMessages(suppressWarnings(
    generateReferenceDistribution2IFC(rdata, iter = 20, ncores = 1, save_rdata = FALSE, ...)
  )))
  norms
}

test_that("the stimulus file records the RNG kind in force when it was generated", {
  for (kind in c("Mersenne-Twister", "L'Ecuyer-CMRG")) {
    local_rng_kind(kind)
    expected <- RNGkind()
    rdata <- make_fixture_rdata(withr::local_tempdir())
    expect_identical(saved_field(rdata, "rng_kind"), expected, info = kind)
  }
})

test_that("a reference is the same whichever kind the scoring session runs", {
  local_rng_kind("L'Ecuyer-CMRG")
  rdata <- make_fixture_rdata(withr::local_tempdir(), n_trials = 8)
  own <- list(default = reference(rdata), seeded = reference(rdata, response_seed = 3),
              images = reference(rdata, reference_method = "images"))

  RNGkind("Mersenne-Twister")
  expect_identical(reference(rdata), own$default)
  expect_identical(reference(rdata, response_seed = 3), own$seeded)
  expect_identical(reference(rdata, reference_method = "images"), own$images)

  # Without the record, the session's kind decides, as it did for every older file.
  mutate_rdata(rdata, .remove = "rng_kind")
  expect_false(identical(reference(rdata), own$default))
  expect_false(identical(reference(rdata, response_seed = 3), own$seeded))
})

test_that("subset and independent-base references replay under the recorded kind too", {
  local_rng_kind("L'Ecuyer-CMRG")
  rdata <- make_fixture_rdata(withr::local_tempdir(), n_trials = 8)
  independent <- make_independent_fixture(withr::local_tempdir())
  subset <- reference(rdata, reference_stimuli = 1:5)
  by_base <- reference(independent, baseimage = "second")

  RNGkind("Mersenne-Twister")
  expect_identical(reference(rdata, reference_stimuli = 1:5), subset)
  expect_identical(reference(independent, baseimage = "second"), by_base)
})

test_that("a cross-kind reference leaves the caller's kind and position alone", {
  local_rng_kind("L'Ecuyer-CMRG")
  rdata <- make_fixture_rdata(withr::local_tempdir())
  RNGkind("Mersenne-Twister")

  set.seed(99)
  expected <- runif(3)
  set.seed(99)
  reference(rdata)
  expect_identical(RNGkind()[1], "Mersenne-Twister")
  expect_identical(runif(3), expected)

  # With no stream to return to, removing the replay's seed is not enough: the
  # session's kind would stay switched.
  rm(".Random.seed", envir = globalenv())
  reference(rdata)
  expect_identical(RNGkind()[1], "Mersenne-Twister")
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))

  # And on an error during the draws.
  set.seed(99)
  local_mocked_bindings(gramNorms = function(gram, iter, ...) {
    runif(5)
    stop("draw failed")
  })
  expect_error(reference(rdata), "draw failed")
  expect_identical(RNGkind()[1], "Mersenne-Twister")
  expect_identical(runif(3), expected)
})

test_that("under the recorded kind, a requested reference still moves the stream", {
  # The documented contract for a same-kind regeneration is unchanged.
  rdata <- make_fixture_rdata(withr::local_tempdir())
  set.seed(99)
  expected <- runif(3)
  set.seed(99)
  reference(rdata)
  expect_false(identical(runif(3), expected))
})

test_that("generateStimuli2IFC() leaves the caller's stream where it was", {
  set.seed(42)
  expected <- runif(3)
  set.seed(42)
  make_fixture_rdata(withr::local_tempdir())
  expect_identical(runif(3), expected)

  rm(".Random.seed", envir = globalenv())
  make_fixture_rdata(withr::local_tempdir())
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))

  # And when the call fails after drawing everything.
  set.seed(42)
  local_mocked_bindings(saveStimulusFile = function(...) stop("save failed"))
  expect_error(make_fixture_rdata(withr::local_tempdir()), "save failed")
  expect_identical(runif(3), expected)
})

test_that("an interrupted generateStimuli2IFC() also gives the stream back", {
  interrupt <- function() signalCondition(structure(class = c("interrupt", "condition"), list()))
  set.seed(42)
  expected <- runif(3)
  set.seed(42)
  local_mocked_bindings(saveStimulusFile = function(...) interrupt())
  expect_identical(tryCatch(make_fixture_rdata(withr::local_tempdir()),
                            interrupt = function(e) "aborted"), "aborted")
  expect_identical(runif(3), expected)
})
