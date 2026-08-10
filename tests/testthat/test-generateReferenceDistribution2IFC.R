test_that("generateReferenceDistribution2IFC saves a reference_norms vector of the requested length", {
  tmp <- withr::local_tempdir()
  rdata_path <- make_fixture_rdata(tmp, img_size = 32, n_trials = 6, nscales = 1, seed = 1)

  # ncores now defaults to parallel::detectCores()-1. Under CRAN-style checks
  # (_R_CHECK_LIMIT_CORES_), parallel::makeCluster() errors if more than 2 cores
  # are requested - which trips on any check runner with more than 3 cores. This
  # test deliberately exercises the default rather than passing ncores, so mock
  # detectCores() to simulate a modest-core machine and keep it deterministic.
  testthat::local_mocked_bindings(detectCores = function(...) 2L, .package = "parallel")

  # iter is kept tiny for test speed; the function warns that iter should be
  # >= 10000 for a reliable InfoVal statistic, which is expected here.
  suppressWarnings(generateReferenceDistribution2IFC(rdata_path, iter = 3))

  e <- new.env()
  load(rdata_path, envir = e)
  expect_true(exists("reference_norms", envir = e))
  expect_length(e$reference_norms, 3)
  expect_false(anyNA(e$reference_norms))
})

test_that("the reference norms are positive and actually vary across iterations", {
  # Length and no-NAs (above) are satisfied by a vector of zeros, or by one
  # where every iteration produced the same number. Both would be wrong in a
  # way that matters: reference_norms is the null distribution InfoVal is
  # scored against, and computeInfoVal2IFC() divides by its mad(), so a
  # degenerate distribution yields Inf or NaN rather than an error.
  tmp <- withr::local_tempdir()
  rdata_path <- make_fixture_rdata(tmp, img_size = 32, n_trials = 6, nscales = 1, seed = 1)

  suppressWarnings(generateReferenceDistribution2IFC(rdata_path, iter = 8, ncores = 1))

  e <- new.env()
  load(rdata_path, envir = e)

  # Each entry is a Frobenius norm of a non-degenerate CI, so strictly positive.
  expect_true(all(e$reference_norms > 0))

  # Responses are redrawn inside the loop, so the iterations must not collapse
  # to one repeated value -- which is what drawing them once outside the loop
  # would produce.
  expect_gt(length(unique(e$reference_norms)), 1)
  expect_gt(mad(e$reference_norms), 0)
})

test_that("an .Rdata file predating noise_type still works, and says so", {
  # Issue #94. Old files have no noise_type, and re-generating the stimuli from
  # one failed outright with "object 'noise_type' not found". The workaround on
  # record was to load the file and assign noise_type by hand.
  tmp <- withr::local_tempdir()
  rdata_path <- make_fixture_rdata(tmp, img_size = 32, n_trials = 6, nscales = 1, seed = 1)

  # Age the fixture by removing the field, which is how a pre-0.3.x file looks.
  e <- new.env()
  load(rdata_path, envir = e)
  rm("noise_type", envir = e)
  save(list = ls(e), file = rdata_path, envir = e)

  # Collected by hand rather than with expect_warning(): the iter < 10000
  # warning fires too, and this needs to assert on one specific warning without
  # swallowing or being tripped by the other.
  collect_warnings <- function(expr) {
    seen <- character()
    withCallingHandlers(
      expr,
      warning = function(cond) {
        seen <<- c(seen, conditionMessage(cond))
        invokeRestart("muffleWarning")
      }
    )
    seen
  }

  warnings_seen <- collect_warnings(
    generateReferenceDistribution2IFC(rdata_path, iter = 3, ncores = 1))

  expect_true(any(grepl("does not contain `noise_type`", warnings_seen, fixed = TRUE)))

  # It has to actually finish, not just warn on the way to the old error.
  after <- new.env()
  load(rdata_path, envir = after)
  expect_length(after$reference_norms, 3)
  expect_false(anyNA(after$reference_norms))

  # And the warning must be specific to the missing field, or it is just noise
  # on every call.
  fresh <- make_fixture_rdata(withr::local_tempdir(), img_size = 32, n_trials = 6,
                              nscales = 1, seed = 1)
  expect_false(any(grepl(
    "does not contain `noise_type`",
    collect_warnings(generateReferenceDistribution2IFC(fresh, iter = 3, ncores = 1)),
    fixed = TRUE)))
})

test_that("the reference distribution is fixed by the stimulus file, not the caller's RNG state", {
  # This is what makes InfoVal reproducible across sessions: the function
  # re-generates the stimuli via generateStimuli2IFC(), which calls
  # set.seed(seed) with the seed stored in the .Rdata file
  # (R/generateStimuli2IFC.R:53). That reset governs the runif() draws in the
  # simulation loop, so the same stimulus file gives the same reference
  # distribution no matter what the ambient RNG was doing beforehand.
  #
  # Verified rather than assumed: seeding 42 vs 99 around the two calls below
  # gives byte-identical norms, and ncores = 1 vs 2 also agree, so the null is
  # portable across machines.
  #
  # This determinism began as a *side effect* of the stimulus rebuild rather
  # than a designed guarantee. That was resolved in 1.2.0: it is now
  # documented as a guarantee under ?generateReferenceDistribution2IFC, the
  # set.seed() call it rests on carries a comment saying so, and the deliberate
  # way to vary the null is the response_seed argument tested below. This test
  # is the guarantee's guard -- removing or moving that set.seed() would
  # silently change every InfoVal ever computed.
  tmp <- withr::local_tempdir()
  rdata_path <- make_fixture_rdata(tmp, img_size = 32, n_trials = 6, nscales = 1, seed = 1)

  read_norms <- function(seed) {
    withr::with_seed(seed, {
      suppressWarnings(generateReferenceDistribution2IFC(rdata_path, iter = 8, ncores = 1))
    })
    e <- new.env()
    load(rdata_path, envir = e)
    e$reference_norms
  }

  expect_equal(read_norms(42), read_norms(99))
})

test_that("response_seed varies the null reproducibly, and NULL leaves the default alone", {
  # The point of the argument is to let a user draw an independent null from the
  # same stimuli -- e.g. to measure how much Monte Carlo error a given iter
  # leaves in their InfoVal. Two claims have to hold together: the same seed
  # must reproduce, and different seeds must actually differ. Without the
  # second, the first is satisfied by an argument that is ignored entirely.
  tmp <- withr::local_tempdir()
  pristine <- make_fixture_rdata(tmp, img_size = 32, n_trials = 6, nscales = 1, seed = 1)

  # Each run gets its own copy, since a run with save_rdata = TRUE rewrites the file.
  norms_for <- function(tag, ...) {
    path <- file.path(tmp, paste0(tag, ".Rdata"))
    file.copy(pristine, path, overwrite = TRUE)
    suppressWarnings(generateReferenceDistribution2IFC(path, iter = 8, ncores = 1, ...))
  }

  expect_equal(norms_for("s99a", response_seed = 99), norms_for("s99b", response_seed = 99))
  expect_false(identical(norms_for("s99a2", response_seed = 99),
                         norms_for("s1234", response_seed = 1234)))

  # And a seeded null is not just the default null relabelled.
  expect_false(identical(norms_for("seeded", response_seed = 99), norms_for("default")))
})

test_that("response_seed reaches the responses, and never the stimulus rebuild", {
  # This is the whole semantic distinction, and the reason the seed is not
  # simply forwarded to generateStimuli2IFC(): handing that function a different
  # seed would rebuild a *different stimulus set*, so the null would describe
  # stimuli the participants never saw.
  #
  # It has to be tested on the call itself. Comparing p/stimuli_params in the
  # .Rdata before and after is vacuous -- the rebuild runs with
  # save_rdata = FALSE and never writes stimuli back, so the file is unchanged
  # whichever seed the rebuild received. A mutant forwarding response_seed to
  # generateStimuli2IFC() passed that version of this test.
  tmp <- withr::local_tempdir()
  rdata_path <- make_fixture_rdata(tmp, img_size = 32, n_trials = 6, nscales = 1, seed = 1)

  # Grab the real function before mocking, so the mock can delegate to it.
  real_generate <- get("generateStimuli2IFC", envir = asNamespace("rcicr"))
  seen_seed <- NULL

  testthat::local_mocked_bindings(
    generateStimuli2IFC = function(..., seed) {
      seen_seed <<- seed
      real_generate(..., seed = seed)
    },
    .package = "rcicr"
  )

  suppressWarnings(generateReferenceDistribution2IFC(
    rdata_path, iter = 8, ncores = 1, response_seed = 99))

  # The stimulus seed stored in the file, not the response seed we passed.
  stored <- new.env()
  load(rdata_path, envir = stored)
  expect_equal(seen_seed, stored$seed)
  expect_false(identical(seen_seed, 99))
})

test_that("save_rdata = FALSE returns the norms without touching the file", {
  # Needed so a one-off check of the null cannot become the stimulus set's
  # permanent reference distribution. The norms are only reachable through the
  # return value in that case, which is why the function returns them at all.
  tmp <- withr::local_tempdir()
  rdata_path <- make_fixture_rdata(tmp, img_size = 32, n_trials = 6, nscales = 1, seed = 1)

  norms <- suppressWarnings(generateReferenceDistribution2IFC(
    rdata_path, iter = 8, ncores = 1, response_seed = 99, save_rdata = FALSE))

  expect_length(norms, 8)
  expect_true(all(norms > 0))

  e <- new.env()
  load(rdata_path, envir = e)
  expect_false(exists("reference_norms", envir = e, inherits = FALSE))

  # Discriminating case: the same call with save_rdata = TRUE does write it,
  # so the assertion above is about the flag and not about the call failing.
  suppressWarnings(generateReferenceDistribution2IFC(
    rdata_path, iter = 8, ncores = 1, response_seed = 99, save_rdata = TRUE))
  e2 <- new.env()
  load(rdata_path, envir = e2)
  expect_equal(e2$reference_norms, norms)
})

test_that("the saved file records which seed produced its norms, and no arguments", {
  # Provenance: a file carrying a deliberately varied null must be
  # distinguishable from one carrying the default, otherwise the two are
  # indistinguishable to anyone who picks the file up later.
  #
  # The second half guards the clobbering trap documented at the top of
  # generateReferenceDistribution2IFC(): an argument left in the frame at save
  # time becomes an *input* to the next call's load(). That is what happened to
  # `rdata` and `ncores`, and is why the seed is recorded under a different name
  # (reference_norms_seed) than the argument that sets it (response_seed).
  tmp <- withr::local_tempdir()
  pristine <- make_fixture_rdata(tmp, img_size = 32, n_trials = 6, nscales = 1, seed = 1)

  read_after <- function(tag, ...) {
    path <- file.path(tmp, paste0(tag, ".Rdata"))
    file.copy(pristine, path, overwrite = TRUE)
    suppressWarnings(generateReferenceDistribution2IFC(path, iter = 8, ncores = 1, ...))
    e <- new.env()
    load(path, envir = e)
    e
  }

  default_run <- read_after("default")
  expect_true(exists("reference_norms_seed", envir = default_run, inherits = FALSE))
  expect_null(default_run$reference_norms_seed)

  seeded_run <- read_after("seeded", response_seed = 99)
  expect_equal(seeded_run$reference_norms_seed, 99)

  expect_false(any(c("response_seed", "save_rdata", "rdata", "ncores", "iter") %in%
                     ls(seeded_run, all.names = TRUE)))
})
