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

test_that("the reference distribution is fixed by the stimulus file, not the caller's RNG state", {
  # This is what makes InfoVal reproducible across sessions: the function
  # re-generates the stimuli via generateStimuli2IFC(), which calls
  # set.seed(seed) with the seed stored in the .Rdata file
  # (R/generateStimuli2IFC.R:53). That reset governs the runif() draws in the
  # simulation loop, so the same stimulus file gives the same reference
  # distribution no matter what the ambient RNG was doing beforehand.
  #
  # Verified rather than assumed: seeding 42 vs 99 around the two calls below
  # gives byte-identical norms.
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
