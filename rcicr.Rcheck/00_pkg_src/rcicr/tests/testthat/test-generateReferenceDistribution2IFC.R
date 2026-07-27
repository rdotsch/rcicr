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
