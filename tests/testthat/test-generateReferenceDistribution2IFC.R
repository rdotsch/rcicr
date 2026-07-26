test_that("generateReferenceDistribution2IFC saves a reference_norms vector of the requested length", {
  tmp <- withr::local_tempdir()
  rdata_path <- make_fixture_rdata(tmp, img_size = 32, n_trials = 6, nscales = 1, seed = 1)

  # generateReferenceDistribution2IFC() hardcodes ncores=parallel::detectCores()-1
  # for its internal generateStimuli2IFC() call with no way to override it. Under
  # CRAN-style checks (_R_CHECK_LIMIT_CORES_), parallel::makeCluster() errors if
  # more than 2 cores are requested - which trips on any check runner with more
  # than 3 cores. Mocking detectCores() simulates a modest-core machine so this
  # test is deterministic regardless of how many cores the real machine has.
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
