test_that("generateReferenceDistribution2IFC saves a reference_norms vector of the requested length", {
  tmp <- withr::local_tempdir()
  rdata_path <- make_fixture_rdata(tmp, img_size = 32, n_trials = 6, nscales = 1, seed = 1)

  # iter is kept tiny for test speed; the function warns that iter should be
  # >= 10000 for a reliable InfoVal statistic, which is expected here.
  suppressWarnings(generateReferenceDistribution2IFC(rdata_path, iter = 3))

  e <- new.env()
  load(rdata_path, envir = e)
  expect_true(exists("reference_norms", envir = e))
  expect_length(e$reference_norms, 3)
  expect_false(anyNA(e$reference_norms))
})
