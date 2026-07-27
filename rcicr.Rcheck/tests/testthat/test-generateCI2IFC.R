test_that("generateCI2IFC is a pure wrapper around generateCI", {
  tmp <- withr::local_tempdir()
  rdata_path <- make_fixture_rdata(tmp, img_size = 32, n_trials = 6, nscales = 1, seed = 1)

  set.seed(1)
  responses <- sample(c(1, -1), 6, replace = TRUE)

  r1 <- generateCI(
    stimuli = 1:6, responses = responses, baseimage = "base", rdata = rdata_path,
    save_as_png = FALSE, scaling = "independent", scaling_constant = 0.1
  )
  r2 <- generateCI2IFC(
    stimuli = 1:6, responses = responses, baseimage = "base", rdata = rdata_path,
    save_as_png = FALSE, scaling = "independent", constant = 0.1
  )

  expect_equal(r1, r2)
})
