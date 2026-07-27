test_that("batchGenerateCI computes one CI per group", {
  # dplyr::progress_estimated() is soft-deprecated (dplyr >= 1.0.0) but still
  # present and functional in the dplyr version this was verified against
  # (1.1.4) - suppressing the deprecation warning, not documenting a failure.
  tmp <- withr::local_tempdir()
  rdata_path <- make_fixture_rdata(tmp, img_size = 32, n_trials = 6, nscales = 1, seed = 1)

  set.seed(1)
  df <- data.frame(pid = rep(1:2, each = 3), stim = 1:6, resp = sample(c(1, -1), 6, replace = TRUE))

  cis <- suppressWarnings(
    batchGenerateCI(
      data = df, by = "pid", stimuli = "stim", responses = "resp",
      baseimage = "base", rdata = rdata_path, save_as_png = FALSE
    )
  )

  expect_length(cis, 2)
  for (ci in cis) {
    expect_named(ci, c("ci", "scaled", "base", "combined"))
    expect_equal(dim(ci$ci), c(32, 32))
  }
})
