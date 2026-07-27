test_that("batchGenerateCI2IFC computes one CI per group", {
  tmp <- withr::local_tempdir()
  rdata_path <- make_fixture_rdata(tmp, img_size = 32, n_trials = 6, nscales = 1, seed = 1)

  set.seed(1)
  df <- data.frame(pid = rep(1:2, each = 3), stim = 1:6, resp = sample(c(1, -1), 6, replace = TRUE))

  cis <- suppressWarnings(
    batchGenerateCI2IFC(
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

test_that("each CI is built from exactly its own group's trials", {
  # Same gap as in test-batchGenerateCI.R: the shape test above never checks
  # which rows fed which CI. batchGenerateCI2IFC() is a separate loop with its
  # own subsetting, so it needs its own check rather than relying on the
  # generateCI2IFC/generateCI wrapper equivalence test.
  tmp <- withr::local_tempdir()
  rdata_path <- make_fixture_rdata(tmp, img_size = 32, n_trials = 6, nscales = 1, seed = 1)

  df <- data.frame(pid = rep(1:2, each = 3), stim = 1:6, resp = c(1, 1, -1, -1, 1, 1))

  cis <- suppressWarnings(
    batchGenerateCI2IFC(
      data = df, by = "pid", stimuli = "stim", responses = "resp",
      baseimage = "base", rdata = rdata_path, save_as_png = FALSE
    )
  )

  expect_named(cis, c("base_pid_1", "base_pid_2"))

  direct1 <- generateCI2IFC(
    stimuli = df$stim[1:3], responses = df$resp[1:3], baseimage = "base",
    rdata = rdata_path, save_as_png = FALSE, scaling = "none"
  )
  direct2 <- generateCI2IFC(
    stimuli = df$stim[4:6], responses = df$resp[4:6], baseimage = "base",
    rdata = rdata_path, save_as_png = FALSE, scaling = "none"
  )

  expect_equal(cis[["base_pid_1"]]$ci, direct1$ci)
  expect_equal(cis[["base_pid_2"]]$ci, direct2$ci)
  expect_false(isTRUE(all.equal(cis[["base_pid_1"]]$ci, cis[["base_pid_2"]]$ci)))
})

test_that("grouping follows the `by` column, not row order", {
  tmp <- withr::local_tempdir()
  rdata_path <- make_fixture_rdata(tmp, img_size = 32, n_trials = 6, nscales = 1, seed = 1)

  df <- data.frame(pid = rep(1:2, times = 3), stim = 1:6, resp = c(1, 1, -1, -1, 1, 1))

  cis <- suppressWarnings(
    batchGenerateCI2IFC(
      data = df, by = "pid", stimuli = "stim", responses = "resp",
      baseimage = "base", rdata = rdata_path, save_as_png = FALSE
    )
  )

  rows <- df$pid == 1
  direct <- generateCI2IFC(
    stimuli = df$stim[rows], responses = df$resp[rows], baseimage = "base",
    rdata = rdata_path, save_as_png = FALSE, scaling = "none"
  )
  expect_equal(cis[["base_pid_1"]]$ci, direct$ci)

  positional <- generateCI2IFC(
    stimuli = df$stim[1:3], responses = df$resp[1:3], baseimage = "base",
    rdata = rdata_path, save_as_png = FALSE, scaling = "none"
  )
  expect_false(isTRUE(all.equal(cis[["base_pid_1"]]$ci, positional$ci)))
})
