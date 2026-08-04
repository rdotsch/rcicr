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

test_that("each CI is built from exactly its own group's trials", {
  # The shape test above passes for any function returning two 32x32 CIs -- it
  # never checks that group 1's CI came from group 1's rows. This does: each
  # element must equal a direct generateCI() on that subset alone. Compare $ci,
  # not $scaled, because the default scaling = 'autoscale' rewrites $scaled
  # across the batch; $ci is the field carrying the published numbers and is
  # untouched by scaling (see test-scaling.R).
  tmp <- withr::local_tempdir()
  rdata_path <- make_fixture_rdata(tmp, img_size = 32, n_trials = 6, nscales = 1, seed = 1)

  # Deliberately asymmetric: the two groups differ in both their stimuli and
  # their response patterns, so a CI built from the wrong group's rows cannot
  # coincide with the right answer.
  df <- data.frame(pid = rep(1:2, each = 3), stim = 1:6, resp = c(1, 1, -1, -1, 1, 1))

  cis <- suppressWarnings(
    batchGenerateCI(
      data = df, by = "pid", stimuli = "stim", responses = "resp",
      baseimage = "base", rdata = rdata_path, save_as_png = FALSE
    )
  )

  # The list is keyed <baseimage>_<by>_<unit>, which is what identifies whose CI
  # is whose once these are written out or passed to autoscale().
  expect_named(cis, c("base_pid_1", "base_pid_2"))

  direct1 <- generateCI(
    stimuli = df$stim[1:3], responses = df$resp[1:3], baseimage = "base",
    rdata = rdata_path, save_as_png = FALSE, scaling = "none"
  )
  direct2 <- generateCI(
    stimuli = df$stim[4:6], responses = df$resp[4:6], baseimage = "base",
    rdata = rdata_path, save_as_png = FALSE, scaling = "none"
  )

  expect_equal(cis[["base_pid_1"]]$ci, direct1$ci)
  expect_equal(cis[["base_pid_2"]]$ci, direct2$ci)

  # Guards against a subsetting bug that hands every group the full data frame:
  # that would make both CIs identical while still satisfying the shape test.
  expect_false(isTRUE(all.equal(cis[["base_pid_1"]]$ci, cis[["base_pid_2"]]$ci)))
})

test_that("grouping follows the `by` column, not row order", {
  # With interleaved group labels, splitting the data positionally (first half /
  # second half) gives a different answer than grouping by `pid`. A contiguous
  # fixture cannot tell those two implementations apart.
  tmp <- withr::local_tempdir()
  rdata_path <- make_fixture_rdata(tmp, img_size = 32, n_trials = 6, nscales = 1, seed = 1)

  df <- data.frame(pid = rep(1:2, times = 3), stim = 1:6, resp = c(1, 1, -1, -1, 1, 1))

  cis <- suppressWarnings(
    batchGenerateCI(
      data = df, by = "pid", stimuli = "stim", responses = "resp",
      baseimage = "base", rdata = rdata_path, save_as_png = FALSE
    )
  )

  rows <- df$pid == 1
  direct <- generateCI(
    stimuli = df$stim[rows], responses = df$resp[rows], baseimage = "base",
    rdata = rdata_path, save_as_png = FALSE, scaling = "none"
  )
  expect_equal(cis[["base_pid_1"]]$ci, direct$ci)

  # ...and the positional split gives a genuinely different CI, so the
  # assertion above is not vacuous.
  positional <- generateCI(
    stimuli = df$stim[1:3], responses = df$resp[1:3], baseimage = "base",
    rdata = rdata_path, save_as_png = FALSE, scaling = "none"
  )
  expect_false(isTRUE(all.equal(cis[["base_pid_1"]]$ci, positional$ci)))
})

test_that("batchGenerateCI ignores rows without a grouping value", {
  tmp <- withr::local_tempdir()
  rdata_path <- make_fixture_rdata(tmp, img_size = 32, n_trials = 6, nscales = 1, seed = 1)

  df <- data.frame(
    pid = c(1, 1, NA, 2, 2, NA),
    stim = 1:6,
    resp = c(1, -1, 1, -1, 1, -1)
  )

  cis <- suppressWarnings(
    batchGenerateCI(
      data = df, by = "pid", stimuli = "stim", responses = "resp",
      baseimage = "base", rdata = rdata_path, save_as_png = FALSE
    )
  )

  expect_named(cis, c("base_pid_1", "base_pid_2"))

  rows <- !is.na(df$pid) & df$pid == 1
  direct <- generateCI(
    stimuli = df$stim[rows], responses = df$resp[rows], baseimage = "base",
    rdata = rdata_path, save_as_png = FALSE, scaling = "none"
  )
  expect_equal(cis[["base_pid_1"]]$ci, direct$ci)
})
