# Every function that writes files takes its destination as a required
# argument. The defaults these replace ('./stimuli', './cis', './zmaps',
# 'zmaps') wrote into whatever the working directory happened to be, which CRAN
# policy does not allow.
#
# Each test runs from a temporary working directory, so a guard that failed to
# fire would create a real directory here rather than in the checkout.

test_that("generateStimuli2IFC requires stimulus_path when it writes", {
  skip_if_not_installed("withr")
  dir <- withr::local_tempdir()
  withr::local_dir(dir)
  png_path <- file.path(dir, "base.png")
  make_square_png(png_path)

  expect_error(
    generateStimuli2IFC(base_face_files = list(base = png_path), n_trials = 2,
                        img_size = 32, seed = 1, ncores = 1, nscales = 1),
    "No stimulus_path"
  )
})

test_that("generateStimuli2IFC needs no path when it writes nothing", {
  skip_if_not_installed("withr")
  dir <- withr::local_tempdir()
  withr::local_dir(dir)
  png_path <- file.path(dir, "base.png")
  make_square_png(png_path)

  noise <- suppressWarnings(
    generateStimuli2IFC(base_face_files = list(base = png_path), n_trials = 2,
                        img_size = 32, seed = 1, ncores = 1, nscales = 1,
                        return_as_dataframe = TRUE, save_as_png = FALSE,
                        save_rdata = FALSE)
  )

  expect_equal(ncol(noise), 2) # one column of noise per trial
  # This path used to create ./stimuli anyway.
  expect_false(dir.exists(file.path(dir, "stimuli")))
})

test_that("generateReferenceDistribution2IFC leaves no stray stimuli directory", {
  skip_if_not_installed("withr")
  dir <- withr::local_tempdir()
  rdata <- make_fixture_rdata(dir)

  run_dir <- withr::local_tempdir()
  withr::local_dir(run_dir)

  suppressWarnings(
    generateReferenceDistribution2IFC(rdata, iter = 3, ncores = 1)
  )

  # It re-generates the stimuli in memory only, so it must write nothing to the
  # working directory.
  expect_equal(list.files(run_dir), character(0))
})

test_that("generateCI requires targetpath when it saves a PNG", {
  skip_if_not_installed("withr")
  dir <- withr::local_tempdir()
  rdata <- make_fixture_rdata(dir)
  withr::local_dir(withr::local_tempdir())

  expect_error(
    generateCI(stimuli = 1:6, responses = rep(c(1, -1), 3),
               baseimage = "base", rdata = rdata),
    "no targetpath"
  )

  expect_error(
    generateCI(stimuli = 1:6, responses = rep(c(1, -1), 3),
               baseimage = "base", rdata = rdata,
               save_as_png = FALSE, save_individual_cis = TRUE,
               participants = rep(c("p1", "p2"), each = 3)),
    "no targetpath"
  )
})

test_that("generateCI requires zmaptargetpath when it draws a z-map", {
  skip_if_not_installed("withr")
  dir <- withr::local_tempdir()
  rdata <- make_fixture_rdata(dir)
  withr::local_dir(withr::local_tempdir())

  expect_error(
    generateCI(stimuli = 1:6, responses = rep(c(1, -1), 3),
               baseimage = "base", rdata = rdata, save_as_png = FALSE,
               zmap = TRUE, n_cores = 1),
    "no zmaptargetpath"
  )
})

test_that("generateCI computes without a path when it writes nothing", {
  skip_if_not_installed("withr")
  dir <- withr::local_tempdir()
  rdata <- make_fixture_rdata(dir)
  run_dir <- withr::local_tempdir()
  withr::local_dir(run_dir)

  ci <- generateCI(stimuli = 1:6, responses = rep(c(1, -1), 3),
                   baseimage = "base", rdata = rdata, save_as_png = FALSE)

  expect_equal(dim(ci$ci), c(32, 32))
  expect_equal(list.files(run_dir), character(0))
})

test_that("generateCI2IFC requires targetpath when it saves a PNG", {
  skip_if_not_installed("withr")
  dir <- withr::local_tempdir()
  rdata <- make_fixture_rdata(dir)
  withr::local_dir(withr::local_tempdir())

  expect_error(
    generateCI2IFC(stimuli = 1:6, responses = rep(c(1, -1), 3),
                   baseimage = "base", rdata = rdata),
    "no targetpath"
  )
})

test_that("the batch functions require targetpath when they save PNGs", {
  skip_if_not_installed("withr")
  dir <- withr::local_tempdir()
  rdata <- make_fixture_rdata(dir)
  withr::local_dir(withr::local_tempdir())

  df <- data.frame(pid = rep(c("p1", "p2"), each = 3), stim = 1:6,
                   resp = rep(c(1, -1), 3), stringsAsFactors = FALSE)

  # The error must name the function the user called, not one it delegates to.
  expect_error(
    batchGenerateCI(data = df, by = "pid", stimuli = "stim",
                    responses = "resp", baseimage = "base", rdata = rdata),
    "no targetpath"
  )
  expect_error(
    batchGenerateCI2IFC(data = df, by = "pid", stimuli = "stim",
                        responses = "resp", baseimage = "base", rdata = rdata),
    "no targetpath"
  )
})

test_that("the batch functions run without a path when they write nothing", {
  skip_if_not_installed("withr")
  dir <- withr::local_tempdir()
  rdata <- make_fixture_rdata(dir)
  run_dir <- withr::local_tempdir()
  withr::local_dir(run_dir)

  df <- data.frame(pid = rep(c("p1", "p2"), each = 3), stim = 1:6,
                   resp = rep(c(1, -1), 3), stringsAsFactors = FALSE)

  # scaling defaults to 'autoscale', so this also exercises autoscale()'s own
  # guard through the pass-through.
  cis <- suppressWarnings(
    batchGenerateCI(data = df, by = "pid", stimuli = "stim",
                    responses = "resp", baseimage = "base", rdata = rdata,
                    save_as_png = FALSE)
  )

  expect_length(cis, 2)
  expect_equal(list.files(run_dir), character(0))
})

test_that("autoscale requires targetpath when it saves PNGs", {
  skip_if_not_installed("withr")
  withr::local_dir(withr::local_tempdir())

  cis <- list(
    a = list(ci = matrix(runif(64, -0.2, 0.2), 8, 8), base = matrix(0.5, 8, 8)),
    b = list(ci = matrix(runif(64, -0.3, 0.3), 8, 8), base = matrix(0.5, 8, 8))
  )

  expect_error(autoscale(cis), "no targetpath")

  # Without save_as_pngs there is nothing to write, so no path is needed.
  # autoscale() reports its constant on stdout, hence capture.output().
  scaled <- NULL
  capture.output(scaled <- autoscale(cis, save_as_pngs = FALSE))
  expect_named(scaled, c("a", "b"))
})

# Backwards compatibility ---------------------------------------------------
#
# Every .Rdata file this package has ever written records the stimulus_path it
# was generated into -- typically the old default, './stimuli', relative to a
# working directory that no longer exists. Removing the path defaults must not
# make those files unreadable, and the recorded path must never be used as a
# destination: load() assigns straight into the reading function's frame.

test_that("an .Rdata file recording a stale stimulus_path still computes a CI", {
  skip_if_not_installed("withr")
  dir <- withr::local_tempdir()
  rdata <- make_fixture_rdata(dir)

  # What an old file looks like: generated into './stimuli' on someone else's
  # machine, years ago.
  mutate_rdata(rdata, stimulus_path = "./stimuli")

  run_dir <- withr::local_tempdir()
  withr::local_dir(run_dir)

  ci <- generateCI(stimuli = 1:6, responses = rep(c(1, -1), 3),
                   baseimage = "base", rdata = rdata, save_as_png = FALSE)

  expect_equal(dim(ci$ci), c(32, 32))
  expect_equal(list.files(run_dir), character(0))
})

test_that("a path stored in the .Rdata never overrides the one passed in", {
  skip_if_not_installed("withr")
  dir <- withr::local_tempdir()
  rdata <- make_fixture_rdata(dir)

  decoy <- withr::local_tempdir()
  wanted <- file.path(withr::local_tempdir(), "cis")
  wanted_zmaps <- file.path(withr::local_tempdir(), "zmaps")

  # No file this package writes contains these names, but load() would let one
  # replace the caller's argument if it did -- the hazard that made every z-map
  # since 1.1.0 use the wrong sigma. captureArgs() is what prevents it.
  mutate_rdata(rdata, stimulus_path = decoy, targetpath = decoy,
               zmaptargetpath = decoy)

  generateCI(stimuli = 1:6, responses = rep(c(1, -1), 3), baseimage = "base",
             rdata = rdata, targetpath = wanted, zmap = TRUE,
             zmaptargetpath = wanted_zmaps, zmapdecoration = FALSE, n_cores = 1)

  expect_length(list.files(wanted, pattern = "\\.png$"), 1)
  expect_length(list.files(wanted_zmaps, pattern = "\\.png$"), 1)
  expect_equal(list.files(decoy), character(0))
})

test_that("generateReferenceDistribution2IFC ignores the recorded stimulus_path", {
  skip_if_not_installed("withr")
  dir <- withr::local_tempdir()
  rdata <- make_fixture_rdata(dir)

  decoy <- withr::local_tempdir()
  mutate_rdata(rdata, stimulus_path = decoy)

  run_dir <- withr::local_tempdir()
  withr::local_dir(run_dir)

  suppressWarnings(
    generateReferenceDistribution2IFC(rdata, iter = 3, ncores = 1)
  )

  expect_equal(list.files(decoy), character(0))
  expect_equal(list.files(run_dir), character(0))
})

test_that("plotZmap requires targetpath", {
  skip_if_not_installed("withr")
  withr::local_dir(withr::local_tempdir())

  zmap <- withr::with_seed(1, matrix(rnorm(64, sd = 5), 8, 8))

  expect_error(
    plotZmap(zmap, sigma = 3, threshold = 3, decoration = FALSE, size = 200),
    "No targetpath"
  )
})
