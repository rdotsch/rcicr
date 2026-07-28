# Regression tests for bugs that HAVE been fixed.
#
# These began life in test-known-bugs.R as deliberately failing tests - the
# executable form of the P0 list in BACKLOG.md. Every one of them now passes,
# so the file was renamed: its job from here on is to stop these bugs coming
# back.
#
# Each test asserts intended behaviour, never the buggy output it replaced.
# Do not "fix" a failure here by asserting the broken result - that locks the
# bug back in. A failure means a regression.

test_that("generateStimuli2IFC saves nscales and sigma in the .Rdata file", {
  # BACKLOG item 2 / issue #81. Because nscales is not saved,
  # generateReferenceDistribution2IFC() silently rebuilds the InfoVal null
  # distribution on the default nscales=5 basis, so InfoVal values reported
  # from non-default nscales are computed against the wrong distribution.
  tmp <- withr::local_tempdir()
  rdata_path <- make_fixture_rdata(tmp, img_size = 32, n_trials = 6, nscales = 1, seed = 1)

  e <- new.env()
  load(rdata_path, envir = e)

  expect_true("nscales" %in% ls(e))
  expect_equal(e$nscales, 1)
  expect_true("sigma" %in% ls(e))
})

test_that("computeInfoVal2IFC(force_gen_ref_dist = TRUE) regenerates the reference distribution", {
  # BACKLOG item 3 / issue #113. The flag short-circuits the lookup-table
  # branch but never reaches the regeneration branch, because reference_norms
  # still exists() after load(). The user is given no indication it was ignored.
  tmp <- withr::local_tempdir()
  rdata_path <- make_fixture_rdata(tmp, img_size = 32, n_trials = 6, nscales = 1, seed = 1)
  seed_reference_norms(rdata_path, n = 5, seed = 1)

  testthat::local_mocked_bindings(detectCores = function(...) 2L, .package = "parallel")

  e <- new.env()
  load(rdata_path, envir = e)
  before <- e$reference_norms

  target_ci <- list(ci = matrix(withr::with_seed(2, rnorm(32 * 32)), 32, 32))
  suppressWarnings(invisible(utils::capture.output(
    computeInfoVal2IFC(target_ci, rdata_path, iter = 3, force_gen_ref_dist = TRUE)
  )))

  e2 <- new.env()
  load(rdata_path, envir = e2)
  expect_false(identical(before, e2$reference_norms))
})

test_that("generateCI accepts tibble columns as well as data.frame columns", {
  # BACKLOG item 4 / issues #70 and #123. tbl[, "col"] stays a 1-column tibble
  # where df[, "col"] drops to a vector, so aggregate() fails with
  # "arguments must have same length". Since readr/dplyr return tibbles by
  # default this is now the normal path for a modern user.
  tmp <- withr::local_tempdir()
  rdata_path <- make_fixture_rdata(tmp, img_size = 32, n_trials = 6, nscales = 1, seed = 1)

  set.seed(1)
  tb <- tibble::tibble(stimulus = 1:6, response = sample(c(1, -1), 6, replace = TRUE))

  ci <- generateCI(
    stimuli = tb[, "stimulus"], responses = tb[, "response"],
    baseimage = "base", rdata = rdata_path, save_as_png = FALSE
  )

  expect_equal(dim(ci$ci), c(32, 32))
  expect_false(anyNA(ci$ci))
})

test_that("generateCI applies a matrix mask", {
  # BACKLOG item 6 (no issue filed). generateCI() branches on
  # `if (!is.na(mask))`, which gets a matrix condition - a hard error since
  # R 4.2. A second bug sits behind it: applyMask() checks a hardcoded 512
  # rather than img_size, so masks fail for any other stimulus size.
  tmp <- withr::local_tempdir()
  rdata_path <- make_fixture_rdata(tmp, img_size = 32, n_trials = 6, nscales = 1, seed = 1)

  mask <- matrix(1, 32, 32)
  mask[1:10, 1:10] <- 0 # black (0) = masked

  set.seed(1)
  ci <- generateCI(
    stimuli = 1:6, responses = sample(c(1, -1), 6, replace = TRUE),
    baseimage = "base", rdata = rdata_path, save_as_png = FALSE, mask = mask
  )

  expect_true(all(is.na(ci$ci[1:10, 1:10])))
  expect_false(anyNA(ci$ci[11:32, 11:32]))
})

test_that("generateStimuli2IFC gives an informative error when the base image is not img_size", {
  # BACKLOG item 7 / issue #124. The base image resize step is commented out
  # (left behind when the biOps dependency was dropped), and only squareness
  # is validated - so a size mismatch surfaces from inside a foreach worker as
  # "non-conformable arrays", naming neither the file nor the sizes.
  tmp <- withr::local_tempdir()
  base_face <- file.path(tmp, "face64.png")
  make_square_png(base_face, size = 64)

  expect_error(
    generateStimuli2IFC(
      base_face_files = list(face = base_face),
      n_trials = 2, img_size = 32, stimulus_path = tmp,
      seed = 1, ncores = 1, nscales = 1, save_as_png = FALSE, save_rdata = FALSE
    ),
    regexp = "64|32|size|dimension",
    info = "the error should name the offending dimensions, not just 'non-conformable arrays'"
  )
})

test_that("generateNoiseImage supports pre-0.3.3 sinusoids/sinIdx noise patterns", {
  # BACKLOG item 8 (no issue filed). The length check reads p$patchIdx before
  # the block that renames sinusoids/sinIdx to patches/patchIdx, so the
  # backward compatibility the code clearly intends never actually works.
  p <- generateNoisePattern(16, nscales = 1)
  p_legacy <- list(sinusoids = p$patches, sinIdx = p$patchIdx)
  params <- rep(1, max(p$patchIdx))

  expect_equal(generateNoiseImage(params, p_legacy), generateNoiseImage(params, p))
})

test_that("simulateNoiseIntensities returns a matrix of noise intensity ranges", {
  # BACKLOG item 5 (no issue filed). Sizes its progress bar with `data[, by]`,
  # but neither is a parameter of the function, so `data` resolves to
  # utils::data and it errors every time. It also ignores img_size, hardcoding
  # generateNoisePattern(img_size = 512).
  # simulateNoiseIntensities() draws a boxplot as a side effect. Without an
  # open device that writes an Rplots.pdf into the working directory, which
  # R CMD check then flags as leftover detritus.
  grDevices::pdf(file = file.path(withr::local_tempdir(), "plot.pdf"))
  withr::defer(grDevices::dev.off())

  result <- simulateNoiseIntensities(nrep = 2, img_size = 32)

  expect_true(is.matrix(result))
  expect_equal(dim(result), c(2, 2))
  expect_true(all(result[, 1] <= result[, 2])) # min <= max per replication
})

test_that("generateReferenceDistribution2IFC does not write its own arguments into the .Rdata file", {
  tmp <- withr::local_tempdir()
  rdata_path <- make_fixture_rdata(tmp, img_size = 32, n_trials = 6, nscales = 1, seed = 1)

  # This function re-generates stimuli without forwarding a stimulus_path, so
  # it always creates ./stimuli relative to the working directory.
  withr::local_dir(tmp)

  suppressWarnings(generateReferenceDistribution2IFC(rdata_path, iter = 3, ncores = 1))

  e <- new.env()
  load(rdata_path, envir = e)

  # load() assigns straight into the caller's frame, so an `rdata` or `ncores`
  # object stored in the file overwrites the argument of the same name on the
  # next call: the caller's ncores would be silently ignored and the save
  # redirected to whatever path the file happened to record.
  expect_false(exists("rdata", envir = e, inherits = FALSE))
  expect_false(exists("ncores", envir = e, inherits = FALSE))
  expect_false(exists("i", envir = e, inherits = FALSE))

  # ...while everything the analysis functions actually need is still there.
  for (nm in c("base_faces", "stimuli_params", "p", "img_size", "seed",
               "n_trials", "nscales", "sigma", "reference_norms")) {
    expect_true(exists(nm, envir = e, inherits = FALSE), info = nm)
  }
})

test_that("computeInfoVal2IFC ignores a stale `rdata` path stored inside the .Rdata file", {
  tmp <- withr::local_tempdir()
  rdata_path <- make_fixture_rdata(tmp, img_size = 32, n_trials = 6, nscales = 1, seed = 1)

  # Simulate a file written by an older generateReferenceDistribution2IFC(),
  # which saved its own `rdata` argument into the file. Here that recorded path
  # no longer exists, so if load() is allowed to overwrite the argument this
  # function was called with, the regeneration below reads a missing file.
  e <- new.env()
  load(rdata_path, envir = e)
  assign("rdata", file.path(tmp, "moved-away.Rdata"), envir = e)
  save(list = ls(e, all.names = TRUE), file = rdata_path, envir = e)

  testthat::local_mocked_bindings(detectCores = function(...) 2L, .package = "parallel")
  withr::local_dir(tmp)

  target_ci <- list(ci = matrix(0.01, 32, 32))
  iv <- suppressWarnings(
    computeInfoVal2IFC(target_ci, rdata_path, iter = 3, force_gen_ref_dist = TRUE)
  )

  expect_type(iv, "double")
  expect_length(iv, 1)
  expect_false(is.na(iv))
})

test_that("generateCI truncates 4096-parameter files for a single trial too", {
  # Pre-0.3.0 rcicr drew 4096 random contrasts per trial while only 4092 patch
  # indices exist (6 orientations x 2 phases x sum(4^0..4^4)); 0.3.0 dropped the
  # 4 redundant draws and generateCI() has truncated old files ever since. The
  # single-trial branch tested `length(params) == 4092` and then truncated to
  # 4092 -- a no-op that could never fire on the 4096-length input it exists
  # for, so this path died in generateNoiseImage(). The multi-trial branch was
  # always correct and is asserted here too, so the fix cannot invert them.
  tmp <- withr::local_tempdir()

  img_size <- 512
  p <- generateNoisePattern()          # defaults: 4092 parameters
  expect_equal(max(p$patchIdx), 4092)

  base_faces <- list(base = matrix(0.5, img_size, img_size))
  stimuli_params <- list(base = matrix(runif(6 * 4096, -1, 1), nrow = 6))
  seed <- 1
  rdata_path <- file.path(tmp, "old4096.Rdata")
  save(base_faces, stimuli_params, p, img_size, seed, file = rdata_path)

  single <- generateCI(
    stimuli = 1, responses = 1, baseimage = "base",
    rdata = rdata_path, save_as_png = FALSE
  )
  multi <- generateCI(
    stimuli = 1:3, responses = c(1, -1, 1), baseimage = "base",
    rdata = rdata_path, save_as_png = FALSE
  )

  expect_equal(dim(single$ci), c(img_size, img_size))
  expect_equal(dim(multi$ci), c(img_size, img_size))
  expect_false(anyNA(single$ci))

  # The single-trial CI must be the noise of trial 1 alone, which is also what
  # the multi-trial route produces from that one row -- so the truncation cannot
  # be silently dropping or shifting parameters.
  expect_equal(single$ci, generateNoiseImage(stimuli_params$base[1, 1:4092], p))
})

test_that("autoscale handles masked CIs instead of returning all NA", {
  # generateCI(mask = ...) sets masked pixels to NA by design, and autoscale()
  # called a bare range() on the result -- so the scaling constant became NA and
  # the function aborted with "missing value where TRUE/FALSE needed".
  # applyScaling() has always guarded its reductions; this path had not.
  ci <- matrix(seq(-0.2, 0.2, length.out = 64), 8, 8)
  ci[1:3, 1:3] <- NA

  cis <- list(p1 = list(ci = ci, base = matrix(0.5, 8, 8)))
  out <- autoscale(cis, save_as_pngs = FALSE)

  # The masked cells stay masked and everything else is scaled into [0, 1].
  expect_true(all(is.na(out$p1$scaled[1:3, 1:3])))
  expect_false(anyNA(out$p1$scaled[4:8, 4:8]))
  expect_gte(min(out$p1$scaled, na.rm = TRUE), 0)
  expect_lte(max(out$p1$scaled, na.rm = TRUE), 1)

  # The constant must come from the unmasked pixels only, so scaling matches
  # what the same CI with those cells simply absent would give.
  constant <- max(abs(range(ci, na.rm = TRUE)))
  expect_equal(out$p1$scaled, (ci + constant) / (2 * constant))
})

test_that("autoscale errors clearly on a fully masked CI", {
  cis <- list(p1 = list(ci = matrix(NA_real_, 8, 8), base = matrix(0.5, 8, 8)))
  expect_error(autoscale(cis, save_as_pngs = FALSE), "entirely NA")
})
