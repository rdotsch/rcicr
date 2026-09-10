# The InfoVal reference is built from the noise the stimulus file saved.
#
# It used to be built by re-generating the stimuli, which reopened every base
# image and rebuilt the noise basis from fields older files do not carry
# (issue #301). These pin that the numbers did not move for files the rebuild
# already reproduced, and that files it could not reach are now scored.

reference_of <- function(rdata, iter = 6, ...) {
  norms <- NULL
  suppressWarnings(utils::capture.output({
    norms <- rcicr::generateReferenceDistribution2IFC(rdata, iter = iter, ncores = 1,
                                                      save_rdata = FALSE, ...)
  }))
  norms
}

test_that("the reference and the stream it leaves behind are unchanged", {
  # Measured on the tree before the change, at full precision.
  tmp <- withr::local_tempdir()
  rdata <- make_fixture_rdata(tmp, img_size = 32, n_trials = 4, nscales = 1, seed = 1)

  expect_equal(reference_of(rdata), c(
    2.33438571356728, 1.85392028318574, 1.74389383373095,
    1.54790538103884, 1.74389383373095, 1.72263882900775
  ))

  # The responses are drawn from whatever the stimulus seed left, so what the
  # caller's stream holds afterwards is part of the contract.
  set.seed(99)
  invisible(reference_of(rdata, iter = 3))
  expect_equal(runif(3), c(0.912875924259424, 0.293603372760117, 0.459065726259723))
})

test_that("a stimulus file whose base image has moved is still scored", {
  tmp <- withr::local_tempdir()
  rdata <- make_fixture_rdata(tmp, img_size = 32, n_trials = 4, nscales = 1, seed = 1)

  intact <- reference_of(rdata)

  e <- new.env()
  load(rdata, envir = e)
  expect_true(file.exists(e$base_face_files$base))
  unlink(e$base_face_files$base)

  # Not merely that it runs: the noise is the same noise, so the answer is too.
  expect_equal(reference_of(rdata), intact)
})

test_that("computeInfoVal2IFC is unaffected by a base image that has moved", {
  tmp <- withr::local_tempdir()
  rdata <- make_fixture_rdata(tmp, img_size = 32, n_trials = 6, nscales = 1, seed = 1)
  ci <- generateCI(1:6, rep(c(1, -1), 3), "base", rdata, save_as_png = FALSE, n_cores = 1)

  suppressWarnings(utils::capture.output({
    before <- computeInfoVal2IFC(ci, rdata, iter = 5)
  }))

  e <- new.env()
  load(rdata, envir = e)
  unlink(e$base_face_files$base)
  # The cached norms would hide the difference, so score from a clean copy.
  clean <- file.path(tmp, "clean.Rdata")
  file.copy(rdata, clean)
  suppressWarnings(utils::capture.output({
    after <- computeInfoVal2IFC(ci, clean, iter = 5, force_gen_ref_dist = TRUE)
  }))

  expect_equal(after, before)
})

test_that("a uniform base with contrast maximization off is scored", {
  tmp <- withr::local_tempdir()
  base_png <- file.path(tmp, "flat.png")
  png::writePNG(matrix(0.5, 32, 32), base_png)
  suppressWarnings(utils::capture.output(generateStimuli2IFC(
    base_face_files = list(flat = base_png), n_trials = 4, img_size = 32,
    stimulus_path = tmp, seed = 1, nscales = 1, ncores = 1,
    maximize_baseimage_contrast = FALSE, save_as_png = FALSE
  )))
  rdata <- list.files(tmp, pattern = "\\.Rdata$", full.names = TRUE)

  norms <- reference_of(rdata)
  expect_length(norms, 6)
  expect_true(all(is.finite(norms)))
})

test_that("each committed legacy fixture is scored", {
  # Their base_face_files point at temp directories from the machines that wrote
  # them, so every one of these failed with "does not exist" before the change.
  tmp <- withr::local_tempdir()
  for (f in list.files(test_path("fixtures"), pattern = "^legacy-rdata.*\\.Rdata$",
                       full.names = TRUE)) {
    copy <- file.path(tmp, basename(f))
    file.copy(f, copy)
    norms <- reference_of(copy, iter = 4)
    expect_length(norms, 4)
    expect_true(all(is.finite(norms)), info = basename(f))
  }
})

test_that("a pre-0.3.0 file spends 4092 draws per trial, not its 4096 columns", {
  # No fixture is pre-0.3.0, so the 4096-column layout is built here: those files
  # hold four parameters per trial that no patch index refers to, and
  # selectStimulusParams() drops them.
  tmp <- withr::local_tempdir()
  p <- generateNoisePattern(16, nscales = 5)
  expect_equal(max(p$patchIdx), 4092)

  n_trials <- 2
  seed <- 7
  stimuli_params <- list(base = withr::with_seed(3, {
    matrix(runif(n_trials * 4096) * 2 - 1, n_trials, 4096)
  }))
  base_faces <- list(base = matrix(0.5, 16, 16))
  base_face_files <- list(base = file.path(tmp, "gone.png"))
  img_size <- 16
  rdata <- file.path(tmp, "pre030.Rdata")
  save(base_face_files, base_faces, img_size, n_trials, p, seed, stimuli_params,
       file = rdata)

  norms <- reference_of(rdata, iter = 4)

  oracle <- function(draws_per_trial) {
    noise <- vapply(seq_len(n_trials), function(i) {
      as.vector(generateNoiseImage(stimuli_params$base[i, 1:4092], p))
    }, numeric(img_size^2))
    set.seed(seed)
    for (trial in seq_len(n_trials)) runif(draws_per_trial)
    vapply(seq_len(4), function(i) {
      responses <- ((runif(n_trials) > 0.5) * 2) - 1
      norm(noise %*% responses / ncol(noise), "f")
    }, numeric(1))
  }

  expect_equal(norms, oracle(4092))
  # Not vacuous: spending the raw 4096 would land somewhere else.
  expect_false(isTRUE(all.equal(norms, oracle(4096))))
})

test_that("a file without nscales uses its saved basis, not the default", {
  # Files written before 1.1.0 do not record nscales. The rebuild assumed 5 and
  # scored them against a basis their stimuli never used; the saved basis is the
  # one participants saw. Documented under Reproducibility impact.
  tmp <- withr::local_tempdir()
  base_png <- make_square_png(file.path(tmp, "base.png"), size = 32, seed = 1)
  suppressWarnings(utils::capture.output(generateStimuli2IFC(
    base_face_files = list(base = base_png), n_trials = 4, img_size = 32,
    stimulus_path = tmp, seed = 1, nscales = 3, ncores = 1, save_as_png = FALSE
  )))
  rdata <- list.files(tmp, pattern = "\\.Rdata$", full.names = TRUE)

  e <- new.env()
  load(rdata, envir = e)
  saved_width <- ncol(e$stimuli_params$base)
  rm("nscales", "sigma", envir = e)
  save(list = ls(e, all.names = TRUE), file = rdata, envir = e)

  norms <- reference_of(rdata, iter = 5)

  e2 <- new.env()
  load(rdata, envir = e2)
  noise <- vapply(seq_len(e2$n_trials), function(i) {
    as.vector(generateNoiseImage(e2$stimuli_params$base[i, ], e2$p))
  }, numeric(e2$img_size^2))
  set.seed(e2$seed)
  for (trial in seq_len(e2$n_trials)) runif(saved_width)
  expected <- vapply(seq_len(5), function(i) {
    responses <- ((runif(e2$n_trials) > 0.5) * 2) - 1
    norm(noise %*% responses / ncol(noise), "f")
  }, numeric(1))

  expect_equal(norms, expected)
  # The values the rebuild produced at the assumed nscales = 5, which these replace.
  expect_false(isTRUE(all.equal(norms, c(
    0.856963968913157, 0.855755848808479, 0.896672059991411,
    0.877211505127452, 0.896672059991411
  ))))
})
