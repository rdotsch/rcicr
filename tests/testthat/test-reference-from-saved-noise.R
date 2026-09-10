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

test_that("a reference cached before the change is regenerated once", {
  # Scoring writes the reference into the file, so an old analysis script rerun
  # on an old stimulus set would otherwise keep returning the value #301 exists
  # to correct: the cache short-circuits the calculation entirely.
  tmp <- withr::local_tempdir()
  base_png <- make_square_png(file.path(tmp, "base.png"), size = 32, seed = 1)
  suppressWarnings(utils::capture.output(generateStimuli2IFC(
    base_face_files = list(base = base_png), n_trials = 4, img_size = 32,
    stimulus_path = tmp, seed = 1, nscales = 3, ncores = 1, save_as_png = FALSE
  )))
  rdata <- list.files(tmp, pattern = "\\.Rdata$", full.names = TRUE)
  ci <- generateCI(1:4, c(1, -1, 1, -1), "base", rdata, save_as_png = FALSE, n_cores = 1)

  # Age the file: no nscales, and a cached distribution with no marker, which is
  # what a pre-1.1.0 file scored by an older rcicr looks like.
  e <- new.env()
  load(rdata, envir = e)
  rm("nscales", "sigma", envir = e)
  e$reference_norms <- rep(0.5, 5)
  e$reference_norms_seed <- NULL
  save(list = ls(e, all.names = TRUE), file = rdata, envir = e)

  warnings_seen <- character()
  withCallingHandlers(
    utils::capture.output(first <- computeInfoVal2IFC(ci, rdata, iter = 5)),
    warning = function(cond) {
      warnings_seen <<- c(warnings_seen, conditionMessage(cond))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("before rcicr built references",
                        warnings_seen, fixed = TRUE)))

  after <- new.env()
  load(rdata, envir = after)
  expect_false(identical(after$reference_norms, rep(0.5, 5)))
  expect_identical(after$reference_norms_source, "saved_noise")

  # Once, not on every call: the marker it wrote settles it.
  second_warnings <- character()
  withCallingHandlers(
    utils::capture.output(second <- computeInfoVal2IFC(ci, rdata, iter = 5)),
    warning = function(cond) {
      second_warnings <<- c(second_warnings, conditionMessage(cond))
      invokeRestart("muffleWarning")
    }
  )
  expect_false(any(grepl("before rcicr built references",
                         second_warnings, fixed = TRUE)))
  expect_equal(second, first)
})

test_that("a deliberately seeded reference is never discarded", {
  # reference_norms_seed records that someone asked for that exact null. It is
  # not a stale default, so the invalidation above must leave it alone.
  tmp <- withr::local_tempdir()
  rdata <- make_fixture_rdata(tmp, img_size = 32, n_trials = 4, nscales = 1, seed = 1)
  ci <- generateCI(1:4, c(1, -1, 1, -1), "base", rdata, save_as_png = FALSE, n_cores = 1)

  e <- new.env()
  load(rdata, envir = e)
  rm("nscales", "sigma", envir = e)
  e$reference_norms <- seq(0.4, 0.8, length.out = 5)
  e$reference_norms_seed <- 99
  save(list = ls(e, all.names = TRUE), file = rdata, envir = e)

  warnings_seen <- character()
  withCallingHandlers(
    utils::capture.output(computeInfoVal2IFC(ci, rdata, iter = 5)),
    warning = function(cond) {
      warnings_seen <<- c(warnings_seen, conditionMessage(cond))
      invokeRestart("muffleWarning")
    }
  )

  expect_false(any(grepl("before rcicr built references",
                         warnings_seen, fixed = TRUE)))
  after <- new.env()
  load(rdata, envir = after)
  expect_equal(after$reference_norms, seq(0.4, 0.8, length.out = 5))
})

test_that("an unmarked cache is refreshed even when the file records nscales", {
  # The marker is the only positive evidence that a cache was built on this
  # file's own noise. A rebuild's correctness depended on the RNG kind of the
  # session that ran it, which set.seed() does not restore and the file does not
  # record, so a modern file's unmarked cache cannot be certified either.
  tmp <- withr::local_tempdir()
  rdata <- make_fixture_rdata(tmp, img_size = 32, n_trials = 4, nscales = 1, seed = 1)
  ci <- generateCI(1:4, c(1, -1, 1, -1), "base", rdata, save_as_png = FALSE, n_cores = 1)

  e <- new.env()
  load(rdata, envir = e)
  expect_true(exists("nscales", envir = e, inherits = FALSE))
  e$reference_norms <- rep(0.5, 5)
  e$reference_norms_seed <- NULL
  save(list = ls(e, all.names = TRUE), file = rdata, envir = e)

  suppressWarnings(utils::capture.output(first <- computeInfoVal2IFC(ci, rdata, iter = 5)))

  after <- new.env()
  load(rdata, envir = after)
  expect_false(identical(after$reference_norms, rep(0.5, 5)))
  expect_identical(after$reference_norms_source, "saved_noise")

  # And once marked it is kept, so the cost is paid once rather than per call.
  marked <- after$reference_norms
  warnings_seen <- character()
  withCallingHandlers(
    utils::capture.output(second <- computeInfoVal2IFC(ci, rdata, iter = 5)),
    warning = function(cond) {
      warnings_seen <<- c(warnings_seen, conditionMessage(cond))
      invokeRestart("muffleWarning")
    }
  )
  expect_false(any(grepl("before rcicr built references", warnings_seen, fixed = TRUE)))
  again <- new.env()
  load(rdata, envir = again)
  expect_identical(again$reference_norms, marked)
  expect_equal(second, first)
})

test_that("a refresh that changes nothing is silent, and survives warn = 2", {
  # Most unmarked caches were already correct, and refreshing reproduces them
  # exactly. Warning on those would put a compatibility notice on stimulus sets
  # nothing happened to, and stop any script running with warnings as errors.
  tmp <- withr::local_tempdir()
  rdata <- make_fixture_rdata(tmp, img_size = 32, n_trials = 4, nscales = 1, seed = 1)
  ci <- generateCI(1:4, c(1, -1, 1, -1), "base", rdata, save_as_png = FALSE, n_cores = 1)

  # A genuine cache from this version, aged by removing only its marker.
  suppressWarnings(utils::capture.output(
    generateReferenceDistribution2IFC(rdata, iter = 20, ncores = 1, save_rdata = TRUE)
  ))
  e <- new.env()
  load(rdata, envir = e)
  genuine <- e$reference_norms
  rm("reference_norms_source", envir = e)
  save(list = ls(e, all.names = TRUE), file = rdata, envir = e)

  warnings_seen <- character()
  withCallingHandlers(
    utils::capture.output(iv <- computeInfoVal2IFC(ci, rdata, iter = 20)),
    warning = function(cond) {
      warnings_seen <<- c(warnings_seen, conditionMessage(cond))
      invokeRestart("muffleWarning")
    }
  )

  # It was refreshed and marked, and came back identical, so nothing is said.
  after <- new.env()
  load(rdata, envir = after)
  expect_identical(after$reference_norms, genuine)
  expect_identical(after$reference_norms_source, "saved_noise")
  expect_false(any(grepl("built references from the saved noise",
                         warnings_seen, fixed = TRUE)))

  # The discriminating case: a cache that really is different does warn.
  e2 <- new.env()
  load(rdata, envir = e2)
  e2$reference_norms <- rep(0.5, 20)
  rm("reference_norms_source", envir = e2)
  save(list = ls(e2, all.names = TRUE), file = rdata, envir = e2)

  moved <- character()
  withCallingHandlers(
    utils::capture.output(computeInfoVal2IFC(ci, rdata, iter = 20)),
    warning = function(cond) {
      moved <<- c(moved, conditionMessage(cond))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("gave different values", moved, fixed = TRUE)))
})

test_that("the warning fires on a real superseded reference, not just a planted one", {
  # The cases above plant an obviously wrong cache to exercise the comparison.
  # This one plants what the old rebuild actually produced for this file --
  # measured on the tree before the change, at nscales = 3 with nscales stripped,
  # where the rebuild assumed the default 5 -- so the warning is shown to fire in
  # the situation it exists for, on values a researcher could really be holding.
  superseded <- c(
    0.856963968913157, 0.855755848808479, 0.896672059991411,
    0.877211505127452, 0.896672059991411
  )

  tmp <- withr::local_tempdir()
  base_png <- make_square_png(file.path(tmp, "base.png"), size = 32, seed = 1)
  suppressWarnings(utils::capture.output(generateStimuli2IFC(
    base_face_files = list(base = base_png), n_trials = 4, img_size = 32,
    stimulus_path = tmp, seed = 1, nscales = 3, ncores = 1, save_as_png = FALSE
  )))
  rdata <- list.files(tmp, pattern = "\\.Rdata$", full.names = TRUE)
  ci <- generateCI(1:4, c(1, -1, 1, -1), "base", rdata, save_as_png = FALSE, n_cores = 1)

  e <- new.env()
  load(rdata, envir = e)
  saved_width <- ncol(e$stimuli_params$base)
  rm("nscales", "sigma", envir = e)
  e$reference_norms <- superseded
  e$reference_norms_seed <- NULL
  save(list = ls(e, all.names = TRUE), file = rdata, envir = e)

  warnings_seen <- character()
  withCallingHandlers(
    utils::capture.output(computeInfoVal2IFC(ci, rdata, iter = length(superseded))),
    warning = function(cond) {
      warnings_seen <<- c(warnings_seen, conditionMessage(cond))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("gave different values", warnings_seen, fixed = TRUE)))

  # And what replaced it is the file's own noise, not merely something else.
  after <- new.env()
  load(rdata, envir = after)
  noise <- vapply(seq_len(after$n_trials), function(i) {
    as.vector(generateNoiseImage(after$stimuli_params$base[i, ], after$p))
  }, numeric(after$img_size^2))
  set.seed(after$seed)
  for (trial in seq_len(after$n_trials)) runif(saved_width)
  expected <- vapply(seq_along(superseded), function(i) {
    responses <- ((runif(after$n_trials) > 0.5) * 2) - 1
    norm(noise %*% responses / ncol(noise), "f")
  }, numeric(1))

  expect_equal(after$reference_norms, expected)
  expect_false(isTRUE(all.equal(after$reference_norms, superseded)))
})

test_that("an automatic refresh keeps the cache's own iteration count", {
  # A file cached at more iterations than the default holds a more precise null.
  # A refresh nobody asked for must not trade that away: the InfoVal would move
  # for a reason unrelated to the noise the null is built on.
  tmp <- withr::local_tempdir()
  rdata <- make_fixture_rdata(tmp, img_size = 32, n_trials = 4, nscales = 1, seed = 1)
  ci <- generateCI(1:4, c(1, -1, 1, -1), "base", rdata, save_as_png = FALSE, n_cores = 1)

  # 40 rather than the 10000 default, so a refresh that ignored it is visible.
  suppressWarnings(utils::capture.output(
    generateReferenceDistribution2IFC(rdata, iter = 40, ncores = 1, save_rdata = TRUE)
  ))
  e <- new.env()
  load(rdata, envir = e)
  genuine <- e$reference_norms
  expect_length(genuine, 40)
  rm("reference_norms_source", envir = e)
  save(list = ls(e, all.names = TRUE), file = rdata, envir = e)

  suppressWarnings(utils::capture.output(computeInfoVal2IFC(ci, rdata)))

  after <- new.env()
  load(rdata, envir = after)
  expect_length(after$reference_norms, 40)
  expect_identical(after$reference_norms, genuine)

  # A caller who names iter still gets it: the rule is about refreshes nobody asked for.
  rm("reference_norms_source", envir = after)
  save(list = ls(after, all.names = TRUE), file = rdata, envir = after)
  suppressWarnings(utils::capture.output(computeInfoVal2IFC(ci, rdata, iter = 12)))
  asked <- new.env()
  load(rdata, envir = asked)
  expect_length(asked$reference_norms, 12)
})
