# Paths a coverage run found no test executing. Each is behaviour a researcher
# depends on, or the error that stands in for a wrong number.

quietly <- function(expr) {
  out <- NULL
  utils::capture.output(out <- suppressWarnings(expr))
  out
}

responses8 <- c(1, -1, -1, 1, 1, 1, -1, 1)

test_that("masked individual CI images hold each participant's masked, scaled CI", {
  rdata <- make_fixture_rdata(withr::local_tempdir(), n_trials = 8)
  out <- withr::local_tempdir()
  participants <- rep(c(3, 1), each = 4)
  mask <- matrix(1, 32, 32)
  mask[1:6, 1:6] <- 0
  quietly(generateCI(1:8, responses8, "base", rdata, participants = participants, mask = mask,
                     save_individual_cis = TRUE, targetpath = out, save_as_png = FALSE,
                     n_cores = 1))
  e <- new.env()
  load(rdata, envir = e)
  for (id in c(1, 3)) {
    rows <- participants == id
    own <- quietly(generateCI((1:8)[rows], responses8[rows], "base", rdata, mask = mask,
                              save_as_png = FALSE))
    expected <- own$combined
    expected[is.na(expected)] <- 0 # writePNG() writes NA as 0
    written <- png::readPNG(file.path(out, "individual_cis", paste0("ci_", id, ".png")))
    expect_equal(written, expected, tolerance = 1 / 255, info = id)
    expect_true(all(written[1:6, 1:6] == 0), info = id)
  }
})

test_that("computeCumulativeCICorrelation() reads a pre-0.3.3 file's s as its p", {
  current <- make_fixture_rdata(withr::local_tempdir(), n_trials = 8)
  renamed <- file.path(withr::local_tempdir(), "old.Rdata")
  file.copy(current, renamed)
  e <- new.env()
  load(renamed, envir = e)
  mutate_rdata(renamed, s = list(sinusoids = e$p$patches, sinIdx = e$p$patchIdx), .remove = "p")
  curve <- function(rdata) quietly(computeCumulativeCICorrelation(1:8, responses8, "base", rdata))
  expected <- curve(current)
  expect_identical(curve(renamed), expected)
  expect_gt(length(unique(expected)), 1)
})

test_that("a stored reference that is not finite stops rather than giving NaN", {
  # A seeded entry is kept as stored, so only the guard stands between it and the InfoVal.
  independent <- make_independent_fixture(withr::local_tempdir())
  ci <- quietly(generateCI(1:12, rep(c(1, -1), 6), "second", independent, save_as_png = FALSE))
  stored <- list(second = list(norms = c(1, NA, 2), response_seed = 5))
  mutate_rdata(independent, reference_norms_by_base = stored)
  expect_error(quietly(computeInfoVal2IFC(ci, independent, baseimage = "second")),
               "Invalid cached reference for baseimage second")
  mutate_rdata(independent, reference_norms_by_base = "not a list")
  expect_error(quietly(computeInfoVal2IFC(ci, independent, baseimage = "second")),
               "reference_norms_by_base must be a list")

  shared <- make_fixture_rdata(withr::local_tempdir(), n_trials = 8)
  sub_ci <- quietly(generateCI(1:6, responses8[1:6], "base", shared, save_as_png = FALSE))
  mutate_rdata(shared, reference_norms_by_stimuli = list(list(
    reference_stimuli = 1:6, baseimage = NULL, norms = c(1, Inf, 2), response_seed = 5
  )))
  expect_error(quietly(computeInfoVal2IFC(sub_ci, shared, reference_stimuli = 1:6)),
               "Invalid cached reference for these reference_stimuli")
  mutate_rdata(shared, reference_norms_by_stimuli = "not a list")
  expect_error(quietly(computeInfoVal2IFC(sub_ci, shared, reference_stimuli = 1:6)),
               "reference_norms_by_stimuli must be a list")
})

test_that("a stimulus file without its noise basis or a valid n_trials stops with that reason", {
  rdata <- make_fixture_rdata(withr::local_tempdir(), n_trials = 8)
  no_basis <- file.path(withr::local_tempdir(), "no-basis.Rdata")
  file.copy(rdata, no_basis)
  mutate_rdata(no_basis, .remove = "p")
  for (method in c("gram", "images")) {
    expect_error(quietly(generateReferenceDistribution2IFC(no_basis, iter = 5, ncores = 1,
                                                           save_rdata = FALSE,
                                                           reference_method = method)),
                 "does not contain its saved noise basis", info = method)
  }
  mutate_rdata(rdata, n_trials = 2.5)
  expect_error(quietly(generateReferenceDistribution2IFC(rdata, iter = 5, ncores = 1, save_rdata = FALSE)),
               "positive integer n_trials")
  expect_error(quietly(generateReferenceDistribution2IFC(rdata, iter = 5, ncores = 1, save_rdata = FALSE,
                                                         reference_stimuli = 1:2)),
               "positive integer n_trials")
  ci <- list(ci = matrix(0, 32, 32))
  expect_error(quietly(batchComputeInfoVal2IFC(list(ci, ci), rdata)), "positive integer n_trials")
})

test_that("a stimulus folder that cannot be written stops before anything is generated", {
  skip_on_os("windows")
  skip_if(unname(Sys.info()[["effective_user"]]) == "root", "root can write to any folder")
  base <- make_square_png(file.path(withr::local_tempdir(), "base.png"), size = 16)
  dir <- withr::local_tempdir()
  Sys.chmod(dir, "555")
  withr::defer(Sys.chmod(dir, "755"))
  generate <- function(...) {
    quietly(generateStimuli2IFC(list(base = base), n_trials = 4, img_size = 16, stimulus_path = dir,
                                seed = 1, ncores = 1, nscales = 1, ...))
  }
  expect_error(generate(), "Could not create .* to reserve .*; check that the folder exists and is writable")
  expect_error(generate(save_rdata = FALSE), "to reserve the stimulus PNGs; check that .* is writable")
  expect_length(list.files(dir, all.files = TRUE, no.. = TRUE), 0)
})

test_that("default_ncores() stays at 2 under R CMD check and leaves one core free otherwise", {
  withr::with_envvar(c("_R_CHECK_LIMIT_CORES_" = "TRUE"), expect_identical(rcicr:::default_ncores(), 2L))
  withr::with_envvar(c("_R_CHECK_LIMIT_CORES_" = ""), {
    expect_equal(rcicr:::default_ncores(), max(1L, parallel::detectCores() - 1))
  })
})

test_that("long lists of offending IDs are cut at five", {
  rdata <- make_fixture_rdata(withr::local_tempdir(), n_trials = 12)
  d <- data.frame(g = letters[1:12], s = 1:12, r = rep(c(1, -1), 6), pid = NA)
  expect_error(batchGenerateCI(d, "g", "s", "r", "base", rdata, save_as_png = FALSE, participants = "pid"),
               "\\(g a, b, c, d, e, \\.\\.\\.\\)")
  expect_error(canonicalReferenceStimuli(rep(1:7, 2), 12), "stimulus 1, 2, 3, 4, 5, \\.\\.\\. more than once")
  cis <- lapply(1:7, function(i) {
    quietly(generateCI(c(i, i), c(1, -1), "base", rdata, save_as_png = FALSE))
  })
  names(cis) <- paste0("c", 1:7)
  seen <- character()
  utils::capture.output(withCallingHandlers(
    suppressWarnings(batchComputeInfoVal2IFC(cis, rdata, iter = 20, response_seed = 1)),
    message = function(m) {
      seen <<- c(seen, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  ))
  expect_match(seen, "7 of the 7 classification images \\(c1, c2, c3, c4, c5, \\.\\.\\.\\) average",
               all = FALSE)
})

test_that("a damaged stimulus file without a backup fails with load()'s own error", {
  path <- file.path(withr::local_tempdir(), "damaged.Rdata")
  writeLines("not an Rdata file", path)
  expect_error(rcicr:::loadRdata(path, new.env()), "bad restore file magic number")
})

test_that("a backup that cannot be copied back is reported, and kept", {
  dir <- withr::local_tempdir()
  file <- file.path(dir, "stimuli.Rdata")
  backup <- paste0(file, ".rcicr-backup")
  writeLines("current", file)
  writeLines("backup", backup)
  local_mocked_bindings(copyInto = function(from, to) invisible(FALSE))
  expect_warning(rcicr:::restoreFromBackup(file, backup), "could not be restored after the failed save")
  expect_true(file.exists(backup))
})

test_that("an all-NA z-map with no scale draws no legend", {
  out <- withr::local_tempdir()
  expect_no_error(plotZmap(matrix(NA_real_, 16, 16), sigma = 1, threshold = 0, decoration = FALSE,
                           targetpath = out, size = 16))
  expect_null(rcicr:::drawZmapLegend(matrix(NA_real_, 4, 4), col = grDevices::grey.colors(3)))
})
