read_reference_fixture <- function(path) {
  e <- new.env(parent = emptyenv())
  load(path, envir = e)
  as.list(e)
}

# Construct trial pixels independently of the reference implementation, then
# apply the historical random-response stream to those pixels.
reference_oracle <- function(saved, baseimage, iter = 24, response_seed = 91,
                             responses = NULL) {
  params <- saved$stimuli_params[[baseimage]]
  if (ncol(params) == 4096) params <- params[, seq_len(4092), drop = FALSE]
  pixels <- vapply(seq_len(saved$n_trials), function(i) {
    as.vector(generateNoiseImage(params[i, ], saved$p))
  }, numeric(saved$img_size ^ 2))
  if (is.null(responses)) {
    responses <- withr::with_seed(response_seed, replicate(iter,
      ((runif(saved$n_trials) > 0.5) * 2) - 1
    ))
  }
  vapply(seq_len(ncol(responses)), function(i) {
    norm(pixels %*% responses[, i] / saved$n_trials, "f")
  }, numeric(1))
}

selected_reference <- function(path, baseimage, iter = 24, response_seed = 91,
                               save_rdata = FALSE, ncores = 1) {
  suppressWarnings(generateReferenceDistribution2IFC(
    path, iter = iter, ncores = ncores, response_seed = response_seed,
    save_rdata = save_rdata, baseimage = baseimage
  ))
}

test_that("each independent base uses its saved trial noise for reference and InfoVal", {
  path <- make_independent_fixture(withr::local_tempdir())
  saved <- read_reference_fixture(path)
  testthat::local_mocked_bindings(default_ncores = function() 1L, .package = "rcicr")
  expected <- lapply(c("first", "second"), function(key) reference_oracle(saved, key))
  expect_false(isTRUE(all.equal(expected[[1]], expected[[2]])))
  for (i in seq_along(expected)) {
    key <- c("first", "second")[i]
    expect_gt(mad(expected[[i]]), 0)
    expect_equal(selected_reference(path, key), expected[[i]], tolerance = 1e-12)
    target <- generateCI(seq_len(saved$n_trials), rep(c(1, -1), 6),
                         baseimage = key, rdata = path, save_as_png = FALSE, scaling = "none")
    info <- suppressWarnings(computeInfoVal2IFC(target, path, iter = 24,
                                                response_seed = 91, baseimage = key))
    expect_equal(info, (norm(target$ci, "f") - median(expected[[i]])) /
                   mad(expected[[i]]), tolerance = 1e-12)
  }
})

test_that("independent defaults retain the old shared-rebuild response stream", {
  path <- make_independent_fixture(withr::local_tempdir())
  saved <- read_reference_fixture(path)
  suppressWarnings(generateStimuli2IFC(saved$base_face_files, saved$n_trials,
                                       saved$img_size, seed = saved$seed, noise_type = saved$noise_type,
                                       nscales = saved$nscales, sigma = saved$sigma, ncores = 1,
                                       return_as_dataframe = TRUE, save_as_png = FALSE, save_rdata = FALSE))
  responses <- replicate(24, ((runif(saved$n_trials) > 0.5) * 2) - 1)
  for (key in c("first", "second")) {
    expected <- reference_oracle(saved, key, responses = responses)
    expect_equal(selected_reference(path, key, response_seed = NULL), expected,
                 tolerance = 1e-12)
  }
  expect_equal(selected_reference(path, "second", ncores = 2),
               selected_reference(path, "second", ncores = 1), tolerance = 1e-12)
})

test_that("selection is required from actual differing matrices and validates keys", {
  path <- make_independent_fixture(withr::local_tempdir())
  mutate_rdata(path, use_same_parameters = TRUE, reference_norms = seq_len(24))
  expect_error(suppressWarnings(generateReferenceDistribution2IFC(
                                                                  path, iter = 24, ncores = 1, save_rdata = FALSE)), "baseimage")
  expect_error(computeInfoVal2IFC(list(ci = matrix(1, 32, 32)), path), "baseimage")
  for (key in list("missing", "", NA_character_, c("first", "second"), 1)) {
    expect_error(selected_reference(path, key), "baseimage")
    expect_error(computeInfoVal2IFC(list(ci = matrix(1, 32, 32)), path,
                                    baseimage = key), "baseimage")
  }
})

test_that("identical matrices preserve the shared route despite independent metadata", {
  path <- make_independent_fixture(withr::local_tempdir())
  saved <- read_reference_fixture(path)
  saved$stimuli_params$second <- saved$stimuli_params$first
  mutate_rdata(path, stimuli_params = saved$stimuli_params)
  actual <- selected_reference(path, NULL)
  expect_equal(actual, reference_oracle(saved, "first"), tolerance = 1e-12)
  expect_equal(selected_reference(path, "second"), actual, tolerance = 1e-12)
})

test_that("base-scoped caches isolate hits and forced refreshes", {
  path <- make_independent_fixture(withr::local_tempdir())
  legacy <- seq(100, 123)
  mutate_rdata(path, reference_norms = legacy, reference_norms_seed = 77)
  first <- selected_reference(path, "first", response_seed = NULL, save_rdata = TRUE)
  second <- selected_reference(path, "second", response_seed = NULL, save_rdata = TRUE)
  saved <- read_reference_fixture(path)
  expect_identical(saved$reference_norms, legacy)
  expect_identical(saved$reference_norms_seed, 77)
  expect_equal(saved$reference_norms_by_base$first$norms, first)
  expect_equal(saved$reference_norms_by_base$second$norms, second)
  expect_true("response_seed" %in% names(saved$reference_norms_by_base$first))
  expect_null(saved$reference_norms_by_base$first$response_seed)
  testthat::local_mocked_bindings(default_ncores = function() 1L, .package = "rcicr")
  target <- list(ci = matrix(1, 32, 32))
  for (key in c("first", "second")) {
    norms <- saved$reference_norms_by_base[[key]]$norms
    expect_equal(computeInfoVal2IFC(target, path, iter = 7, baseimage = key),
                 (norm(target$ci, "f") - median(norms)) / mad(norms))
  }
  suppressWarnings(computeInfoVal2IFC(target, path, iter = 31,
                                      force_gen_ref_dist = TRUE, baseimage = "first"))
  updated <- read_reference_fixture(path)
  expect_length(updated$reference_norms_by_base$first$norms, 31)
  expect_identical(updated$reference_norms_by_base$second,
                   saved$reference_norms_by_base$second)
  expect_identical(updated$reference_norms, legacy)
})

test_that("unscoped references never supply an independent base cache miss", {
  path <- make_independent_fixture(withr::local_tempdir())
  mutate_rdata(path, reference_norms = seq(100, 123))
  testthat::local_mocked_bindings(default_ncores = function() 1L, .package = "rcicr")
  target <- list(ci = matrix(1, 32, 32))
  norms <- selected_reference(path, "second", response_seed = NULL)
  actual <- suppressWarnings(computeInfoVal2IFC(target, path, iter = 24,
                                                baseimage = "second"))
  expect_equal(actual, (norm(target$ci, "f") - median(norms)) / mad(norms))
  expect_equal(read_reference_fixture(path)$reference_norms_by_base$second$norms,
               norms)
})

test_that("seeded InfoVal is read-only and loaded collisions cannot replace arguments", {
  path <- make_independent_fixture(withr::local_tempdir())
  decoy <- file.path(dirname(path), "absent.Rdata")
  mutate_rdata(path, baseimage = "first", rdata = decoy, iter = 2,
               response_seed = 2, save_rdata = TRUE, ncores = 99,
               target_ci = list(ci = matrix(100, 32, 32)), force_gen_ref_dist = FALSE)
  expected <- reference_oracle(read_reference_fixture(path), "second")
  before <- tools::md5sum(path)
  expect_equal(selected_reference(path, "second"), expected, tolerance = 1e-12)
  testthat::local_mocked_bindings(default_ncores = function() 1L, .package = "rcicr")
  target <- list(ci = matrix(1, 32, 32))
  actual <- suppressWarnings(computeInfoVal2IFC(target, path, iter = 24,
                                                response_seed = 91, baseimage = "second"))
  expect_equal(actual, (norm(target$ci, "f") - median(expected)) / mad(expected))
  expect_identical(tools::md5sum(path), before)
  expect_false(file.exists(decoy))
})

test_that("independent references work without original image files", {
  path <- make_independent_fixture(withr::local_tempdir())
  saved <- read_reference_fixture(path)
  expected <- reference_oracle(saved, "second")
  unlink(unlist(saved$base_face_files))
  expect_equal(selected_reference(path, "second"), expected, tolerance = 1e-12)
})

test_that("legacy basis names and unused parameter columns retain selected noise", {
  path <- make_independent_fixture(withr::local_tempdir(), nscales = 5)
  saved <- read_reference_fixture(path)
  expected <- reference_oracle(saved, "second")
  params <- lapply(saved$stimuli_params, function(x) cbind(x, matrix(999, nrow(x), 4)))
  old_p <- list(sinusoids = saved$p$patches, sinIdx = saved$p$patchIdx)
  mutate_rdata(path, stimuli_params = params, s = old_p, .remove = "p")
  expect_equal(selected_reference(path, "second"), expected, tolerance = 1e-12)
})
