# reference_method = "gram" (the default) computes the reference through the
# stimulus Gram matrix; "images" is the calculation of rcicr 1.5.0 and earlier
# (#354). They differ by rounding only, measured in
# analyses/gram-reference-accuracy.md; these tests hold "gram" to that on every
# basis layout a stimulus file can carry, and "images" to bit-identity.

reference_of <- function(rdata, iter = 20, response_seed = 5, ...) {
  norms <- NULL
  run <- function() {
    generateReferenceDistribution2IFC(rdata, iter = iter, ncores = 1, response_seed = response_seed,
                                      save_rdata = FALSE, ...)
  }
  suppressWarnings(utils::capture.output(norms <- run()))
  norms
}

# The arithmetic of the reference before the Gram route existed, spelled out.
rendered_arithmetic <- function(rdata, iter = 20, response_seed = 5, base = NULL, ids = NULL) {
  e <- new.env()
  load(rdata, envir = e)
  p <- if (exists("s", envir = e, inherits = FALSE)) e$s else e$p
  params <- rcicr:::savedReferenceParams(e, if (is.null(base)) names(e$stimuli_params)[1] else base)
  if (!is.null(ids)) params <- params[ids, , drop = FALSE]
  noise <- suppressWarnings(vapply(seq_len(nrow(params)), function(i) {
    as.vector(generateNoiseImage(params[i, ], p))
  }, numeric(e$img_size^2)))
  set.seed(response_seed)
  vapply(seq_len(iter), function(i) {
    r <- ((runif(ncol(noise)) > 0.5) * 2) - 1
    norm((noise %*% as.matrix(r)) / ncol(noise), "f")
  }, numeric(1))
}

legacy_copy <- function(version) {
  dest <- file.path(withr::local_tempdir(.local_envir = parent.frame()), "legacy.Rdata")
  file.copy(test_path("fixtures", paste0("legacy-rdata-", version, ".Rdata")), dest) # nolint: object_usage_linter.
  dest
}

test_that("the Gram route matches the rendered arithmetic on every basis layout", {
  current <- make_fixture_rdata(withr::local_tempdir(), img_size = 32, n_trials = 6, nscales = 3)
  expect_equal(reference_of(current), rendered_arithmetic(current), tolerance = 1e-12)

  for (version in c("1.0.1", "1.0.1-gabor", "1.1.0")) {
    rdata <- legacy_copy(version)
    expect_equal(reference_of(rdata), rendered_arithmetic(rdata), tolerance = 1e-12,
                 info = version)
  }

  # Pre-0.3.3 files name the basis s = list(sinusoids, sinIdx).
  renamed <- make_fixture_rdata(withr::local_tempdir(), img_size = 32, n_trials = 6, nscales = 3)
  e <- new.env()
  load(renamed, envir = e)
  mutate_rdata(renamed, s = list(sinusoids = e$p$patches, sinIdx = e$p$patchIdx), .remove = "p")
  expect_equal(reference_of(renamed), rendered_arithmetic(renamed), tolerance = 1e-12)
})

test_that("a 0-based basis matches the rendered arithmetic and warns once", {
  rdata <- make_fixture_rdata(withr::local_tempdir(), img_size = 16, n_trials = 6, nscales = 2)
  mutate_rdata(rdata, p = generateNoisePattern(16, nscales = 2, pre_0.3.0 = TRUE))

  warned <- character()
  norms <- NULL
  record <- function(w) {
    warned <<- c(warned, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
  run <- function() {
    generateReferenceDistribution2IFC(rdata, iter = 20, ncores = 1, response_seed = 5,
                                      save_rdata = FALSE)
  }
  utils::capture.output(norms <- withCallingHandlers(run(), warning = record))
  expect_equal(sum(grepl("patch indices start at 0", warned)), 1)
  expect_equal(norms, rendered_arithmetic(rdata), tolerance = 1e-12)
})

test_that("both ways of building G give the Gram matrix of the rendered noise", {
  gram_of_rendered <- function(params, p) {
    crossprod(vapply(seq_len(nrow(params)), function(i) as.vector(generateNoiseImage(params[i, ], p)),
                     numeric(length(p$patches[, , 1]))))
  }
  # nscales = 1 has 12 parameters, fewer than pixels x trials: the cross-product
  # build. nscales = 5 has 4,092: the rendered build. At 32px each is one block
  # of image columns; at 256px each is summed over four.
  for (img_size in c(32, 256)) {
    for (nscales in c(1, 5)) {
      rdata <- make_fixture_rdata(withr::local_tempdir(), img_size = img_size, n_trials = 6,
                                  nscales = nscales)
      e <- new.env()
      load(rdata, envir = e)
      params <- e$stimuli_params$base
      expect_equal(rcicr:::stimulusGram(params, e$p), gram_of_rendered(params, e$p),
                   tolerance = 1e-12, info = paste(img_size, nscales))
    }
  }
})

test_that("images reproduces the rendered arithmetic exactly, and gram matches it for any size", {
  # 8 x 8 = 64 pixels: 100 stimuli outnumber them, 64 equal them.
  for (n_trials in c(100, 64, 20)) {
    rdata <- make_fixture_rdata(withr::local_tempdir(), img_size = 8, n_trials = n_trials, nscales = 1)
    expect_identical(reference_of(rdata, reference_method = "images"), rendered_arithmetic(rdata),
                     info = n_trials)
    expect_equal(reference_of(rdata), rendered_arithmetic(rdata), tolerance = 1e-12, info = n_trials)
  }
})

test_that("G is built from whichever of the rendered noise and the cross-product is smaller", {
  # 4,092 parameters at five scales: the cross-product is 4,092^2 elements.
  expect_true(rcicr:::gramFromRenderedNoise(npix = 4092^2 / 4, n_trials = 4, nparams = 4092))
  expect_false(rcicr:::gramFromRenderedNoise(npix = 4092^2 / 4, n_trials = 5, nparams = 4092))
  # The package default: 512 x 512 pixels and 770 trials build through the cross-product.
  expect_false(rcicr:::gramFromRenderedNoise(npix = 512^2, n_trials = 770, nparams = 4092))
  expect_true(rcicr:::gramFromRenderedNoise(npix = 64^2, n_trials = 100, nparams = 4092))
})

test_that("images reproduces the rendered arithmetic on the subset and independent-base paths", {
  rdata <- make_fixture_rdata(withr::local_tempdir(), img_size = 8, n_trials = 100, nscales = 1)
  expect_identical(reference_of(rdata, reference_stimuli = 1:70, reference_method = "images"),
                   rendered_arithmetic(rdata, ids = 1:70))

  dir <- withr::local_tempdir()
  base <- make_square_png(file.path(dir, "base.png"), size = 8)
  suppressWarnings(
    generateStimuli2IFC(
      base_face_files = list(first = base, second = base), n_trials = 64, img_size = 8,
      stimulus_path = dir, seed = 17, use_same_parameters = FALSE, nscales = 1, ncores = 1,
      save_as_png = FALSE
    )
  )
  independent <- list.files(dir, pattern = "\\.Rdata$", full.names = TRUE)[1]
  expect_identical(reference_of(independent, baseimage = "second", reference_method = "images"),
                   rendered_arithmetic(independent, base = "second"))
})

test_that("computing with gram gives the notice, and images gives none", {
  rdata <- make_fixture_rdata(withr::local_tempdir(), img_size = 32, n_trials = 6, nscales = 1)
  messages_for <- function(method) {
    seen <- character()
    record <- function(m) {
      seen <<- c(seen, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
    run <- function() {
      generateReferenceDistribution2IFC(rdata, iter = 20, ncores = 1, save_rdata = FALSE,
                                        reference_method = method)
    }
    utils::capture.output(withCallingHandlers(suppressWarnings(run()), message = record))
    seen
  }
  gram <- messages_for("gram")
  expect_length(gram, 1)
  expect_match(gram, 'computed with reference_method = "gram"')
  expect_match(gram, 'pass reference_method = "images" to reproduce them bit for bit')
  expect_length(messages_for("images"), 0)
  expect_error(messages_for("pixels"), "should be one of")
})

test_that("a stored reference is reused under either method, and a new one records its method", {
  quiet_messages <- function(expr) {
    seen <- character()
    out <- utils::capture.output(withCallingHandlers(suppressWarnings(expr), message = function(m) {
      seen <<- c(seen, conditionMessage(m))
      invokeRestart("muffleMessage")
    }))
    list(messages = seen, output = out)
  }
  rdata <- make_fixture_rdata(withr::local_tempdir(), img_size = 32, n_trials = 6, nscales = 1)
  mutate_rdata(rdata, reference_method = "planted")
  computed <- quiet_messages(generateReferenceDistribution2IFC(rdata, iter = 20, ncores = 1))
  expect_match(computed$messages, 'computed with reference_method = "gram"')
  saved <- new.env()
  load(rdata, envir = saved)
  expect_identical(saved$reference_norms_method, "gram")
  expect_identical(saved$reference_method, "planted")

  ci <- NULL
  utils::capture.output(ci <- generateCI(1:6, rep(c(1, -1), 3), "base", rdata, save_as_png = FALSE))
  for (method in c("images", "gram")) {
    reused <- quiet_messages(computeInfoVal2IFC(ci, rdata, reference_method = method))
    expect_true(any(grepl("Using reference distribution found", reused$output)), info = method)
    expect_false(any(grepl("computed with reference_method", reused$messages)), info = method)
  }

  independent <- make_independent_fixture(withr::local_tempdir())
  quiet_messages(generateReferenceDistribution2IFC(independent, iter = 20, ncores = 1,
                                                   baseimage = "second", reference_method = "images"))
  quiet_messages(generateReferenceDistribution2IFC(independent, iter = 20, ncores = 1,
                                                   baseimage = "second", reference_stimuli = 1:3))
  saved <- new.env()
  load(independent, envir = saved)
  expect_identical(saved$reference_norms_by_base$second$method, "images")
  expect_identical(saved$reference_norms_by_stimuli[[1]]$method, "gram")
})

test_that("the Gram route holds its tolerance at a realistic size", {
  # CI runs this on each platform's BLAS; the analysis measured only the
  # reference BLAS.
  skip_on_cran()
  rdata <- make_fixture_rdata(withr::local_tempdir(), img_size = 128, n_trials = 300, nscales = 5)
  expect_equal(reference_of(rdata, iter = 200), rendered_arithmetic(rdata, iter = 200),
               tolerance = 1e-12)
})
