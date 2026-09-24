test_that("generateNoisePattern returns the expected list structure", {
  p <- generateNoisePattern(img_size = 32, nscales = 2)
  expect_named(p, c("patches", "patchIdx", "noise_type", "generator_version"))
})

test_that("number of patch layers scales with nscales", {
  for (nscales in c(1, 2, 3)) {
    p <- generateNoisePattern(img_size = 16, nscales = nscales)
    expect_equal(dim(p$patches)[3], nscales * 6 * 2)
  }
})

test_that("max patch index matches the documented nparams formula", {
  expected <- list("1" = 12, "2" = 60, "5" = 4092)
  for (nscales in c(1, 2, 5)) {
    p <- generateNoisePattern(img_size = 16, nscales = nscales)
    expect_equal(max(p$patchIdx), expected[[as.character(nscales)]])
    expect_equal(max(p$patchIdx), sum(6 * 2 * (2^(0:(nscales - 1)))^2))
  }
})

test_that("noise_type is passed through", {
  p <- generateNoisePattern(16, nscales = 1, noise_type = "gabor")
  expect_equal(p$noise_type, "gabor")
})

test_that("generator_version matches the installed package version", {
  p <- generateNoisePattern(16, nscales = 1)
  expect_equal(as.character(p$generator_version), as.character(utils::packageVersion("rcicr")))
})

test_that("pre_0.3.0 controls whether patch indices start at 0 or 1", {
  p_legacy <- generateNoisePattern(16, nscales = 1, pre_0.3.0 = TRUE)
  expect_equal(min(p_legacy$patchIdx), 0)

  p_current <- generateNoisePattern(16, nscales = 1, pre_0.3.0 = FALSE)
  expect_equal(min(p_current$patchIdx), 1)
})

test_that("patches and patchIdx contain no NAs", {
  p <- generateNoisePattern(16, nscales = 1)
  expect_false(anyNA(p$patches))
  expect_false(anyNA(p$patchIdx))
})

test_that("an img_size the finest scale cannot tile stops with the constraint and the ways out", {
  expect_error(generateNoisePattern(100, nscales = 5),
               "divisible by 2^(nscales - 1) = 16 for nscales = 5", fixed = TRUE)
  expect_error(generateNoisePattern(100, nscales = 5),
               "Use img_size 96 or 112, or nscales = 3 or fewer.", fixed = TRUE)
  # Below one tile only the larger size is offered.
  expect_error(generateNoisePattern(10, nscales = 5), "Use img_size 16, or nscales = 2 or fewer.",
               fixed = TRUE)
  expect_error(generateNoisePattern(7, nscales = 2), "or nscales = 1.", fixed = TRUE)
})

test_that("divisibility is exactly the condition: every tileable size still works", {
  grid <- expand.grid(img_size = c(16, 24, 30, 32, 48, 50, 64, 96, 100, 120, 128), nscales = 1:5,
                      noise_type = c("sinusoid", "gabor"), stringsAsFactors = FALSE)
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    call <- quote(generateNoisePattern(g$img_size, nscales = g$nscales, noise_type = g$noise_type))
    if (g$img_size %% 2^(g$nscales - 1) == 0) {
      p <- eval(call)
      expect_equal(dim(p$patches)[1:2], c(g$img_size, g$img_size))
    } else {
      expect_error(eval(call), "must be divisible by")
    }
  }
})

test_that("generateStimuli2IFC stops on such a size before writing anything", {
  dir <- file.path(withr::local_tempdir(), "stimuli")
  base <- tempfile(fileext = ".png")
  make_square_png(base, size = 100)
  expect_error(
    suppressWarnings(generateStimuli2IFC(list(face = base), n_trials = 2, img_size = 100,
                                         stimulus_path = dir, ncores = 1, nscales = 5)),
    "Use img_size 96 or 112", fixed = TRUE
  )
  expect_false(dir.exists(dir))
})
