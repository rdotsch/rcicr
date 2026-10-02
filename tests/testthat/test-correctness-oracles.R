# Correctness against oracles built from the documented definitions, not from
# the code under test or from values this repository computed earlier.

quietly <- function(expr) {
  out <- NULL
  utils::capture.output(out <- suppressMessages(suppressWarnings(expr)))
  out
}

test_that("a sinusoid follows its formula, endpoints included", {
  # sin(2 pi (x cos a + y sin a) + phase), x and y running from 0 to `cycles`
  # in img_size equal steps, both ends included (MATLAB's linspace).
  for (angle in c(0, 30, 90, 150)) {
    for (phase in c(0, pi / 2)) {
      n <- 12
      x <- seq(0, 2, length.out = n)
      a <- angle * pi / 180
      expected <- outer(x, x, function(row, col) sin(2 * pi * (col * cos(a) + row * sin(a)) + phase))
      expect_equal(generateSinusoid(n, 2, angle, phase, 0.7), 0.7 * expected, tolerance = 1e-12,
                   info = paste(angle, phase))
    }
  }
})

test_that("the noise basis lays out its patches as documented", {
  # Layers run scale (1, 2, 4, ...) x orientation (0 to 150 by 30) x phase (0,
  # pi/2). A layer at scale s tiles an s x s grid of patches of side
  # img_size / s, each two cycles of its sinusoid, and numbers them down the
  # columns, continuing from the previous layer.
  img_size <- 16
  nscales <- 3
  p <- generateNoisePattern(img_size, nscales = nscales)
  layer <- 0
  index <- 0
  for (s in 2^(0:(nscales - 1))) {
    side <- img_size / s
    for (orientation in seq(0, 150, by = 30)) {
      for (phase in c(0, pi / 2)) {
        layer <- layer + 1
        x <- seq(0, 2, length.out = side)
        a <- orientation * pi / 180
        patch <- outer(x, x, function(row, col) sin(2 * pi * (col * cos(a) + row * sin(a)) + phase))
        expect_equal(p$patches[, , layer], kronecker(matrix(1, s, s), patch), tolerance = 1e-12,
                     info = layer)
        for (col in seq_len(s)) {
          for (row in seq_len(s)) {
            index <- index + 1
            block <- p$patchIdx[(row - 1) * side + seq_len(side), (col - 1) * side + seq_len(side), layer]
            expect_true(all(block == index), info = paste(layer, row, col))
          }
        }
      }
    }
  }
  expect_equal(layer, dim(p$patches)[3])
  expect_equal(index, max(p$patchIdx))
})

test_that("each stimulus image is the base averaged with its noise, or the noise inverted", {
  # ori = ((noise + 0.3) / 0.6 + base) / 2, inv the same with -noise, where
  # noise is the trial's saved parameters rendered through the saved basis.
  dir <- withr::local_tempdir()
  base_png <- make_square_png(file.path(withr::local_tempdir(), "base.png"), size = 16)
  quietly(generateStimuli2IFC(list(face = base_png), n_trials = 4, img_size = 16,
                              stimulus_path = dir, seed = 3, ncores = 1, nscales = 2))
  e <- new.env()
  load(list.files(dir, pattern = "\\.Rdata$", full.names = TRUE), envir = e)
  base <- e$base_faces$face
  for (trial in 1:4) {
    weights <- e$stimuli_params$face[trial, ]
    noise <- apply(array(weights[e$p$patchIdx], dim(e$p$patches)) * e$p$patches, 1:2, mean)
    for (side in c("ori", "inv")) {
      signed <- if (side == "ori") noise else -noise
      expected <- pmin(pmax(((signed + 0.3) / 0.6 + base) / 2, 0), 1)
      written <- png::readPNG(file.path(dir, sprintf("rcic_face_3_%05d_%s.png", trial, side)))
      expect_equal(written, expected, tolerance = 1 / 255, ignore_attr = TRUE, info = paste(trial, side))
    }
    ori <- png::readPNG(file.path(dir, sprintf("rcic_face_3_%05d_ori.png", trial)))
    inv <- png::readPNG(file.path(dir, sprintf("rcic_face_3_%05d_inv.png", trial)))
    expect_gt(max(abs(ori - inv)), 0.05)
  }
})

# Random responses give a CI with no signal, so their InfoVals, scored against
# a reference from random responses to the same stimuli, have median 0 and
# MAD 1 up to sampling error. That only holds if generateCI() and the
# reference scale a CI the same way; a factor between them shifts the median.
null_infovals <- function(rdata, base, stimuli, n_ci = 150, ...) {
  set.seed(11)
  cis <- lapply(seq_len(n_ci), function(i) {
    quietly(generateCI(stimuli, sample(c(1, -1), length(stimuli), TRUE), base, rdata,
                       save_as_png = FALSE))
  })
  quietly(batchComputeInfoVal2IFC(cis, rdata, iter = 4000, response_seed = 2, ...))
}

test_that("InfoVal is calibrated under the null: median 0 and MAD 1", {
  skip_on_cran()
  rdata <- make_fixture_rdata(withr::local_tempdir(), img_size = 32, n_trials = 40, nscales = 2)
  full <- null_infovals(rdata, "base", 1:40)
  expect_lt(abs(median(full)), 0.25)
  expect_lt(abs(mad(full) - 1), 0.25)
  # The two reference methods agree to rounding (#354).
  expect_equal(null_infovals(rdata, "base", 1:40, reference_method = "images"), full,
               tolerance = 1e-10)

  # The subset reference matches a CI built from those stimuli...
  subset <- null_infovals(rdata, "base", 1:25, reference_stimuli = 1:25)
  expect_lt(abs(median(subset)), 0.25)
  expect_lt(abs(mad(subset) - 1), 0.25)
  # ...and the full reference does not: fewer stimuli give larger norms.
  mismatched <- null_infovals(rdata, "base", 1:25)
  expect_gt(median(mismatched), 1)

  independent <- make_independent_fixture(withr::local_tempdir())
  second <- null_infovals(independent, "second", 1:12, baseimage = "second")
  expect_lt(abs(median(second)), 0.35)
  expect_lt(abs(mad(second) - 1), 0.35)
})

test_that("the stored base face is the colour channels' mean, stretched to 0..1", {
  dir <- withr::local_tempdir()
  set.seed(4)
  img <- array(stats::runif(16 * 16 * 3, 0.2, 0.7), dim = c(16, 16, 3))
  path <- file.path(dir, "rgb.png")
  png::writePNG(img, path)
  out <- withr::local_tempdir()
  quietly(generateStimuli2IFC(list(face = path), n_trials = 2, img_size = 16, stimulus_path = out,
                              seed = 1, ncores = 1, nscales = 1, save_as_png = FALSE))
  e <- new.env()
  load(list.files(out, pattern = "\\.Rdata$", full.names = TRUE), envir = e)
  grey <- (png::readPNG(path)[, , 1] + png::readPNG(path)[, , 2] + png::readPNG(path)[, , 3]) / 3
  expect_equal(e$base_faces$face, (grey - min(grey)) / (max(grey) - min(grey)))
  expect_identical(range(e$base_faces$face), c(0, 1))
})

test_that("each point of the cumulative curve correlates the CI so far with the final CI", {
  rdata <- make_fixture_rdata(withr::local_tempdir(), n_trials = 12)
  e <- new.env()
  load(rdata, envir = e)
  stimuli <- c(4, 9, 1, 12, 7, 3, 10, 2)
  responses <- c(1, -1, 1, 1, -1, -1, 1, -1)
  noise <- lapply(stimuli, function(s) generateNoiseImage(e$stimuli_params$base[s, ], e$p))
  running <- function(k) Reduce(`+`, Map(`*`, noise[seq_len(k)], responses[seq_len(k)])) / k
  expected <- vapply(seq(1, 8, by = 3), function(k) stats::cor(as.vector(running(k)), as.vector(running(8))),
                     numeric(1))
  curve <- quietly(computeCumulativeCICorrelation(stimuli, responses, "base", rdata, step = 3))
  expect_equal(as.vector(curve), expected, tolerance = 1e-12)
  expect_gt(length(unique(round(expected, 6))), 2)
})
