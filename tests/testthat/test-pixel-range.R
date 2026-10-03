# Values outside [0, 1] are clipped where an image is written or drawn: png's
# writePNG() wraps a value above 1 to a dark pixel (#371), and rasterImage()
# stops on either side (#373). Returned matrices are never clamped.

quietly <- function(expr) {
  out <- NULL
  utils::capture.output(out <- suppressMessages(suppressWarnings(expr)))
  out
}

# White but for one dark corner, so contrast maximization leaves it white.
white_base <- function(dir, size = 32) {
  base <- matrix(1, size, size)
  base[1:4, 1:4] <- 0
  path <- file.path(dir, "white.png")
  png::writePNG(base, path)
  path
}

test_that("clampUnit clips both ends and keeps NA", {
  expect_identical(rcicr:::clampUnit(c(-0.5, 0, 0.3, 1, 1.004, NA)), c(0, 0, 0.3, 1, 1, NA))
  m <- matrix(c(-1, 2, 0.5, NA), 2)
  expect_identical(dim(rcicr:::clampUnit(m)), dim(m))
})

test_that("stimulus pixels above 1 are written white, not wrapped to dark", {
  dir <- withr::local_tempdir()
  quietly(generateStimuli2IFC(list(face = white_base(dir)), n_trials = 6, img_size = 32,
                              stimulus_path = dir, seed = 1, ncores = 1, nscales = 1))
  e <- new.env()
  load(list.files(dir, "Rdata$", full.names = TRUE), envir = e)
  over <- 0
  for (trial in 1:6) {
    noise <- generateNoiseImage(e$stimuli_params$face[trial, ], e$p)
    for (side in c("ori", "inv")) {
      sign <- if (side == "ori") 1 else -1
      computed <- ((sign * noise + 0.3) / 0.6 + e$base_faces$face) / 2
      written <- png::readPNG(rcicr:::stimulusPngPath(dir, "rcic", "face", 1, trial, side))
      expect_equal(written, pmin(pmax(computed, 0), 1), tolerance = 1 / 255, info = side)
      over <- over + sum(computed > 1)
    }
  }
  # Without pixels above 1 this test could not fail.
  expect_gt(over, 0)
})

test_that("CI pixels above 1 are written white under 'constant' and weighted 'none'", {
  dir <- withr::local_tempdir()
  stim <- file.path(dir, "stimuli")
  quietly(generateStimuli2IFC(list(face = white_base(dir)), n_trials = 8, img_size = 32,
                              stimulus_path = stim, seed = 1, ncores = 1, nscales = 2,
                              save_as_png = FALSE))
  rdata <- list.files(stim, "Rdata$", full.names = TRUE)
  responses <- c(1, -1, -1, 1, 1, 1, -1, 1)
  for (case in list(list(scaling = "constant", scaling_constant = 0.01, responses = responses),
                    list(scaling = "none", scaling_constant = 0.1, responses = 40 * responses))) {
    out <- file.path(dir, case$scaling)
    ci <- quietly(generateCI(1:8, case$responses, "face", rdata, targetpath = out, n_cores = 1,
                             scaling = case$scaling, scaling_constant = case$scaling_constant))
    written <- png::readPNG(file.path(out, "ci_face.png"))
    expect_gt(sum(ci$combined > 1), 0)
    expect_true(all(written[ci$combined > 1] == 1), info = case$scaling)
    expect_equal(written, pmin(pmax(ci$combined, 0), 1), tolerance = 1 / 255)
    # The returned image keeps its values.
    expect_gt(max(ci$combined), 1)
  }
})

test_that("a z-map is drawn over a background outside [0, 1]", {
  rdata <- make_fixture_rdata(withr::local_tempdir(), n_trials = 8)
  out <- withr::local_tempdir()
  ci <- quietly(generateCI(1:8, c(1, -1, -1, 1, 1, 1, -1, 1), "base", rdata, save_as_png = FALSE,
                           scaling = "none", zmap = TRUE, zmapdecoration = FALSE,
                           zmaptargetpath = out, n_cores = 1))
  expect_lt(min(ci$combined), 0)
  expect_true(file.exists(file.path(out, "base.png")))

  zmap <- matrix(seq(-5, 5, length.out = 200 * 200), 200, 200)
  bg <- matrix(seq(-0.2, 1.2, length.out = 200 * 200), 200, 200)
  for (decoration in c(TRUE, FALSE)) {
    expect_no_error(plotZmap(zmap, bgimage = bg, sigma = 3, decoration = decoration,
                             targetpath = out, filename = paste0("z", decoration), size = 200))
  }
})

test_that("a nativeRaster background is drawn as it is, not clamped", {
  out <- withr::local_tempdir()
  path <- file.path(out, "bg.png")
  png::writePNG(matrix(seq(0, 1, length.out = 200 * 200), 200, 200), path)
  native <- png::readPNG(path, native = TRUE)
  zmap <- matrix(NA_real_, 200, 200)
  zmap[1:10, 1:10] <- 5
  draw <- function(bg, name) {
    plotZmap(zmap, bgimage = bg, sigma = 3, decoration = FALSE, targetpath = out,
             filename = name, size = 200)
    png::readPNG(file.path(out, paste0(name, ".png")))
  }
  expect_identical(draw(native, "native"), draw(png::readPNG(path), "numeric"))
})
