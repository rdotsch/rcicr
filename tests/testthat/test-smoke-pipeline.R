# Skipped on CRAN: spawns parallel workers and writes PNGs. Keeps running in CI.
skip_on_cran()

test_that("generateStimuli2IFC produces stimuli, PNGs, and a loadable .Rdata end-to-end", {
  tmp <- withr::local_tempdir()
  input_dir <- file.path(tmp, "input")
  output_dir <- file.path(tmp, "output")
  dir.create(input_dir)
  dir.create(output_dir)

  png_path <- file.path(input_dir, "base.png")
  make_square_png(png_path, size = 32)

  generateStimuli2IFC(
    base_face_files = list(base = png_path),
    n_trials = 4,
    img_size = 32,
    stimulus_path = output_dir,
    seed = 1,
    ncores = 1,
    nscales = 1,
    save_as_png = TRUE,
    save_rdata = TRUE
  )

  pngs <- list.files(output_dir, pattern = "\\.png$")
  expect_equal(length(pngs), 4 * 2) # n_trials * (orig + inv), single base face

  rdata_files <- list.files(output_dir, pattern = "\\.Rdata$", full.names = TRUE)
  expect_length(rdata_files, 1)

  e <- new.env()
  load(rdata_files[1], envir = e)
  expect_true(all(c("p", "base_faces", "stimuli_params", "img_size", "seed") %in% ls(e)))
  expect_equal(e$img_size, 32)
})

test_that("the .Rdata records the rcicr version that actually wrote it", {
  # generator_version was a hardcoded '0.4.0' literal from 2016 until
  # 1.1.0.9000, so it described no file written after 0.4.0 and was useless for
  # the provenance purpose it exists for (issue #29). It must now agree with the
  # installed version, and with p$generator_version, which was already correct --
  # the two disagreeing is what made the field untrustworthy in the first place.
  tmp <- withr::local_tempdir()
  rdata_path <- make_fixture_rdata(tmp, img_size = 32, n_trials = 4, nscales = 1, seed = 1)

  e <- new.env()
  load(rdata_path, envir = e)

  installed <- as.character(utils::packageVersion("rcicr"))
  expect_equal(as.character(e$generator_version), installed)
  expect_equal(as.character(e$p$generator_version), installed)

  # The specific failure being guarded against: any run reporting the frozen
  # literal. Kept as its own assertion so a future re-freeze reads as such
  # rather than as a stale-installed-package problem.
  expect_false(identical(as.character(e$generator_version), "0.4.0"))
})

test_that("full pipeline: generateStimuli2IFC -> generateCI produces a classification image", {
  tmp <- withr::local_tempdir()
  input_dir <- file.path(tmp, "input")
  output_dir <- file.path(tmp, "output")
  dir.create(input_dir)
  dir.create(output_dir)

  png_path <- file.path(input_dir, "base.png")
  make_square_png(png_path, size = 32)

  generateStimuli2IFC(
    base_face_files = list(base = png_path),
    n_trials = 8, img_size = 32, stimulus_path = output_dir,
    seed = 1, ncores = 1, nscales = 1,
    save_as_png = FALSE, save_rdata = TRUE
  )
  rdata_path <- list.files(output_dir, pattern = "\\.Rdata$", full.names = TRUE)[1]

  set.seed(1)
  responses <- sample(c(1, -1), 8, replace = TRUE)

  ci <- generateCI(
    stimuli = 1:8, responses = responses, baseimage = "base",
    rdata = rdata_path, save_as_png = FALSE
  )

  expect_type(ci, "list")
  expect_equal(dim(ci$ci), c(32, 32))
  expect_false(anyNA(ci$ci))
})
