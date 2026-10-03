
test_that("an n_trials that is not a positive whole number stops before anything is written", {
  base <- tempfile(fileext = ".png")
  make_square_png(base, size = 16) # nolint: object_usage_linter.
  for (n in list(0, -1, 2.5, NA, c(2, 3), "4")) {
    for (png in c(TRUE, FALSE)) {
      dir <- file.path(withr::local_tempdir(), "stimuli")
      expect_error(generateStimuli2IFC(list(face = base), n_trials = n, img_size = 16,
                                       stimulus_path = dir, ncores = 1, nscales = 1,
                                       save_as_png = png),
                   "n_trials must be a positive whole number", info = paste(format(n), png))
      expect_false(dir.exists(dir))
    }
  }
})
