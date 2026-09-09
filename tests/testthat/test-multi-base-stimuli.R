# Stimulus writing across several base images.
#
# generateStimuli2IFC() loops over base faces inside its per-trial foreach body;
# these pin what that loop must produce for every base image, and what it must
# not recompute when nothing is written.

test_that("generateStimuli2IFC writes every base image's stimuli when returning a data frame", {
  # Issue #302. The returned noise used to leave the foreach body from inside
  # the loop over base faces, so every base after the first was skipped and the
  # stimulus set was silently incomplete.
  for (shared in c(TRUE, FALSE)) {
    tmp <- withr::local_tempdir()
    bases <- list(
      one = make_square_png(file.path(tmp, "one.png"), size = 32, seed = 1),
      two = make_square_png(file.path(tmp, "two.png"), size = 32, seed = 2)
    )
    out <- file.path(tmp, "stimuli")

    frame <- generateStimuli2IFC(
      base_face_files = bases, n_trials = 2, img_size = 32, stimulus_path = out,
      seed = 31, nscales = 1, use_same_parameters = shared, ncores = 1,
      return_as_dataframe = TRUE, save_as_png = TRUE, save_rdata = TRUE
    )

    files <- list.files(out, pattern = "\\.png$")
    expect_length(files, 8)
    expect_length(grep("_one_", files), 4)
    expect_length(grep("_two_", files), 4)

    e <- new.env()
    load(list.files(out, pattern = "\\.Rdata$", full.names = TRUE), envir = e)

    # Every base is checked against its own parameter matrix. Comparing base two
    # only against base one's stimuli would pass even if it were written from
    # base one's noise, because the base faces differ.
    if (!shared) {
      expect_gt(max(abs(e$stimuli_params$one - e$stimuli_params$two)), 0.01)
    }

    signs <- c(ori = 1, inv = -1)
    for (base in c("one", "two")) {
      for (trial in 1:2) {
        noise <- generateNoiseImage(e$stimuli_params[[base]][trial, ], e$p)
        for (suffix in names(signs)) {
          expected <- ((signs[[suffix]] * noise + 0.3) / 0.6 + e$base_faces[[base]]) / 2
          stimulus <- sprintf("rcic_%s_31_%05d_%s.png", base, trial, suffix)
          written <- png::readPNG(file.path(out, stimulus))
          # Compared through the same writer rather than against the matrix:
          # 8-bit quantisation, and what png does with the few pixels that noise
          # beyond the nominal +/-0.3 puts outside [0, 1], then apply to both.
          reference <- tempfile(fileext = ".png")
          png::writePNG(expected, reference)
          expect_equal(written, png::readPNG(reference),
                       info = paste(base, trial, suffix, "shared:", shared))
        }
      }
    }

    # The returned frame is unchanged: one noise image per trial, the first
    # base's, exactly as before the fix.
    first_noise <- vapply(1:2, function(i) {
      as.vector(generateNoiseImage(e$stimuli_params$one[i, ], e$p))
    }, numeric(1024))
    expect_equal(dim(frame), c(1024L, 2L))
    expect_equal(as.matrix(frame), first_noise, ignore_attr = TRUE)
  }
})

test_that("generateStimuli2IFC computes noise once per trial when it writes no PNGs", {
  # Issue #302's fix must not make reference generation pay for base images
  # whose noise it never uses: it asks for the data frame and no PNGs, and the
  # frame holds only the first base image's noise.
  tmp <- withr::local_tempdir()
  bases <- lapply(setNames(1:3, c("one", "two", "three")), function(i) {
    make_square_png(file.path(tmp, paste0(i, ".png")), size = 32, seed = i)
  })

  calls <- 0L
  trace(generateNoiseImage, tracer = function() calls <<- calls + 1L,
        print = FALSE, where = asNamespace("rcicr"))
  on.exit(untrace(generateNoiseImage, where = asNamespace("rcicr")), add = TRUE)

  invisible(generateStimuli2IFC(
    base_face_files = bases, n_trials = 4, img_size = 32, seed = 31, nscales = 1,
    use_same_parameters = FALSE, ncores = 1, return_as_dataframe = TRUE,
    save_as_png = FALSE, save_rdata = FALSE
  ))

  expect_identical(calls, 4L)
})
