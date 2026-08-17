# The ncores == 1 path runs the foreach loops in the current process via
# registerDoSEQ() instead of building a one-worker cluster. That is a change of
# execution path, so the property that has to hold is that it changes nothing
# about the result.
#
# It is safe for a specific reason worth stating: neither parallel loop in this
# package draws random numbers. Stimulus parameters are drawn under set.seed()
# *before* the loop begins, so there is no per-worker RNG stream that could
# diverge from the sequential one. If a future change introduces a random draw
# inside one of those loops, this test is what should catch it.

test_that("ncores = 1 and ncores = 2 produce identical stimuli and CIs", {
  skip_on_cran()

  tmp <- withr::local_tempdir()
  base <- file.path(tmp, "base.png")
  make_square_png(base, size = 32)

  generate_with <- function(ncores) {
    dir <- file.path(tmp, paste0("run", ncores))
    dir.create(dir)
    suppressWarnings(
      generateStimuli2IFC(
        base_face_files = list(base = base),
        n_trials = 20, img_size = 32, stimulus_path = dir,
        seed = 1, ncores = ncores, nscales = 2,
        save_as_png = FALSE, save_rdata = TRUE
      )
    )
    list.files(dir, pattern = "\\.Rdata$", full.names = TRUE)[1]
  }

  serial_rdata <- generate_with(1)
  parallel_rdata <- generate_with(2)

  serial <- new.env()
  parallel_env <- new.env()
  load(serial_rdata, envir = serial)
  load(parallel_rdata, envir = parallel_env)

  # The stimulus parameters are the entire contract with later analysis: if
  # these match, every classification image computed from them matches too.
  expect_identical(serial$stimuli_params, parallel_env$stimuli_params)
  expect_identical(serial$p, parallel_env$p)

  # ...and check that a CI computed from either file matches. Note this does
  # *not* exercise generateCI()'s own parallel loop: with no `participants` it
  # takes the single-shot generateCINoise() path, which has no foreach in it.
  # That loop is covered by test-participant-ci-equivalence.R.
  responses <- rep(c(1, -1), 10)
  ci_serial <- generateCI(1:20, responses, "base", serial_rdata, save_as_png = FALSE)
  ci_parallel <- generateCI(1:20, responses, "base", parallel_rdata, save_as_png = FALSE)

  expect_identical(ci_serial$ci, ci_parallel$ci)
})
