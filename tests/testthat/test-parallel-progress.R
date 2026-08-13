# The progress bars are created in the parent but were ticked inside the foreach
# body, so under parallel execution -- the default -- each worker updated its own
# private copy and the user watched 0% until the job finished (issue #178).
#
# What makes this observable is that capture.output() sinks the parent's output
# while the workers' (makeCluster(outfile = "")) goes straight to the terminal,
# past the sink. So the percentages captured here are exactly the ones a user
# would have seen, which is the property that was broken.
#
# The assertion is "more than one distinct percentage", not an exact sequence:
# tasks complete in whatever order the cluster returns them, and pinning the
# ticks would buy nothing and fail intermittently.

progress_pcts <- function(expr) {
  out <- capture.output(suppressWarnings(expr), type = "output")
  unique(unlist(regmatches(out, gregexpr("[0-9]+%", out))))
}

expect_progress <- function(pcts, label) {
  expect_gt(length(pcts), 1)
  expect_true("100%" %in% pcts, info = label)
}

test_that("generateStimuli2IFC reports progress serially and in parallel", {
  tmp <- withr::local_tempdir()
  base <- file.path(tmp, "base.png")
  make_square_png(base, size = 32)

  generate_with <- function(ncores) {
    dir <- file.path(tmp, paste0("run", ncores))
    dir.create(dir)
    progress_pcts(generateStimuli2IFC(
      base_face_files = list(base = base),
      n_trials = 6, img_size = 32, stimulus_path = dir,
      seed = 1, ncores = ncores, nscales = 1,
      save_as_png = FALSE, save_rdata = TRUE
    ))
  }

  expect_progress(generate_with(1), "ncores = 1")

  skip_on_cran()
  expect_progress(generate_with(2), "ncores = 2")
})

test_that("generateCI reports progress from both of its parallel loops", {
  tmp <- withr::local_tempdir()
  rdata <- make_fixture_rdata(tmp, img_size = 32, n_trials = 6)
  responses <- rep(c(1, -1), 3)

  # The participant-CI loop and the z-map loop are separate foreach calls with
  # separate bars; only `participants` reaches the first one.
  participant_ci <- function(n_cores) {
    progress_pcts(generateCI(
      1:6, responses, "base", rdata, save_as_png = FALSE,
      participants = c(1, 1, 2, 2, 3, 3), n_cores = n_cores
    ))
  }

  # Only the 't.test' method builds a z-map in a foreach loop; the default
  # 'quick' method has no loop and so no progress bar to fix.
  zmap_ci <- function(n_cores) {
    progress_pcts(generateCI(
      1:6, responses, "base", rdata, save_as_png = FALSE,
      zmap = TRUE, zmapmethod = "t.test", zmapdecoration = FALSE,
      zmaptargetpath = file.path(tmp, paste0("z", n_cores)), n_cores = n_cores
    ))
  }

  expect_progress(participant_ci(1), "participants, n_cores = 1")
  expect_progress(zmap_ci(1), "zmap, n_cores = 1")

  skip_on_cran()
  expect_progress(participant_ci(2), "participants, n_cores = 2")
  expect_progress(zmap_ci(2), "zmap, n_cores = 2")
})
