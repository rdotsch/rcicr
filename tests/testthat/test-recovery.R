# Does the pipeline actually work? Every other test checks a part in isolation:
# that a function returns the right shape, errors when it should, or matches a
# pinned number. None of them checks the thing the package exists to do, which
# is to recover an observer's internal template from their 2IFC choices.
#
# This test simulates an observer whose template is known exactly, feeds their
# choices through generateCI(), and asks whether the classification image looks
# like the template they were responding to. If a future refactor breaks the
# link between responses and parameters -- a sign flip, an off-by-one in the
# stimulus index, responses paired with the wrong rows -- the shape assertions
# elsewhere would all still pass and this is the test that would fail.
#
# What it does NOT catch, deliberately: any transformation applied consistently
# by generateNoiseImage() itself, such as a transposition, because the template
# is built through that same function and the error cancels. The oracle test in
# test-generateNoiseImage.R covers that; this one covers the wiring around it.
#
# Kept small on purpose: 32px, nscales = 2, 300 trials, about 6 seconds.

simulate_observer <- function(rdata, template_seed, response_seed, internal_noise = 1) {
  e <- new.env()
  load(rdata, envir = e)
  p <- e$p
  params <- e$stimuli_params[["base"]]
  n_trials <- nrow(params)

  # The observer's template is itself a noise image, so it is exactly
  # representable in the basis the stimuli are drawn from. That isolates what
  # is being tested (does the weighting recover it?) from a separate question
  # (can this basis express an arbitrary picture?).
  set.seed(template_seed)
  template <- generateNoiseImage(rnorm(max(p$patchIdx)), p)

  # On each trial the observer picks whichever of the pair looks more like the
  # template. The inverted stimulus carries exactly the negated noise, so the
  # sign of the match decides the choice; internal_noise makes them imperfect.
  stack <- vapply(seq_len(n_trials),
                  function(i) generateNoiseImage(params[i, ], p),
                  matrix(0, e$img_size, e$img_size))
  evidence <- apply(stack, 3, function(z) base::sum(z * template))
  evidence <- evidence / sd(evidence)

  set.seed(response_seed)
  responses <- ifelse(evidence + rnorm(n_trials, 0, internal_noise) > 0, 1, -1)

  list(template = template, responses = responses, n_trials = n_trials)
}

# One stimulus set for the whole file: generating it is the slow part, and all
# three tests want the same one. Torn down when the file finishes.
recovery_rdata <- local({
  dir <- withr::local_tempdir(.local_envir = testthat::teardown_env())
  png_path <- file.path(dir, "base.png")
  make_square_png(png_path, size = 32)
  suppressWarnings(
    generateStimuli2IFC(
      base_face_files = list(base = png_path),
      n_trials = 300, img_size = 32, stimulus_path = dir,
      seed = 1, ncores = 1, nscales = 2,
      save_as_png = FALSE, save_rdata = TRUE
    )
  )
  list.files(dir, pattern = "\\.Rdata$", full.names = TRUE)[1]
})

test_that("the CI recovers a simulated observer's template", {
  rdata <- recovery_rdata
  obs <- simulate_observer(rdata, template_seed = 1, response_seed = 1)

  ci <- generateCI(
    stimuli = seq_len(obs$n_trials), responses = obs$responses,
    baseimage = "base", rdata = rdata, save_as_png = FALSE
  )
  observed <- cor(as.vector(ci$ci), as.vector(obs$template))

  # Measured across 10 different templates the correlation ranged 0.71-0.84;
  # 0.5 leaves room for the run-to-run variation without being vacuous.
  expect_gt(observed, 0.5)

  # And it must beat chance. The null here is the same responses in a shuffled
  # order: that keeps the number of 1s and -1s identical and destroys only the
  # pairing between a response and the stimulus that produced it, which is
  # exactly the relationship the CI is supposed to exploit. A fixed threshold
  # would not do -- with 60 patches and 300 trials a random response vector can
  # reach |r| = 0.45 by chance, so the null is calibrated rather than assumed.
  set.seed(42)
  null_cors <- vapply(1:20, function(i) {
    shuffled <- generateCI(
      stimuli = seq_len(obs$n_trials), responses = sample(obs$responses),
      baseimage = "base", rdata = rdata, save_as_png = FALSE
    )
    cor(as.vector(shuffled$ci), as.vector(obs$template))
  }, numeric(1))

  expect_gt(observed, max(null_cors))
})

test_that("a more reliable observer yields a better CI than a noisier one", {
  rdata <- recovery_rdata
  clean <- simulate_observer(rdata, 2, 1, internal_noise = 0.2)
  noisy <- simulate_observer(rdata, 2, 1, internal_noise = 4.0)

  ci_of <- function(o) {
    ci <- generateCI(seq_len(o$n_trials), o$responses, "base", rdata,
                     save_as_png = FALSE)
    cor(as.vector(ci$ci), as.vector(o$template))
  }
  expect_gt(ci_of(clean), ci_of(noisy))
})

test_that("antiCI is the negation of the CI", {
  rdata <- recovery_rdata
  obs <- simulate_observer(rdata, template_seed = 3, response_seed = 1)

  args <- list(stimuli = seq_len(obs$n_trials), responses = obs$responses,
               baseimage = "base", rdata = rdata, save_as_png = FALSE)
  ci <- do.call(generateCI, args)
  anti <- do.call(generateCI, c(args, list(antiCI = TRUE)))

  expect_equal(anti$ci, -ci$ci)
  expect_lt(cor(as.vector(anti$ci), as.vector(obs$template)), -0.5)
})
