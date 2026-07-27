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

  # And it must beat chance. This is a one-sided Monte Carlo permutation test:
  # the statistic is the correlation above, and the null distribution comes
  # from permuting the response labels across trials. That keeps the noise
  # images and the count of 1s and -1s exactly as they are, and destroys only
  # the pairing between a response and the stimulus that produced it -- which
  # is the relationship generateCI() exploits, so it is the right thing to
  # break. Exchangeability holds by construction, since the per-trial
  # parameter vectors are drawn i.i.d.
  #
  # WHY A PERMUTATION TEST AND NOT A PARAMETRIC ONE ON r
  #
  # A t-test on Pearson r would take df = n_pixels - 2, here 1022, and be
  # wildly anticonservative, because the pixels are nowhere near independent.
  # The noise basis has only 60 patches, so the effective dimensionality of a
  # 32x32 CI is about 60, not 1024, and neighbouring pixels are strongly
  # autocorrelated by construction -- that is what a sinusoid basis *is*.
  #
  # The size of the error is not subtle. Permuting responses here produces
  # null correlations reaching |r| = 0.45, which a parametric test on 1022 df
  # would report as overwhelmingly significant. Any fixed threshold picked
  # without measuring this would be either flaky or vacuous.
  #
  # Every null CI is built from the same basis as the real one, so the null
  # inherits the identical spatial autocorrelation. That is what makes it the
  # correct reference distribution rather than merely a convenient one.
  #
  # B = 100 puts the smallest achievable p-value at 1/(B+1) = 0.0099. Note the
  # seed: this is a deterministic regression guard, not live inference, so
  # that figure is the resolution of the procedure and not a result. What
  # actually makes the test safe is the margin, measured across 10 different
  # templates at this B: the real CI beat the largest of 100 permutations by
  # 0.26 to 0.48, never less.
  set.seed(42)
  null_cors <- vapply(1:100, function(i) {
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
