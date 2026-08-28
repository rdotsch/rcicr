# Does the method recover a representation that is actually there?
#
# Everything else in the latent module tests machinery. This tests the claim:
# an observer with a hidden target latent produces responses from which the
# direction to that target can be read back.

skip_on_cran()

test_that("the recovered direction points at the observer's target", {
  withr::local_options(rcicr.experimental = TRUE)
  out <- withr::local_tempdir()

  generator <- recovery_generator()
  stimuli <- generateStimuliLatent2IFC(generator, n_trials = 400, stimulus_path = out,
                                       seed = 1, save_as_png = FALSE)

  target <- stimuli$base_latent + generator$latent_sd * c(2, -1, 1.5, 0, -0.5, 1)
  responses <- simulate_latent_observer(stimuli$latent_params, stimuli$base_latent, target)

  ci <- generateLatentCI(stimuli = seq_len(400), responses = responses,
                         rdata = stimuli$rdata, save_as_png = FALSE)

  expected <- expected_latent_direction(generator, stimuli$base_latent, target)
  expect_gt(stats::cor(ci$direction, expected), 0.9)
})

test_that("permuting the responses destroys the recovery", {
  withr::local_options(rcicr.experimental = TRUE)
  out <- withr::local_tempdir()

  generator <- recovery_generator()
  stimuli <- generateStimuliLatent2IFC(generator, n_trials = 400, stimulus_path = out,
                                       seed = 2, save_as_png = FALSE)

  target <- stimuli$base_latent + generator$latent_sd * c(1, 1, -1, -1, 1, -1)
  responses <- simulate_latent_observer(stimuli$latent_params, stimuli$base_latent, target)
  expected <- expected_latent_direction(generator, stimuli$base_latent, target)

  real <- generateLatentCI(stimuli = seq_len(400), responses = responses,
                           rdata = stimuli$rdata, save_as_png = FALSE)

  withr::with_seed(9, {
    shuffled <- generateLatentCI(stimuli = seq_len(400), responses = sample(responses),
                                 rdata = stimuli$rdata, save_as_png = FALSE)
  })

  expect_gt(stats::cor(real$direction, expected), 0.9)
  expect_lt(abs(stats::cor(shuffled$direction, expected)), 0.6)
})

test_that("more trials recover the target more closely", {
  withr::local_options(rcicr.experimental = TRUE)
  out <- withr::local_tempdir()

  generator <- recovery_generator()
  stimuli <- generateStimuliLatent2IFC(generator, n_trials = 800, stimulus_path = out,
                                       seed = 3, save_as_png = FALSE)

  target <- stimuli$base_latent + generator$latent_sd * c(1.5, -2, 0.5, 1, -1, 0)
  responses <- simulate_latent_observer(stimuli$latent_params, stimuli$base_latent, target)
  expected <- expected_latent_direction(generator, stimuli$base_latent, target)

  correlation_at <- function(n) {
    ci <- generateLatentCI(stimuli = seq_len(n), responses = responses[seq_len(n)],
                           rdata = stimuli$rdata, save_as_png = FALSE)
    stats::cor(ci$direction, expected)
  }

  expect_gt(correlation_at(800), correlation_at(50))
})

test_that("informational value separates a real observer from a random one", {
  withr::local_options(rcicr.experimental = TRUE)
  out <- withr::local_tempdir()

  generator <- recovery_generator()
  stimuli <- generateStimuliLatent2IFC(generator, n_trials = 300, stimulus_path = out,
                                       seed = 4, save_as_png = FALSE)

  target <- stimuli$base_latent + generator$latent_sd * c(2, 1, -1, 1, -2, 0.5)
  responses <- simulate_latent_observer(stimuli$latent_params, stimuli$base_latent, target)

  real <- generateLatentCI(stimuli = seq_len(300), responses = responses,
                           rdata = stimuli$rdata, save_as_png = FALSE)
  random <- withr::with_seed(11, {
    generateLatentCI(stimuli = seq_len(300),
                     responses = sample(c(-1, 1), 300, replace = TRUE),
                     rdata = stimuli$rdata, save_as_png = FALSE)
  })

  real_iv <- computeLatentInfoVal2IFC(real, stimuli$rdata, seq_len(300),
                                      iter = 500, response_seed = 1)
  random_iv <- computeLatentInfoVal2IFC(random, stimuli$rdata, seq_len(300),
                                        iter = 500, response_seed = 1)

  expect_gt(real_iv$infoVal, random_iv$infoVal)
  expect_lt(real_iv$p, 0.01)
  expect_gt(random_iv$p, 0.01)
})

test_that("the rendered classification image moves towards the target's face", {
  withr::local_options(rcicr.experimental = TRUE)
  out <- withr::local_tempdir()

  generator <- recovery_generator()
  stimuli <- generateStimuliLatent2IFC(generator, n_trials = 400, stimulus_path = out,
                                       seed = 5, save_as_png = FALSE)

  target <- stimuli$base_latent + generator$latent_sd * c(2, -2, 1, -1, 0.5, 0)
  responses <- simulate_latent_observer(stimuli$latent_params, stimuli$base_latent, target)

  ci <- generateLatentCI(stimuli = seq_len(400), responses = responses,
                         rdata = stimuli$rdata, save_as_png = FALSE)

  # The point of working in a face space rather than in pixels: the answer is
  # judged as an image, so the rendered classification image has to be closer to
  # the target's face than the base face is.
  target_image <- renderLatent(generator, target)[1, , ]
  expect_lt(sum((ci$ci_image - target_image)^2),
            sum((ci$base_image - target_image)^2))
})

test_that("an anti-classification image points the other way", {
  withr::local_options(rcicr.experimental = TRUE)
  out <- withr::local_tempdir()

  generator <- recovery_generator()
  stimuli <- generateStimuliLatent2IFC(generator, n_trials = 200, stimulus_path = out,
                                       seed = 6, save_as_png = FALSE)

  target <- stimuli$base_latent + generator$latent_sd * c(1, -1, 1, -1, 1, -1)
  responses <- simulate_latent_observer(stimuli$latent_params, stimuli$base_latent, target)

  ci <- generateLatentCI(stimuli = seq_len(200), responses = responses,
                         rdata = stimuli$rdata, save_as_png = FALSE)
  anti <- generateLatentCI(stimuli = seq_len(200), responses = responses,
                           rdata = stimuli$rdata, save_as_png = FALSE, antiCI = TRUE)

  expect_equal(anti$direction, -ci$direction)
  expect_equal(stats::cor(anti$scaled_direction, ci$scaled_direction), -1)
})
