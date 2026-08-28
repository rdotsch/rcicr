# The adaptive search. The first test is the one that decides whether this mode
# should exist at all; the rest are its machinery.

test_that("the search finds where a representation is, not just which way it lies", {
  skip_on_cran()
  withr::local_options(rcicr.experimental = TRUE)

  generator <- recovery_generator(n_components = 6, n_faces = 14)
  out <- withr::local_tempdir()
  budget <- 200

  # A single round of trials carries no information about distance: the constant
  # that turns the estimated direction into a face is the researcher's choice.
  # The claim for the adaptive mode is that re-centring supplies the distance.
  # If this does not hold the mode should not ship, so it is measured here
  # rather than asserted in the documentation.
  closer <- 0
  closer_than_oracle <- 0
  for (rep in 1:5) {
    withr::with_seed(rep, {
      base <- generator$latent_mean
      target <- base + generator$latent_sd * stats::rnorm(6, sd = 2)
    })

    stimuli <- generateStimuliLatent2IFC(generator, n_trials = budget, stimulus_path = out,
                                         seed = 100 + rep, save_as_png = FALSE)
    responses <- simulate_latent_observer(stimuli$latent_params, base, target)
    one_shot <- generateLatentCI(stimuli = seq_len(budget), responses = responses,
                                 rdata = stimuli$rdata, save_as_png = FALSE)

    adaptive <- searchLatent2IFC(generator, n_generations = 10, n_trials = budget / 10,
                                 respond = latent_observer_callback(target),
                                 save_as_png = FALSE, seed = 200 + rep)

    reached <- latent_sd_distance(generator, adaptive$latent, target)

    if (reached < latent_sd_distance(generator, one_shot$latent_ci, target)) {
      closer <- closer + 1
    }

    # The stronger comparison: the best distance any scaling constant could have
    # given the one-shot method, chosen with the answer already in hand. A
    # researcher cannot pick that constant, so this is a bound the one-shot
    # method cannot reach in practice, and the search still has to beat it.
    oracle <- min(vapply(seq(0.5, 8, by = 0.5), function(constant) {
      ci <- generateLatentCI(stimuli = seq_len(budget), responses = responses,
                             rdata = stimuli$rdata, save_as_png = FALSE,
                             latent_scaling = 'sd', scaling_constant = constant)
      latent_sd_distance(generator, ci$latent_ci, target)
    }, numeric(1)))

    if (reached < oracle) {
      closer_than_oracle <- closer_than_oracle + 1
    }
  }

  expect_equal(closer, 5)
  expect_equal(closer_than_oracle, 5)
})

test_that("the search moves towards the target and the step shrinks as it arrives", {
  skip_on_cran()
  withr::local_options(rcicr.experimental = TRUE)

  generator <- recovery_generator(n_components = 6, n_faces = 14)
  withr::with_seed(3, target <- generator$latent_mean + generator$latent_sd * stats::rnorm(6, sd = 2))

  search <- searchLatent2IFC(generator, n_generations = 8, n_trials = 40,
                             respond = latent_observer_callback(target),
                             save_as_png = FALSE, seed = 1)

  start <- latent_sd_distance(generator, search$base_latent, target)
  finish <- latent_sd_distance(generator, search$latent, target)
  expect_lt(finish, start / 3)

  # The step is not on a schedule: it grows while generations agree and shrinks
  # when they disagree. Having arrived, the last steps have to be smaller than
  # the largest one taken on the way.
  expect_lt(mean(utils::tail(search$generations$step_size, 3)),
            max(search$generations$step_size))
})

test_that("the trajectory and diagnostics describe every generation", {
  withr::local_options(rcicr.experimental = TRUE)

  generator <- recovery_generator(n_components = 4, n_faces = 8)
  withr::with_seed(5, target <- generator$latent_mean + generator$latent_sd * stats::rnorm(4))

  search <- searchLatent2IFC(generator, n_generations = 4, n_trials = 12,
                             respond = latent_observer_callback(target),
                             save_as_png = FALSE, seed = 1)

  expect_equal(dim(search$trajectory), c(5, 4))
  expect_equal(search$trajectory[1, ], generator$latent_mean)
  expect_equal(search$trajectory[5, ], search$latent)

  expect_equal(nrow(search$generations), 4)
  expect_named(search$generations, c('generation', 'latent_sigma', 'step_size',
                                     'direction_norm', 'cosine_with_previous'))
  expect_true(is.na(search$generations$cosine_with_previous[1]))
  expect_false(any(is.na(search$generations$cosine_with_previous[-1])))
  expect_equal(dim(search$latent_image), c(16, 16))
})

test_that("sigma_decay anneals the perturbation size and 1 holds it", {
  withr::local_options(rcicr.experimental = TRUE)

  generator <- recovery_generator(n_components = 4, n_faces = 8)
  respond <- function(trial) 1

  annealed <- searchLatent2IFC(generator, n_generations = 4, n_trials = 6,
                               respond = respond, save_as_png = FALSE,
                               sigma_decay = 0.5, latent_sigma = 2)
  expect_equal(annealed$generations$latent_sigma, c(2, 1, 0.5, 0.25))

  held <- searchLatent2IFC(generator, n_generations = 3, n_trials = 6,
                           respond = respond, save_as_png = FALSE)
  expect_equal(held$generations$latent_sigma, c(1, 1, 1))
})

test_that("stimuli are written per generation when asked", {
  withr::local_options(rcicr.experimental = TRUE)

  generator <- recovery_generator(n_components = 3, n_faces = 6)
  out <- withr::local_tempdir()

  searchLatent2IFC(generator, n_generations = 2, n_trials = 3, stimulus_path = out,
                   respond = function(trial) 1, seed = 1)

  expect_equal(length(list.files(out, pattern = 'gen001.*_ori\\.png$')), 3)
  expect_equal(length(list.files(out, pattern = 'gen002.*_inv\\.png$')), 3)
})

test_that("a respond function that answers with anything else is rejected", {
  withr::local_options(rcicr.experimental = TRUE)

  generator <- recovery_generator(n_components = 3, n_faces = 6)

  expect_error(
    searchLatent2IFC(generator, n_generations = 1, n_trials = 3,
                     respond = function(trial) 0, save_as_png = FALSE),
    'respond must return 1'
  )
})

test_that("resumable mode writes a generation and takes its responses back", {
  withr::local_options(rcicr.experimental = TRUE)

  generator <- recovery_generator(n_components = 4, n_faces = 8)
  out <- withr::local_tempdir()
  withr::with_seed(7, target <- generator$latent_mean + generator$latent_sd * stats::rnorm(4, sd = 2))

  step <- searchLatent2IFC(generator, n_generations = 3, n_trials = 20,
                           stimulus_path = out, save_as_png = FALSE, seed = 1)
  expect_equal(step$generation, 1)
  expect_equal(step$remaining, 2)
  expect_true(file.exists(step$state))

  centres <- list(step$centre)
  while (step$remaining > 0) {
    answers <- simulate_latent_observer(step$latent_params, step$centre, target)
    step <- searchLatent2IFC(generator, n_generations = 3, n_trials = 20,
                             stimulus_path = out, save_as_png = FALSE, seed = 1,
                             state = step$state, responses = answers)
    if (step$remaining > 0) {
      centres[[length(centres) + 1]] <- step$centre
    }
  }

  # The last call, with nothing left to run, returns the finished search.
  answers <- simulate_latent_observer(step$latent_params, step$centre, target)
  final <- searchLatent2IFC(generator, n_generations = 3, n_trials = 20,
                            stimulus_path = out, save_as_png = FALSE, seed = 1,
                            state = step$state, responses = answers)

  expect_equal(final$remaining, 0)
  expect_null(final$state)
  expect_equal(nrow(final$trajectory), 4)
  expect_equal(nrow(final$generations), 3)
  expect_lt(latent_sd_distance(generator, final$latent, target),
            latent_sd_distance(generator, generator$latent_mean, target))
})

test_that("resumable mode says what it is missing", {
  withr::local_options(rcicr.experimental = TRUE)

  generator <- recovery_generator(n_components = 3, n_faces = 6)
  out <- withr::local_tempdir()

  step <- searchLatent2IFC(generator, n_generations = 2, n_trials = 5,
                           stimulus_path = out, save_as_png = FALSE)

  err <- expect_error(searchLatent2IFC(generator, n_generations = 2, n_trials = 5,
                                       stimulus_path = out, save_as_png = FALSE,
                                       state = step$state))
  expect_match(conditionMessage(err), 'no responses')
  expect_match(conditionMessage(err), 'the 5 answers', fixed = TRUE)

  expect_error(searchLatent2IFC(generator, n_generations = 2, n_trials = 5,
                                stimulus_path = out, save_as_png = FALSE,
                                state = step$state, responses = rep(1, 3)),
               'had 5 trials but 3 responses')

  expect_error(searchLatent2IFC(generator, n_generations = 2, n_trials = 5,
                                stimulus_path = out, save_as_png = FALSE,
                                responses = rep(1, 5)),
               'without a state file')
})

test_that("a state file from another generator is refused", {
  withr::local_options(rcicr.experimental = TRUE)

  generator <- recovery_generator(n_components = 3, n_faces = 6)
  other <- recovery_generator(n_components = 3, n_faces = 9)
  out <- withr::local_tempdir()

  step <- searchLatent2IFC(generator, n_generations = 2, n_trials = 5,
                           stimulus_path = out, save_as_png = FALSE)

  expect_error(searchLatent2IFC(other, n_generations = 2, n_trials = 5,
                                stimulus_path = out, save_as_png = FALSE,
                                state = step$state, responses = rep(1, 5)),
               'not the one that made these stimuli')
})

test_that("a file that is not a search state says so", {
  withr::local_options(rcicr.experimental = TRUE)

  generator <- recovery_generator(n_components = 3, n_faces = 6)
  fx <- latent_fixture()

  expect_error(searchLatent2IFC(generator, n_generations = 2, n_trials = 5,
                                stimulus_path = withr::local_tempdir(),
                                save_as_png = FALSE, state = fx$stimuli$rdata,
                                responses = rep(1, 5)),
               'does not hold a search state')
})

test_that("respond and state cannot both be given", {
  withr::local_options(rcicr.experimental = TRUE)

  generator <- recovery_generator(n_components = 3, n_faces = 6)

  expect_error(searchLatent2IFC(generator, n_generations = 1, n_trials = 3,
                                save_as_png = FALSE, respond = function(trial) 1,
                                state = 'somewhere.Rdata'),
               'either respond or state')
})

test_that("stimulus_path is required unless nothing is written", {
  withr::local_options(rcicr.experimental = TRUE)

  generator <- recovery_generator(n_components = 3, n_faces = 6)

  expect_error(searchLatent2IFC(generator, n_generations = 1, n_trials = 3,
                                respond = function(trial) 1),
               'No stimulus_path')
  expect_error(searchLatent2IFC(generator, n_generations = 1, n_trials = 3),
               'No stimulus_path')
})

test_that("a base latent of the wrong width is rejected", {
  withr::local_options(rcicr.experimental = TRUE)

  generator <- recovery_generator(n_components = 3, n_faces = 6)

  expect_error(searchLatent2IFC(generator, n_generations = 1, n_trials = 3,
                                respond = function(trial) 1, save_as_png = FALSE,
                                base_latent = c(1, 2)),
               'latent_dim is 3')
})

test_that("the step grows while generations agree and shrinks when they do not", {
  withr::local_options(rcicr.experimental = TRUE)

  generator <- recovery_generator(n_components = 4, n_faces = 8)

  withr::with_seed(4, target <- generator$latent_mean + generator$latent_sd * stats::rnorm(4, sd = 3))
  run <- searchLatent2IFC(generator, n_generations = 8, n_trials = 20,
                          respond = latent_observer_callback(target),
                          save_as_png = FALSE, seed = 1,
                          alpha = 1, step_grow = 2, step_shrink = 0.5)

  # Every step after the first follows the sign of that generation's agreement
  # with the one before it, and nothing else.
  sizes <- run$generations$step_size
  cosines <- run$generations$cosine_with_previous
  expect_equal(sizes[1], 1)
  expect_equal(sizes[-1],
               sizes[-length(sizes)] * ifelse(cosines[-1] > 0, 2, 0.5))

  # A target close to the start makes the search step past it and reverse, so
  # the step has to come back down rather than run away.
  withr::with_seed(2, target <- generator$latent_mean + generator$latent_sd * stats::rnorm(4, sd = 0.2))
  closing <- searchLatent2IFC(generator, n_generations = 12, n_trials = 30,
                              respond = latent_observer_callback(target),
                              save_as_png = FALSE, seed = 1, alpha = 4)

  expect_lt(min(closing$generations$step_size), 4)
  expect_lt(utils::tail(closing$generations$step_size, 1), max(closing$generations$step_size))
})

test_that("the search reaches a distant representation without alpha being tuned for it", {
  skip_on_cran()
  withr::local_options(rcicr.experimental = TRUE)

  generator <- recovery_generator(n_components = 6, n_faces = 14)

  # The step size cannot be chosen in advance, because how far the
  # representation lies from the base face is exactly what is unknown. The
  # adaptive step has to cope with both without being told which it is facing.
  reached <- function(spread) {
    withr::with_seed(round(spread * 10),
                     target <- generator$latent_mean + generator$latent_sd * stats::rnorm(6, sd = spread))
    search <- searchLatent2IFC(generator, n_generations = 12, n_trials = 30,
                               respond = latent_observer_callback(target),
                               save_as_png = FALSE, seed = 1)
    latent_sd_distance(generator, search$latent, target) /
      latent_sd_distance(generator, generator$latent_mean, target)
  }

  # Closes at least 80 percent of the distance whether the target is near or far.
  expect_lt(reached(0.5), 0.2)
  expect_lt(reached(5), 0.2)
})
