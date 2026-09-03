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

test_that("a resumed search keeps the design it was started with", {
  withr::local_options(rcicr.experimental = TRUE)

  generator <- recovery_generator(n_components = 4, n_faces = 8)
  out <- withr::local_tempdir()

  # Started with non-default settings, then resumed naming only the state and
  # the responses, which is how the documented call reads. The second generation
  # must not fall back to this call's defaults: sigma has to be 0.5 * 0.8, not 1.
  step <- searchLatent2IFC(generator, n_generations = 3, n_trials = 8,
                           stimulus_path = out, save_as_png = FALSE, seed = 1,
                           latent_sigma = 0.5, sigma_decay = 0.8)

  step <- searchLatent2IFC(generator, stimulus_path = out, save_as_png = FALSE,
                           state = step$state, responses = rep(c(1, -1), 4))

  expect_equal(nrow(step$latent_params), 8)
  expect_equal(step$remaining, 1)

  answers <- rep(c(1, -1), 4)
  final <- searchLatent2IFC(generator, stimulus_path = out, save_as_png = FALSE,
                            state = step$state, responses = answers)
  final <- searchLatent2IFC(generator, stimulus_path = out, save_as_png = FALSE,
                            state = final$state, responses = answers)

  expect_equal(final$generations$latent_sigma, c(0.5, 0.4, 0.32))
})

test_that("every setting survives a second resume, not just the first", {
  withr::local_options(rcicr.experimental = TRUE)

  generator <- recovery_generator(n_components = 4, n_faces = 8)
  out <- withr::local_tempdir()

  # The settings are written back into each new state, so one that is restored
  # for the current generation but not carried forward reverts to this call's
  # default a resume later. step_grow and step_shrink are the ones that do not
  # affect generation 2 at all and so only show up here.
  step <- searchLatent2IFC(generator, n_generations = 4, n_trials = 8,
                           stimulus_path = out, save_as_png = FALSE, seed = 1,
                           alpha = 3, step_grow = 2, step_shrink = 0.5)

  answers <- rep(c(1, -1), 4)
  for (i in 1:3) {
    step <- searchLatent2IFC(generator, stimulus_path = out, save_as_png = FALSE,
                             state = step$state, responses = answers)
  }
  final <- searchLatent2IFC(generator, stimulus_path = out, save_as_png = FALSE,
                            state = step$state, responses = answers)

  sizes <- final$generations$step_size
  cosines <- final$generations$cosine_with_previous
  expect_equal(sizes[1], 3)
  expect_equal(sizes[-1], sizes[-length(sizes)] * ifelse(cosines[-1] > 0, 2, 0.5))
})

test_that("a resumed search reproduces itself from the seed across sessions", {
  withr::local_options(rcicr.experimental = TRUE)

  generator <- recovery_generator(n_components = 4, n_faces = 8)
  answers <- rep(c(1, -1, 1, -1, -1, 1, 1, -1, -1, 1), 1)

  # set.seed(seed) runs only when the first state is written, so without the RNG
  # stream travelling in that state a search resumed in another session draws its
  # later generations from whatever that session's stream happened to be. Two
  # runs separated by an unrelated draw stand in for two sessions.
  run <- function(disturb) {
    out <- withr::local_tempdir(.local_envir = parent.frame())
    step <- searchLatent2IFC(generator, n_generations = 3, n_trials = 10,
                             stimulus_path = out, save_as_png = FALSE, seed = 7)
    if (disturb) {
      stats::rnorm(1000)
    }
    step <- searchLatent2IFC(generator, stimulus_path = out, save_as_png = FALSE,
                             state = step$state, responses = answers)
    step$latent_params
  }

  expect_equal(run(FALSE), run(TRUE))
})

test_that("a resumed generation rejects responses that are not 1 or -1", {
  withr::local_options(rcicr.experimental = TRUE)

  generator <- recovery_generator(n_components = 3, n_faces = 6)
  out <- withr::local_tempdir()

  step <- searchLatent2IFC(generator, n_generations = 2, n_trials = 4,
                           stimulus_path = out, save_as_png = FALSE)

  # Responses collected in another program arrive as whatever that program
  # wrote. A 0/1 coding would otherwise pass into the weighted mean and move the
  # search to the wrong centre with no error at all.
  expect_error(searchLatent2IFC(generator, stimulus_path = out, save_as_png = FALSE,
                                state = step$state, responses = c(0, 1, 1, 0)),
               'takes responses coded 1')
  expect_error(searchLatent2IFC(generator, stimulus_path = out, save_as_png = FALSE,
                                state = step$state, responses = c(1, -1, NA, 1)),
               'takes responses coded 1')
})

test_that("a resumed search keeps the seed it was started with", {
  withr::local_options(rcicr.experimental = TRUE)

  generator <- recovery_generator(n_components = 3, n_faces = 6)
  out <- withr::local_tempdir()

  # The seed names the stimulus files and the state file. Falling back to the
  # current call's default, two searches begun under different seeds and resumed
  # into one directory would overwrite each other's generations, and whichever
  # files survived would name a seed that did not produce them.
  step <- searchLatent2IFC(generator, n_generations = 3, n_trials = 4,
                           stimulus_path = out, save_as_png = FALSE, seed = 77)
  expect_true(grepl('seed_77', basename(step$state)))

  step <- searchLatent2IFC(generator, stimulus_path = out, save_as_png = FALSE,
                           state = step$state, responses = rep(c(1, -1), 2))
  expect_true(grepl('seed_77', basename(step$state)))
})

test_that("callback mode renders each stimulus once, in batches", {
  withr::local_options(rcicr.experimental = TRUE)
  base <- recovery_generator(n_components = 4, n_faces = 8, size = 8)

  # A generator that records every call, so the number of renders and the size
  # of each is observable rather than inferred.
  calls <- new.env(parent = emptyenv())
  calls$sizes <- integer(0)
  counting <- base
  counting$render <- function(latents) {
    calls$sizes <- c(calls$sizes, nrow(latents))
    base$render(latents)
  }

  out <- withr::local_tempdir()
  searchLatent2IFC(counting, n_generations = 1, n_trials = 10, stimulus_path = out,
                   respond = function(trial) 1, batch_size = 4, seed = 1)

  # No call anywhere exceeds batch_size. Rendering a whole generation at once
  # gave calls of 10, which is what a renderer that has to batch cannot survive.
  expect_true(all(calls$sizes <= 4))

  # Single-latent calls are validateGenerator()'s contract probe, which renders
  # latent_mean before the search begins. The generation's own renders are the
  # rest: ten trials at batch_size 4 is three batches, each rendering the
  # originals and then the inversions, so 4, 4, 4, 4, 2, 2.
  batches <- calls$sizes[calls$sizes > 1]
  expect_identical(sort(batches), c(2L, 2L, 4L, 4L, 4L, 4L))

  # Twenty images for ten trials: each original and each inversion rendered once.
  # Rendering separately for the archive and for the callback doubled this.
  expect_identical(sum(batches), 20L)
})

test_that("the callback sees the pixels that were archived", {
  withr::local_options(rcicr.experimental = TRUE)
  base <- recovery_generator(n_components = 4, n_faces = 8, size = 8)

  # A generator whose output changes on every call. Rendering once for the PNG
  # and again for the callback would show a participant an image the saved file
  # does not match, which nothing downstream could detect.
  drift <- new.env(parent = emptyenv())
  drift$n <- 0
  shifting <- base
  shifting$render <- function(latents) {
    drift$n <- drift$n + 1
    rendered <- base$render(latents)
    pmin(pmax(rendered + drift$n / 1000, 0), 1)
  }

  out <- withr::local_tempdir()
  seen <- new.env(parent = emptyenv())
  seen$originals <- list()

  searchLatent2IFC(shifting, n_generations = 1, n_trials = 4, stimulus_path = out,
                   respond = function(trial) {
                     seen$originals[[trial$trial]] <- trial$original
                     1
                   },
                   batch_size = 2, seed = 1)

  files <- sort(list.files(out, pattern = '_ori[.]png$', full.names = TRUE))
  expect_length(files, 4)
  for (i in seq_along(files)) {
    expect_equal(png::readPNG(files[i]), seen$originals[[i]], tolerance = 1 / 255)
  }
})

test_that("two searches in one directory do not overwrite each other's state", {
  withr::local_options(rcicr.experimental = TRUE)
  generator <- recovery_generator(n_components = 4, n_faces = 8, size = 8)
  out <- withr::local_tempdir()

  # Same stimulus_path and the same seed, which is the default, so the seed
  # cannot tell the two apart. Every generation used to write the same state
  # filename: the second search overwrote the first, and resuming the first then
  # applied its responses to the second's centre with nothing to say so.
  first <- searchLatent2IFC(generator, n_generations = 3, n_trials = 4,
                            stimulus_path = out, save_as_png = FALSE)
  second <- searchLatent2IFC(generator, n_generations = 3, n_trials = 4,
                             stimulus_path = out, save_as_png = FALSE,
                             base_latent = generator$latent_mean + generator$latent_sd)

  expect_false(first$state == second$state)
  expect_true(file.exists(first$state))
  expect_true(file.exists(second$state))

  # And the first search's own centre survives, rather than the second's.
  centre_of <- function(file) {
    env <- new.env()
    load(file, envir = env)
    env$search_state$centre
  }
  expect_equal(centre_of(first$state), generator$latent_mean)
  expect_false(isTRUE(all.equal(centre_of(second$state), generator$latent_mean)))

  # A resume keeps its own run, so generation 2 does not collide either.
  resumed <- searchLatent2IFC(generator, stimulus_path = out, save_as_png = FALSE,
                              state = first$state, responses = rep(c(1, -1), 2))
  expect_false(resumed$state == second$state)
  expect_true(grepl(sub('^search_([^_]+)_.*', '\\1', basename(first$state)),
                    basename(resumed$state), fixed = TRUE))
})

test_that("rendering settings survive a resume", {
  withr::local_options(rcicr.experimental = TRUE)
  generator <- recovery_generator(n_components = 4, n_faces = 8, size = 8)
  out <- withr::local_tempdir()

  # batch_size and save_as_png have defaults, which is what makes leaving them
  # out of the state dangerous: the documented resume names only the state and
  # the responses, so both would quietly revert.
  calls <- new.env(parent = emptyenv())
  calls$sizes <- integer(0)
  counting <- generator
  counting$render <- function(latents) {
    calls$sizes <- c(calls$sizes, nrow(latents))
    generator$render(latents)
  }

  step <- searchLatent2IFC(counting, n_generations = 3, n_trials = 8,
                           stimulus_path = out, batch_size = 2, save_as_png = FALSE,
                           seed = 1)
  expect_length(list.files(out, pattern = '_ori[.]png$'), 0)

  calls$sizes <- integer(0)
  step <- searchLatent2IFC(counting, stimulus_path = out, state = step$state,
                           responses = rep(c(1, -1), 4))

  # Still no images, and still rendered two at a time rather than eight.
  expect_length(list.files(out, pattern = '_ori[.]png$'), 0)
  expect_true(all(calls$sizes <= 2))
})

test_that("resuming one generation twice does not overwrite the first result", {
  withr::local_options(rcicr.experimental = TRUE)
  generator <- recovery_generator(n_components = 4, n_faces = 8, size = 8)
  out <- withr::local_tempdir()

  first <- searchLatent2IFC(generator, n_generations = 3, n_trials = 4,
                            stimulus_path = out, save_as_png = FALSE)

  # Two resumes of the same state with different answers: comparing two codings,
  # or redoing a mistaken batch. run_id, generation and seed are all unchanged,
  # so the second used to overwrite the first and the path the first returned
  # then described the second's responses.
  one <- searchLatent2IFC(generator, stimulus_path = out, state = first$state,
                          responses = c(1, 1, -1, -1))
  two <- searchLatent2IFC(generator, stimulus_path = out, state = first$state,
                          responses = c(-1, -1, 1, 1))

  expect_false(one$state == two$state)
  expect_true(file.exists(one$state))

  centre_of <- function(file) {
    env <- new.env()
    load(file, envir = env)
    env$search_state$centre
  }

  # Opposite answers move the centre opposite ways, so the two states must differ
  # and the first must still hold its own.
  expect_false(isTRUE(all.equal(centre_of(one$state), centre_of(two$state))))
})

test_that("a second resume refuses to replace the images of the first", {
  withr::local_options(rcicr.experimental = TRUE)
  generator <- recovery_generator(n_components = 4, n_faces = 8, size = 8)
  out <- withr::local_tempdir()

  first <- searchLatent2IFC(generator, n_generations = 3, n_trials = 2,
                            stimulus_path = out, save_as_png = TRUE)

  searchLatent2IFC(generator, stimulus_path = out, state = first$state,
                   responses = c(1, -1))

  # Unlike the state file, generation images cannot be renamed to stay unique:
  # their names are what a task script builds to show generation g, trial 42.
  expect_error(
    searchLatent2IFC(generator, stimulus_path = out, state = first$state,
                     responses = c(-1, 1)),
    'would overwrite'
  )
})
