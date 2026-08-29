# generateLatentCI()'s machinery: scaling, participants, files and refusals.

test_that("the direction is the response-weighted mean of the perturbations", {
  withr::local_options(rcicr.experimental = TRUE)
  fx <- latent_fixture()
  responses <- rep(c(1, -1), 10)

  ci <- generateLatentCI(stimuli = 1:20, responses = responses,
                         rdata = fx$stimuli$rdata, save_as_png = FALSE)

  expect_equal(ci$direction, colMeans(fx$stimuli$latent_params * responses))
  expect_equal(ci$latent_ci, fx$stimuli$base_latent + ci$scaled_direction)
})

test_that("the classification latent renders to the returned image", {
  withr::local_options(rcicr.experimental = TRUE)
  fx <- latent_fixture()

  ci <- generateLatentCI(stimuli = 1:20, responses = rep(c(1, -1), 10),
                         rdata = fx$stimuli$rdata, save_as_png = FALSE)

  expect_equal(ci$ci_image, renderLatent(fx$generator, ci$latent_ci)[1, , ])
  expect_equal(ci$base_image, renderLatent(fx$generator, ci$base_latent)[1, , ])
  expect_equal(dim(ci$ci_image), c(16, 16))
})

test_that("sd scaling moves the face a stated number of standard deviations", {
  withr::local_options(rcicr.experimental = TRUE)
  fx <- latent_fixture()

  for (constant in c(1, 2, 4)) {
    ci <- generateLatentCI(stimuli = 1:20, responses = rep(c(1, -1), 10),
                           rdata = fx$stimuli$rdata, save_as_png = FALSE,
                           latent_scaling = 'sd', scaling_constant = constant)
    in_sd <- sqrt(mean((ci$scaled_direction / fx$generator$latent_sd)^2))
    expect_equal(in_sd, constant)
  }
})

test_that("sd scaling makes magnitudes incomparable and constant scaling keeps them", {
  withr::local_options(rcicr.experimental = TRUE)
  fx <- latent_fixture(n_trials = 200, n_components = 4)

  consistent <- simulate_latent_observer(
    fx$stimuli$latent_params, fx$stimuli$base_latent,
    fx$stimuli$base_latent + fx$generator$latent_sd * c(2, -1, 1, 0)
  )
  noisy <- withr::with_seed(4, sample(c(-1, 1), 200, replace = TRUE))

  scaled_norm <- function(responses, scaling, constant) {
    ci <- generateLatentCI(stimuli = seq_len(200), responses = responses,
                           rdata = fx$stimuli$rdata, save_as_png = FALSE,
                           latent_scaling = scaling, scaling_constant = constant)
    sqrt(sum(ci$scaled_direction^2))
  }

  # Under sd scaling both come out the same size, which is why the picture
  # cannot be read as evidence of consistency.
  expect_equal(scaled_norm(consistent, 'sd', 2), scaled_norm(noisy, 'sd', 2),
               tolerance = 0.3)

  # Under constant scaling the consistent observer's direction stays longer.
  expect_gt(scaled_norm(consistent, 'constant', 0.1), scaled_norm(noisy, 'constant', 0.1))
})

test_that("none and constant scaling do what they say", {
  withr::local_options(rcicr.experimental = TRUE)
  fx <- latent_fixture()
  responses <- rep(c(1, -1), 10)

  raw <- generateLatentCI(stimuli = 1:20, responses = responses, rdata = fx$stimuli$rdata,
                          save_as_png = FALSE, latent_scaling = 'none')
  divided <- generateLatentCI(stimuli = 1:20, responses = responses, rdata = fx$stimuli$rdata,
                              save_as_png = FALSE, latent_scaling = 'constant',
                              scaling_constant = 0.1)

  expect_equal(raw$scaled_direction, raw$direction)
  expect_equal(divided$scaled_direction, raw$direction / 0.1)
})

test_that("an unknown scaling warns and leaves the direction alone", {
  withr::local_options(rcicr.experimental = TRUE)
  fx <- latent_fixture()

  expect_warning(
    ci <- generateLatentCI(stimuli = 1:20, responses = rep(c(1, -1), 10),
                           rdata = fx$stimuli$rdata, save_as_png = FALSE,
                           latent_scaling = 'independent'),
    'Unknown latent_scaling'
  )
  expect_equal(ci$scaled_direction, ci$direction)
})

test_that("participants are averaged rather than pooled", {
  withr::local_options(rcicr.experimental = TRUE)
  fx <- latent_fixture(n_trials = 30)

  # One participant ran 25 trials, the other 5. Pooling would let the first
  # decide the answer almost alone.
  participants <- c(rep('a', 25), rep('b', 5))
  responses <- c(rep(1, 25), rep(-1, 5))

  ci <- generateLatentCI(stimuli = 1:30, responses = responses,
                         participants = participants,
                         rdata = fx$stimuli$rdata, save_as_png = FALSE)

  expect_equal(dim(ci$pid_directions), c(2, 3))
  expect_equal(rownames(ci$pid_directions), c('a', 'b'))
  expect_equal(ci$direction, colMeans(ci$pid_directions))

  a <- colMeans(fx$stimuli$latent_params[1:25, ])
  b <- -colMeans(fx$stimuli$latent_params[26:30, ])
  expect_equal(ci$pid_directions[1, ], a)
  expect_equal(ci$pid_directions[2, ], b)
})

test_that("repeated presentations of a stimulus are averaged when pooling", {
  withr::local_options(rcicr.experimental = TRUE)
  fx <- latent_fixture()

  # Stimulus 1 seen twice with opposite answers cancels, leaving stimulus 2.
  ci <- generateLatentCI(stimuli = c(1, 1, 2), responses = c(1, -1, 1),
                         rdata = fx$stimuli$rdata, save_as_png = FALSE,
                         latent_scaling = 'none')

  expected <- colMeans(rbind(fx$stimuli$latent_params[1, ] * 0,
                             fx$stimuli$latent_params[2, ] * 1))
  expect_equal(ci$direction, expected)
})

test_that("a PNG is written, named, and prefixed for an anti-CI", {
  withr::local_options(rcicr.experimental = TRUE)
  fx <- latent_fixture()
  target <- withr::local_tempdir()

  generateLatentCI(stimuli = 1:20, responses = rep(c(1, -1), 10),
                   rdata = fx$stimuli$rdata, targetpath = target, filename = 'demo')
  generateLatentCI(stimuli = 1:20, responses = rep(c(1, -1), 10),
                   rdata = fx$stimuli$rdata, targetpath = target, filename = 'demo',
                   antiCI = TRUE)

  expect_true(file.exists(file.path(target, 'ci_demo.png')))
  expect_true(file.exists(file.path(target, 'antici_demo.png')))
  expect_equal(dim(png::readPNG(file.path(target, 'ci_demo.png'))), c(16, 16))
})

test_that("targetpath is required unless nothing is written", {
  withr::local_options(rcicr.experimental = TRUE)
  fx <- latent_fixture()

  expect_error(generateLatentCI(stimuli = 1:20, responses = rep(c(1, -1), 10),
                                rdata = fx$stimuli$rdata),
               'No targetpath')
})

test_that("a supplied generator must be the one that made the stimuli", {
  withr::local_options(rcicr.experimental = TRUE)
  fx <- latent_fixture()
  other_dir <- withr::local_tempdir()
  other <- latentGeneratorPCA(make_face_set(other_dir, n = 7), n_components = 3,
                              img_size = 16)

  expect_silent(generateLatentCI(stimuli = 1:20, responses = rep(c(1, -1), 10),
                                 rdata = fx$stimuli$rdata, generator = fx$generator,
                                 save_as_png = FALSE))
  expect_error(generateLatentCI(stimuli = 1:20, responses = rep(c(1, -1), 10),
                                rdata = fx$stimuli$rdata, generator = other,
                                save_as_png = FALSE),
               'not the one that made these stimuli')
})

test_that("a file whose generator cannot travel asks for the generator back", {
  withr::local_options(rcicr.experimental = TRUE)
  fx <- latent_fixture()

  spec <- generatorSpec(fx$generator)
  spec$portable <- FALSE
  spec$state <- NULL
  spec$kind <- 'stylegan'
  mutate_rdata(fx$stimuli$rdata, generator_spec = spec)

  err <- expect_error(generateLatentCI(stimuli = 1:20, responses = rep(c(1, -1), 10),
                                       rdata = fx$stimuli$rdata, save_as_png = FALSE))
  expect_match(conditionMessage(err), 'too ')
  expect_match(conditionMessage(err), 'generator = <your generator>', fixed = TRUE)
})

test_that("stimulus numbers outside the file are rejected", {
  withr::local_options(rcicr.experimental = TRUE)
  fx <- latent_fixture()

  expect_error(generateLatentCI(stimuli = c(1, 99), responses = c(1, -1),
                                rdata = fx$stimuli$rdata, save_as_png = FALSE),
               'between 1 and 20')
})

test_that("mismatched stimuli and responses are rejected", {
  withr::local_options(rcicr.experimental = TRUE)
  fx <- latent_fixture()

  # Valid codes, wrong count, so this tests the length check rather than
  # tripping the response-coding one on the way past.
  expect_error(generateLatentCI(stimuli = 1:20, responses = rep(1, 5),
                                rdata = fx$stimuli$rdata, save_as_png = FALSE),
               'same length')
})

test_that("informational value refuses a classification image from other stimuli", {
  withr::local_options(rcicr.experimental = TRUE)
  fx <- latent_fixture()
  # A different face set, so the second generator is genuinely a different one:
  # the same images would rebuild the same generator and match legitimately.
  other <- latent_fixture(n_faces = 9)

  ci <- generateLatentCI(stimuli = 1:20, responses = rep(c(1, -1), 10),
                         rdata = fx$stimuli$rdata, save_as_png = FALSE)

  expect_error(computeLatentInfoVal2IFC(ci, other$stimuli$rdata, 1:20,
                                        rep(c(1, -1), 10), iter = 10),
               'not computed from these inputs')
})

test_that("the null is reproducible under a response_seed and free without one", {
  withr::local_options(rcicr.experimental = TRUE)
  fx <- latent_fixture()

  ci <- generateLatentCI(stimuli = 1:20, responses = rep(c(1, -1), 10),
                         rdata = fx$stimuli$rdata, save_as_png = FALSE)

  first <- computeLatentInfoVal2IFC(ci, fx$stimuli$rdata, 1:20, rep(c(1, -1), 10),
                                    iter = 200, response_seed = 5)
  again <- computeLatentInfoVal2IFC(ci, fx$stimuli$rdata, 1:20, rep(c(1, -1), 10),
                                    iter = 200, response_seed = 5)
  expect_equal(first$infoVal, again$infoVal)
  expect_equal(first$reference_median, again$reference_median)

  expect_gte(first$p, 0)
  expect_lte(first$p, 1)
  expect_equal(first$iter, 200)
})

test_that("the null preserves the responses the participant actually gave", {
  withr::local_options(rcicr.experimental = TRUE)
  fx <- latent_fixture(n_trials = 60)

  # A participant who pressed one key throughout carries no information about
  # the stimuli at all: every permutation of their responses is the same vector,
  # so the observed direction has to sit in the middle of its own null. Drawing
  # fresh balanced responses instead builds a null this participant could never
  # have produced, and reports their key-press bias as signal.
  biased <- rep(1, 60)
  ci <- generateLatentCI(stimuli = 1:60, responses = biased,
                         rdata = fx$stimuli$rdata, save_as_png = FALSE)

  result <- computeLatentInfoVal2IFC(ci, fx$stimuli$rdata, 1:60, biased,
                                     iter = 300, response_seed = 1)

  expect_equal(result$p, 1)
  expect_equal(result$observed_norm, result$reference_median)
})

test_that("a lopsided but informative participant is still scored fairly", {
  withr::local_options(rcicr.experimental = TRUE)
  fx <- latent_fixture(n_trials = 200, n_components = 4)

  # Three quarters of one answer, but the answers still track the stimuli. The
  # null holds that imbalance fixed, so what is left to measure is the pairing.
  target <- fx$stimuli$base_latent + fx$generator$latent_sd * c(2, -1, 1, 0)
  honest <- simulate_latent_observer(fx$stimuli$latent_params,
                                     fx$stimuli$base_latent, target)
  lopsided <- honest
  lopsided[seq(2, 200, by = 4)] <- 1

  ci <- generateLatentCI(stimuli = seq_len(200), responses = lopsided,
                         rdata = fx$stimuli$rdata, save_as_png = FALSE)
  result <- computeLatentInfoVal2IFC(ci, fx$stimuli$rdata, seq_len(200), lopsided,
                                     iter = 300, response_seed = 1)

  expect_lt(result$p, 0.05)
})

test_that("informational value is tied to the stimulus set, not just the generator", {
  withr::local_options(rcicr.experimental = TRUE)
  fx <- latent_fixture()
  out <- withr::local_tempdir()

  # The same generator, a different seed. The generator fingerprints match by
  # design, so only a stimulus-set fingerprint can tell these apart, and without
  # one the direction gets scored against unrelated perturbations.
  other <- generateStimuliLatent2IFC(fx$generator, n_trials = 20,
                                     stimulus_path = out, seed = 99,
                                     save_as_png = FALSE)

  ci <- generateLatentCI(stimuli = 1:20, responses = rep(c(1, -1), 10),
                         rdata = fx$stimuli$rdata, save_as_png = FALSE)

  expect_error(computeLatentInfoVal2IFC(ci, other$rdata, 1:20, rep(c(1, -1), 10),
                                        iter = 10),
               'not computed from these inputs')
})

test_that("the null shuffles trials rather than collapsed means", {
  withr::local_options(rcicr.experimental = TRUE)
  fx <- latent_fixture()

  # Codex's case. Stimulus 1 is seen twice, so pooling collapses it to a mean,
  # and permuting the collapsed means is a different null from permuting the
  # answers and collapsing afterwards. Here the two do not even overlap.
  stimuli <- c(1, 1, 2)
  responses <- c(1, 1, -1)

  ci <- generateLatentCI(stimuli = stimuli, responses = responses,
                         rdata = fx$stimuli$rdata, save_as_png = FALSE)
  result <- computeLatentInfoVal2IFC(ci, fx$stimuli$rdata, stimuli, responses,
                                     iter = 100, response_seed = 1)

  params <- fx$stimuli$latent_params
  latent_sd <- fx$generator$latent_sd
  norm_of <- function(weights) {
    sqrt(mean((colSums(params[c(1, 2), ] * weights) / length(weights) / latent_sd)^2))
  }

  # Every arrangement of the three answers, collapsed the way pooling collapses
  # them: (1,1,-1) gives stimulus 1 a mean of 1, the other two give it 0.
  trial_level <- unique(vapply(
    list(c(1, 1, -1), c(1, -1, 1), c(-1, 1, 1)),
    function(r) norm_of(c(mean(r[1:2]), r[3])), numeric(1)
  ))

  # Permuting the collapsed means instead: (1, -1) and (-1, 1).
  mean_level <- vapply(list(c(1, -1), c(-1, 1)), norm_of, numeric(1))

  observed <- unique(round(result$reference_median, 10))
  expect_true(all(round(observed, 10) %in% round(trial_level, 10)))
  expect_false(any(round(observed, 10) %in% round(setdiff(mean_level, trial_level), 10)))
})

test_that("a per-participant image is scored against a per-participant null", {
  withr::local_options(rcicr.experimental = TRUE)
  fx <- latent_fixture(n_trials = 40)

  # Unequal trial counts. Pooling the null while the classification image
  # averaged per participant compares the observed direction against one built
  # by a different estimator, and lets the participant who ran more trials
  # decide the answer.
  participants <- c(rep('a', 32), rep('b', 8))
  responses <- rep(c(1, -1), 20)

  ci <- generateLatentCI(stimuli = 1:40, responses = responses,
                         participants = participants,
                         rdata = fx$stimuli$rdata, save_as_png = FALSE)

  result <- computeLatentInfoVal2IFC(ci, fx$stimuli$rdata, 1:40, responses,
                                     participants, iter = 200, response_seed = 1)
  expect_true(is.finite(result$infoVal))

  # The same trials analysed pooled are a different analysis, and must not be
  # scored against this image.
  expect_error(computeLatentInfoVal2IFC(ci, fx$stimuli$rdata, 1:40, responses,
                                        iter = 10),
               'not computed from these inputs')
})

test_that("informational value is tied to the trials and answers analysed", {
  withr::local_options(rcicr.experimental = TRUE)
  fx <- latent_fixture()

  ci <- generateLatentCI(stimuli = 1:20, responses = rep(c(1, -1), 10),
                         rdata = fx$stimuli$rdata, save_as_png = FALSE)

  # Same file, same generator, different answers.
  expect_error(computeLatentInfoVal2IFC(ci, fx$stimuli$rdata, 1:20,
                                        rep(c(-1, 1), 10), iter = 10),
               'not computed from these inputs')

  # Same file, same answers, a different subset of its trials.
  expect_error(computeLatentInfoVal2IFC(ci, fx$stimuli$rdata, 1:10,
                                        rep(c(1, -1), 5), iter = 10),
               'not computed from these inputs')
})

test_that("the null shuffles within each participant, not across them", {
  withr::local_options(rcicr.experimental = TRUE)
  fx <- latent_fixture(n_trials = 40)

  # One participant pressed one key throughout, the other pressed the other.
  # Every within-participant permutation is the identity, so the observed
  # direction has to sit exactly at its own null. A global shuffle would mix
  # their answers into vectors neither of them could have given and report the
  # difference between their key biases as signal.
  participants <- c(rep('a', 20), rep('b', 20))
  responses <- c(rep(1, 20), rep(-1, 20))

  ci <- generateLatentCI(stimuli = 1:40, responses = responses,
                         participants = participants,
                         rdata = fx$stimuli$rdata, save_as_png = FALSE)
  result <- computeLatentInfoVal2IFC(ci, fx$stimuli$rdata, 1:40, responses,
                                     participants, iter = 200, response_seed = 1)

  expect_equal(result$p, 1)
  expect_equal(result$observed_norm, result$reference_median)
})

test_that("a degenerate null gives a readable informational value", {
  withr::local_options(rcicr.experimental = TRUE)
  fx <- latent_fixture(n_trials = 30)

  # One arrangement of the responses means the null has no spread, so the
  # standardised score is 0 / 0. The limit is reported instead, because NaN in
  # the headline number is worse than a defined zero for exactly the case the
  # permutation null was rewritten to handle.
  biased <- rep(1, 30)
  ci <- generateLatentCI(stimuli = 1:30, responses = biased,
                         rdata = fx$stimuli$rdata, save_as_png = FALSE)
  result <- computeLatentInfoVal2IFC(ci, fx$stimuli$rdata, 1:30, biased,
                                     iter = 100, response_seed = 1)

  expect_equal(result$reference_mad, 0)
  expect_false(is.nan(result$infoVal))
  expect_equal(result$infoVal, 0)
})

test_that("responses that are not 1 or -1 are refused", {
  withr::local_options(rcicr.experimental = TRUE)
  fx <- latent_fixture()

  # 0/1 is what plenty of experiment software writes, and it errors nowhere
  # downstream: it silently turns the estimator into something meaningless and
  # still renders a plausible face.
  expect_error(generateLatentCI(stimuli = 1:20, responses = rep(c(0, 1), 10),
                                rdata = fx$stimuli$rdata, save_as_png = FALSE),
               'takes responses coded 1')
  expect_error(generateLatentCI(stimuli = 1:20,
                                responses = c(rep(c(1, -1), 9), NA, 1),
                                rdata = fx$stimuli$rdata, save_as_png = FALSE),
               'takes responses coded 1')
})

test_that("responses that agree on every summary statistic are still told apart", {
  withr::local_options(rcicr.experimental = TRUE)
  fixture <- latent_fixture(n_trials = 6)

  # Six statistics over a vector of 1s and -1s are close to no information: the
  # count of each is fixed by the sum, so only the index-weighted sum separates
  # the arrangements. These two agree on all six, so a summary-statistic
  # fingerprint accepts one classification image against the other's responses
  # and scores it against a null it has nothing to do with.
  first <- c(1, -1, -1, 1, -1, -1)
  second <- c(-1, 1, 1, -1, -1, -1)

  summarise <- function(x) {
    c(length(x), sum(x), sum(x^2), sum(x * seq_along(x)), min(x), max(x))
  }
  expect_equal(summarise(first), summarise(second))
  expect_false(identical(first, second))

  ci <- generateLatentCI(1:6, first, fixture$stimuli$rdata, save_as_png = FALSE)

  expect_error(
    computeLatentInfoVal2IFC(ci, fixture$stimuli$rdata, 1:6, second, iter = 5),
    'not computed from these inputs'
  )
  expect_no_error(
    computeLatentInfoVal2IFC(ci, fixture$stimuli$rdata, 1:6, first, iter = 5)
  )
})

test_that("a participant vector shorter than the trials is refused", {
  withr::local_options(rcicr.experimental = TRUE)
  fixture <- latent_fixture(n_trials = 20)
  responses <- rep(c(1, -1), 10)

  # Two identifiers for twenty trials recycles through `participants == pid`,
  # dealing the trials out alternately and producing a group classification
  # image over participants nobody ran.
  expect_error(
    generateLatentCI(1:20, responses, fixture$stimuli$rdata,
                     participants = c('a', 'b'), save_as_png = FALSE),
    'one participant identifier per trial'
  )
  expect_error(
    computeLatentInfoVal2IFC(list(), fixture$stimuli$rdata, 1:20, responses,
                             participants = c('a', 'b'), iter = 5),
    'one participant identifier per trial'
  )

  # One per trial is still accepted, and so is the pooling default.
  expect_no_error(
    generateLatentCI(1:20, responses, fixture$stimuli$rdata,
                     participants = rep(c('a', 'b'), 10), save_as_png = FALSE)
  )
  expect_no_error(
    generateLatentCI(1:20, responses, fixture$stimuli$rdata, save_as_png = FALSE)
  )
})
