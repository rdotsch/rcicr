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

  expect_error(generateLatentCI(stimuli = 1:20, responses = 1:5,
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

  expect_error(computeLatentInfoVal2IFC(ci, other$stimuli$rdata, 1:20, iter = 10),
               'not computed from these stimuli')
})

test_that("the null is reproducible under a response_seed and free without one", {
  withr::local_options(rcicr.experimental = TRUE)
  fx <- latent_fixture()

  ci <- generateLatentCI(stimuli = 1:20, responses = rep(c(1, -1), 10),
                         rdata = fx$stimuli$rdata, save_as_png = FALSE)

  first <- computeLatentInfoVal2IFC(ci, fx$stimuli$rdata, 1:20, iter = 200, response_seed = 5)
  again <- computeLatentInfoVal2IFC(ci, fx$stimuli$rdata, 1:20, iter = 200, response_seed = 5)
  expect_equal(first$infoVal, again$infoVal)
  expect_equal(first$reference_median, again$reference_median)

  expect_gte(first$p, 0)
  expect_lte(first$p, 1)
  expect_equal(first$iter, 200)
})
