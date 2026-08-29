# Latent stimulus generation, and the .Rdata contract it writes.

test_that("stimuli and an .Rdata file are written, one pair per trial", {
  withr::local_options(rcicr.experimental = TRUE)
  dir <- withr::local_tempdir()
  out <- withr::local_tempdir()

  generator <- latentGeneratorPCA(make_face_set(dir), n_components = 3, img_size = 16)
  result <- generateStimuliLatent2IFC(generator, n_trials = 4, stimulus_path = out, seed = 1)

  expect_equal(length(list.files(out, pattern = '_ori\\.png$')), 4)
  expect_equal(length(list.files(out, pattern = '_inv\\.png$')), 4)
  expect_true(file.exists(result$rdata))
  expect_equal(dim(result$latent_params), c(4, 3))
})

test_that("the trial number in a filename indexes latent_params", {
  withr::local_options(rcicr.experimental = TRUE)
  dir <- withr::local_tempdir()
  out <- withr::local_tempdir()

  generator <- latentGeneratorPCA(make_face_set(dir), n_components = 3, img_size = 16)
  result <- generateStimuliLatent2IFC(generator, n_trials = 3, stimulus_path = out,
                                      label = 'demo', seed = 7)

  # Trial 2's file has to be the render of base_latent + latent_params[2, ], or
  # every classification image computed from these stimuli is built on the
  # wrong rows.
  written <- png::readPNG(file.path(out, 'demo_7_00002_ori.png'))
  expected <- renderLatent(generator, result$base_latent + result$latent_params[2, ])[1, , ]
  expect_equal(written, expected, tolerance = 1 / 255)

  inverted <- png::readPNG(file.path(out, 'demo_7_00002_inv.png'))
  expect_equal(inverted,
               renderLatent(generator, result$base_latent - result$latent_params[2, ])[1, , ],
               tolerance = 1 / 255)
})

test_that("the same seed and generator reproduce the same stimuli", {
  withr::local_options(rcicr.experimental = TRUE)
  dir <- withr::local_tempdir()

  generator <- latentGeneratorPCA(make_face_set(dir), n_components = 3, img_size = 16)

  first <- generateStimuliLatent2IFC(generator, n_trials = 5, stimulus_path = withr::local_tempdir(),
                                     seed = 42, save_as_png = FALSE)
  again <- generateStimuliLatent2IFC(generator, n_trials = 5, stimulus_path = withr::local_tempdir(),
                                     seed = 42, save_as_png = FALSE)
  expect_identical(first$latent_params, again$latent_params)

  different <- generateStimuliLatent2IFC(generator, n_trials = 5, stimulus_path = withr::local_tempdir(),
                                         seed = 43, save_as_png = FALSE)
  expect_false(identical(first$latent_params, different$latent_params))
})

test_that("latent_sigma scales the perturbations in units of latent_sd", {
  withr::local_options(rcicr.experimental = TRUE)
  dir <- withr::local_tempdir()
  out <- withr::local_tempdir()

  generator <- latentGeneratorPCA(make_face_set(dir, n = 10), n_components = 4, img_size = 16)

  one <- generateStimuliLatent2IFC(generator, n_trials = 200, stimulus_path = out,
                                   latent_sigma = 1, seed = 3, save_as_png = FALSE,
                                   save_rdata = FALSE)
  half <- generateStimuliLatent2IFC(generator, n_trials = 200, stimulus_path = out,
                                    latent_sigma = 0.5, seed = 3, save_as_png = FALSE,
                                    save_rdata = FALSE)

  expect_equal(half$latent_params, one$latent_params / 2)
  # sigma = 1 means "about as far as the training faces vary", so the sampled
  # spread has to track latent_sd rather than being unit-scaled.
  expect_equal(apply(one$latent_params, 2, stats::sd), generator$latent_sd,
               tolerance = 0.25)
})

test_that("the base latent defaults to the centre of the space and can be moved", {
  withr::local_options(rcicr.experimental = TRUE)
  dir <- withr::local_tempdir()
  out <- withr::local_tempdir()

  generator <- latentGeneratorPCA(make_face_set(dir), n_components = 3, img_size = 16)

  centred <- generateStimuliLatent2IFC(generator, n_trials = 2, stimulus_path = out,
                                       save_as_png = FALSE, save_rdata = FALSE)
  expect_equal(centred$base_latent, generator$latent_mean)

  moved <- generateStimuliLatent2IFC(generator, n_trials = 2, stimulus_path = out,
                                     base_latent = c(1, 2, 3), save_as_png = FALSE,
                                     save_rdata = FALSE)
  expect_equal(moved$base_latent, c(1, 2, 3))
})

test_that("a base latent of the wrong width is rejected", {
  withr::local_options(rcicr.experimental = TRUE)
  dir <- withr::local_tempdir()

  generator <- latentGeneratorPCA(make_face_set(dir), n_components = 3, img_size = 16)
  expect_error(
    generateStimuliLatent2IFC(generator, n_trials = 2, stimulus_path = withr::local_tempdir(),
                              base_latent = c(1, 2), save_as_png = FALSE, save_rdata = FALSE),
    'latent_dim is 3'
  )
})

test_that("stimulus_path is required when anything is written", {
  withr::local_options(rcicr.experimental = TRUE)
  dir <- withr::local_tempdir()

  generator <- latentGeneratorPCA(make_face_set(dir), n_components = 3, img_size = 16)

  expect_error(generateStimuliLatent2IFC(generator, n_trials = 2), 'No stimulus_path')
  expect_silent(generateStimuliLatent2IFC(generator, n_trials = 2, save_as_png = FALSE,
                                          save_rdata = FALSE))
})

test_that("batch_size changes nothing about the stimuli", {
  withr::local_options(rcicr.experimental = TRUE)
  dir <- withr::local_tempdir()

  generator <- latentGeneratorPCA(make_face_set(dir), n_components = 3, img_size = 16)

  one_at_a_time <- withr::local_tempdir()
  all_at_once <- withr::local_tempdir()
  generateStimuliLatent2IFC(generator, n_trials = 5, stimulus_path = one_at_a_time,
                            seed = 2, batch_size = 1, save_rdata = FALSE)
  generateStimuliLatent2IFC(generator, n_trials = 5, stimulus_path = all_at_once,
                            seed = 2, batch_size = 100, save_rdata = FALSE)

  for (f in list.files(one_at_a_time)) {
    expect_equal(png::readPNG(file.path(one_at_a_time, f)),
                 png::readPNG(file.path(all_at_once, f)), info = f)
  }
})

test_that("the .Rdata carries the fields the analysis half reads", {
  withr::local_options(rcicr.experimental = TRUE)
  dir <- withr::local_tempdir()
  out <- withr::local_tempdir()

  generator <- latentGeneratorPCA(make_face_set(dir), n_components = 3, img_size = 16)
  result <- generateStimuliLatent2IFC(generator, n_trials = 4, stimulus_path = out,
                                      latent_sigma = 0.5, seed = 5, save_as_png = FALSE)

  env <- new.env()
  load(result$rdata, envir = env)
  expect_setequal(ls(env), c('latent_params', 'base_latent', 'latent_sigma',
                             'generator_spec', 'n_trials', 'img_size', 'seed',
                             'label', 'stimulus_path', 'generator_version'))
  expect_equal(env$latent_sigma, 0.5)
  expect_equal(env$img_size, 16L)
})

test_that("a PCA generator travels inside the file and rebuilds identically", {
  withr::local_options(rcicr.experimental = TRUE)
  dir <- withr::local_tempdir()
  out <- withr::local_tempdir()

  generator <- latentGeneratorPCA(make_face_set(dir), n_components = 3, img_size = 16)
  result <- generateStimuliLatent2IFC(generator, n_trials = 2, stimulus_path = out,
                                      save_as_png = FALSE)

  stored <- loadLatentStimulusParams(result$rdata)
  rebuilt <- generatorFromSpec(stored$generator_spec)

  expect_identical(rebuilt$fingerprint, generator$fingerprint)
  probe <- matrix(c(1, -1, 0.5), nrow = 1)
  expect_equal(renderLatent(rebuilt, probe), renderLatent(generator, probe))
})

test_that("a generator that cannot travel leaves a fingerprint and no renderer", {
  withr::local_options(rcicr.experimental = TRUE)

  external <- rcicrGenerator(
    kind = 'stub', latent_dim = 2, img_size = 4, space = 'w',
    latent_mean = c(0, 0), latent_sd = c(1, 1),
    render = function(latents) array(0.5, dim = c(nrow(latents), 4, 4)),
    fingerprint = 'stub:abc'
  )
  spec <- generatorSpec(external)

  expect_false(spec$portable)
  expect_null(spec$state)
  expect_identical(spec$fingerprint, 'stub:abc')
  expect_null(generatorFromSpec(spec))
})

test_that("a mismatched generator is refused rather than rendering another face", {
  withr::local_options(rcicr.experimental = TRUE)
  dir <- withr::local_tempdir()
  other_dir <- withr::local_tempdir()

  generator <- latentGeneratorPCA(make_face_set(dir), n_components = 3, img_size = 16)
  other <- latentGeneratorPCA(make_face_set(other_dir, n = 7), n_components = 3, img_size = 16)

  spec <- generatorSpec(generator)
  expect_true(matchGenerator(generator, spec, 'generateLatentCI'))

  err <- expect_error(matchGenerator(other, spec, 'generateLatentCI'))
  expect_match(conditionMessage(err), 'not the one that made these stimuli')
})

test_that("the loader names the field a truncated file is missing", {
  withr::local_options(rcicr.experimental = TRUE)
  dir <- withr::local_tempdir()
  out <- withr::local_tempdir()

  generator <- latentGeneratorPCA(make_face_set(dir), n_components = 3, img_size = 16)
  result <- generateStimuliLatent2IFC(generator, n_trials = 2, stimulus_path = out,
                                      save_as_png = FALSE)

  mutate_rdata(result$rdata, .remove = 'base_latent')
  expect_error(loadLatentStimulusParams(result$rdata), 'did not contain base_latent')
})

test_that("a pixel-noise file loaded as a latent one says which function to use", {
  dir <- withr::local_tempdir()
  rdata <- make_fixture_rdata(dir)

  expect_error(loadLatentStimulusParams(rdata), 'pixel-noise stimulus file')
})

test_that("loading a latent file cannot overwrite the caller's variables", {
  withr::local_options(rcicr.experimental = TRUE)
  dir <- withr::local_tempdir()
  out <- withr::local_tempdir()

  generator <- latentGeneratorPCA(make_face_set(dir), n_components = 3, img_size = 16)
  result <- generateStimuliLatent2IFC(generator, n_trials = 2, stimulus_path = out,
                                      label = 'kept', seed = 9, save_as_png = FALSE)

  # load() assigns into the calling frame, and every saved name below is also an
  # argument name somewhere in the package. loadLatentStimulusParams() loads
  # into an isolated environment instead, so these bindings have to survive the
  # call unchanged. This is the guard that three separate bugs came from
  # lacking on the pixel-noise side.
  img_size <- 'mine'
  seed <- 'mine'
  label <- 'mine'
  n_trials <- 'mine'
  stimulus_path <- 'mine'
  base_latent <- 'mine'
  latent_sigma <- 'mine'
  latent_params <- 'mine'
  generator_spec <- 'mine'
  generator_version <- 'mine'

  stored <- loadLatentStimulusParams(result$rdata)

  expect_equal(dim(stored$latent_params), c(2, 3))
  for (nm in c('img_size', 'seed', 'label', 'n_trials', 'stimulus_path',
               'base_latent', 'latent_sigma', 'latent_params', 'generator_spec',
               'generator_version')) {
    expect_identical(get(nm), 'mine', info = nm)
  }
})

test_that("two stimulus sets in one directory do not overwrite each other", {
  withr::local_options(rcicr.experimental = TRUE)
  dir <- withr::local_tempdir()
  generator <- latentGeneratorPCA(make_face_set(dir, n = 6), n_components = 3,
                                  img_size = 16)
  out <- withr::local_tempdir()

  # Same label and seed, made well within one minute of each other, differing
  # only in latent_sigma. A timestamp resolved to the minute gave both the same
  # filename, so the first call's returned path pointed at the second call's
  # perturbations and responses to the first set would be analysed against them.
  first <- generateStimuliLatent2IFC(generator, n_trials = 4, stimulus_path = out,
                                     latent_sigma = 1, seed = 1, save_as_png = FALSE)
  second <- generateStimuliLatent2IFC(generator, n_trials = 4, stimulus_path = out,
                                      latent_sigma = 3, seed = 1, save_as_png = FALSE)

  expect_false(first$rdata == second$rdata)
  expect_true(file.exists(first$rdata))
  expect_true(file.exists(second$rdata))

  # And the first file still holds the first call's own perturbations.
  stored <- loadLatentStimulusParams(first$rdata)
  expect_equal(stored$latent_params, first$latent_params)
  expect_equal(stored$latent_sigma, 1)
})
